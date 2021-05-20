//
// Created by mlaupichler on 17.05.21.
//

#ifndef KAHYPAR_LOCAL_GAIN_POLICY_H
#define KAHYPAR_LOCAL_GAIN_POLICY_H


#include <vector>

#include "tbb/enumerable_thread_specific.h"

#include "kahypar/meta/policy_registry.h"
#include "kahypar/meta/typelist.h"
#include "kahypar/partition/metrics.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

    /// Same as GainPolicy, except made for only one thread at a time, i.e. without thread local storage
    template <class Derived = Mandatory,
            class HyperGraph = Mandatory>
    class LocalGainPolicy : public kahypar::meta::PolicyBase {
        using DeltaGain = Gain;
        using TmpScores = parallel::scalable_vector<Gain>;

    public:
        LocalGainPolicy(const Context& context) :
                _context(context),
                _delta(0),
                _tmp_scores(context.partition.k, 0) { }

        Move computeMaxGainMove(HyperGraph& hypergraph,
                                const HypernodeID hn,
                                const bool rebalance = false) {
            return static_cast<Derived*>(this)->computeMaxGainMoveImpl(hypergraph, hn, rebalance);
        }

        inline void computeDeltaForHyperedge(const HyperedgeID he,
                                             const HyperedgeWeight edge_weight,
                                             const HypernodeID edge_size,
                                             const HypernodeID pin_count_in_from_part_after,
                                             const HypernodeID pin_count_in_to_part_after) {
            static_cast<Derived*>(this)->computeDeltaForHyperedgeImpl(
                    he, edge_weight, edge_size,
                    pin_count_in_from_part_after, pin_count_in_to_part_after);
        }

        // ! Returns the delta in the objective function for all moves
        // ! performed relative to the last call to reset()
        Gain localDelta() const {
            return _delta;
        }

        // ! Returns the delta of all moves performed relative to the last call of reset()
        Gain delta() const {
            return localDelta();
        }

        void reset() {
            _delta = 0;
        }

    protected:
        const Context& _context;
        DeltaGain _delta;
        TmpScores _tmp_scores;
    };

    template <class HyperGraph = Mandatory>
    class LocalKm1Policy : public LocalGainPolicy<LocalKm1Policy<HyperGraph>, HyperGraph> {
        using Base = LocalGainPolicy<LocalKm1Policy<HyperGraph>, HyperGraph>;

        static constexpr bool enable_heavy_assert = false;

    public:
        LocalKm1Policy(const Context& context,
                  bool disable_randomization = false) :
                Base(context),
                _disable_randomization(disable_randomization) { }

        Move computeMaxGainMoveImpl(HyperGraph& hypergraph,
                                    const HypernodeID hn,
                                    const bool rebalance) {
            HEAVY_REFINEMENT_ASSERT([&] {
                for (PartitionID k = 0; k < _context.partition.k; ++k) {
                    if (_tmp_scores[k] != 0) {
                        return false;
                    }
                }
                return true;
            } (), "Scores and valid parts not correctly reset");

            PartitionID from = hypergraph.partID(hn);
            Gain internal_weight = 0;
            for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
                HypernodeID pin_count_in_from_part = hypergraph.pinCountInPart(he, from);
                HyperedgeWeight he_weight = hypergraph.edgeWeight(he);

                // In case, there is more one than one pin left in from part, we would
                // increase the connectivity, if we would move the pin to one block
                // no contained in the connectivity set. In such cases, we can only
                // increase the connectivity of a hyperedge and therefore gather
                // the edge weight of all those edges and add it later to move gain
                // to all other blocks.
                if ( pin_count_in_from_part > 1 ) {
                    internal_weight += he_weight;
                }

                // Substract edge weight of all incident blocks.
                // Note, in case the pin count in from part is greater than one
                // we will later add that edge weight to the gain (see internal_weight).
                for (const PartitionID& to : hypergraph.connectivitySet(he)) {
                    if (from != to) {
                        _tmp_scores[to] -= he_weight;
                    }
                }
            }

            Move best_move { from, from, hn, rebalance ? std::numeric_limits<Gain>::max() : 0 };
            HypernodeWeight hn_weight = hypergraph.nodeWeight(hn);
            int cpu_id = sched_getcpu();
            utils::Randomize& rand = utils::Randomize::instance();
            for (PartitionID to = 0; to < _context.partition.k; ++to) {
                if (from != to) {
                    Gain score = _tmp_scores[to] + internal_weight;
                    bool new_best_gain = (score < best_move.gain) ||
                                         (score == best_move.gain &&
                                          !_disable_randomization &&
                                          rand.flipCoin(cpu_id));
                    if (new_best_gain &&
                        hypergraph.partWeight(to) + hn_weight <= _context.partition.max_part_weights[to]) {
                        best_move.to = to;
                        best_move.gain = score;
                    }
                }
                _tmp_scores[to] = 0;
            }
            return best_move;
        }

        inline void computeDeltaForHyperedgeImpl(const HyperedgeID he,
                                                 const HyperedgeWeight edge_weight,
                                                 const HypernodeID edge_size,
                                                 const HypernodeID pin_count_in_from_part_after,
                                                 const HypernodeID pin_count_in_to_part_after) {
            _delta += km1Delta(he, edge_weight, edge_size,
                                        pin_count_in_from_part_after,
                                        pin_count_in_to_part_after);
        }

        using Base::_context;
        using Base::_delta;
        using Base::_tmp_scores;
        bool _disable_randomization;
    };

    template <class HyperGraph = Mandatory>
    class LocalCutPolicy : public LocalGainPolicy<LocalCutPolicy<HyperGraph>, HyperGraph> {
        using Base = LocalGainPolicy<LocalCutPolicy<HyperGraph>, HyperGraph>;

        static constexpr bool enable_heavy_assert = false;

    public:
        LocalCutPolicy(const Context& context,
                  bool disable_randomization = false) :
                Base(context),
                _disable_randomization(disable_randomization) { }

        Move computeMaxGainMoveImpl(HyperGraph& hypergraph,
                                    const HypernodeID hn,
                                    const bool rebalance) {
            HEAVY_REFINEMENT_ASSERT([&] {
                for (PartitionID k = 0; k < _context.partition.k; ++k) {
                    if (_tmp_scores[k] != 0) {
                        return false;
                    }
                }
                return true;
            } (), "Scores and valid parts not correctly reset");

            PartitionID from = hypergraph.partID(hn);
            Gain internal_weight = 0;
            for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
                PartitionID connectivity = hypergraph.connectivity(he);
                HypernodeID pin_count_in_from_part = hypergraph.pinCountInPart(he, from);
                HyperedgeWeight weight = hypergraph.edgeWeight(he);
                if (connectivity == 1) {
                    ASSERT(hypergraph.edgeSize(he) > 1);
                    // In case, the hyperedge is a non-cut hyperedge, we would increase
                    // the cut, if we move vertex hn to an other block.
                    internal_weight += weight;
                } else if (connectivity == 2 && pin_count_in_from_part == 1) {
                    for (const PartitionID& to : hypergraph.connectivitySet(he)) {
                        // In case there are only two blocks contained in the current
                        // hyperedge and only one pin left in the from part of the hyperedge,
                        // we would make the current hyperedge a non-cut hyperedge when moving
                        // vertex hn to the other block.
                        if (from != to) {
                            _tmp_scores[to] -= weight;
                        }
                    }
                }
            }

            Move best_move { from, from, hn, rebalance ? std::numeric_limits<Gain>::max() : 0 };
            HypernodeWeight hn_weight = hypergraph.nodeWeight(hn);
            int cpu_id = sched_getcpu();
            utils::Randomize& rand = utils::Randomize::instance();
            for (PartitionID to = 0; to < _context.partition.k; ++to) {
                if (from != to) {
                    Gain score = _tmp_scores[to] + internal_weight;
                    bool new_best_gain = (score < best_move.gain) ||
                                         (score == best_move.gain &&
                                          !_disable_randomization &&
                                          rand.flipCoin(cpu_id));
                    if (new_best_gain && hypergraph.partWeight(to) + hn_weight <=
                                         _context.partition.max_part_weights[to]) {
                        best_move.to = to;
                        best_move.gain = score;
                    }
                }
                _tmp_scores[to] = 0;
            }
            return best_move;
        }

        inline void computeDeltaForHyperedgeImpl(const HyperedgeID he,
                                                 const HyperedgeWeight edge_weight,
                                                 const HypernodeID edge_size,
                                                 const HypernodeID pin_count_in_from_part_after,
                                                 const HypernodeID pin_count_in_to_part_after) {
            _delta += cutDelta(he, edge_weight, edge_size,
                                        pin_count_in_from_part_after,
                                        pin_count_in_to_part_after);
        }

        using Base::_context;
        using Base::_delta;
        using Base::_tmp_scores;
        bool _disable_randomization;
    };
}  // namespace mt_kahypar


#endif //KAHYPAR_LOCAL_GAIN_POLICY_H

//
// Created by mlaupichler on 04.05.21.
//

#ifndef KAHYPAR_ASYNC_LP_REFINER_H
#define KAHYPAR_ASYNC_LP_REFINER_H

#include <mt-kahypar/partition/refinement/i_refiner.h>
#include <mt-kahypar/partition/refinement/policies/local_gain_policy.h>
#include <mt-kahypar/datastructures/async/array_lock_manager.h>
#include <mt-kahypar/datastructures/thread_safe_fast_reset_flag_array.h>
#include <mt-kahypar/datastructures/async/thread_safe_flag_array.h>

namespace mt_kahypar {

    /// Label Propagation Refiner to be used in asynchronous uncoarsening for localized refinement. Uses global locking
    /// datastructure to allow multiple concurrent localized label propagations as well as concurrent uncontractions.
    /// Can not be used for global refinement or (only) rebalancing as it always requires seed nodes for refining!
    template <template <typename> class LocalGainPolicy> class AsyncLPRefiner : public IAsyncRefiner {

    private:
        using GainCalculator = LocalGainPolicy<PartitionedHypergraph>;
        using ActiveNodes = std::vector<HypernodeID>;
        using NextActiveNodes = std::vector<HypernodeID>;
        using VisitedEdges = std::vector<HyperedgeID>;
        using MovedNodes = ds::StreamingVector<HypernodeID>;

        static constexpr bool debug = false;
        static constexpr bool enable_heavy_assert = true;

    public:
        explicit AsyncLPRefiner(Hypergraph &hypergraph, const Context &context, const TaskGroupID task_group_id,
                                ds::GroupLockManager *lockManager,
                                ds::ThreadSafeFlagArray <HypernodeID> &node_anti_duplicator,
                                ds::ThreadSafeFlagArray <HyperedgeID> &edge_anti_duplicator) :
        _context(context),
        _task_group_id(task_group_id),
        _gain(context),
        _active_nodes(),
        _contraction_group_id(ds::invalidGroupID),
        _lock_manager(lockManager),
        _rng(),
        _next_active(node_anti_duplicator),
        _visited_he(edge_anti_duplicator) {
            ASSERT(_next_active.size() == hypergraph.initialNumNodes());
            ASSERT(_visited_he.size() == hypergraph.initialNumEdges());
        }

        AsyncLPRefiner(const AsyncLPRefiner&) = delete;
        AsyncLPRefiner(AsyncLPRefiner&&) = delete;

        AsyncLPRefiner & operator= (const AsyncLPRefiner &) = delete;
        AsyncLPRefiner & operator= (AsyncLPRefiner &&) = delete;

        ~AsyncLPRefiner() override {
            // Expect all AsyncLPRefiners between Refinement tasks (relevant for initial partitioning) to be built and destroyed
            // simultaneously, so when one gets destroyed, all should be done with any refinement calls and therefore, all
            // flags on the shared anti-duplicator arrays should have been reset

            HEAVY_REFINEMENT_ASSERT(_next_active.checkAllFalse());
            HEAVY_REFINEMENT_ASSERT(_visited_he.checkAllFalse());
        }

    private:

        bool refineImpl(PartitionedHypergraph& hypergraph,
                        const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                        metrics::ThreadSafeMetrics& best_metrics,
                        double,
                        ds::StreamingVector<HypernodeID>& moved_nodes) final ;

        void resetForGroup(ds::ContractionGroupID groupID) override;

        void labelPropagation(PartitionedHypergraph &hypergraph, MovedNodes &moved_nodes);

        bool labelPropagationRound(PartitionedHypergraph &hypergraph,
                                   NextActiveNodes &next_active_nodes,
                                   VisitedEdges &visited_edges, MovedNodes &moved_nodes);

        template<typename F>
        bool moveVertex(PartitionedHypergraph& hypergraph,
                        const HypernodeID hn,
                        NextActiveNodes& next_active_nodes,
                        VisitedEdges& visited_edges,
                        const F& objective_delta) {
            bool is_moved = false;
            ASSERT(hn != kInvalidHypernode);
            if ( hypergraph.isBorderNode(hn) ) {
                ASSERT(hypergraph.nodeIsEnabled(hn));

                Move best_move = _gain.computeMaxGainMove(hypergraph, hn);
                // We perform a move if it either improves the solution quality or, in case of a
                // zero gain move, the balance of the solution.
                const bool positive_gain = best_move.gain < 0;
                const bool zero_gain_move = (_context.refinement.label_propagation.rebalancing &&
                                             best_move.gain == 0 &&
                                             hypergraph.partWeight(best_move.from) - 1 >
                                             hypergraph.partWeight(best_move.to) + 1 &&
                                             hypergraph.partWeight(best_move.to) <
                                             _context.partition.perfect_balance_part_weights[best_move.to]);
                const bool perform_move = positive_gain || zero_gain_move;
                if (best_move.from != best_move.to && perform_move) {
                    PartitionID from = best_move.from;
                    PartitionID to = best_move.to;

                    Gain delta_before = _gain.localDelta();
                    bool changed_part = changeNodePart(hypergraph, hn, from, to, objective_delta);
                    if (changed_part) {
                        // In case the move to block 'to' was successful, we verify that the "real" gain
                        // of the move is either equal to our computed gain or if not, still improves
                        // the solution quality.
                        Gain move_delta = _gain.localDelta() - delta_before;
                        bool accept_move = (move_delta == best_move.gain || move_delta <= 0);
                        if (accept_move) {
                            DBG << "Move hypernode" << hn << "from block" << from << "to block" << to
                                << "with gain" << best_move.gain << "( Real Gain: " << move_delta << ")";

                            // Set all neighbors of the vertex to active
                            for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
                                if ( hypergraph.edgeSize(he) <=
                                     ID(_context.refinement.label_propagation.hyperedge_size_activation_threshold) ) {
                                    if ( _visited_he.compare_and_set_to_true(he) ) {
                                        visited_edges.push_back(he);
                                        for (const HypernodeID& pin : hypergraph.pins(he)) {
                                            //Make sure that pin is not in the intermediate state between being
                                            // reactivated as pin and being enabled that occurs during uncontraction
                                            if (!hypergraph.nodeIsEnabled(pin)) continue;
                                            if (_next_active.compare_and_set_to_true(pin)) {
                                                next_active_nodes.push_back(pin);
                                            }
                                        }
                                    }
                                }
                            }
                            if ( _next_active.compare_and_set_to_true(hn) ) {
                                next_active_nodes.push_back(hn);
                            }
                            is_moved = true;
                        } else {
                            DBG << "Revert move of hypernode" << hn << "from block" << from << "to block" << to
                                << "( Expected Gain:" << best_move.gain << ", Real Gain:" << move_delta << ")";
                            // In case, the real gain is not equal to the computed gain and
                            // worsens the solution quality we revert the move.
                            ASSERT(hypergraph.partID(hn) == to);
                            changeNodePart(hypergraph, hn, to, from, objective_delta);
                        }
                    }
                }
            }

            return is_moved;
        }

        template<typename F>
        bool changeNodePart(PartitionedHypergraph& phg,
                            const HypernodeID hn,
                            const PartitionID from,
                            const PartitionID to,
                            const F& objective_delta) {
            ASSERT(phg.nodeIsEnabled(hn));
            bool success = false;
            if ( _context.partition.paradigm == Paradigm::nlevel && phg.isGainCacheInitialized()) {
                success = phg.changeNodePartWithGainCacheUpdate(hn, from, to,
                                                                _context.partition.max_part_weights[to], [] { }, objective_delta);
            } else {
                success = phg.changeNodePart(hn, from, to,
                                             _context.partition.max_part_weights[to], []{}, objective_delta);
            }
            return success;
        }

        const Context& _context;
        const TaskGroupID _task_group_id;
        GainCalculator _gain;
        ActiveNodes _active_nodes;
        // todo mlaupichler The AsyncLPRefiner should not be bound so tightly to the idea of ContractionGroups. Use a template argument for the OwnerID (also in the _lock_manager pointer) in order to allow any OwnerID for the locking.
        // ! The ID of the contraction group that the seed refinement nodes are from. Used as an identifier for locking
        // nodes so uncontracting and refinement of a contraction group use the same locks.
        ds::ContractionGroupID _contraction_group_id;

        // ! A pointer to the LockManager that is used by uncontraction and refinement operations during the entire
        // uncoarsening.
        ds::GroupLockManager* _lock_manager;

        // ! Mersenne-twister random number generator for shuffling vectors
        std::mt19937 _rng;

        // ! Anti-duplicator flag-arrays that are shared with other refiners
        ds::ThreadSafeFlagArray<HypernodeID>& _next_active;
        ds::ThreadSafeFlagArray<HyperedgeID>& _visited_he;

    };

    using AsyncLPKm1Refiner = AsyncLPRefiner<LocalKm1Policy>;
    using AsyncLPCutRefiner = AsyncLPRefiner<LocalCutPolicy>;
}  // namespace kahypar



#endif //KAHYPAR_ASYNC_LP_REFINER_H

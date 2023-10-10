/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Nikolai Maas <nikolai.maas@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include "mt-kahypar/partition/refinement/deterministic/deterministic_jet_refiner.h"

#include "tbb/parallel_for.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/utilities.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/datastructures/streaming_vector.h"

#include <tbb/parallel_reduce.h>

namespace mt_kahypar {

template<typename GraphAndGainTypes>
bool DeterministicJetRefiner<GraphAndGainTypes>::refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
    const vec<HypernodeID>&,
    Metrics& best_metrics, const double time_limit) {
    Metrics current_metrics = best_metrics;
    const HyperedgeWeight input_quality = best_metrics.quality;
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    resizeDataStructuresForCurrentK();

    auto afterburner = [&](const HypernodeID hn, auto add_node_fn) {
        const PartitionID from = phg.partID(hn);
        const auto [gain, to] = _gains_and_target[hn];
        Gain total_gain = 0;
        for (const HyperedgeID& he : phg.incidentEdges(hn)) {
            HypernodeID pin_count_in_from_part_after = 0;
            HypernodeID pin_count_in_to_part_after = 1;
            for (const HypernodeID& pin : phg.pins(he)) {
                if (pin != hn) {
                    // Jet uses an order based on the precomputed gain values:
                    // If the precomputed gain of another node is better than for the current node
                    // (or the gain is equal and the id is smaller), we assume the node is already
                    // moved to its target part.
                    auto [gain_p, to_p] = _gains_and_target[pin];
                    PartitionID part = (gain_p < gain || (gain_p == gain && pin < hn)) ? to_p : phg.partID(pin);
                    if (part == from) {
                        pin_count_in_from_part_after++;
                    } else if (part == to) {
                        pin_count_in_to_part_after++;
                    }
                }
            }
            SynchronizedEdgeUpdate sync_update;
            sync_update.he = he;
            sync_update.edge_weight = phg.edgeWeight(he);
            sync_update.edge_size = phg.edgeSize(he);
            sync_update.pin_count_in_from_part_after = pin_count_in_from_part_after;
            sync_update.pin_count_in_to_part_after = pin_count_in_to_part_after;
            total_gain += AttributedGains::gain(sync_update);
        }

        if (total_gain <= 0) {
            add_node_fn();
            _locks.set(hn);
        }
    };
    _current_partition_is_best = true;
    size_t rounds_without_improvement = 0;
    const size_t max_rounds = _context.refinement.deterministic_refinement.jet.fixed_n_iterations;
    const size_t max_rounds_without_improvement = _context.refinement.deterministic_refinement.jet.num_iterations;
    for (size_t i = 0; rounds_without_improvement < max_rounds_without_improvement && (max_rounds == 0 || i < max_rounds); ++i) {

        if (_current_partition_is_best) {
            storeCurrentPartition(phg, _best_partition);
        } else {
            storeCurrentPartition(phg, _current_partition);
        }

        HEAVY_REFINEMENT_ASSERT(noInvalidPartitions(phg, _best_partition));

        computeActiveNodesFromGraph(phg, true);

        HEAVY_REFINEMENT_ASSERT(arePotentialMovesToOtherParts(phg, _active_nodes), "active nodes");

        // label prop round
        _locks.reset();
        ds::StreamingVector<HypernodeID> tmp_final_moves;
        tbb::parallel_for(UL(0), _active_nodes.size(), [&](size_t j) {
            const auto n = _active_nodes[j];
            afterburner(n, [&] {tmp_final_moves.stream(n);});
        });
        _moves = tmp_final_moves.copy_parallel();
        HEAVY_REFINEMENT_ASSERT(arePotentialMovesToOtherParts(phg, _moves), "moves");

        // Apply all moves
        auto range = tbb::blocked_range<size_t>(UL(0), _moves.size());
        auto accum = [&](const tbb::blocked_range<size_t>& r, const Gain& init) -> Gain {
            Gain my_gain = init;
            for (size_t i = r.begin(); i < r.end(); ++i) {
                my_gain += performMoveWithAttributedGain(phg, _moves[i]);
            }
            return my_gain;
        };
        Gain gain = tbb::parallel_reduce(range, 0, accum, std::plus<>());

        current_metrics.quality -= gain;
        current_metrics.imbalance = metrics::imbalance(phg, _context);
        HEAVY_REFINEMENT_ASSERT(current_metrics.quality == metrics::quality(phg, _context, false),
            V(current_metrics.quality) << V(metrics::quality(phg, _context, false)));

        // rebalance
        // TODO: This is not deterministic yet
        recomputePenalties(phg, false);
        if (!metrics::isBalanced(phg, _context)) {
            DBG << "[JET] starting rebalancing with quality " << current_metrics.quality << " and imbalance " << current_metrics.imbalance;
            mt_kahypar_partitioned_hypergraph_t part_hg = utils::partitioned_hg_cast(phg);
            _rebalancer.refine(part_hg, {}, current_metrics, time_limit);
            current_metrics.imbalance = metrics::imbalance(phg, _context);
            DBG << "[JET] finished rebalancing with quality " << current_metrics.quality << " and imbalance " << current_metrics.imbalance;
            recomputePenalties(phg, true);
            HEAVY_REFINEMENT_ASSERT(phg.checkTrackedPartitionInformation(_gain_cache));
        }
        ASSERT(current_metrics.quality == metrics::quality(phg, _context, false));
        ++rounds_without_improvement;
        if (current_metrics.quality < best_metrics.quality && metrics::isBalanced(phg, _context)) {
            if (best_metrics.quality - current_metrics.quality > _context.refinement.deterministic_refinement.jet.relative_improvement_threshold * best_metrics.quality) {
                rounds_without_improvement = 0;
            }
            best_metrics = current_metrics;
            _current_partition_is_best = true;
        } else {
            _current_partition_is_best = false;
        }
        DBG << "[JET] Finished iteration " << i << " with quality " << current_metrics.quality << " and imbalance " << current_metrics.imbalance;
    }

    if (!_current_partition_is_best) {
        DBG << "[JET] Rollback to best partition with value " << best_metrics.quality;
        rollbackToBestPartition(phg);
        recomputePenalties(phg, true);
        HEAVY_REFINEMENT_ASSERT(phg.checkTrackedPartitionInformation(_gain_cache));
    }
    HEAVY_REFINEMENT_ASSERT(best_metrics.quality == metrics::quality(phg, _context, false),
        V(best_metrics.quality) << V(metrics::quality(phg, _context, false)));
    return best_metrics.quality < input_quality;
}


template<typename GraphAndGainTypes>
void DeterministicJetRefiner<GraphAndGainTypes>::computeActiveNodesFromGraph(const PartitionedHypergraph& phg, bool first_round) {
    const bool top_level = (phg.initialNumNodes() == _top_level_num_nodes);
    _active_nodes.clear();
    auto process_node = [&](const HypernodeID hn, auto add_node_fn) {
        const bool is_border = phg.isBorderNode(hn);
        const bool is_locked = !first_round || _locks[hn];
        if (!is_border || is_locked) {
            _gains_and_target[hn] = { 0, phg.partID(hn) };
        } else {
            const double gain_factor = top_level ? _context.refinement.deterministic_refinement.jet.negative_gain_factor_fine :
                _context.refinement.deterministic_refinement.jet.negative_gain_factor_coarse;
            RatingMap& tmp_scores = _gain_computation.localScores();
            Gain isolated_block_gain = 0;
            _gain_computation.precomputeGains(phg, hn, tmp_scores, isolated_block_gain, true);
            // Note: rebalance=true is important here to allow negative gain moves
            Move best_move = _gain_computation.computeMaxGainMoveForScores(phg, tmp_scores, isolated_block_gain, hn,
                /*rebalance=*/true,
                /*consider_non_adjacent_blocks=*/false,
                /*allow_imbalance=*/true);
            tmp_scores.clear();
            bool accept_node = (best_move.gain <= 0 || best_move.gain < std::floor(gain_factor * isolated_block_gain))
                && best_move.to != phg.partID(hn);
            if (accept_node) {
                _gains_and_target[hn] = { best_move.gain, best_move.to };
                add_node_fn();
            } else {
                _gains_and_target[hn] = { 0, phg.partID(hn) };
            }
        }
    };
    ds::StreamingVector<HypernodeID> tmp_active_nodes;
    // compute gain for every node 
    phg.doParallelForAllNodes([&](const HypernodeID& hn) {
        process_node(hn, [&] {tmp_active_nodes.stream(hn);});
    });
    _active_nodes = tmp_active_nodes.copy_parallel();
}

template <typename GraphAndGainTypes>
void DeterministicJetRefiner<GraphAndGainTypes>::initializeImpl(mt_kahypar_partitioned_hypergraph_t& phg) {
    _rebalancer.initialize(phg);
}

template<typename GraphAndGainTypes>
Gain DeterministicJetRefiner<GraphAndGainTypes>::performMoveWithAttributedGain(
    PartitionedHypergraph& phg, const HypernodeID hn) {
    const auto from = phg.partID(hn);
    const auto [gain, to] = _gains_and_target[hn];
    Gain attributed_gain = 0;
    auto objective_delta = [&](const SynchronizedEdgeUpdate& sync_update) {
        attributed_gain -= AttributedGains::gain(sync_update);
    };
    ASSERT(to >= 0 && to < _current_k);
    changeNodePart(phg, hn, from, to, objective_delta);
    return attributed_gain;
}

template <typename GraphAndGainTypes>
void DeterministicJetRefiner<GraphAndGainTypes>::storeCurrentPartition(const PartitionedHypergraph& phg,
    parallel::scalable_vector<PartitionID>& parts) {
    phg.doParallelForAllNodes([&](const HypernodeID hn) {
        parts[hn] = phg.partID(hn);
    });
}

template <typename GraphAndGainTypes>
void DeterministicJetRefiner<GraphAndGainTypes>::rollbackToBestPartition(PartitionedHypergraph& phg) {
    // This function is passed as lambda to the changeNodePart function and used
    // to calculate the "real" delta of a move (in terms of the used objective function).
    auto objective_delta = [&](const SynchronizedEdgeUpdate& sync_update) {
        _gain_computation.computeDeltaForHyperedge(sync_update);
    };

    auto reset_node = [&](const HypernodeID hn) {
        const PartitionID part_id = phg.partID(hn);
        if (part_id != _best_partition[hn]) {
            ASSERT(_best_partition[hn] != kInvalidPartition);
            ASSERT(_best_partition[hn] >= 0 && _best_partition[hn] < _current_k);
            changeNodePart(phg, hn, part_id, _best_partition[hn], objective_delta);
        }
    };
    phg.doParallelForAllNodes(reset_node);
    _current_partition_is_best = true;
}

template <typename GraphAndGainTypes>
void DeterministicJetRefiner<GraphAndGainTypes>::recomputePenalties(const PartitionedHypergraph& phg,
    bool did_rebalance) {
    parallel::scalable_vector<PartitionID>& current_parts = _current_partition_is_best ? _best_partition : _current_partition;
    auto recompute = [&](const HypernodeID hn) {
        const bool node_was_moved = (phg.partID(hn) != current_parts[hn]);
        if (node_was_moved) {
            _gain_cache.recomputeInvalidTerms(phg, hn);
        } else {
            ASSERT(_gain_cache.penaltyTerm(hn, phg.partID(hn)) == _gain_cache.recomputePenaltyTerm(phg, hn));
        }
    };

    // TODO this isn't needed for graphs --> skip
    if (_gain_cache.isInitialized()) {
        if (did_rebalance) {
            phg.doParallelForAllNodes(recompute);
        } else {
            tbb::parallel_for(UL(0), _active_nodes.size(), [&](const size_t j) {
                recompute(_active_nodes[j]);
            });
        }
    }
}

template <typename GraphAndGainTypes>
bool DeterministicJetRefiner<GraphAndGainTypes>::arePotentialMovesToOtherParts(const PartitionedHypergraph& phg, const parallel::scalable_vector<HypernodeID>& moves) {
    for (auto hn : moves) {
        const auto [gain, to] = _gains_and_target[hn];
        if (to == phg.partID(hn)) {
            DBG << "Trying to move node " << hn << " to own part expecting gain " << gain;
            return false;
        }
    }
    return true;
}

template <typename GraphAndGainTypes>
bool DeterministicJetRefiner<GraphAndGainTypes>::noInvalidPartitions(const PartitionedHypergraph& phg, const parallel::scalable_vector<PartitionID>& parts) {
    phg.doParallelForAllNodes([&](const HypernodeID hn) {
        ASSERT(parts[hn] != kInvalidPartition);
        ASSERT(parts[hn] < _current_k);
        unused(hn);
    });
    unused(parts);
    return true;
}

namespace {
#define DETERMINISTIC_JET_REFINER(X) DeterministicJetRefiner<X>
}


INSTANTIATE_CLASS_WITH_VALID_TRAITS(DETERMINISTIC_JET_REFINER)
}; // namespace
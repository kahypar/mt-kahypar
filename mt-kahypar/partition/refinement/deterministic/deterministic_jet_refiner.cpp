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
    Metrics& best_metrics, const double) {
    utils::Timer& timer = utils::Utilities::instance().getTimer(_context.utility_id);
    Metrics current_metrics = best_metrics;
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    const HyperedgeWeight input_quality = best_metrics.quality;
    bool was_already_balanced = metrics::isBalanced(phg, _context);
    const bool top_level = (phg.initialNumNodes() == _top_level_num_nodes);

    resizeDataStructuresForCurrentK();

    auto afterburner = [&](const HypernodeID hn, auto add_node_fn) {
        Gain total_gain = 0;
        const PartitionID from = phg.partID(hn);
        const auto [gain, to] = _gains_and_target[hn];
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
        timer.start_timer("active_nodes", "Active Nodes");
        computeActiveNodesFromGraph(phg);
        timer.stop_timer("active_nodes");
        HEAVY_REFINEMENT_ASSERT(arePotentialMovesToOtherParts(phg, _active_nodes), "active nodes");
        timer.start_timer("afterburner", "Afterburner");
        // label prop round
        _locks.reset();
        tmp_active_nodes.clear_parallel();

        if (phg.is_graph) {
            tbb::parallel_for(UL(0), _active_nodes.size(), [&](size_t j) {
                const auto n = _active_nodes[j];
                afterburner(n, [&] {tmp_active_nodes.stream(n);});
            });
            _moves = tmp_active_nodes.copy_parallel();
        } else {
            hypergraphAfterburner(phg);
        }
        HEAVY_REFINEMENT_ASSERT(arePotentialMovesToOtherParts(phg, _moves), "moves");
        timer.stop_timer("afterburner");

        // Apply all moves
        timer.start_timer("apply_moves", "Apply Moves");
        if (phg.is_graph) {
            tbb::parallel_for(0UL, _moves.size(), [&](const size_t i) {
                performMoveWithAttributedGain<true>(phg, _moves[i]);
            });
        } else {
            auto range = tbb::blocked_range<size_t>(UL(0), _active_nodes.size());
            auto accum = [&](const tbb::blocked_range<size_t>& r, const Gain& init) -> Gain {
                Gain my_gain = init;
                for (size_t i = r.begin(); i < r.end(); ++i) {
                    const HypernodeID hn = _active_nodes[i];
                    if (_afterburner_gain[hn] <= 0) {
                        _locks.set(hn);
                        my_gain += performMoveWithAttributedGain<false>(phg, hn);
                    }
                }
                return my_gain;
            };
            Gain gain = tbb::parallel_reduce(range, 0, accum, std::plus<>());
            //DBG << V(i) << V(current_metrics.quality) << V(gain);
            current_metrics.quality -= gain;
        }
        timer.stop_timer("apply_moves");
        // rebalance
        if (!metrics::isBalanced(phg, _context)) {
            DBG << "[JET] starting rebalancing with quality " << current_metrics.quality << " and imbalance " << metrics::imbalance(phg, _context);
            const bool run_until_balanced = rounds_without_improvement == max_rounds_without_improvement - 1 && was_already_balanced;
            if (top_level) {
                timer.start_timer("top_level_rebalance", "Top Level Rebalance");
            } else {
                timer.start_timer("rebalance", "Rebalance");
            }
            mt_kahypar_partitioned_hypergraph_t part_hg = utils::partitioned_hg_cast(phg);
            _rebalancer.jetRebalance(part_hg, current_metrics, run_until_balanced);
            if (top_level) {
                timer.stop_timer("top_level_rebalance");
            } else {
                timer.stop_timer("rebalance");
            }
            DBG << "[JET] finished rebalancing with quality " << current_metrics.quality << " and imbalance " << metrics::imbalance(phg, _context);
        }
        timer.start_timer("reb_quality", "Quality after Rebalancing");
        if (phg.is_graph) {
            current_metrics.quality += calculateGainDelta(phg);
        }
        current_metrics.imbalance = metrics::imbalance(phg, _context);
        timer.stop_timer("reb_quality");
        ASSERT(current_metrics.quality == metrics::quality(phg, _context, false), V(current_metrics.quality) << V(metrics::quality(phg, _context, false)));
        ++rounds_without_improvement;
        // if the parition was ever balanced => look for balanced partition with better quality
        // if the partition was never balanced => look for less imbalanced partition regardless of quality
        const bool is_balanced = metrics::isBalanced(phg, _context);
        if ((current_metrics.quality < best_metrics.quality && was_already_balanced && is_balanced) || (!was_already_balanced && current_metrics.imbalance < best_metrics.imbalance)) {
            if (best_metrics.quality - current_metrics.quality > _context.refinement.deterministic_refinement.jet.relative_improvement_threshold * best_metrics.quality) {
                rounds_without_improvement = 0;
            }
            best_metrics = current_metrics;
            _current_partition_is_best = true;
            was_already_balanced |= is_balanced;
        } else {
            _current_partition_is_best = false;
        }
        DBG << "[JET] Finished iteration " << i << " with quality " << current_metrics.quality << " and imbalance " << current_metrics.imbalance;
    }

    phg.resetEdgeSynchronization();
    if (!_current_partition_is_best) {
        DBG << "[JET] Rollback to best partition with value " << best_metrics.quality;
        if (phg.is_graph) {
            rollbackToBestPartition<true>(phg);
        } else {
            rollbackToBestPartition<false>(phg);
        }
        timer.start_timer("reb_quality", "Quality after Rebalancing");
        timer.stop_timer("reb_quality");
    }
    HEAVY_REFINEMENT_ASSERT(best_metrics.quality == metrics::quality(phg, _context, false),
        V(best_metrics.quality) << V(metrics::quality(phg, _context, false)));

    return best_metrics.quality < input_quality;
}


template<typename GraphAndGainTypes>
void DeterministicJetRefiner<GraphAndGainTypes>::computeActiveNodesFromGraph(const PartitionedHypergraph& phg) {
    const bool top_level = (phg.initialNumNodes() == _top_level_num_nodes);
    _active_nodes.clear();
    auto process_node = [&](const HypernodeID hn, auto add_node_fn) {
        _part_before_round[hn] = phg.partID(hn);
        const bool is_border = phg.isBorderNode(hn);
        const bool is_locked = _locks[hn];
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
    tmp_active_nodes.clear_sequential();
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

template <typename GraphAndGainTypes>
void DeterministicJetRefiner<GraphAndGainTypes>::storeCurrentPartition(const PartitionedHypergraph& phg,
    parallel::scalable_vector<PartitionID>& parts) {
    phg.doParallelForAllNodes([&](const HypernodeID hn) {
        parts[hn] = phg.partID(hn);
    });
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
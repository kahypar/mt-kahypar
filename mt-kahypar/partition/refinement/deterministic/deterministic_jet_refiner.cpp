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
    const vec<HypernodeID>& refinement_nodes,
    Metrics& best_metrics, const double time_limit) {
    const HyperedgeWeight input_quality = best_metrics.quality;
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(phg);
    Gain overall_improvement = 0;

    // resize data structures for current k
    if (_context.partition.k != _current_k) {
        _current_k = _context.partition.k;
        _gain_computation.changeNumberOfBlocks(_current_k);
    }
    computeActiveNodesFromGraph(phg, true);

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
            total_gain += AttributedGains::gain(he, phg.edgeWeight(he), phg.edgeSize(he),
                pin_count_in_from_part_after, pin_count_in_to_part_after);
        }

        if (total_gain <= 0) {
            add_node_fn();
            _locks.set(hn);
        }
    };

    ds::StreamingVector<HypernodeID> tmp_final_moves; // TODO: Use actual moves to be able to revert
    for (size_t i = 0; i < _context.refinement.deterministic_refinement.jet.fixed_n_iterations; ++i) {
        // label prop round
        _locks.reset();
        tbb::parallel_for(UL(0), _active_nodes.size(), [&](size_t j) {
            afterburner(_active_nodes[j], tmp_final_moves.stream(_active_nodes[j]));
        });

        _moves = tmp_final_moves.copy_parallel();

        // TODO: Apply all moves
        // are corresponding values updated properly?
        auto range = tbb::blocked_range<size_t>(UL(0), _moves.size());
        auto accum = [&](const tbb::blocked_range<size_t>& r, const Gain& init) -> Gain {
            Gain my_gain = init;
            for (size_t i = r.begin(); i < r.end(); ++i) {
                my_gain += performMoveWithAttributedGain(phg, _moves[i], true);
            }
            return my_gain;
        };
        Gain gain = tbb::parallel_reduce(range, 0, accum, std::plus<>());
        if (gain < 0) {
            // revert?
        } else {
            best_metrics.quality -= gain;
            best_metrics.imbalance = metrics::imbalance(phg, _context);
        }
        // rebalance
        if (!metrics::isBalanced(phg, _context)) {
            _rebalancer.refine(phg, {}, best_metrics, time_limit);
        }
        // TODO: Probably rollback here
        // TODO: Find out whether saving partitions or reverting moves is the better option

    }
    return best_metrics.quality < input_quality;
}


template<typename GraphAndGainTypes>
void DeterministicJetRefiner<GraphAndGainTypes>::computeActiveNodesFromGraph(const PartitionedHypergraph& phg, bool first_round) {
    const bool top_level = (phg.initialNumNodes() == _top_level_num_nodes);

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
            _gain_computation.precomputeGains(phg, hn, tmp_scores, isolated_block_gain);
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
        process_node(hn, tmp_active_nodes.stream(hn));
    });
    _active_nodes = tmp_active_nodes.copy_parallel();
}

template <typename GraphAndGainTypes>
void DeterministicJetRefiner<GraphAndGainTypes>::initializeImpl(mt_kahypar_partitioned_hypergraph_t& phg) {
    _rebalancer.initialize(phg);
}

/*
 * for configs where we don't know exact gains --> have to trace the overall improvement with attributed gains
 */
template<typename GraphAndGainTypes>
Gain DeterministicJetRefiner<GraphAndGainTypes>::performMoveWithAttributedGain(
    PartitionedHypergraph& phg, const Move& m, bool activate_neighbors) {
    Gain attributed_gain = 0;
    auto objective_delta = [&](const SynchronizedEdgeUpdate& sync_update) {
        attributed_gain -= AttributedGains::gain(sync_update);
    };
    return attributed_gain;
}
};
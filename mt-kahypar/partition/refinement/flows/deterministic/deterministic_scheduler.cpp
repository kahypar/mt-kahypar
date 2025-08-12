/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Robert Krause <robert.krause@student.kit.edu>
 * Copyright (C) 2025 Nikolai Maas <nikolai.maas@kit.edu>
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

#include "mt-kahypar/partition/refinement/flows/deterministic/deterministic_scheduler.h"

namespace mt_kahypar {


template<typename GraphAndGainTypes>
bool DeterministicFlowRefinementScheduler<GraphAndGainTypes>::refineImpl(
    mt_kahypar_partitioned_hypergraph_t& hypergraph,
    const parallel::scalable_vector<HypernodeID>&,
    Metrics& best_metrics,
    const double) {
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    if (num_hypernodes == phg.initialNumNodes()) {
        _schedule.setTopLevelFlag();
    }
    utils::Timer& timer = utils::Utilities::instance().getTimer(_context.utility_id);

    Metrics current_metrics = best_metrics;
    _schedule.initialize(_quotient_graph);

    HyperedgeWeight overall_delta = 0;
    HyperedgeWeight minImprovement = 0;
    std::atomic<HyperedgeWeight> round_delta(std::numeric_limits<HyperedgeWeight>::max());
    minImprovement = _context.refinement.flows.min_relative_improvement_per_round * current_metrics.quality;
    while (round_delta >= minImprovement && _schedule.hasActiveBlocks()) {
        size_t numScheduledBlocks = _schedule.getNextMatching(_scheduled_blocks, _quotient_graph);
        round_delta = 0;
        minImprovement = _context.refinement.flows.min_relative_improvement_per_round * current_metrics.quality;
        while (numScheduledBlocks > 0) {
            tbb::parallel_for(0UL, _refiner.size(), [&](const size_t refiner_idx) {
                ScheduledPair sp;
                while (_scheduled_blocks.try_pop(sp)) {
                    timer.start_timer("region_growing", "Grow Region", true);
                    const Subhypergraph sub_hg = _constructor.construct(sp.bp, _quotient_graph, phg, /*deterministic=*/true);
                    timer.stop_timer("region_growing");

                    HyperedgeWeight improvement = 0;
                    auto start = std::chrono::high_resolution_clock::now();
                    if ( sub_hg.numNodes() > 0 ) {
                        auto partitioned_hg = utils::partitioned_hg_const_cast(phg);
                        _refiner[refiner_idx]->initialize(partitioned_hg);
                        MoveSequence moves = _refiner[refiner_idx]->refine(partitioned_hg, sub_hg, start);

                        timer.start_timer("apply_moves", "Apply Moves", true);
                        improvement = applyMoves(moves, phg);
                        timer.stop_timer("apply_moves");

                        reportResults(sp.bp.i, sp.bp.j, moves);
                    }

                    round_delta += improvement;
                    _quotient_graph.edge(sp.bp).total_improvement += improvement;
                    _quotient_graph.edge(sp.bp).markAsPreviouslyScheduled();
                    // TODO: num improvements?
                }
            });
            _new_cut_hes.clear();
            numScheduledBlocks = _schedule.getNextMatching(_scheduled_blocks, _quotient_graph);
            ASSERT(metrics::isBalanced(phg, _context));
        }
        overall_delta += round_delta;
        current_metrics.quality -= round_delta;
        _schedule.resetForNewRound(_quotient_graph);
    }

    // Update metrics statistics
    HEAVY_REFINEMENT_ASSERT(best_metrics.quality - overall_delta == metrics::quality(phg, _context),
        V(best_metrics.quality) << V(overall_delta) << V(metrics::quality(phg, _context)));
    best_metrics.quality -= overall_delta;
    best_metrics.imbalance = metrics::imbalance(phg, _context);

    // Update Gain Cache
    if (_context.forceGainCacheUpdates() && _gain_cache.isInitialized()) {
        phg.doParallelForAllNodes([&](const HypernodeID& hn) {
            if (_was_moved[hn]) {
                _gain_cache.recomputeInvalidTerms(phg, hn);
                _was_moved[hn] = uint8_t(false);
            }
        });
    }
    return overall_delta > 0;
}

namespace {
#define DETERMINISTIC_FLOW_REFINEMENT_SCHEDULER(X) DeterministicFlowRefinementScheduler<X>
}

INSTANTIATE_CLASS_WITH_VALID_TRAITS(DETERMINISTIC_FLOW_REFINEMENT_SCHEDULER)

}
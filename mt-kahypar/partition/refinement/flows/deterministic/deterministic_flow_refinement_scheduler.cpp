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

#include "mt-kahypar/partition/refinement/flows/deterministic/deterministic_flow_refinement_scheduler.h"

#include <tbb/concurrent_queue.h>

#include "mt-kahypar/macros.h"

namespace mt_kahypar {

template<typename GraphAndGainTypes>
DeterministicFlowRefinementScheduler<GraphAndGainTypes>::DeterministicFlowRefinementScheduler(const HypernodeID num_hypernodes,
                                                                                              const HyperedgeID num_hyperedges,
                                                                                              const Context& context,
                                                                                              GainCache& gain_cache) :
  _context(context),
  _gain_cache(gain_cache),
  _current_k(context.partition.k),
  _num_hyperedges(num_hyperedges),
  _quotient_graph(num_hyperedges, context),
  _schedule(context, _quotient_graph),
  _refiner(),
  _constructor(num_hypernodes, num_hyperedges, context),
  _was_moved(num_hypernodes,  static_cast<uint8_t>(false)) {
  for (size_t i = 0; i < _context.refinement.flows.num_parallel_searches; ++i) {
    _refiner.emplace_back(nullptr);
  }
}

template<typename GraphAndGainTypes>
DeterministicFlowRefinementScheduler<GraphAndGainTypes>::DeterministicFlowRefinementScheduler(const HypernodeID num_hypernodes,
                                                                                              const HyperedgeID num_hyperedges,
                                                                                              const Context& context,
                                                                                              gain_cache_t gain_cache) :
  DeterministicFlowRefinementScheduler(num_hypernodes, num_hyperedges, context,
    GainCachePtr::cast<GainCache>(gain_cache)) {}

template<typename GraphAndGainTypes>
bool DeterministicFlowRefinementScheduler<GraphAndGainTypes>::refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                                         const parallel::scalable_vector<HypernodeID>&,
                                                                         Metrics& best_metrics,
                                                                         const double) {
  PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
  utils::Timer& timer = utils::Utilities::instance().getTimer(_context.utility_id);
  Metrics current_metrics = best_metrics;

  HyperedgeWeight overall_delta = 0;
  HyperedgeWeight min_improvement = _context.refinement.flows.min_relative_improvement_per_round * current_metrics.quality;
  std::atomic<HyperedgeWeight> round_delta(std::numeric_limits<HyperedgeWeight>::max());
  tbb::concurrent_queue<BlockPair> scheduled_blocks;
  for (size_t round = 0; round_delta >= min_improvement; ++round) {
    size_t num_scheduled_blocks = _schedule.getNextMatching(scheduled_blocks);
    if (num_scheduled_blocks == 0) break;

    round_delta = 0;
    while (num_scheduled_blocks > 0) {
      DBG << "Next matching:" << num_scheduled_blocks << "blocks";
      _new_cut_hes.clear_sequential();

      tbb::parallel_for(UL(0), _context.refinement.flows.num_parallel_searches, [&](const size_t refiner_idx) {
        BlockPair blocks;
        while (scheduled_blocks.try_pop(blocks)) {
          timer.start_timer("region_growing", "Grow Region", true);
          const Subhypergraph sub_hg = _constructor.construct(blocks, _quotient_graph, phg, /*deterministic=*/true);
          timer.stop_timer("region_growing");

          HyperedgeWeight improvement = 0;
          auto start = std::chrono::high_resolution_clock::now();
          if ( sub_hg.numNodes() > 0 ) {
            auto partitioned_hg = utils::partitioned_hg_const_cast(phg);
            if (_refiner[refiner_idx] == nullptr) {
              _refiner[refiner_idx] = constructFlowRefiner();
            }
            _refiner[refiner_idx]->initialize(partitioned_hg);
            MoveSequence moves = _refiner[refiner_idx]->refine(partitioned_hg, sub_hg, start);

            timer.start_timer("apply_moves", "Apply Moves", true);
            improvement = applyMoves(blocks, moves, phg);
            timer.stop_timer("apply_moves");

            _schedule.reportResults(blocks, improvement);
          }

          round_delta.fetch_add(improvement, std::memory_order_relaxed);
          _quotient_graph.edge(blocks).reportImprovement(improvement);
        }
      });

      // adding HEs to the quotient graph must happen after all moves are applied,
      // since the result depends on the current connectivity set of the HE
      auto hes_to_add = _new_cut_hes.copy_parallel();
      tbb::parallel_for(UL(0), hes_to_add.size(), [&](const size_t i) {
        NewCutHyperedge hyperedge = hes_to_add[i];
        _quotient_graph.addNewCutHyperedge(phg, hyperedge.he, hyperedge.block);
      });
      num_scheduled_blocks = _schedule.getNextMatching(scheduled_blocks);
      ASSERT(metrics::isBalanced(phg, _context));
    }
    overall_delta += round_delta;
    current_metrics.quality -= round_delta;
    _schedule.resetForNewRound();
    DBG << "Round" << round << "finished, delta =" << round_delta
        << "(current quality=" << current_metrics.quality << ")";
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
        _was_moved[hn] = static_cast<uint8_t>(false);
      }
    });
  }
  return overall_delta > 0;
}

template<typename GraphAndGainTypes>
void DeterministicFlowRefinementScheduler<GraphAndGainTypes>::initializeImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph) {
  PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
  resizeDataStructuresForCurrentK();

  utils::Timer& timer = utils::Utilities::instance().getTimer(_context.utility_id);
  timer.start_timer("initialize_quotient_graph", "Initialize Quotient Graph");
  _quotient_graph.initialize(phg);
  _schedule.initialize(_quotient_graph.isInputHypergraph());
  timer.stop_timer("initialize_quotient_graph");
}

template<typename GraphAndGainTypes>
void DeterministicFlowRefinementScheduler<GraphAndGainTypes>::resizeDataStructuresForCurrentK() {
  if ( _current_k != _context.partition.k ) {
    _current_k = _context.partition.k;
    _quotient_graph.changeNumberOfBlocks(_current_k);
    _constructor.changeNumberOfBlocks(_current_k);
  }
}

template<typename GraphAndGainTypes>
std::unique_ptr<IFlowRefiner> DeterministicFlowRefinementScheduler<GraphAndGainTypes>::constructFlowRefiner() {
  return std::make_unique<FlowRefiner<GraphAndGainTypes>>(
    _num_hyperedges, _context, true /*deterministic*/);
}

namespace {

template<typename PartitionedHypergraph, typename GainCache, typename F>
bool changeNodePart(PartitionedHypergraph& phg,
                    GainCache& gain_cache,
                    const HypernodeID hn,
                    const PartitionID from,
                    const PartitionID to,
                    const F& objective_delta,
                    const bool gain_cache_update) {
  bool success = false;
  if (gain_cache_update && gain_cache.isInitialized()) {
    success = phg.changeNodePart(gain_cache, hn, from, to,
      std::numeric_limits<HypernodeWeight>::max(), [] {}, objective_delta);
  } else {
    success = phg.changeNodePart(hn, from, to,
      std::numeric_limits<HypernodeWeight>::max(), [] {}, objective_delta);
  }
  ASSERT(success);
  return success;
}

}

template<typename GraphAndGainTypes>
HyperedgeWeight DeterministicFlowRefinementScheduler<GraphAndGainTypes>::applyMoves(const BlockPair& blocks,
                                                                                    MoveSequence& sequence,
                                                                                    PartitionedHypergraph& phg) {
  const bool gain_cache_update = _context.forceGainCacheUpdates();
  HyperedgeWeight improvement = 0;

  for (const Move& move : sequence.moves) {
    ASSERT(move.from == phg.partID(move.node));
    if (move.from != move.to) {
      auto delta_func = [&](const SynchronizedEdgeUpdate& sync_update) {
        // Collect hyperedges with new blocks in its connectivity set
        if (sync_update.pin_count_in_to_part_after == 1) {
          _new_cut_hes.stream(NewCutHyperedge {sync_update.he, move.to});
        }
        improvement -= AttributedGains::gain(sync_update);
      };

      changeNodePart(phg, _gain_cache, move.node, move.from,
        move.to, delta_func, gain_cache_update);
      _was_moved[move.node] = static_cast<uint8_t>(true);
    }
  }

  ASSERT(!FlowNetworkConstruction::is_exact_model || improvement == sequence.expected_improvement);
  DBG << ( improvement > 0 ? GREEN : "" )
      << "Moved Nodes =" << sequence.moves.size()
      << ", Improvement =" << improvement
      << ", Block Pair =" << blocks.i << "," << blocks.j << ( improvement > 0 ? END : "" );
  return improvement;
}

namespace {
#define DETERMINISTIC_FLOW_REFINEMENT_SCHEDULER(X) DeterministicFlowRefinementScheduler<X>
}

INSTANTIATE_CLASS_WITH_VALID_TRAITS(DETERMINISTIC_FLOW_REFINEMENT_SCHEDULER)

}

#include "deterministic_flow_scheduler.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h" 


namespace mt_kahypar {

template<typename GraphAndGainTypes>
bool DeterministicFlowScheduler<GraphAndGainTypes>::refineImpl(
  mt_kahypar_partitioned_hypergraph_t& hypergraph,
  const vec<HypernodeID>&,
  Metrics& best_metrics,
  const double) {
  constexpr size_t maxIterations = std::numeric_limits<size_t>::max();
  PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);

  bool terminate = false;
  for (size_t round = 0; round < 20 && !terminate; ++round) {
    HyperedgeWeight improvement_this_round = 0;
    size_t processed_block_pairs = 0;
    while (true) {
      SearchID search_id = _quotient_graph.requestNewSearchDeterministic();
      if (search_id == QuotientGraph<TypeTraits>::INVALID_SEARCH_ID) {
        terminate = processed_block_pairs == 0;
        break;
      }
      ++processed_block_pairs;

      const Subhypergraph sub_hg = _constructor.construct(search_id, _quotient_graph, phg);
      _quotient_graph.finalizeConstruction(search_id);
      HyperedgeWeight delta = 0;
      bool improved_solution = false;
      if (sub_hg.numNodes() > 0) {
        MoveSequence sequence = _refiner.refine(search_id, phg, sub_hg);
        if (!sequence.moves.empty() && sequence.expected_improvement >= 0) {
          //timer.start_timer("apply_moves", "Apply Moves", true);
          improvement_this_round += sequence.expected_improvement;
          delta = applyMoves(search_id, sequence);
          assert(delta == sequence.expected_improvement);
          improved_solution = sequence.state == MoveSequenceState::SUCCESS && delta > 0;
         // timer.stop_timer("apply_moves");
        }
      }
      _quotient_graph.finalizeSearchDeterministic(search_id, improved_solution ? delta : 0);
    }
  }

  // for (whfc::TimeReporter& local_timer : timers_thread_specific) {
  //   timer.merge(local_timer, "Refinement", "total");
  // }

  return true;//iterations_done.load();
}

template<typename GraphAndGainTypes>
void DeterministicFlowScheduler<GraphAndGainTypes>::initializeImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph) {
  PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
  _phg = &phg;
  resizeDataStructuresForCurrentK();

  // Initialize Part Weights
  for (PartitionID i = 0; i < _context.partition.k; ++i) {
    _part_weights[i] = phg.partWeight(i);
    _max_part_weights[i] = std::max(
      phg.partWeight(i), _context.partition.max_part_weights[i]);
  }

  utils::Timer& timer = utils::Utilities::instance().getTimer(_context.utility_id);
  timer.start_timer("initialize_quotient_graph", "Initialize Quotient Graph");
  _quotient_graph.initialize(phg);
  timer.stop_timer("initialize_quotient_graph");

  const size_t max_parallism = _context.refinement.flows.num_parallel_searches;
  DBG << "Initial Active Block Pairs =" << _quotient_graph.numActiveBlockPairs()
    << ", Initial Num Threads =" << max_parallism;
  _refiner.initialize(phg);
}

template<typename GraphAndGainTypes>
void DeterministicFlowScheduler<GraphAndGainTypes>::resizeDataStructuresForCurrentK() {
  if (_current_k != _context.partition.k) {
    _current_k = _context.partition.k;
    // Note that in general changing the number of blocks should not resize
    // any data structure as we initialize the scheduler with the final
    // number of blocks. This is just a fallback if someone changes this in the future.
    if (static_cast<size_t>(_current_k) > _part_weights.size()) {
      _part_weights.resize(_current_k);
      _max_part_weights.resize(_current_k);
    }
    _quotient_graph.changeNumberOfBlocks(_current_k);
    _constructor.changeNumberOfBlocks(_current_k);
  }
}

namespace {

struct NewCutHyperedge {
  HyperedgeID he;
  PartitionID block;
};

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

template<typename PartitionedHypergraph, typename GainCache, typename F>
void applyMoveSequence(PartitionedHypergraph& phg,
  GainCache& gain_cache,
  const MoveSequence& sequence,
  const F& objective_delta,
  const bool gain_cache_update,
  vec<uint8_t>& was_moved,
  vec<NewCutHyperedge>& new_cut_hes) {
  for (const Move& move : sequence.moves) {
    ASSERT(move.from == phg.partID(move.node));
    if (move.from != move.to) {
      changeNodePart(phg, gain_cache, move.node, move.from,
        move.to, objective_delta, gain_cache_update);
      was_moved[move.node] = uint8_t(true);
      // If move increases the pin count of some hyperedges in block 'move.to' to one 1
      // we set the corresponding block here.
      int i = new_cut_hes.size() - 1;
      while (i >= 0 && new_cut_hes[i].block == kInvalidPartition) {
        new_cut_hes[i].block = move.to;
        --i;
      }
    }
  }
}

template<typename PartitionedHypergraph, typename GainCache, typename F>
void revertMoveSequence(PartitionedHypergraph& phg,
  GainCache& gain_cache,
  const MoveSequence& sequence,
  const F& objective_delta,
  const bool gain_cache_update) {
  for (const Move& move : sequence.moves) {
    if (move.from != move.to) {
      ASSERT(phg.partID(move.node) == move.to);
      changeNodePart(phg, gain_cache, move.node, move.to,
        move.from, objective_delta, gain_cache_update);
    }
  }
}

template<typename TypeTraits>
void addCutHyperedgesToQuotientGraph(QuotientGraph<TypeTraits>& quotient_graph,
  const vec<NewCutHyperedge>& new_cut_hes) {
  for (const NewCutHyperedge& new_cut_he : new_cut_hes) {
    ASSERT(new_cut_he.block != kInvalidPartition);
    quotient_graph.addNewCutHyperedge(new_cut_he.he, new_cut_he.block);
  }
}

} // namespace

template<typename GraphAndGainTypes>
HyperedgeWeight DeterministicFlowScheduler<GraphAndGainTypes>::applyMoves(MoveSequence& sequence) {
  ASSERT(_phg);

  // TODO: currently we lock the applyMoves method
  // => find something smarter here
  _apply_moves_lock.lock();

  // Compute Part Weight Deltas
  vec<HypernodeWeight> part_weight_deltas(_context.partition.k, 0);
  for (Move& move : sequence.moves) {
    move.from = _phg->partID(move.node);
    if (move.from != move.to) {
      const HypernodeWeight node_weight = _phg->nodeWeight(move.node);
      part_weight_deltas[move.from] -= node_weight;
      part_weight_deltas[move.to] += node_weight;
    }
  }

  HyperedgeWeight improvement = 0;
  vec<NewCutHyperedge> new_cut_hes;
  auto delta_func = [&](const SynchronizedEdgeUpdate& sync_update) {
    improvement -= AttributedGains::gain(sync_update);

    // Collect hyperedges with new blocks in its connectivity set
    if (sync_update.pin_count_in_to_part_after == 1) {
      // the corresponding block will be set in applyMoveSequence(...) function
      new_cut_hes.emplace_back(NewCutHyperedge{ sync_update.he, kInvalidPartition });
    }
  };

  // Update part weights atomically
  PartWeightUpdateResult update_res = partWeightUpdate(part_weight_deltas, false);
  if (update_res.is_balanced) {
    // Apply move sequence to partition
    applyMoveSequence(*_phg, _gain_cache, sequence, delta_func,
      _context.forceGainCacheUpdates(), _was_moved, new_cut_hes);


    sequence.state = MoveSequenceState::SUCCESS;
    DBG << (improvement > 0 ? GREEN : "") << "SUCCESS -"
      << "Moved Nodes =" << sequence.moves.size()
      << ", Expected Improvement =" << sequence.expected_improvement
      << ", Real Improvement =" << improvement
      << ", " << (improvement > 0 ? END : "");

  } else {
    sequence.state = MoveSequenceState::VIOLATES_BALANCE_CONSTRAINT;
    DBG << RED << "Move sequence violated balance constraint ( Moved Nodes ="
      << sequence.moves.size()
      << ", Expected Improvement =" << sequence.expected_improvement
      << ")" << END;
  }

  _apply_moves_lock.unlock();

  if (sequence.state == MoveSequenceState::SUCCESS && improvement > 0) {
    addCutHyperedgesToQuotientGraph(_quotient_graph, new_cut_hes);
  }

  return improvement;
}

template<typename GraphAndGainTypes>
typename DeterministicFlowScheduler<GraphAndGainTypes>::PartWeightUpdateResult
DeterministicFlowScheduler<GraphAndGainTypes>::partWeightUpdate(const vec<HypernodeWeight>& part_weight_deltas,
  const bool rollback) {
  const HypernodeWeight multiplier = rollback ? -1 : 1;
  PartWeightUpdateResult res;
  _part_weights_lock.lock();
  PartitionID i = 0;
  for (; i < _context.partition.k; ++i) {
    if (_part_weights[i] + multiplier * part_weight_deltas[i] > _max_part_weights[i]) {
      DBG << "Move sequence violated balance constraint of block" << i
        << "(Max =" << _max_part_weights[i]
        << ", Actual =" << (_part_weights[i] + multiplier * part_weight_deltas[i]) << ")";
      res.is_balanced = false;
      res.overloaded_block = i;
      res.overload_weight = (_part_weights[i] + multiplier *
        part_weight_deltas[i]) - _max_part_weights[i];
      // Move Sequence Violates Balance Constraint => Rollback
      --i;
      for (; i >= 0; --i) {
        _part_weights[i] -= multiplier * part_weight_deltas[i];
      }
      break;
    }
    _part_weights[i] += multiplier * part_weight_deltas[i];
  }
  _part_weights_lock.unlock();
  return res;
}
}
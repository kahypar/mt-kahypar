/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/partition/refinement/deterministic/rebalancing/deterministic_rebalancer.h"


#include <boost/dynamic_bitset.hpp>

#include <tbb/parallel_for_each.h>
#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_for.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/cast.h"

#include "external_tools/parlaylib/include/parlay/primitives.h"

namespace mt_kahypar {
static constexpr size_t ABSOLUTE_MAX_ROUNDS = 30;

float transformGain(Gain gain_, HypernodeWeight wu) {
  float gain = gain_;
  if (gain > 0) {
    gain /= wu;
  } else if (gain < 0) {
    gain *= wu;
  }
  return gain;
}

template <typename GraphAndGainTypes>
bool DeterministicRebalancer<GraphAndGainTypes>::refineInternal(mt_kahypar_partitioned_hypergraph_t& hypergraph,
  Metrics& best_metrics,
  bool run_until_balanced) {
  PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
  resizeDataStructuresForCurrentK();
  if (_max_part_weights == nullptr) {
    _max_part_weights = &_context.partition.max_part_weights[0];
  }
  _gain_computation.reset();
  initializeDataStructures(phg);
  size_t iteration = 0;
  while (_num_imbalanced_parts > 0 && iteration < ABSOLUTE_MAX_ROUNDS && (run_until_balanced || _context.refinement.deterministic_refinement.jet.max_rebalancing_rounds == 0 || iteration < _context.refinement.deterministic_refinement.jet.max_rebalancing_rounds)) {
    weakRebalancingRound(phg);
    HEAVY_REFINEMENT_ASSERT(checkPreviouslyOverweightParts(phg)); // NOTE this is a light assertion, you can use the normal ASSERT macro
    updateImbalance(phg);
    ++iteration;
  }
  if (!phg.is_graph) {
    Gain delta = _gain_computation.delta();
    HEAVY_REFINEMENT_ASSERT(best_metrics.quality + delta == metrics::quality(phg, _context),
      V(best_metrics.quality) << V(delta) << V(metrics::quality(phg, _context)));
    best_metrics.quality += delta;
    DBG << "[REBALANCE] " << "delta=" << delta;
  }
  DBG << "[REBALANCE] " << "  imbalance=" << metrics::imbalance(phg, _context);
  _max_part_weights = nullptr;
  return true;
}

template <typename  GraphAndGainTypes>
void DeterministicRebalancer< GraphAndGainTypes>::initializeDataStructures(const PartitionedHypergraph& phg) {
  updateImbalance(phg);
}


template <typename  GraphAndGainTypes>
void DeterministicRebalancer<GraphAndGainTypes>::updateImbalance(const PartitionedHypergraph& phg) {
  _num_imbalanced_parts = 0;
  _num_valid_targets = 0;
  for (PartitionID part = 0; part < _context.partition.k; ++part) { // TODO: Not worth parallelizing?!
    if (imbalance(phg, part) > 0) {
      ++_num_imbalanced_parts;
    } else if (isValidTarget(phg, part, 0)) {
      ++_num_valid_targets;
    }
  }
}

template <typename GraphAndGainTypes>
rebalancer::RebalancingMove DeterministicRebalancer<GraphAndGainTypes>::computeGainAndTargetPart(const PartitionedHypergraph& phg,
  const HypernodeID hn,
  bool non_adjacent_blocks) {
  const HypernodeWeight hn_weight = phg.nodeWeight(hn);
  RatingMap& tmp_scores = _gain_computation.localScores();
  Gain isolated_block_gain = 0;
  _gain_computation.precomputeGains(phg, hn, tmp_scores, isolated_block_gain, true);

  Gain best_gain = isolated_block_gain;
  PartitionID best_target = kInvalidPartition;
  for (const auto& entry : tmp_scores) {
    const PartitionID to = entry.key;
    const Gain gain = _gain_computation.gain(entry.value, isolated_block_gain);
    if (gain <= best_gain && isValidTarget(phg, to, hn_weight)) {
      best_gain = gain;
      best_target = to;
    }
  }

  // if no adjacent block with free capacity exists, we need to consider non-adjacent blocks
  if (non_adjacent_blocks && best_target == kInvalidPartition) {
    // we start with a block that is chosen by random, to ensure a reasonable distribution of nodes
    // to target blocks (note: this does not always result in a uniform distribution since some blocks
    // are not an acceptable target, but it should be good enough)
    const PartitionID start = hn % _current_k;
    PartitionID to = start;
    do {
      if (isValidTarget(phg, to, hn_weight)
        && !tmp_scores.contains(to)) {
        best_target = to;
        break;
      }

      ++to;
      if (to == _context.partition.k) {
        to = 0;
      }
    } while (to != start);
    // assertion does not always hold with tight balance constraint or large node weights
    // ASSERT(best_target != kInvalidPartition);
  }
  tmp_scores.clear();
  return { hn, best_target, transformGain(best_gain, hn_weight) };
}

template <typename  GraphAndGainTypes>
void DeterministicRebalancer<GraphAndGainTypes>::weakRebalancingRound(PartitionedHypergraph& phg) {
  for (auto& moves : tmp_potential_moves){
    moves.clear_sequential();
  }

  // calculate gain and target for each node in a overweight part
  // group moves by source part
  phg.doParallelForAllNodes([&](const HypernodeID hn) {
    const PartitionID from = phg.partID(hn);
    const HypernodeWeight weight = phg.nodeWeight(hn);
    if (imbalance(phg, from) > 0 && mayMoveNode(phg, from, weight)) {
      const auto& triple = computeGainAndTargetPart(phg, hn, true); // NOTE How can this be a reference?
      if (from != triple.to && triple.to != kInvalidPartition) {
        tmp_potential_moves[from].stream(triple);
      }
    }
  });
  tbb::parallel_for(0UL, _moves.size(), [&](const size_t i) {   // NOTE use a variable that indicates that this is a block?
    if (tmp_potential_moves[i].size() > 0) {
      _moves[i] = tmp_potential_moves[i].copy_parallel();
      const size_t move_size = _moves[i].size();

      if (move_size > 0) {
        // sort the moves from each overweight part by priority
        parlay::sort_inplace(_moves[i], [&](const rebalancer::RebalancingMove& a, const rebalancer::RebalancingMove& b) {
          return a.priority < b.priority || (a.priority == b.priority && a.hn > b.hn);
        });
        // calculate perfix sum for each source-part to know which moves to execute (prefix_sum > current_weight - max_weight)
        size_t last_move_idx = 0;
        if (move_size > _context.refinement.deterministic_refinement.jet.seq_find_rebalancing_moves) {
          if (move_size > _move_weights[i].size()) {
            _move_weights[i].resize(move_size);
          }
          tbb::parallel_for(0UL, move_size, [&](const size_t j) {
            _move_weights[i][j] = phg.nodeWeight(_moves[i][j].hn);
          });
          parallel_prefix_sum(_move_weights[i].begin(), _move_weights[i].begin() + move_size, _move_weights[i].begin(), std::plus<HypernodeWeight>(), 0);
          // NOTE did you test if this gives the same result as the sequential implementation?
          last_move_idx = std::upper_bound(_move_weights[i].begin(), _move_weights[i].begin() + move_size, phg.partWeight(i) - _max_part_weights[i] - 1) - _move_weights[i].begin();
          ++last_move_idx;
        } else {
          HypernodeWeight sum = 0;
          for (; last_move_idx < move_size && sum <= phg.partWeight(i) - _max_part_weights[i] - 1; ++last_move_idx) {
            sum += phg.nodeWeight(_moves[i][last_move_idx].hn);
          }
        }
        if (phg.is_graph) {
          tbb::parallel_for(0UL, last_move_idx, [&](const size_t j) {
            const auto move = _moves[i][j];   // NOTE make this a const ref
            changeNodePart<true>(phg, move.hn, i, move.to, false);
          });
        } else {
          tbb::parallel_for(0UL, last_move_idx, [&](const size_t j) {
            const auto move = _moves[i][j];   // NOTE make this a const ref
            changeNodePart<false>(phg, move.hn, i, move.to, false);
          });
        }
      }
    }
  });
}

// explicitly instantiate so the compiler can generate them when compiling this cpp file
namespace {
#define DETERMINISTIC_REBALANCER(X) DeterministicRebalancer<X>
}

// explicitly instantiate so the compiler can generate them when compiling this cpp file
INSTANTIATE_CLASS_WITH_VALID_TRAITS(DETERMINISTIC_REBALANCER)
}

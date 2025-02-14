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

#include <tbb/parallel_sort.h>

namespace mt_kahypar {

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
bool DeterministicRebalancer<GraphAndGainTypes>::refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
  const vec<HypernodeID>&,
  Metrics& best_metrics,
  double) {
  PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
  resizeDataStructuresForCurrentK();
  if (_max_part_weights == nullptr) {
    _max_part_weights = &_context.partition.max_part_weights[0];
  }
  _gain_computation.reset();
  initializeDataStructures(phg);

  while (_num_imbalanced_parts > 0) {
    weakRebalancingRound(phg);
    HEAVY_REFINEMENT_ASSERT(checkPreviouslyOverweightParts(phg));
    updateImbalance(phg);
  }

  Gain delta = _gain_computation.delta();
  HEAVY_REFINEMENT_ASSERT(best_metrics.quality + delta == metrics::quality(phg, _context),
    V(best_metrics.quality) << V(delta) << V(metrics::quality(phg, _context)));
  best_metrics.quality += delta;
  best_metrics.imbalance = metrics::imbalance(phg, _context);
  DBG << "[REBALANCE] " << V(delta) << "  imbalance=" << best_metrics.imbalance;
  _max_part_weights = nullptr;
  return delta < 0;
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
  timer.start_timer("gain_computation", "Gain Computation");
  phg.doParallelForAllNodes([&](const HypernodeID hn) {
    const PartitionID from = phg.partID(hn);
    const HypernodeWeight weight = phg.nodeWeight(hn);
    if (imbalance(phg, from) > 0 && mayMoveNode(phg, from, weight)) {
      tmp_potential_moves[from].stream(computeGainAndTargetPart(phg, hn, true));
    }
  });
  timer.stop_timer("gain_computation");
  tbb::parallel_for(0UL, _moves.size(), [&](const size_t i) {
    _moves[i] = tmp_potential_moves[i].copy_parallel();
    if (_moves[i].size() > 0) {
      tbb::parallel_sort(_moves[i].begin(), _moves[i].end(), [&](const rebalancer::RebalancingMove& a, const rebalancer::RebalancingMove& b) {
        return a.priority < b.priority || (a.priority == b.priority && a.hn > b.hn);
      });
      // calculate perfix sum for each source-part to know which moves to execute (prefix_sum > current_weight - max_weight)
      _move_weights[i].resize(_moves[i].size());
      tbb::parallel_for(0UL, _moves[i].size(), [&](const size_t j) {
        _move_weights[i][j] = phg.nodeWeight(_moves[i][j].hn);
      });
      parallel_prefix_sum(_move_weights[i].begin(), _move_weights[i].end(), _move_weights[i].begin(), std::plus<HypernodeWeight>(), 0);
      const size_t last_move_idx = std::upper_bound(_move_weights[i].begin(), _move_weights[i].end(), phg.partWeight(i) - _max_part_weights[i] - 1) - _move_weights[i].begin();
      tbb::parallel_for(0UL, last_move_idx + 1, [&](const size_t j) {
        const auto move = _moves[i][j];
        changeNodePart(phg, move.hn, i, move.to, false);
      });
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

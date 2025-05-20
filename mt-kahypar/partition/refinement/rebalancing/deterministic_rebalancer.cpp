/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Robert Krause <robert.krause@student.kit.edu>
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

#include "mt-kahypar/partition/refinement/rebalancing/deterministic_rebalancer.h"

#include <tbb/parallel_for_each.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/cast.h"

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

  size_t iteration = 0;
  const size_t max_rounds = _context.refinement.rebalancing.det_max_rounds == 0 ? ABSOLUTE_MAX_ROUNDS : std::min(ABSOLUTE_MAX_ROUNDS, _context.refinement.rebalancing.det_max_rounds);
  while (_num_imbalanced_parts > 0 && iteration < max_rounds) {
    weakRebalancingRound(phg);
    ASSERT(checkPreviouslyOverweightParts(phg));
    updateImbalance(phg);
    ++iteration;
  }

  if constexpr (!PartitionedHypergraph::is_graph) {
    Gain delta = _gain_computation.delta();
    HEAVY_REFINEMENT_ASSERT(best_metrics.quality + delta == metrics::quality(phg, _context),
      V(best_metrics.quality) << V(delta) << V(metrics::quality(phg, _context)));
    best_metrics.quality += delta;
  }
  best_metrics.imbalance = metrics::imbalance(phg, _context);
  DBG << "[REBALANCE] " << "  imbalance=" << best_metrics.imbalance;
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
  for (PartitionID part = 0; part < _context.partition.k; ++part) {
    if (imbalance(phg, part) > 0) {
      ++_num_imbalanced_parts;
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
      if (to == _current_k) {
        to = 0;
      }
    } while (to != start);
    ASSERT(best_target == kInvalidPartition || (best_target >= 0 && best_target < _current_k));
  }
  tmp_scores.clear();
  return { hn, best_target, transformGain(best_gain, hn_weight) };
}

template <typename GraphAndGainTypes>
void DeterministicRebalancer<GraphAndGainTypes>::weakRebalancingRound(PartitionedHypergraph& phg) {
  ASSERT(_current_k == _context.partition.k && _current_k == phg.k());
  for (auto& moves : _tmp_potential_moves) {
    moves.clear_sequential();
  }
  for (PartitionID part = 0; part < _current_k; ++part) {
    _previous_block_imbalance[part] = imbalance(phg, part);
  }

  // calculate gain and target for each node in a overweight part
  // group moves by source part
  phg.doParallelForAllNodes([&](const HypernodeID hn) {
    const PartitionID from = phg.partID(hn);
    const HypernodeWeight weight = phg.nodeWeight(hn);
    if (_previous_block_imbalance[from] > 0 && mayMoveNode(phg, from, weight)) {
      const auto triple = computeGainAndTargetPart(phg, hn, true);
      if (from != triple.to && triple.to != kInvalidPartition) {
        _tmp_potential_moves[from].stream(triple);
      }
    }
  });

  tbb::parallel_for(0, _current_k, [&](const PartitionID part) {
    if (_tmp_potential_moves[part].size() > 0) {
      _moves[part] = _tmp_potential_moves[part].copy_parallel();
      const size_t move_size = _moves[part].size();
      ASSERT(move_size > 0);

      // sort the moves from each overweight part by priority
      tbb::parallel_sort(_moves[part], [&](const rebalancer::RebalancingMove& a, const rebalancer::RebalancingMove& b) {
        return a.priority < b.priority || (a.priority == b.priority && a.hn > b.hn);
      });

      // calculate perfix sum for each source-part to know which moves to execute (prefix_sum > current_weight - max_weight)
      const HypernodeWeight imbalance_of_block = _previous_block_imbalance[part];
      size_t last_move_idx = 0;
      bool success = false;  // whether the block is balanced afterwards
      if (move_size > _context.refinement.rebalancing.det_moves_sequential) {
        if (move_size > _move_weights[part].size()) {
          _move_weights[part].resize(move_size);
        }
        tbb::parallel_for(0UL, move_size, [&](const size_t j) {
          _move_weights[part][j] = phg.nodeWeight(_moves[part][j].hn);
        });
        parallel_prefix_sum(_move_weights[part].begin(), _move_weights[part].begin() + move_size, _move_weights[part].begin(), std::plus<HypernodeWeight>(), 0);
        last_move_idx = std::upper_bound(_move_weights[part].begin(), _move_weights[part].begin() + move_size, imbalance_of_block - 1) - _move_weights[part].begin();
        if (last_move_idx < move_size) {
          success = true;
          ++last_move_idx;
        }
      } else {
        HypernodeWeight sum = 0;
        for (; last_move_idx < move_size && sum < imbalance_of_block; ++last_move_idx) {
          sum += phg.nodeWeight(_moves[part][last_move_idx].hn);
        }
        success = (sum >= imbalance_of_block);
      }

      tbb::parallel_for(0UL, last_move_idx, [&](const size_t j) {
        const auto& move = _moves[part][j];
        ASSERT(move.to == kInvalidPartition || (move.to >= 0 && move.to < _current_k));
        changeNodePart(phg, move.hn, part, move.to, false);
      });
      _block_has_only_heavy_vertices[part] = static_cast<uint8_t>(!success);
    } else if (_previous_block_imbalance[part] > 0) {
      _block_has_only_heavy_vertices[part] = static_cast<uint8_t>(true);
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

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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/scalable_sort.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar {
static constexpr size_t ABSOLUTE_MAX_ROUNDS = 30;

static float transformGain(Gain gain_, HNWeightConstRef wu, HNWeightConstRef current_imbalance, HNWeightConstRef max_part_weight) {
  float gain = gain_;

  if (wu.dimension() == 1) {
    if (gain > 0) {
      gain /= wu.at(0);
    } else if (gain < 0) {
      gain *= wu.at(0);
    }
  } else {
    float relevant_weight_fraction = 0;
    for (Dimension d = 0; d < wu.dimension(); ++d) {
      if (current_imbalance.at(d) > 0) {
        relevant_weight_fraction += wu.at(d) / static_cast<float>(max_part_weight.at(d));
      }
    }
    if (gain > 0 && relevant_weight_fraction == 0) {
      gain = std::numeric_limits<float>::max();
    } else if (gain > 0) {
      gain /= relevant_weight_fraction;
    } else if (gain < 0) {
      gain *= relevant_weight_fraction;
    }
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
  updateImbalance(phg);
  _gain_computation.reset();

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
  phg.resetEdgeSynchronization();
  best_metrics.imbalance = metrics::imbalance(phg, _context);
  return true;
}

template <typename  GraphAndGainTypes>
void DeterministicRebalancer<GraphAndGainTypes>::updateImbalance(const PartitionedHypergraph& phg) {
  _num_imbalanced_parts = 0;
  for (PartitionID part = 0; part < _context.partition.k; ++part) {
    _current_imbalance[part] = phg.partWeight(part) - _context.partition.max_part_weights[part];
    if ( !(_current_imbalance[part] <= weight::broadcast(0, phg.dimension())) ) {
      ++_num_imbalanced_parts;
    }
  }
}

template <typename GraphAndGainTypes>
rebalancer::RebalancingMove DeterministicRebalancer<GraphAndGainTypes>::computeGainAndTargetPart(const PartitionedHypergraph& phg,
                                                                                                 const HypernodeID hn,
                                                                                                 bool non_adjacent_blocks) {
  const HNWeightConstRef hn_weight = phg.nodeWeight(hn);
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

  PartitionID from = phg.partID(hn);
  return { hn, best_target, transformGain(best_gain, hn_weight, _current_imbalance[from], _context.partition.max_part_weights[from])};
}

template <typename GraphAndGainTypes>
void DeterministicRebalancer<GraphAndGainTypes>::weakRebalancingRound(PartitionedHypergraph& phg) {
  ASSERT(_current_k == _context.partition.k && _current_k == phg.k());
  for (auto& moves : _tmp_potential_moves) {
    moves.clear_sequential();
  }

  // calculate gain and target for each node in a overweight part
  // group moves by source part
  phg.doParallelForAllNodes([&](const HypernodeID hn) {
    const PartitionID from = phg.partID(hn);
    const HNWeightConstRef weight = phg.nodeWeight(hn);
    if ( !(_current_imbalance[from] <= weight::broadcast(0, phg.dimension())) && mayMoveNode(phg, from, weight) && !phg.isFixed(hn)) {
      const auto triple = computeGainAndTargetPart(phg, hn, true);
      if (from != triple.to && triple.to != kInvalidPartition) {
        _tmp_potential_moves[from].stream(triple);
      }
    }
  });

  tbb::parallel_for(0, _current_k, [&](const PartitionID part) {
    if (_tmp_potential_moves[part].size() > 0) {
      _moves[part] = _tmp_potential_moves[part].copy_parallel();
      ASSERT(_moves[part].size() > 0);

      // sort the moves from each overweight part by priority
      parallel::scalable_sort(_moves[part], [](const rebalancer::RebalancingMove& a, const rebalancer::RebalancingMove& b) {
        return a.priority < b.priority || (a.priority == b.priority && a.hn > b.hn);
      });

      // determine which moves to execute
      size_t last_move_idx = 0;
      AllocatedHNWeight sum(phg.dimension(), 0);

      // TODO: this is not quite optimal for multiple dimensions, we should sort anew (or smth) when one dimension is balanced
      for (; last_move_idx < _moves[part].size() && !(sum >= _current_imbalance[part]); ++last_move_idx) {
        sum += phg.nodeWeight(_moves[part][last_move_idx].hn);
      }

      bool success = (sum >= _current_imbalance[part]);
      _block_has_only_heavy_vertices[part] = static_cast<uint8_t>(!success);

      // execute moves in parallel
      tbb::parallel_for(UL(0), last_move_idx, [&](const size_t j) {
        const auto& move = _moves[part][j];
        ASSERT(move.to == kInvalidPartition || (move.to >= 0 && move.to < _current_k));
        changeNodePart(phg, move.hn, part, move.to);
      });
    } else if ( !(_current_imbalance[part] <= weight::broadcast(0, phg.dimension())) ) {
      _block_has_only_heavy_vertices[part] = static_cast<uint8_t>(true);
    }
  });
}


template <typename GraphAndGainTypes>
void DeterministicRebalancer<GraphAndGainTypes>::resizeDataStructuresForCurrentK() {
  const Dimension dimension = _current_imbalance.dimension();
  // If the number of blocks changes, we resize data structures
  // (can happen during deep multilevel partitioning)
  if (_current_k != _context.partition.k) {
    _current_k = _context.partition.k;
    _gain_computation.changeNumberOfBlocks(_current_k);
    _moves.resize(_current_k);
    _tmp_potential_moves.resize(_current_k);
    _current_imbalance.replaceWith(_current_k, dimension, 0, false);
    _block_has_only_heavy_vertices.resize(_current_k);
  }
}

template <typename GraphAndGainTypes>
bool DeterministicRebalancer<GraphAndGainTypes>::mayMoveNode(const PartitionedHypergraph& phg, PartitionID part, HNWeightConstRef hn_weight) const {
  const auto weight_difference = weight::toNonAtomic(phg.partWeight(part)) - _context.partition.perfect_balance_part_weights[part];
  const auto allowed_weight = _context.refinement.rebalancing.det_heavy_vertex_exclusion_factor * weight_difference;
  return hn_weight <= allowed_weight || _block_has_only_heavy_vertices[part];
}

template <typename GraphAndGainTypes>
bool DeterministicRebalancer<GraphAndGainTypes>::isValidTarget(const PartitionedHypergraph& hypergraph, PartitionID part, HNWeightConstRef hn_weight) const {
  const HNWeightConstRef part_weight =  weight::toNonAtomic(hypergraph.partWeight(part));
  return (part_weight < deadzoneForPart(part)) &&
    part_weight + hn_weight <= _context.partition.max_part_weights[part];
}

template <typename GraphAndGainTypes>
bool DeterministicRebalancer<GraphAndGainTypes>::changeNodePart(PartitionedHypergraph& phg,
                                                                const HypernodeID hn,
                                                                const PartitionID from,
                                                                const PartitionID to) {
  // This function is passed as lambda to the changeNodePart function and used
  // to calculate the "real" delta of a move (in terms of the used objective function).
  auto objective_delta = [&](const SynchronizedEdgeUpdate& sync_update) {
      _gain_computation.computeDeltaForHyperedge(sync_update);
  };
  const auto inf_weight = weight::broadcast(std::numeric_limits<HNWeightScalar>::max(), phg.dimension());

  bool success = false;
  success = PartitionedHypergraph::is_graph ? phg.changeNodePartNoSync(hn, from, to) : phg.changeNodePart(hn, from, to, inf_weight, [] {}, objective_delta);
  ASSERT(success);
  return success;
}

template <typename GraphAndGainTypes>
bool DeterministicRebalancer<GraphAndGainTypes>::checkPreviouslyOverweightParts(const PartitionedHypergraph& phg) const {
  for (PartitionID i = 0; i < _current_k; ++i) {
    const auto partWeight = phg.partWeight(i);
    unused(partWeight);
    if ( !(_current_imbalance[i] <= weight::broadcast(0, phg.dimension())) ) {
        const auto& max_part_weights = _context.partition.max_part_weights;
        ASSERT(partWeight <= max_part_weights[i] || _block_has_only_heavy_vertices[i], V(partWeight) << V(max_part_weights[i]));
        unused(max_part_weights);
    }
  }
  return true;
}

// explicitly instantiate so the compiler can generate them when compiling this cpp file
namespace {
#define DETERMINISTIC_REBALANCER(X) DeterministicRebalancer<X>
}

// explicitly instantiate so the compiler can generate them when compiling this cpp file
INSTANTIATE_CLASS_WITH_VALID_TRAITS(DETERMINISTIC_REBALANCER)
}

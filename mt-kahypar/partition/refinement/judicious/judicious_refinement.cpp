/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@student.kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "mt-kahypar/partition/refinement/judicious/judicious_refinement.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/partition/metrics.h"
#include <algorithm>
#include <boost/range/irange.hpp>
#include <cmath>
#include <mt-kahypar/datastructures/priority_queue.h>

namespace mt_kahypar {
bool JudiciousRefiner::refineImpl(PartitionedHypergraph &phg,
                                  const vec<HypernodeID> &refinement_nodes,
                                  kahypar::Metrics &metrics, double) {

  unused(refinement_nodes);
  if (_reached_lower_bound)
    return false;
  if (!_is_initialized)
    throw std::runtime_error(
        "Call initialize on judicious refinement before calling refine");
  DBG << "Initial judicious load:" << V(metrics::judiciousLoad(phg));
  ASSERT(_last_load == metrics::judiciousLoad(phg) || _last_load == 0 || _rebalancing);
  _part_loads.clear();
  for (PartitionID i = 0; i < _context.partition.k; ++i) {
    _part_loads.insert(i, phg.partLoad(i));
  }
  const HyperedgeWeight initial_max_load = _part_loads.topKey();
  _num_bad_refinements = 0;
  bool done = false;
  do {
    const PartitionID from_block = _part_loads.top();
    calculateRefinementNodes(phg, from_block);
    const HyperedgeWeight max_load_before_refinement = phg.partLoad(from_block);
    HighResClockTimepoint refinement_start =
        std::chrono::high_resolution_clock::now();
    doRefinement(phg, from_block);
    HighResClockTimepoint refinement_stop =
        std::chrono::high_resolution_clock::now();
    double refinement_time =
        std::chrono::duration<double>(refinement_stop - refinement_start)
            .count();
    finalizeRefinementRound(phg, refinement_time, from_block,
                            max_load_before_refinement);
    done = shouldRefinementStop(phg, from_block, max_load_before_refinement);
  } while (!done);
  finalizeRefinement(phg, initial_max_load);
  const HyperedgeWeight judicious_load = _part_loads.topKey();
  HyperedgeWeight min_load = judicious_load;
  for (PartitionID i = 0; i < _context.partition.k; ++i) {
    if (phg.partLoad(i) > 0) {
      min_load = std::min(min_load, phg.partLoad(i));
    }
  }
  const double load_ratio =
      static_cast<double>(judicious_load) / min_load;
  if (load_ratio <= _context.refinement.judicious.min_load_ratio && !_rebalancing) {
    _rebalancing = true;
    rebalancePartition(phg, metrics);
  }
  _rebalancing = false;
  metrics.imbalance = metrics::imbalance(phg, _context);
  return initial_max_load - _part_loads.topKey() > 0;
}

void JudiciousRefiner::finalizeRefinementRound(
    const PartitionedHypergraph &phg, const double refinement_time,
    const PartitionID block, const HyperedgeWeight load_before) {
  _context.refinement.judicious.max_block_time =
      std::max(_context.refinement.judicious.max_block_time, refinement_time);
  DBG << "Reduced load of block" << block << "by"
      << (load_before - phg.partLoad(block)) << "[Time:" << refinement_time
      << "s]";
  for (PartitionID i = 0; i < _context.partition.k; ++i) {
    _part_loads.adjustKey(i, phg.partLoad(i));
  }
  if (block == _part_loads.top()) {
    DBG << RED << "Heaviest part has not changed" << END;
  }
}

bool JudiciousRefiner::shouldRefinementStop(const PartitionedHypergraph &phg,
                                            const PartitionID block,
                                            const HyperedgeWeight load_before) {
  const HyperedgeWeight current_max_load = _part_loads.topKey();
  HyperedgeWeight min_part_load = current_max_load;
  for (PartitionID i = 0; i < _context.partition.k; ++i) {
    if (phg.partLoad(i) > 0) {
      min_part_load = std::min(min_part_load, phg.partLoad(i));
    }
  }
  const double load_ratio =
      static_cast<double>(current_max_load) / min_part_load;
  DBG << V(load_ratio);
  HyperedgeWeight delta = load_before - current_max_load;
  if (delta <= 0) {
    _num_bad_refinements++;
  }
  return load_ratio < _context.refinement.judicious.min_load_ratio ||
         block == _part_loads.top() ||
         // NOTE: last case should never happen but I am paranoid
         _num_bad_refinements >= static_cast<size_t>(_context.partition.k);
}

void JudiciousRefiner::finalizeRefinement(
    const PartitionedHypergraph &phg, const HyperedgeWeight initial_max_load) {
  HyperedgeWeight max_load = phg.partLoad(0);
  HyperedgeWeight min_load = phg.partLoad(0);
  for (PartitionID i = 1; i < _context.partition.k; ++i) {
    max_load = std::max(max_load, phg.partLoad(i));
    min_load = std::min(min_load, phg.partLoad(i));
  }
  ASSERT(initial_max_load >= max_load);
  DBG << "improved judicious load by" << initial_max_load - max_load;
  DBG << V(metrics::judiciousLoad(phg));
  _last_load = max_load;
  if (static_cast<HyperedgeID>(max_load) ==
          _context.refinement.judicious.max_degree ||
      max_load == min_load) {
    _reached_lower_bound = true;
    DBG << "Reached lower bound" << V(max_load) << V(min_load)
        << V(_context.refinement.judicious.max_degree);
  }
}

void JudiciousRefiner::initializeImpl(PartitionedHypergraph &phg) {
  if (!phg.isGainCacheInitialized()) {
    phg.initializeGainCache();
  }
  _is_initialized = true;
}

void JudiciousRefiner::calculateRefinementNodes(
    const PartitionedHypergraph &phg, const PartitionID p) {
  _refinement_nodes.clear();
  tbb::enumerable_thread_specific<vec<HypernodeID>> ets_refinement_nodes;

  // thread local refinement node calculation
  tbb::parallel_for(tbb::blocked_range<HypernodeID>(0, phg.initialNumNodes()),
                    [&](const tbb::blocked_range<HypernodeID> &r) {
                      auto &tl_refinement_nodes = ets_refinement_nodes.local();
                      // tl_refinement_nodes.reserve(phg.partWeight(p));
                      for (HypernodeID u = r.begin(); u < r.end(); ++u) {
                        if (phg.nodeIsEnabled(u) && phg.partID(u) == p) {
                          tl_refinement_nodes.push_back(u);
                        }
                      }
                    });

  for (const auto &tl_refinement_nodes : ets_refinement_nodes) {
    _refinement_nodes.insert(_refinement_nodes.end(),
                             tl_refinement_nodes.begin(),
                             tl_refinement_nodes.end());
  }
}

void JudiciousRefiner::doRefinement(PartitionedHypergraph &phg,
                                    const PartitionID part_id) {
  DBG << V(_refinement_nodes.size());
  _pq.init(phg, _refinement_nodes, part_id);
  // do not need load of current from-block in the PQ as it only induces updates
  // with info that is not needed
  _part_loads.deleteTop();
  const HyperedgeWeight initial_from_load = phg.partLoad(part_id);
  HyperedgeWeight from_load = initial_from_load;
  auto delta_func = [&](const HyperedgeID he, const HyperedgeWeight,
                        const HypernodeID,
                        const HypernodeID pin_count_in_from_part_after,
                        const HypernodeID pin_count_in_to_part_after) {
    // Gains of the pins of a hyperedge can only change in the following
    // situations.
    if (pin_count_in_from_part_after == 1 || pin_count_in_to_part_after == 1) {
      _edges_with_gain_changes.push_back(he);
    }
  };
  Move move;
  vec<HypernodeID> moved_nodes;
  while (_pq.findNextMove(phg, move)) {
    ASSERT(move.to != kInvalidPartition);
    ASSERT(move.from == part_id);
    phg.changeNodePartWithGainCacheUpdate(
        move.node, move.from, move.to,
        std::numeric_limits<HypernodeWeight>::max(), [] {}, delta_func);

    const HyperedgeWeight to_load = phg.partLoad(move.to);
    from_load = phg.partLoad(move.from);
    _part_loads.adjustKey(move.to, to_load);
    // want to accept all moves which do not increase the judicious load, i.e
    // all moves s.t. to's load will not exceed from's load before the move.
    // Moves that increase the judicious load should not be done, because they
    // cannot help escape local minima, because if a node is moved from node A
    // to B to C, it could be moved to C directly and its neighbors which
    // decrease the load of B are moved when B becomes the heaviest
    // TLDR:
    // escaping local minima would require some extra tricks, moving nodes to
    // heavy blocks and hoping for improvements down the line is not feasible;
    moved_nodes.push_back(move.node);
    if (move.gain < 0) {
      phg.changeNodePartWithGainCacheUpdate(move.node, move.to, move.from);
      DBG << "Negative gain move which would increase judicious load";
      break;
    }
    ASSERT(phg.partLoad(move.to) <= initial_from_load);

    if (_part_loads.topKey() >=
        from_load * _context.refinement.judicious.part_load_margin)
      break;
    updateNeighbors(phg, move);
  }
  _edges_with_gain_changes.clear();
  _pq.reset();
  _part_loads.insert(part_id, from_load);
  tbb::parallel_for(0UL, moved_nodes.size(), [&](const MoveID i) {
    phg.recomputeMoveFromBenefit(moved_nodes[i]);
  });
}

void JudiciousRefiner::updateNeighbors(PartitionedHypergraph &phg,
                                       const Move &move) {
  for (const HyperedgeID &he : _edges_with_gain_changes) {
    for (const HypernodeID &v : phg.pins(he)) {
      if (phg.partID(v) == move.from &&
          _gain_update_state[v] != _gain_update_time) {
        _pq.updateGain(phg, v);
        _gain_update_state[v] = _gain_update_time;
      }
    }
  }
  _edges_with_gain_changes.clear();
  if (++_gain_update_time == 0) {
    _gain_update_state.assign(_gain_update_state.size(), 0);
    _gain_update_time = 1;
  }
}

void JudiciousRefiner::rebalancePartition(PartitionedHypergraph &phg, kahypar::Metrics &metrics) {
  ASSERT(mt_kahypar::metrics::judiciousLoad(phg) == _part_loads.topKey());
  const HyperedgeWeight initial_judicious_load = _part_loads.topKey();
  DBG << "Rebalancing" << V(initial_judicious_load);

  vec<PartitionID> part_ids(phg.initialNumNodes());
  for (size_t i = 0; i < phg.initialNumNodes(); ++i) {
    part_ids[i] = phg.partID(i);
  }

  ds::ExclusiveHandleHeap<ds::Heap<HyperedgeWeight, HypernodeID>> parts(_context.partition.k);
  for (PartitionID p = 0; p < _context.partition.k; ++p) {
    parts.insert(p, phg.partLoad(p));
  }

  auto bestToBlock = [&](const HypernodeID v) {
    auto targets = boost::irange(0, _context.partition.k);
    return *std::min_element(targets.begin(), targets.end(), [&phg, v] (const PartitionID a, const PartitionID b) {
      if (phg.partID(v) == a) return false;
      if (phg.partID(v) == b) return true;
      return phg.moveToPenalty(v, a) < phg.moveToPenalty(v, b);
    });
  };

  ds::ExclusiveHandleHeap<ds::MaxHeap<HyperedgeWeight, HypernodeID>> benefit_pq(
      phg.initialNumNodes());
  const size_t num_rounds = parts.size() / 2;
  for (PartitionID i = 0; i < static_cast<PartitionID>(num_rounds); ++i) {
    PartitionID p = parts.top();
    parts.deleteTop();
    calculateRefinementNodes(phg, p);
    benefit_pq.clear();
    _edges_with_gain_changes.clear();
    for (const auto v : _refinement_nodes) {
      benefit_pq.insert(v, phg.moveFromBenefit(v));
    }

    auto delta_func = [&](const HyperedgeID he, const HyperedgeWeight,
                          const HypernodeID,
                          const HypernodeID pin_count_in_from_part_after,
                          const HypernodeID) {
      if (pin_count_in_from_part_after == 1) {
        _edges_with_gain_changes.push_back(he);
      }
    };

    auto update_neighbors = [&]() {
      for (const HyperedgeID &he : _edges_with_gain_changes) {
        for (const HypernodeID &v : phg.pins(he)) {
          if (_gain_update_state[v] != _gain_update_time && benefit_pq.contains(v)) {
            benefit_pq.adjustKey(v, phg.moveFromBenefit(v));
            _gain_update_state[v] = _gain_update_time;
          }
        }
      }
      _edges_with_gain_changes.clear();
      if (++_gain_update_time == 0) {
        _gain_update_state.assign(_gain_update_state.size(), 0);
        _gain_update_time = 1;
      }
    };

    size_t j = 0;
    while (j < _refinement_nodes.size() * 0.1) {
      if (benefit_pq.empty()) break;
      const HypernodeID v = benefit_pq.top();
      const HyperedgeWeight benefit = benefit_pq.topKey();
      benefit_pq.deleteTop();
      PartitionID to = bestToBlock(v);
      if (benefit - phg.moveToPenalty(v, to) >= 0) {
        phg.changeNodePartWithGainCacheUpdate(v, p, to,
                                              std::numeric_limits<HypernodeWeight>::max(),
                                              [] {},
                                              delta_func);
        ++j;
        if (parts.contains(to)) {
          parts.adjustKey(to, phg.partLoad(to));
        }
        update_neighbors();
      }
    }
    DBG << "Moved" << V(j) << "nodes during rebalancing";
  }
  refine(phg, {}, metrics, 0);
  if (initial_judicious_load < _part_loads.topKey()) {
    DBG << RED << "Did not improve judicious load after rebalancing! Reverting changes!" << END;
    phg.resetPartition();
    for (size_t i = 0; i < phg.initialNumNodes(); ++i) {
      phg.setOnlyNodePart(i, part_ids[i]);
    }
    phg.initializePartition();
    _last_load = initial_judicious_load;
  } else {
    DBG << "Rebalancing improved judicious load by" << initial_judicious_load - _part_loads.topKey();
  }
  _edges_with_gain_changes.clear();
}
} // namespace mt_kahypar

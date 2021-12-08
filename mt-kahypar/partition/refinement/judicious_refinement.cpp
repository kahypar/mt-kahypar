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

#include "mt-kahypar/partition/refinement/judicious_refinement.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include <cmath>
#include <mt-kahypar/datastructures/priority_queue.h>
#include "mt-kahypar/partition/metrics.h"

namespace mt_kahypar {
  bool JudiciousRefiner::refineImpl(
              PartitionedHypergraph& phg,
              const vec<HypernodeID>& refinement_nodes,
              kahypar::Metrics& metrics,
              double) {

    unused(refinement_nodes);
    if (!_is_initialized) throw std::runtime_error("Call initialize on judicious refinement before calling refine");
    if (debug) {
      LOG << "Initial judicious load: " << V(metrics::judiciousLoad(phg));
      ASSERT(_last_load == metrics::judiciousLoad(phg) || _last_load == 0);
    }
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      _part_loads.insert(i, phg.partLoad(i));
    }
    const HyperedgeWeight initial_max_load = _part_loads.topKey();
    HyperedgeWeight current_max_load = initial_max_load;
    size_t num_bad_refinements = 0;
    bool done = false;
    while (!done) {
      calculateRefinementNodes(phg);
      const PartitionID heaviest_part = _part_loads.top();
      const Gain last_best_improvement = _best_improvement;
      doRefinement(phg, heaviest_part);
      if (debug) {
        LOG << "Improved best state by " << (_best_improvement - last_best_improvement);
      }
      for (PartitionID i = 0; i < _context.partition.k; ++i) {
        _part_loads.adjustKey(i, phg.partLoad(i));
      }
      if (debug && heaviest_part == _part_loads.top()) {
        LOG << RED << "Heaviest part has not changed" << END;
      }
      current_max_load = _part_loads.topKey();
      _total_improvement = initial_max_load - current_max_load;
      HyperedgeWeight min_part_load = current_max_load;
      for (PartitionID i = 0; i < _context.partition.k; ++i) {
        min_part_load = std::min(min_part_load, phg.partLoad(i));
      }
      const double load_ratio = static_cast<double>(current_max_load) / min_part_load;
      HyperedgeWeight delta = _best_improvement - last_best_improvement;
      if (delta <= 0 || heaviest_part == _part_loads.top()) {
        num_bad_refinements++;
      } else {
        num_bad_refinements = 0;
      }
      if (load_ratio < _context.refinement.judicious.min_load_ratio || num_bad_refinements >= 2) {   // (Review Note) This alone will not suffice as stopping criterion. must also include whether heaviest block yielded improvement
        done = true;
      }
    }
    revertToBestLocalPrefix(phg, 0);
    current_max_load = phg.partLoad(0);
    for (PartitionID i = 1; i < _context.partition.k; ++i) {
      current_max_load = std::max(current_max_load, phg.partLoad(i));
    }
    ASSERT(initial_max_load >= current_max_load);
    ASSERT(_best_improvement == initial_max_load - current_max_load);
    if (debug) {
      LOG << "improved judicious load by " << initial_max_load - current_max_load;
      _last_load = metrics::judiciousLoad(phg);
      LOG << V(metrics::judiciousLoad(phg));
    }
    metrics.imbalance = metrics::imbalance(phg, _context);
    reset();
    return initial_max_load - current_max_load > 0;
  }

  void JudiciousRefiner::initializeImpl(PartitionedHypergraph& phg) {
    if ( !phg.isGainCacheInitialized()) {
      phg.initializeGainCache();
    }
    _is_initialized = true;

  }

  void JudiciousRefiner::reset() {
    _best_improvement = 0;
    _total_improvement = 0;
    _moves.clear();
    _part_loads.clear();
  }

  void JudiciousRefiner::calculateRefinementNodes(PartitionedHypergraph& phg) {
    for (auto& b : _refinement_nodes) {
      b.clear();
    }
    tbb::enumerable_thread_specific<vec<vec<HypernodeID>>> ets_refinement_nodes;

    // thread local refinement node calculation
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(0, phg.initialNumNodes()),
                      [&](const tbb::blocked_range<HypernodeID> &r) {
                        auto &tl_refinement_nodes = ets_refinement_nodes.local();
                        tl_refinement_nodes.resize(_context.partition.k);
                        for (HypernodeID u = r.begin(); u < r.end(); ++u) {
                          if (phg.nodeIsEnabled(u)) {
                            tl_refinement_nodes[phg.partID(u)].push_back(u);
                          }
                        }
                      });

    for (const auto &tl_refinement_nodes : ets_refinement_nodes) {
      tbb::parallel_for(PartitionID(0), _context.partition.k, [&](const auto i) {
        _refinement_nodes[i].insert(_refinement_nodes[i].end(), tl_refinement_nodes[i].begin(),
                               tl_refinement_nodes[i].end());
      });
    }
  }

  void JudiciousRefiner::doRefinement(PartitionedHypergraph& phg, PartitionID part_id) {
    auto& refinement_nodes = _refinement_nodes[part_id];
    _gain_cache.setActivePart(part_id);
    for (HypernodeID v : refinement_nodes) {
      _gain_cache.insert(phg, v);
    }
    _part_loads.deleteTop();
    const HyperedgeWeight initial_from_load = phg.partLoad(part_id);
    HyperedgeWeight from_load = initial_from_load;
    // disable to-Blocks that are too large
    _gain_cache.initBlockPQ(phg, initial_from_load);
    auto delta_func = [&](const HyperedgeID he,
                          const HyperedgeWeight,
                          const HypernodeID,
                          const HypernodeID pin_count_in_from_part_after,
                          const HypernodeID pin_count_in_to_part_after) {
      // Gains of the pins of a hyperedge can only change in the following situations.
      if (pin_count_in_from_part_after == 0 || pin_count_in_from_part_after == 1 ||
          pin_count_in_to_part_after == 1 || pin_count_in_to_part_after == 2) {
        _edgesWithGainChanges.push_back(he);
      }
    };
    Move move;
    bool done = false;
    size_t initial_num_moves = _moves.size();
    vec<HypernodeID> move_nodes;
    while (!done) {
      JudiciousGainCache::pqStatus status = _gain_cache.findNextMove(phg, move);
      if (status == JudiciousGainCache::pqStatus::empty) {
        if (debug) {
          LOG << "Abort due to empty PQ";
        }
        break;
      } else if (status == JudiciousGainCache::pqStatus::rollback) {
        if (debug) {
          LOG << "Did rollback";
        }
        revertToBestLocalPrefix(phg, initial_num_moves, true);
        continue;
      }
      if (move.to == kInvalidPartition) {
        continue;
      }
      ASSERT(move.from == part_id);
      phg.changeNodePartWithGainCacheUpdate(move.node, move.from, move.to,
                                            std::numeric_limits<HypernodeWeight>::max(),
                                            []{}, delta_func);

      _moves.push_back(move);
      const HyperedgeWeight to_load = phg.partLoad(move.to);
      from_load = phg.partLoad(move.from);
      _gain_cache.updateEnabledBlocks(move.to, from_load, to_load);
      _part_loads.adjustKey(move.to, to_load);
      Gain gain = initial_from_load - std::max(_part_loads.topKey(), from_load);
      if (_total_improvement + gain >= _best_improvement) {
        _best_improvement = _total_improvement + gain;
        for (size_t i = initial_num_moves; i < _moves.size(); ++i) {
          move_nodes.push_back(_moves[i].node);
        }
        _moves.clear();
        initial_num_moves = 0;
      }
      if (_moves.size() > refinement_nodes.size() * _context.refinement.judicious.abort_factor) {
        if (debug) {
          LOG << "Abort due to too many negative gain moves";
        }
        revertToBestLocalPrefix(phg, initial_num_moves);
        done = true;
      } else if (_part_loads.topKey() >= std::max(from_load, _part_loads.keyOfSecond()) * _context.refinement.judicious.part_load_margin) {
        done = true;
      } else {
        updateNeighbors(phg, move);
      }
    }
    _edgesWithGainChanges.clear();
    _gain_cache.resetGainCache();
    _part_loads.insert(part_id, from_load);
    for (size_t i = initial_num_moves; i < _moves.size(); ++i) {
      move_nodes.push_back(_moves[i].node);
    }
    tbb::parallel_for(0UL, move_nodes.size(), [&](const MoveID i) {
      phg.recomputeMoveFromBenefit(move_nodes[i]);
    });
  }

  void JudiciousRefiner::updateNeighbors(PartitionedHypergraph& phg, const Move& move) {
    for (HyperedgeID e : _edgesWithGainChanges) {
      if (phg.edgeSize(e) < _context.partition.ignore_hyperedge_size_threshold) {
        for (HypernodeID v : phg.pins(e)) {
          if (_neighbor_deduplicator[v] != _deduplication_time) {
            if (phg.partID(v) == move.from) {
              _gain_cache.updateGain(phg, v, move);
            }
            _neighbor_deduplicator[v] = _deduplication_time;
          }
        }
      }
    }
    _edgesWithGainChanges.clear();
    if (++_deduplication_time == 0) {
      _neighbor_deduplicator.assign(_neighbor_deduplicator.size(), 0);
      _deduplication_time = 1;
    }
  }

  void JudiciousRefiner::revertToBestLocalPrefix(PartitionedHypergraph& phg, size_t bestGainIndex, bool update_gain_cache) {
    auto delta_func = [&](const HyperedgeID he,
                          const HyperedgeWeight,
                          const HypernodeID,
                          const HypernodeID pin_count_in_from_part_after,
                          const HypernodeID pin_count_in_to_part_after) {
      // Gains of the pins of a hyperedge can only change in the following situations.
      if (pin_count_in_from_part_after == 0 || pin_count_in_from_part_after == 1 ||
          pin_count_in_to_part_after == 1 || pin_count_in_to_part_after == 2) {
        _edgesWithGainChanges.push_back(he);
      }
    };
    // (Review Note) bestGainIndex is still for km1 gain. At this point we want to look at whether we actually improve the max load
    if (debug) {
      LOG << "reverting" << (_moves.size() - bestGainIndex) << "moves";
    }
    while (_moves.size() > bestGainIndex) {
      Move& m = _moves.back();
      if (update_gain_cache) {
        phg.changeNodePartWithGainCacheUpdate(m.node, m.to, m.from, std::numeric_limits<HypernodeWeight>::max(), []{}, delta_func);
        _gain_cache.insert(phg, m.node);
        updateNeighbors(phg, m);
      } else {
        phg.changeNodePartWithGainCacheUpdate(m.node, m.to, m.from);
      }
      m.invalidate();
      _moves.pop_back();
    }
  }
}

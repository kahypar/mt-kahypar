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
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      _part_loads.insert(i, phg.partLoad(i));
    }
    const HyperedgeWeight initial_max_load = _part_loads.topKey();
    HyperedgeWeight current_max_load = initial_max_load;
    _best_improvement_index = 0;
    _best_improvement = 0;
    _estimated_improvement = 0;
    size_t num_negative_refinements = 0;
    bool done = false;
    while (!done) {
      calculateRefinementNodes(phg);
      const PartitionID heaviest_part = _part_loads.top();
      const Gain last_best_improvement = _best_improvement;
      doRefinement(phg, heaviest_part);
      for (PartitionID i = 0; i < _context.partition.k; ++i) {
        _part_loads.adjustKey(i, phg.partLoad(i));
      }
      current_max_load = _part_loads.topKey();
      HyperedgeWeight min_part_load = current_max_load;
      for (PartitionID i = 0; i < _context.partition.k; ++i) {
        min_part_load = std::min(min_part_load, phg.partLoad(i));
      }
      const double load_ratio = static_cast<double>(current_max_load) / min_part_load;
      /*! TODO: maybe only abort the second time this happens
       *  \todo maybe only abort the second time this happens
       */
      HyperedgeWeight delta = _best_improvement - last_best_improvement;
      if (delta <= 0) {
        num_negative_refinements++;
      } else {
        num_negative_refinements = 0;
      }
      if (load_ratio < _min_load_ratio || num_negative_refinements >= 2) {   // (Review Note) This alone will not suffice as stopping criterion. must also include whether heaviest block yielded improvement
        done = true;
      }
    }
    revertToBestLocalPrefix(phg, _best_improvement_index);
    _moves.clear();
    if (debug) {
      LOG << "improved judicious load by " << _best_improvement;
      LOG << V(metrics::judiciousLoad(phg));
    }
    _part_loads.clear();
    metrics.imbalance = metrics::imbalance(phg, _context);
    return current_max_load < initial_max_load;
  }

  void JudiciousRefiner::initializeImpl(PartitionedHypergraph& phg) {
    if ( !phg.isGainCacheInitialized()) {
      phg.initializeGainCache();
    }
    _is_initialized = true;

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
    for (HypernodeID v : refinement_nodes) {
      /*! TODO: maybe cut of at specific #nodes
       *  \todo maybe cut of at specific #nodes
       */
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
    while (!done) {
      JudiciousGainCache::pqStatus status = _gain_cache.findNextMove(phg, move);
      if (status == JudiciousGainCache::pqStatus::empty) done = true;
      else if (status == JudiciousGainCache::pqStatus::rollback) {
        if (debug) {
          LOG << "Did rollback";
        }
        revertToBestLocalPrefix(phg, _best_improvement_index, true);
        _moves.clear();
        _best_improvement_index = 0;
        _estimated_improvement = 0;
        _best_improvement = 0;
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
      const HyperedgeWeight new_to_load = phg.partLoad(move.to);
      from_load = phg.partLoad(move.from);
      _gain_cache.updateEnabledBlocks(move.to, from_load, new_to_load);
      //_part_loads.adjustKey(move.from, from_load);
      _part_loads.adjustKey(move.to, new_to_load);
      _estimated_improvement = initial_from_load - std::max(_part_loads.topKey(), from_load);
      if (_estimated_improvement >= _best_improvement) {
        _best_improvement = _estimated_improvement;
        _best_improvement_index = _moves.size();
      }
      if (_part_loads.topKey() >= from_load * _part_load_margin) {
        done = true;
      } else {
        updateNeighbors(phg, move);
      }
    }
    _gain_cache.resetGainCache();
    _part_loads.insert(part_id, from_load);
    tbb::parallel_for(0UL, _moves.size(), [&](const MoveID i) {
      phg.recomputeMoveFromBenefit(_moves[i].node);
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
    while (_moves.size() > bestGainIndex) {
      Move& m = _moves.back();
      if (update_gain_cache) {
        phg.changeNodePartWithGainCacheUpdate(m.node, m.to, m.from, std::numeric_limits<HypernodeWeight>::max(), []{}, delta_func);
        updateNeighbors(phg, m);
      } else {
        phg.changeNodePartWithGainCacheUpdate(m.node, m.to, m.from);
      }
      m.invalidate();
      _moves.pop_back();
    }
  }
}

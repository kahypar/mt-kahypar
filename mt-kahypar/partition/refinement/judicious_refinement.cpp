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
    /*! TODO: keep moved nodes and distibute them onto buckets again
     *  \todo keep moved nodes and distibute them onto buckets again
     */
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      _part_loads.insert(i, phg.partLoad(i));
    }
    bool done = false;
    Gain overall_improvement = 0;
    while (!done) {
      calculateRefinementNodes(phg);
      const PartitionID heaviest_part = _part_loads.top();
      overall_improvement += doRefinement(phg, heaviest_part);
      const HyperedgeWeight max_part_load = _part_loads.topKey();
      HyperedgeWeight min_part_load = max_part_load;
      for (PartitionID i = 0; i < _context.partition.k; ++i) {
        min_part_load = std::min(min_part_load, phg.partLoad(i));
      }
      const double load_ratio = static_cast<double>(max_part_load) / min_part_load;
      if (load_ratio < _min_load_ratio) {
        done = true;
      }
    }
    _part_loads.clear();
    _move_status.assign(_hypergraph.initialNumNodes(), false);
    metrics.km1 -= overall_improvement;
    metrics.imbalance = metrics::imbalance(phg, _context);
    ASSERT(metrics.km1 == metrics::km1(phg), V(metrics.km1) << V(metrics::km1(phg)));
    return overall_improvement > 0;
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

  Gain JudiciousRefiner::doRefinement(PartitionedHypergraph& phg, PartitionID part_id) {
    auto& refinement_nodes = _refinement_nodes[part_id];
    for (HypernodeID v : refinement_nodes) {
      /*! TODO: maybe cut of at specific #nodes
       *  \todo maybe cut of at specific #nodes
       */
      if (!_move_status[v]) {
        _gain_cache.insert(phg, v);
      }
    }
    // disable to-Blocks that are too large
    const HyperedgeWeight from_load = _part_loads.topKey();
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      _gain_cache.updateEnabledBlocks(part_id, from_load, phg.partLoad(i));
    }
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
    size_t bestImprovementIndex = 0;
    Gain estimatedImprovement = 0;
    Gain bestImprovement = 0;
    // set to true by gain cache if enabled PQs were exhausted and we need to perform a rollback before enabling all other PQs
    bool should_refiner_perform_rollback = false;
    while (!done && _gain_cache.findNextMove(phg, move, should_refiner_perform_rollback)) {
      if (should_refiner_perform_rollback) {
        revertToBestLocalPrefix(phg, bestImprovementIndex, true);
        should_refiner_perform_rollback = false;
        _moves.clear();
        bestImprovementIndex = 0;
        continue;
      }
      bool moved = false;
      if (move.to != kInvalidPartition) {
        ASSERT(!_move_status[move.node]);
        ASSERT(move.from == part_id);
        moved = phg.changeNodePartWithGainCacheUpdate(move.node, move.from, move.to,
                                                      std::numeric_limits<HypernodeWeight>::max(),
                                                      []{}, delta_func);

      }
      if (moved) {
        _move_status[move.node] = true;
        estimatedImprovement += move.gain;
        _moves.push_back(move);
        const HyperedgeWeight new_to_load = phg.partLoad(move.to);
        const HyperedgeWeight new_from_load = phg.partLoad(move.from);
        _gain_cache.updateEnabledBlocks(move.to, new_from_load, new_to_load);
        if (new_to_load >= new_from_load * _part_load_margin) {
          done = true;
        }
        if (estimatedImprovement > bestImprovement) {
          bestImprovement = estimatedImprovement;
          bestImprovementIndex = _moves.size();
        }
        updateNeighbors(phg, move);
      }
    }
    revertToBestLocalPrefix(phg, bestImprovementIndex);
    _gain_cache.resetGainCache();
    tbb::parallel_for(0UL, _moves.size(), [&](const MoveID i) {
      phg.recomputeMoveFromBenefit(_moves[i].node);
    });
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      _part_loads.adjustKey(i, phg.partLoad(i));
    }
    _moves.clear();
    return bestImprovement;
  }

  void JudiciousRefiner::updateNeighbors(PartitionedHypergraph& phg, const Move& move) {
    for (HyperedgeID e : _edgesWithGainChanges) {
      if (phg.edgeSize(e) < _context.partition.ignore_hyperedge_size_threshold) {
        for (HypernodeID v : phg.pins(e)) {
          if (_neighbor_deduplicator[v] != _deduplication_time) {
            if (phg.partID(v) == move.from && !_move_status[v]) {
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

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
    Gain overall_improvement = 0;
    /*! TODO: keep moved nodes and distibute them onto buckets again
     *  \todo keep moved nodes and distibute them onto buckets again
     */
    calculateRefinementNodes(phg);
    for (size_t i = 0; i < 1; i++) {
      for (PartitionID i = 0; i < _context.partition.k; ++i) {
        _part_weights.insert(i, phg.partWeight(i));
      }
      bool done = false;
      Gain improvement = 0;
      Gain last_improvement = 0;
      size_t consecutive_runs_with_too_little_improvement = 0;
      while (!done) {
        const PartitionID heaviest_part = _part_weights.top();
        last_improvement = improvement;
        improvement += doRefinement(phg, heaviest_part);
        const double improvement_fraction = last_improvement == 0 ? 1.0 : 1.f * improvement / last_improvement;
        if (improvement_fraction < _min_improvement) {
          consecutive_runs_with_too_little_improvement++;
        } else {
          consecutive_runs_with_too_little_improvement = 0;
        }
        done = consecutive_runs_with_too_little_improvement > 2;
      }
      overall_improvement += improvement;
      _part_weights.clear();
      _move_status.assign(_hypergraph.initialNumNodes(), false);
    }
    metrics.km1 -= overall_improvement;
    /*LOG << V(overall_improvement);*/
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
    const HypernodeWeight from_weight = _part_weights.topKey();
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      _gain_cache.updateEnabledBlocks(part_id, from_weight, phg.partWeight(i));
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
    size_t bestImprovementIndex = 0;
    Gain estimatedImprovement = 0;
    Gain bestImprovement = 0;
    bool done = false;
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
        const HypernodeWeight new_to_weight = phg.partWeight(move.to);
        const HypernodeWeight new_from_weight = phg.partWeight(move.from);
        _gain_cache.updateEnabledBlocks(move.to, new_from_weight, new_to_weight);
        if (new_to_weight >= new_from_weight * _part_weight_margin) {
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
      _part_weights.adjustKey(i, phg.partWeight(i));
    }
    _moves.clear();
    return bestImprovement;
  }

  void JudiciousRefiner::updateNeighbors(PartitionedHypergraph& phg, const Move& move) {
    for (HyperedgeID e : _edgesWithGainChanges) {
      /*! TODO: need to ignore large edges when doing this sequential?
       *  \todo need to ignore large edges when doing this sequential?
       */
      if (phg.edgeSize(e) < _context.partition.ignore_hyperedge_size_threshold) {
        for (HypernodeID v : phg.pins(e)) {
          if (phg.partID(v) == move.from && !_move_status[v]) {
            _gain_cache.updateGain(phg, v, move);
          }
        }
      }
    }
    _edgesWithGainChanges.clear();
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

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
namespace mt_kahypar {
  bool JudiciousRefiner::refineImpl(
              PartitionedHypergraph& phg,
              const vec<HypernodeID>& refinement_nodes,
              kahypar::Metrics& metrics,
              double) {

    unused(refinement_nodes);
    unused(metrics);
    if (!_is_initialized) throw std::runtime_error("Call initialize on judicious refinement before calling refine");
    calculateBorderNodes(phg);
    ds::ExclusiveHandleHeap<ds::MaxHeap<HypernodeWeight, PartitionID>> part_weights(static_cast<size_t>(_context.partition.k));
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      part_weights.insert(i, phg.partWeight(i));
    }
    bool done = false;
    while (!done) {
      const PartitionID heaviest_part = part_weights.top();
      auto result = doRefinement(phg, heaviest_part);
      done = !result.first || 0;
    }
    return true;
  }

  void JudiciousRefiner::initializeImpl(PartitionedHypergraph& phg) {
    if ( !phg.isGainCacheInitialized()) {
      phg.initializeGainCache();
    }
    _is_initialized = true;

  }

  void JudiciousRefiner::calculateBorderNodes(PartitionedHypergraph& phg) {
    tbb::enumerable_thread_specific<vec<vec<HypernodeID>>> ets_border_nodes;

    // thread local border node calculation
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(0, phg.initialNumNodes()),
                      [&](const tbb::blocked_range<HypernodeID> &r) {
                        auto &tl_border_nodes = ets_border_nodes.local();
                        tl_border_nodes.resize(_context.partition.k);
                        for (HypernodeID u = r.begin(); u < r.end(); ++u) {
                          if (phg.nodeIsEnabled(u) && phg.isBorderNode(u)) {
                            tl_border_nodes[phg.partID(u)].push_back(u);
                          }
                        }
                      });

    _border_nodes.resize(_context.partition.k);

    for (const auto &tl_border_nodes : ets_border_nodes) {
      tbb::parallel_for(PartitionID(0), _context.partition.k, [&](const auto i) {
        _border_nodes[i].insert(_border_nodes[i].end(), tl_border_nodes[i].begin(),
                               tl_border_nodes[i].end());
      });
    }
  }

  std::pair<bool, Gain> JudiciousRefiner::doRefinement(PartitionedHypergraph& phg, PartitionID part_id) {
    auto& refinement_nodes = _border_nodes[part_id];
    for (HypernodeID v : refinement_nodes) {
      _gain_cache.insert(phg, v);
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
    while (!done && _gain_cache.findNextMove(phg, move)) {
      bool moved = false;
      if (move.to != kInvalidPartition) {
        moved = phg.changeNodePartWithGainCacheUpdate(move.node, move.from, move.to,
                                                      std::numeric_limits<HypernodeWeight>::max(),
                                                      []{}, delta_func);

      }
      if (moved) {
        estimatedImprovement += move.gain;
        _moves.push_back(move);
        if (estimatedImprovement > bestImprovement) {
          bestImprovement = estimatedImprovement;
          bestImprovementIndex = _moves.size();
        }
        updateNeighbors(phg, move);
      }
    }
    revertToBestLocalPrefix(phg, bestImprovementIndex);
    _gain_cache.clearPQs();
    return std::make_pair(true, bestImprovement);
  }

  void JudiciousRefiner::updateNeighbors(PartitionedHypergraph& phg, const Move& move) {
    for (HyperedgeID e : _edgesWithGainChanges) {
      /*! TODO: need to ignore large edges when doing this sequential?
       *  \todo need to ignore large edges when doing this sequential?
       */
      if (phg.edgeSize(e) < _context.partition.ignore_hyperedge_size_threshold) {
        for (HypernodeID v : phg.pins(e)) {
          _gain_cache.updateGain(phg, v, move);
        }
      }
    }
    _edgesWithGainChanges.clear();
  }

  void JudiciousRefiner::revertToBestLocalPrefix(PartitionedHypergraph& phg, size_t bestGainIndex) {
    while (_moves.size() > bestGainIndex) {
      Move& m = _moves.back();
      phg.changeNodePartWithGainCacheUpdate(m.node, m.to, m.from);
      m.invalidate();
      _moves.pop_back();
    }
  }
}

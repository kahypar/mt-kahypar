/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2022 Nikolai Maas <nikolai.maas@student.kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "snodes_sync_coarsening.h"

#include "mt-kahypar/datastructures/separated_nodes.h"
#include "mt-kahypar/partition/coarsening/separated_nodes/snodes_coarsening_pass.h"

namespace mt_kahypar {
namespace star_partitioning {
using ds::SepNodesStack;
using ds::SeparatedNodes;

void coarsenSynchronized(SepNodesStack& stack, const Hypergraph& original_hg, const vec<Level>& levels,
                         const Context& context, const HypernodeID& start_num_nodes, const HypernodeID& target_num_nodes) {
  size_t j = 0;
  HypernodeID current_num_nodes = start_num_nodes;
  for (size_t i = 0; i < levels.size(); ++i) {
    const double reduction_factor = std::pow(static_cast<double>(target_num_nodes) / static_cast<double>(current_num_nodes),
                                             1 / static_cast<double>(levels.size() - i));
    HypernodeID current_step_start_nodes = stack.coarsest().numVisibleNodes();
    bool reveal_batch = true;
    const Hypergraph* hg;
    do {
      if (j < levels.size()) {
        hg = &levels[j].contractedHypergraph();
        const Hypergraph& old_hg = (j == 0) ? original_hg : levels[j - 1].contractedHypergraph();
        current_step_start_nodes += old_hg.replaySeparated(levels[j].communities(), stack.coarsest());
        stack.coarsest().contract(levels[j].communities(), hg->initialNumNodes());
        if (i == j) {
          stack.coarsest().revealNextBatch();
          reveal_batch = false;
        }
        ++j;
      }

      stack.coarsest().initializeOutwardEdges();
      const HypernodeID tmp_target = std::max(0.65 * stack.coarsest().numVisibleNodes(),
                                              reduction_factor * current_step_start_nodes);
      SNodesCoarseningStage stage = (j == levels.size()) ? SNodesCoarseningStage::D1_TWINS : SNodesCoarseningStage::ON_LARGE_GRAPH;
      SNodesCoarseningPass c_pass(stack.coarsest(), *hg, context, tmp_target, stage);
      vec<HypernodeID> communities;
      current_num_nodes -= c_pass.run(communities, false);
      stack.coarsen(std::move(communities));
    } while (static_cast<double>(stack.coarsest().numVisibleNodes()) / current_step_start_nodes > reduction_factor
             && stack.coarsest().numVisibleNodes() > target_num_nodes);

    if (reveal_batch) {
      stack.coarsest().revealNextBatch();
    }
  }
}

} // namepace star_partitioning
} // namespace mt_kahypar

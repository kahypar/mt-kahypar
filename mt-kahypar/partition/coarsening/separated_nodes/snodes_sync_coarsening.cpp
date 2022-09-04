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

void replayToSynchronizeLevels(SepNodesStack& stack, const Hypergraph& original_hg, const vec<Level>& levels) {
  ASSERT(stack.numLevels() == levels.size() + 1);
  SepNodesStack new_stack(stack.finest().createCopyFromSavepoint());
  const Hypergraph* hg = &original_hg;
  for (size_t i = 0; i < levels.size(); ++i) {
    hg->replaySeparated(levels[i].communities(), new_stack.coarsest());
    new_stack.coarsest().revealNextBatch();
    hg = &levels[i].contractedHypergraph();
    new_stack.coarsen(std::move(stack.move_mapping(i, false)));
    new_stack.coarsest().contract(levels[i].communities(), hg->initialNumNodes());
  }
  std::swap(stack, new_stack);
}

void coarsenSynchronized(SepNodesStack& stack, const Hypergraph& original_hg, const vec<Level>& levels,
                         const Context& context, const HypernodeID& start_num_nodes, const HypernodeID& target_num_nodes,
                         Array<PartitionID>* part_ids) {
  ASSERT(part_ids == nullptr || part_ids->size() == start_num_nodes);

  stack.onliest().setSavepoint();
  size_t j = 0;
  HypernodeID current_num_nodes = start_num_nodes;
  SNodesCoarseningStage stage = SNodesCoarseningStage::ON_LARGE_GRAPH;
  for (size_t i = 0; i < levels.size(); ++i) {
    const double reduction_factor = std::pow(static_cast<double>(target_num_nodes) / static_cast<double>(current_num_nodes),
                                             1 / static_cast<double>(levels.size() - i));
    HypernodeID current_step_start_nodes = stack.coarsest().numVisibleNodes();
    bool reveal_batch = true;
    bool any_contracted = true;
    const Hypergraph* hg;
    do {
      if (j < levels.size() && (!context.coarsening.sep_nodes_coarsening_levelwise || i == j)) {
        hg = &levels[j].contractedHypergraph();
        const Hypergraph& old_hg = (j == 0) ? original_hg : levels[j - 1].contractedHypergraph();
        current_step_start_nodes += old_hg.replaySeparated(levels[j].communities(), stack.coarsest());
        stack.coarsest().contract(levels[j].communities(), hg->initialNumNodes());
        ++j;
      }
      if (reveal_batch) {
        stack.coarsest().revealNextBatch();
        reveal_batch = false;
      }
      ASSERT(!context.coarsening.sep_nodes_coarsening_levelwise || j == i + 1);

      double tmp_factor = reduction_factor * current_step_start_nodes / stack.coarsest().numVisibleNodes();
      if (tmp_factor < 0.65) {
        tmp_factor = std::max(0.65, std::sqrt(tmp_factor));
      }
      const HypernodeID tmp_target = tmp_factor * stack.coarsest().numVisibleNodes();
      if ((j == levels.size() && stage == SNodesCoarseningStage::ON_LARGE_GRAPH)
          || (context.coarsening.sep_nodes_coarsening_levelwise && j == stack.numLevels())) {
        stage = SNodesCoarseningStage::D1_TWINS;
      }
      stack.coarsest().initializeOutwardEdges();
      SNodesCoarseningPass c_pass(stack.coarsest(), *hg, context, tmp_target, stage);
      if (part_ids != nullptr) {
        c_pass.setPartIDs(*part_ids);
      }
      vec<HypernodeID> communities;
      const HypernodeID n_contracted = c_pass.run(communities);
      any_contracted = n_contracted > 0;
      current_num_nodes -= n_contracted;
      stage = c_pass.stage();
      stage = previous(stage);
      stack.coarsen(std::move(communities));

      if (part_ids != nullptr) {
        const size_t diff = part_ids->size() - current_num_nodes;
        Array<PartitionID> new_part_ids(current_num_nodes, kInvalidPartition);
        const vec<HypernodeID>& communities = stack.mapping(0);
        tbb::parallel_for(0UL, part_ids->size(), [&] (const size_t& pos) {
          if (pos < communities.size()) {
            const HypernodeID mapped_node = communities[pos];
            ASSERT(new_part_ids[mapped_node] == kInvalidPartition || new_part_ids[mapped_node] == (*part_ids)[pos]);
            new_part_ids[mapped_node] = (*part_ids)[pos];
          } else {
            new_part_ids[pos - diff] = (*part_ids)[pos];
          }
        });
        std::swap(*part_ids, new_part_ids);
      }
    } while (static_cast<double>(stack.coarsest().numVisibleNodes()) / current_step_start_nodes > reduction_factor
             && stack.coarsest().numVisibleNodes() > target_num_nodes && any_contracted);

    stack.contractToNLevels(i + 2);
  }

  replayToSynchronizeLevels(stack, original_hg, levels);

  ASSERT(stack.numLevels() == levels.size() + 1);
  ASSERT([&] {
    for (size_t i = 0; i < stack.numLevels(); ++i) {
      const Hypergraph& hg = (i == 0) ? original_hg : levels[i - 1].contractedHypergraph();
      if (stack.atLevel(i, false).numGraphNodes() !=  hg.initialNumNodes()) {
        return false;
      }
    }
    return true;
  }(), "Constructed separated nodes stack does not match graph levels!");
}

} // namepace star_partitioning
} // namespace mt_kahypar

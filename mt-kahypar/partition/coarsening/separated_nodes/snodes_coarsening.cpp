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

#include "snodes_coarsening.h"

#include "mt-kahypar/datastructures/separated_nodes.h"
#include "mt-kahypar/partition/coarsening/separated_nodes/snodes_coarsening_pass.h"

namespace mt_kahypar {
namespace star_partitioning {
using ds::SepNodesStack;
using ds::SeparatedNodes;

void coarsen(SepNodesStack& stack, const Hypergraph& coarsened_hg,
             const Context& context, const HypernodeID& target_num_nodes) {
  stack.coarsest().revealAll();
  HypernodeID previous_num_nodes = stack.coarsest().numNodes();
  SNodesCoarseningStage stage = SNodesCoarseningStage::D1_TWINS;
  while (stack.coarsest().numNodes() > target_num_nodes) {
    stack.coarsest().initializeOutwardEdges();
    const HypernodeID tmp_target = std::max(static_cast<HypernodeID>(0.65 * stack.coarsest().numNodes()), target_num_nodes);
    SNodesCoarseningPass c_pass(coarsened_hg, context, tmp_target, stage);
    vec<HypernodeID> communities;
    c_pass.run(communities);
    stack.coarsen(std::move(communities));
    stage = c_pass.stage();
    stage = previous(stage);

    if (stack.coarsest().numNodes() < static_cast<HypernodeID>(0.8 * previous_num_nodes)) {
      previous_num_nodes = stack.coarsest().numNodes();
    } else {
      break;
    }
  }
}

} // namepace star_partitioning
} // namespace mt_kahypar

/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "gmock/gmock.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/refinement/ilp/ilp_model.h"

using ::testing::Test;

namespace mt_kahypar {

class AILPModel : public Test {
public:

  AILPModel() :
    hg(HypergraphFactory::construct(TBBNumaArena::GLOBAL_TASK_GROUP,
        7, 4, {{0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6}})),
    phg(),
    context(),
    env() {
    phg = PartitionedHypergraph(2, TBBNumaArena::GLOBAL_TASK_GROUP, hg);
    phg.setOnlyNodePart(0, 0);
    phg.setOnlyNodePart(1, 0);
    phg.setOnlyNodePart(2, 0);
    phg.setOnlyNodePart(3, 0);
    phg.setOnlyNodePart(4, 1);
    phg.setOnlyNodePart(5, 1);
    phg.setOnlyNodePart(6, 1);
    phg.initializePartition(TBBNumaArena::GLOBAL_TASK_GROUP);
    context.partition.max_part_weights.assign(2, 4);
  }

  Hypergraph hg;
  PartitionedHypergraph phg;
  Context context;
  GRBEnv env;
};

TEST_F(AILPModel, TEST) {
  ILPHypergraph ilp_hg(phg, { 1, 3, 4 });
  ILPModel model(ilp_hg, context, env);
  model.construct();
  model.solve();
}


} // namespace mt_kahypar
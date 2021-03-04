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
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/ilp/ilp_model.h"

using ::testing::Test;

namespace mt_kahypar {

void assignVerticesToBlocks(PartitionedHypergraph& phg,
                            const vec<PartitionID>& partition) {
  ASSERT(phg.initialNumNodes() == ID(partition.size()));
  for ( const HypernodeID& hn : phg.nodes() ) {
    phg.setOnlyNodePart(hn, partition[hn]);
  }
  phg.initializePartition(TBBNumaArena::GLOBAL_TASK_GROUP);
}

void applyILPSolution(PartitionedHypergraph& phg,
                      ILPModel& model,
                      const vec<HypernodeID>& nodes) {
  for ( const HypernodeID& hn : nodes ) {
    const PartitionID from = phg.partID(hn);
    const PartitionID to = model.partID(hn);
    if ( from != to ) {
      phg.changeNodePart(hn, from, to);
    }
  }
}

TEST(AILPModel, OptimizesConnectivityMetric1) {
  // Setup Hypergraph
  Hypergraph hg = HypergraphFactory::construct(
    TBBNumaArena::GLOBAL_TASK_GROUP, 7, 4,
    {{0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6}});
  const PartitionID k = 2;
  PartitionedHypergraph phg(k, TBBNumaArena::GLOBAL_TASK_GROUP, hg);
  assignVerticesToBlocks(phg, { 0, 0, 0, 0, 1, 1, 1 });

  // Setup Context
  Context context;
  context.partition.max_part_weights.assign(2, 4);

  // Setup ILP Problem
  const vec<HypernodeID> nodes = {1, 3, 4};
  ILPHypergraph ilp_hg(phg, nodes);
  GRBEnv env;
  ILPModel model(ilp_hg, context, env);
  model.construct();

  // Solve ILP Problem
  const HyperedgeWeight model_objective_before = model.getObjective();
  const HyperedgeWeight connectivity_before = metrics::km1(phg);
  model.solve();
  const HyperedgeWeight model_objective_after = model.getObjective();
  const HyperedgeWeight delta = model_objective_before - model_objective_after;

  // Apply ILP Solution
  applyILPSolution(phg, model, nodes);

  // Verify Correctness of Solution
  const HyperedgeWeight connectivity_after = metrics::km1(phg);
  ASSERT_EQ(delta, 1);
  ASSERT_EQ(connectivity_before - delta, connectivity_after);
  for ( PartitionID i = 0; i < k; ++i ) {
    ASSERT_LE(phg.partWeight(i), context.partition.max_part_weights[i]);
  }
}


} // namespace mt_kahypar
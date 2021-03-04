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

void solve(PartitionedHypergraph& phg,
           const Context& context,
           const vec<HypernodeID>& nodes,
           const HyperedgeWeight& expected_delta,
           bool verify_delta = true,
           bool is_parallel = false) {
  // Setup ILP Problem
  ILPHypergraph ilp_hg(phg);
  ilp_hg.initialize(nodes);
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
  if ( verify_delta ) {
    ASSERT_EQ(expected_delta, delta);
  }
  if ( !is_parallel ) {
    ASSERT_EQ(connectivity_before - delta, connectivity_after);
  }
  for ( PartitionID i = 0; i < phg.k(); ++i ) {
    ASSERT_LE(phg.partWeight(i), context.partition.max_part_weights[i]);
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

  // Solve ILP Problem
  const vec<HypernodeID> nodes = {1, 3, 4};
  solve(phg, context, nodes, 1);
}

TEST(AILPModel, OptimizesConnectivityMetric2) {
  // Setup Hypergraph
  Hypergraph hg = HypergraphFactory::construct(
    TBBNumaArena::GLOBAL_TASK_GROUP, 9, 5,
    {{0, 1, 3, 6}, {1, 2, 5}, {5, 8}, {4, 7, 8}, {6, 7}});
  const PartitionID k = 3;
  PartitionedHypergraph phg(k, TBBNumaArena::GLOBAL_TASK_GROUP, hg);
  assignVerticesToBlocks(phg, { 0, 0, 0, 1, 1, 2, 1, 2, 2 });

  // Setup Context
  Context context;
  context.partition.max_part_weights.assign(3, 4);

  // Solve ILP Problem
  const vec<HypernodeID> nodes = {1, 2, 4, 5};
  solve(phg, context, nodes, 1);
}

TEST(AILPModel, OptimizesConnectivityMetric3) {
  // Setup Hypergraph
  Hypergraph hg = HypergraphFactory::construct(
    TBBNumaArena::GLOBAL_TASK_GROUP, 9, 4,
    {{0, 1, 3, 6}, {1, 2, 5}, {4, 7, 8}, {6, 7}});
  const PartitionID k = 3;
  PartitionedHypergraph phg(k, TBBNumaArena::GLOBAL_TASK_GROUP, hg);
  assignVerticesToBlocks(phg, { 0, 0, 0, 1, 1, 2, 1, 2, 2 });

  // Setup Context
  Context context;
  context.partition.max_part_weights.assign(3, 4);

  // Solve ILP Problem
  const vec<HypernodeID> nodes = {1, 2, 4, 5};
  solve(phg, context, nodes, 2);
}


TEST(AILPModel, OptimizesConnectivityMetric4) {
  // Setup Hypergraph
  Hypergraph hg = HypergraphFactory::construct(
    TBBNumaArena::GLOBAL_TASK_GROUP, 9, 5,
    {{0, 3, 6}, {1, 2, 5}, {5, 8}, {4, 7, 8}, {6, 7}});
  const PartitionID k = 3;
  PartitionedHypergraph phg(k, TBBNumaArena::GLOBAL_TASK_GROUP, hg);
  assignVerticesToBlocks(phg, { 0, 0, 0, 1, 1, 2, 1, 2, 2 });

  // Setup Context
  Context context;
  context.partition.max_part_weights.assign(3, 6);

  // Solve ILP Problem
  const vec<HypernodeID> nodes = {1, 2, 4, 5};
  solve(phg, context, nodes, 2);
}

TEST(AILPModel, OptimizesConnectivityMetric5) {
  // Setup Hypergraph
  Hypergraph hg = HypergraphFactory::construct(
    TBBNumaArena::GLOBAL_TASK_GROUP, 9, 5,
    {{0, 3, 6}, {1, 2, 5}, {5, 8}, {4, 7, 8}, {6, 7}});
  const PartitionID k = 3;
  PartitionedHypergraph phg(k, TBBNumaArena::GLOBAL_TASK_GROUP, hg);
  assignVerticesToBlocks(phg, { 0, 0, 0, 1, 1, 2, 1, 2, 2 });

  // Setup Context
  Context context;
  context.partition.max_part_weights.assign(3, 6);

  // Solve ILP Problem
  const vec<HypernodeID> nodes = {1, 2, 5};
  solve(phg, context, nodes, 1);
}


TEST(AILPModel, SolvesTwoILPsInParallel) {
  // Setup Hypergraph
  Hypergraph hg = HypergraphFactory::construct(
    TBBNumaArena::GLOBAL_TASK_GROUP, 9, 5,
    {{0, 1, 3, 6}, {2, 5}, {5, 8}, {4, 7, 8}, {6, 7}});
  const PartitionID k = 3;
  PartitionedHypergraph phg(k, TBBNumaArena::GLOBAL_TASK_GROUP, hg);
  assignVerticesToBlocks(phg, { 0, 0, 0, 1, 1, 2, 1, 2, 2 });

  // Setup Context
  Context context;
  context.partition.max_part_weights.assign(3, 5);

  // Solve ILP Problem
  tbb::parallel_invoke([&] {
    const vec<HypernodeID> nodes = {0, 1};
    solve(phg, context, nodes, 1, true, true);
  }, [&] {
    const vec<HypernodeID> nodes = {2, 4};
    solve(phg, context, nodes, 2, true, true);
  });

  ASSERT_EQ(1, metrics::km1(phg));
}

Hypergraph generateRandomHypergraph(const HypernodeID num_hypernodes,
                                    const HyperedgeID num_hyperedges,
                                    const HypernodeID max_edge_size) {
  parallel::scalable_vector<parallel::scalable_vector<HypernodeID>> hyperedges;
  utils::Randomize& rand = utils::Randomize::instance();
  for ( size_t i = 0; i < num_hyperedges; ++i ) {
    parallel::scalable_vector<HypernodeID> net;
    const size_t edge_size = rand.getRandomInt(2, max_edge_size, sched_getcpu());
    for ( size_t i = 0; i < edge_size; ++i ) {
      const HypernodeID pin = rand.getRandomInt(0, num_hypernodes - 1, sched_getcpu());
      if ( std::find(net.begin(), net.end(), pin) == net.end() ) {
        net.push_back(pin);
      }
    }
    hyperedges.emplace_back(std::move(net));
  }
  return HypergraphFactory::construct(
    TBBNumaArena::GLOBAL_TASK_GROUP, num_hypernodes, num_hyperedges, hyperedges);
}

TEST(AILPModel, StressTest) {
  // Setup Hypergraph
  const HypernodeID num_hypernodes = 50;
  const HypernodeID num_hyperedges = 50;
  const HypernodeID max_edge_size = 10;
  Hypergraph hg = generateRandomHypergraph(num_hypernodes, num_hyperedges, max_edge_size);
  const PartitionID k = 4;
  PartitionedHypergraph phg(k, TBBNumaArena::GLOBAL_TASK_GROUP, hg);
  for ( const HypernodeID& hn : phg.nodes() ) {
    phg.setOnlyNodePart(hn, hn % k);
  }
  phg.initializePartition(TBBNumaArena::GLOBAL_TASK_GROUP);

  // Setup Context
  Context context;
  context.partition.max_part_weights.assign(k, num_hypernodes / k + 5);

  // Solve ILP Problem
  const HypernodeID num_nodes_in_ilp = 20;
  vec<HypernodeID> nodes;
  for ( HypernodeID hn = 0; hn < num_nodes_in_ilp; ++hn ) {
    nodes.push_back(hn);
  }
  solve(phg, context, nodes, 0, false, false);
}

} // namespace mt_kahypar
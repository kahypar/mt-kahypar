/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/datastructures/graph.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/preprocessing/community_detection/parallel_louvain.h"
#include "mt-kahypar/io/hypergraph_io.h"

using ::testing::Test;

namespace mt_kahypar::community_detection {

class ALouvain : public ds::HypergraphFixture {

  using Base = ds::HypergraphFixture;

 public:

  ALouvain() :
    Base(),
    graph(nullptr),
    context(),
    karate_club_hg(),
    karate_club_graph(nullptr) {
    context.partition.graph_filename = "../tests/instances/karate_club.graph.hgr";
    context.preprocessing.community_detection.edge_weight_function = LouvainEdgeWeight::uniform;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;

    graph = std::make_unique<Graph>(hypergraph, LouvainEdgeWeight::uniform);
    karate_club_hg = io::readHypergraphFile(
      context.partition.graph_filename, TBBNumaArena::GLOBAL_TASK_GROUP);
    karate_club_graph = std::make_unique<Graph>(karate_club_hg, LouvainEdgeWeight::uniform);
  }

  static void SetUpTestSuite() {
    TBBNumaArena::instance(1);
  }

  static void TearDownTestSuite() {
    TBBNumaArena::instance().terminate();
  }

  using Base::hypergraph;
  std::unique_ptr<Graph> graph;
  Context context;
  ds::StaticHypergraph karate_club_hg;
  std::unique_ptr<Graph> karate_club_graph;
};

ds::Clustering clustering(const std::vector<PartitionID>& communities) {
  ds::Clustering c(communities.size());
  for ( size_t i = 0; i < communities.size(); ++i ) {
    c[i] = communities[i];
  }
  return c;
}

TEST_F(ALouvain, ComputesMaxGainMove1) {
  ParallelLocalMovingModularity plm(context, graph->numNodes());
  ds::Clustering communities = clustering( { 0, 1, 0, 2, 3, 4, 5, 1, 2, 3, 4 } );
  plm.initializeClusterVolumes(*graph, communities);
  PartitionID to = plm.computeMaxGainCluster(
    *graph, communities, 7, plm._local_large_incident_cluster_weight.local());
  ASSERT_EQ(0, to);
}

TEST_F(ALouvain, ComputesMaxGainMove2) {
  ParallelLocalMovingModularity plm(context, graph->numNodes());
  ds::Clustering communities = clustering( { 0, 1, 0, 3, 3, 4, 5, 1, 2, 3, 4 } );
  plm.initializeClusterVolumes(*graph, communities);
  PartitionID to = plm.computeMaxGainCluster(
    *graph, communities, 8, plm._local_large_incident_cluster_weight.local());
  ASSERT_EQ(3, to);
}

TEST_F(ALouvain, ComputesMaxGainMove3) {
  ParallelLocalMovingModularity plm(context, graph->numNodes());
  ds::Clustering communities = clustering( { 0, 1, 0, 2, 3, 4, 5, 1, 2, 3, 4 } );
  plm.initializeClusterVolumes(*graph, communities);
  PartitionID to = plm.computeMaxGainCluster(
    *graph, communities, 8, plm._local_large_incident_cluster_weight.local());
  ASSERT_EQ(2, to);
}

TEST_F(ALouvain, ComputesMaxGainMove4) {
  ParallelLocalMovingModularity plm(context, graph->numNodes());
  ds::Clustering communities = clustering( { 0, 1, 0, 2, 3, 4, 5, 1, 2, 3, 4 } );
  plm.initializeClusterVolumes(*graph, communities);
  PartitionID to = plm.computeMaxGainCluster(
    *graph, communities, 9, plm._local_large_incident_cluster_weight.local());
  ASSERT_EQ(3, to);
}

TEST_F(ALouvain, ComputesMaxGainMove5) {
  ParallelLocalMovingModularity plm(context, graph->numNodes());
  ds::Clustering communities = clustering( { 0, 1, 0, 2, 2, 4, 5, 1, 2, 3, 4 } );
  plm.initializeClusterVolumes(*graph, communities);
  PartitionID to = plm.computeMaxGainCluster(
    *graph, communities, 9, plm._local_large_incident_cluster_weight.local());
  ASSERT_EQ(2, to);
}

TEST_F(ALouvain, ComputesMaxGainMove6) {
  ParallelLocalMovingModularity plm(context, graph->numNodes());
  ds::Clustering communities = clustering( { 0, 1, 0, 2, 2, 4, 5, 1, 2, 3, 4 } );
  plm.initializeClusterVolumes(*graph, communities);
  PartitionID to = plm.computeMaxGainCluster(
    *graph, communities, 10, plm._local_large_incident_cluster_weight.local());
  ASSERT_EQ(4, to);
}

TEST_F(ALouvain, ComputesMaxGainMove7) {
  ParallelLocalMovingModularity plm(context, graph->numNodes());
  ds::Clustering communities = clustering( { 0, 1, 0, 2, 2, 4, 0, 1, 2, 3, 4 } );
  plm.initializeClusterVolumes(*graph, communities);
  PartitionID to = plm.computeMaxGainCluster(
    *graph, communities, 10, plm._local_large_incident_cluster_weight.local());
  ASSERT_EQ(0, to);
}

TEST_F(ALouvain, ComputesMaxGainMove8) {
  ParallelLocalMovingModularity plm(context, graph->numNodes());
  ds::Clustering communities = clustering( { 0, 1, 0, 2, 2, 4, 0, 1, 1, 3, 4 } );
  plm.initializeClusterVolumes(*graph, communities);
  PartitionID to = plm.computeMaxGainCluster(
    *graph, communities, 0, plm._local_large_incident_cluster_weight.local());
  ASSERT_EQ(1, to);
}

TEST_F(ALouvain, ComputesMaxGainMove9) {
  ParallelLocalMovingModularity plm(context, graph->numNodes());
  ds::Clustering communities = clustering( { 0, 1, 0, 2, 2, 4, 0, 1, 3, 3, 4 } );
  plm.initializeClusterVolumes(*graph, communities);
  PartitionID to = plm.computeMaxGainCluster(
    *graph, communities, 4, plm._local_large_incident_cluster_weight.local());
  ASSERT_EQ(3, to);
}

TEST_F(ALouvain, ComputesMaxGainMove10) {
  ParallelLocalMovingModularity plm(context, graph->numNodes());
  ds::Clustering communities = clustering( { 0, 1, 0, 2, 2, 0, 4, 1, 3, 3, 4 } );
  plm.initializeClusterVolumes(*graph, communities);
  PartitionID to = plm.computeMaxGainCluster(
    *graph, communities, 6, plm._local_large_incident_cluster_weight.local());
  ASSERT_EQ(4, to);
}

TEST_F(ALouvain, KarateClubTest) {
  ds::Clustering communities = run_parallel_louvain(*karate_club_graph, context);
  std::vector<PartitionID> expected_comm = { 1, 1, 1, 1, 0, 0, 0, 1, 3, 1, 0, 1, 1, 1, 3, 3, 0, 1,
                                             3, 1, 3, 1, 3, 2, 2, 2, 3, 2, 2, 3, 3, 2, 3, 3 };
  ds::Clustering exp(expected_comm.size());
  for (size_t i = 0; i < expected_comm.size(); ++i)
    exp[i] = expected_comm[i];

  LOG << V(communities.size()) << V(exp.size()) << V(expected_comm.size()) << V(karate_club_graph->numNodes());
  LOG << V(metrics::modularity(*karate_club_graph, communities));
  LOG << V(metrics::modularity(*karate_club_graph, exp));

  for (const NodeID u : karate_club_graph->nodes()) {
    LOG << V(u) << V(expected_comm[u]) << V(communities[u]);
  }
  for ( const NodeID u : karate_club_graph->nodes() ) {
    ASSERT_EQ(expected_comm[u], communities[u]);
  }
}

}  // namespace mt_kahypar

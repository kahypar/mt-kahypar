/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/io/hypergraph_io.h"

#include "mt-kahypar/partition/initial_partitioning/flat/bfs_initial_partitioner.h"

#include "mt-kahypar/partition/coarsening/deterministic_multilevel_coarsener.h"
#include "mt-kahypar/partition/refinement/deterministic/deterministic_label_propagation.h"
#include "mt-kahypar/partition/preprocessing/community_detection/parallel_louvain.h"

using ::testing::Test;

namespace mt_kahypar {

  class DeterminismTest : public Test {

  public:
    DeterminismTest() :
            hypergraph(),
            partitioned_hypergraph(),
            context(),
            metrics() {
      context.partition.graph_filename = "../tests/instances/powersim.mtx.hgr";
      context.partition.mode = kahypar::Mode::direct_kway;
      context.partition.epsilon = 0.25;
      context.partition.verbose_output = false;
      context.partition.k = 2;

      // Shared Memory
      context.shared_memory.num_threads = std::thread::hardware_concurrency();

      // Initial Partitioning
      context.initial_partitioning.mode = InitialPartitioningMode::recursive_bisection;
      context.initial_partitioning.runs = 1;

      context.partition.deterministic = true;

      // preprocessing
      context.preprocessing.community_detection.num_sub_rounds_deterministic = 16;
      context.preprocessing.community_detection.max_pass_iterations = 5;
      context.preprocessing.community_detection.min_vertex_move_fraction = 0.01;

      // coarsening
      context.coarsening.num_sub_rounds_deterministic = 3;
      context.coarsening.contraction_limit = 320;
      context.coarsening.max_allowed_node_weight = std::numeric_limits<HypernodeWeight>::max();
      context.coarsening.minimum_shrink_factor = 1.0;
      context.coarsening.maximum_shrink_factor = 4.0;

      context.partition.objective = kahypar::Objective::km1;

      // Read hypergraph
      hypergraph = io::readHypergraphFile(context.partition.graph_filename, TBBNumaArena::GLOBAL_TASK_GROUP, true);
      partitioned_hypergraph = PartitionedHypergraph(
              context.partition.k, TBBNumaArena::GLOBAL_TASK_GROUP, hypergraph);
      context.setupPartWeights(hypergraph.totalWeight());
    }

    void initialPartition() {
      Context ip_context(context);
      ip_context.refinement.label_propagation.algorithm = LabelPropagationAlgorithm::do_nothing;
      InitialPartitioningDataContainer ip_data(partitioned_hypergraph, ip_context, TBBNumaArena::GLOBAL_TASK_GROUP);
      BFSInitialPartitioner& initial_partitioner = *new(tbb::task::allocate_root())
      BFSInitialPartitioner(InitialPartitioningAlgorithm::bfs, ip_data, ip_context, 420, 0);
      tbb::task::spawn_root_and_wait(initial_partitioner);
      ip_data.apply();
      metrics.km1 = metrics::km1(partitioned_hypergraph);
      metrics.cut = metrics::hyperedgeCut(partitioned_hypergraph);
      metrics.imbalance = metrics::imbalance(partitioned_hypergraph, context);
    }

    Hypergraph hypergraph;
    PartitionedHypergraph partitioned_hypergraph;
    Context context;
    kahypar::Metrics metrics;
    static constexpr size_t num_repetitions = 5;
  };

  TEST_F(DeterminismTest, Preprocessing) {
    context.preprocessing.community_detection.low_memory_contraction = true;

    LouvainEdgeWeight edge_weight_type;
    if (static_cast<double>(hypergraph.initialNumEdges()) /
        static_cast<double>(hypergraph.initialNumNodes()) < 0.75) {
      edge_weight_type = LouvainEdgeWeight::degree;
    } else {
      edge_weight_type = LouvainEdgeWeight::uniform;
    }

    Graph graph(hypergraph, edge_weight_type);
    ds::Clustering first;
    for (size_t i = 0; i < num_repetitions; ++i) {
      ds::Clustering communities = community_detection::run_parallel_louvain(graph, context);
      if (i == 0) {
        first = std::move(communities);
      } else {
        ASSERT_EQ(first, communities);
      }
    }
  }

  TEST_F(DeterminismTest, Coarsening) {
    Hypergraph first;
    for (size_t i = 0; i < num_repetitions; ++i) {
      DeterministicMultilevelCoarsener coarsener(hypergraph, context, 0, true);
      coarsener.coarsen();
      if (i == 0) {
        first = coarsener.coarsestHypergraph().copy();
      } else {
        const Hypergraph& other = coarsener.coarsestHypergraph();
        ASSERT_EQ(other.initialNumNodes(), first.initialNumNodes());
        ASSERT_EQ(other.initialNumEdges(), first.initialNumEdges());
        ASSERT_EQ(other.initialNumPins(), first.initialNumPins());
        vec<HyperedgeID> inets_first, inets_other;
        for (HypernodeID u : first.nodes()) {
          for (HyperedgeID e : first.incidentEdges(u)) inets_first.push_back(e);
          for (HyperedgeID e : other.incidentEdges(u)) inets_other.push_back(e);
          ASSERT_EQ(inets_first, inets_other);
          inets_first.clear(); inets_other.clear();
        }

        vec<HypernodeID> pins_first, pins_other;
        for (HyperedgeID e : first.edges()) {
          for (HypernodeID v : first.pins(e)) pins_first.push_back(v);
          for (HypernodeID v : other.pins(e)) pins_other.push_back(v);
          ASSERT_EQ(pins_first, pins_other);
          pins_first.clear(); pins_other.clear();
        }
      }
    }
  }

  TEST_F(DeterminismTest, Refinement) {
    initialPartition();
    vec<PartitionID> initial_partition(hypergraph.initialNumNodes());
    for (HypernodeID u : hypergraph.nodes()) {
      initial_partition[u] = partitioned_hypergraph.partID(u);
    }

    vec<PartitionID> first(hypergraph.initialNumNodes());
    for (size_t i = 0; i < num_repetitions; ++i) {
      partitioned_hypergraph.resetPartition();
      for (HypernodeID u : hypergraph.nodes()) {
        partitioned_hypergraph.setNodePart(u, initial_partition[u]);
      }

      DeterministicLabelPropagationRefiner refiner(hypergraph, context, 0);
      refiner.initialize(partitioned_hypergraph);
      vec<HypernodeID> dummy_refinement_nodes;
      kahypar::Metrics my_metrics = metrics;
      refiner.refine(partitioned_hypergraph, dummy_refinement_nodes, my_metrics, 0.0);

      if (i == 0) {
        for (HypernodeID u : hypergraph.nodes()) {
          first[u] = partitioned_hypergraph.partID(u);
        }
      } else {
        for (HypernodeID u : hypergraph.nodes()) {
          ASSERT_EQ(first[u], partitioned_hypergraph.partID(u));
        }
      }

    }


  }

}  // namespace mt_kahypar

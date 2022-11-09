/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/


#include "gmock/gmock.h"

#include <thread>

#include "tbb/parallel_invoke.h"

#include "libmtkahypar.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/partition/context.h"

using ::testing::Test;

namespace mt_kahypar {

  static constexpr bool debug = false;

  TEST(MtKaHyPar, ReadHypergraphFile) {
    mt_kahypar_context_t* context = mt_kahypar_context_new();
    mt_kahypar_load_preset(context, SPEED);

    mt_kahypar_hypergraph_t* hypergraph =
      mt_kahypar_read_hypergraph_from_file("test_instances/ibm01.hgr", context, HMETIS);

    ASSERT_EQ(12752, mt_kahypar_num_hypernodes(hypergraph));
    ASSERT_EQ(14111, mt_kahypar_num_hyperedges(hypergraph));
    ASSERT_EQ(50566, mt_kahypar_num_pins(hypergraph));
    ASSERT_EQ(12752, mt_kahypar_hypergraph_weight(hypergraph));

    mt_kahypar_free_context(context);
    mt_kahypar_free_hypergraph(hypergraph);
  }

  TEST(MtKaHyPar, ReadGraphFile) {
    mt_kahypar_context_t* context = mt_kahypar_context_new();
    mt_kahypar_load_preset(context, SPEED);

    mt_kahypar_graph_t* graph =
      mt_kahypar_read_graph_from_file("test_instances/delaunay_n15.graph", context, METIS);

    ASSERT_EQ(32768, mt_kahypar_num_nodes(graph));
    ASSERT_EQ(98274, mt_kahypar_num_edges(graph));
    ASSERT_EQ(32768, mt_kahypar_graph_weight(graph));

    mt_kahypar_free_context(context);
    mt_kahypar_free_graph(graph);
  }

  TEST(MtKaHyPar, ReadGraphFileAsHypergraph) {
    mt_kahypar_context_t* context = mt_kahypar_context_new();
    mt_kahypar_load_preset(context, SPEED);

    mt_kahypar_hypergraph_t* hypergraph =
      mt_kahypar_read_hypergraph_from_file("test_instances/delaunay_n15.graph", context, METIS);

    ASSERT_EQ(32768, mt_kahypar_num_hypernodes(hypergraph));
    ASSERT_EQ(98274, mt_kahypar_num_hyperedges(hypergraph));
    ASSERT_EQ(196548, mt_kahypar_num_pins(hypergraph));
    ASSERT_EQ(32768, mt_kahypar_hypergraph_weight(hypergraph));

    mt_kahypar_free_context(context);
    mt_kahypar_free_hypergraph(hypergraph);
  }

  TEST(MtKaHyPar, ConstructUnweightedHypergraph) {
    const mt_kahypar_hypernode_id_t num_vertices = 7;
    const mt_kahypar_hyperedge_id_t num_hyperedges = 4;

    std::unique_ptr<size_t[]> hyperedge_indices = std::make_unique<size_t[]>(5);
    hyperedge_indices[0] = 0; hyperedge_indices[1] = 2; hyperedge_indices[2] = 6;
    hyperedge_indices[3] = 9; hyperedge_indices[4] = 12;

    std::unique_ptr<mt_kahypar_hyperedge_id_t[]> hyperedges = std::make_unique<mt_kahypar_hyperedge_id_t[]>(12);
    hyperedges[0] = 0;  hyperedges[1] = 2;                                        // Hyperedge 0
    hyperedges[2] = 0;  hyperedges[3] = 1; hyperedges[4] = 3;  hyperedges[5] = 4; // Hyperedge 1
    hyperedges[6] = 3;  hyperedges[7] = 4; hyperedges[8] = 6;                     // Hyperedge 2
    hyperedges[9] = 2; hyperedges[10] = 5; hyperedges[11] = 6;                    // Hyperedge 3

    mt_kahypar_hypergraph_t* hypergraph = mt_kahypar_create_hypergraph(
      num_vertices, num_hyperedges, hyperedge_indices.get(), hyperedges.get(), nullptr, nullptr);

    ASSERT_EQ(7, mt_kahypar_num_hypernodes(hypergraph));
    ASSERT_EQ(4, mt_kahypar_num_hyperedges(hypergraph));
    ASSERT_EQ(12, mt_kahypar_num_pins(hypergraph));
    ASSERT_EQ(7, mt_kahypar_hypergraph_weight(hypergraph));

    mt_kahypar_free_hypergraph(hypergraph);
  }

  TEST(MtKaHyPar, ConstructUnweightedGraph) {
    const mt_kahypar_hypernode_id_t num_vertices = 5;
    const mt_kahypar_hyperedge_id_t num_hyperedges = 6;

    std::unique_ptr<mt_kahypar_hypernode_id_t[]> edges =
      std::make_unique<mt_kahypar_hypernode_id_t[]>(12);
    edges[0] = 0;  edges[1] = 1;
    edges[2] = 0;  edges[3] = 2;
    edges[4] = 1;  edges[5] = 2;
    edges[6] = 1;  edges[7] = 3;
    edges[8] = 2;  edges[9] = 3;
    edges[10] = 3; edges[11] = 4;

    mt_kahypar_graph_t* graph = mt_kahypar_create_graph(
      num_vertices, num_hyperedges, edges.get(), nullptr, nullptr);

    ASSERT_EQ(5, mt_kahypar_num_nodes(graph));
    ASSERT_EQ(6, mt_kahypar_num_edges(graph));
    ASSERT_EQ(5, mt_kahypar_graph_weight(graph));

    mt_kahypar_free_graph(graph);
  }

  TEST(MtKaHyPar, ConstructHypergraphWithNodeWeights) {
    const mt_kahypar_hypernode_id_t num_vertices = 7;
    const mt_kahypar_hyperedge_id_t num_hyperedges = 4;

    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> vertex_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(7);
    vertex_weights[0] = 1; vertex_weights[1] = 2; vertex_weights[2] = 3; vertex_weights[3] = 4;
    vertex_weights[4] = 5; vertex_weights[5] = 6; vertex_weights[6] = 7;

    std::unique_ptr<size_t[]> hyperedge_indices = std::make_unique<size_t[]>(5);
    hyperedge_indices[0] = 0; hyperedge_indices[1] = 2; hyperedge_indices[2] = 6;
    hyperedge_indices[3] = 9; hyperedge_indices[4] = 12;

    std::unique_ptr<mt_kahypar_hyperedge_id_t[]> hyperedges = std::make_unique<mt_kahypar_hyperedge_id_t[]>(12);
    hyperedges[0] = 0;  hyperedges[1] = 2;                                        // Hyperedge 0
    hyperedges[2] = 0;  hyperedges[3] = 1; hyperedges[4] = 3;  hyperedges[5] = 4; // Hyperedge 1
    hyperedges[6] = 3;  hyperedges[7] = 4; hyperedges[8] = 6;                     // Hyperedge 2
    hyperedges[9] = 2; hyperedges[10] = 5; hyperedges[11] = 6;                    // Hyperedge 3

    mt_kahypar_hypergraph_t* hypergraph = mt_kahypar_create_hypergraph(
      num_vertices, num_hyperedges, hyperedge_indices.get(), hyperedges.get(), nullptr, vertex_weights.get());

    ASSERT_EQ(7, mt_kahypar_num_hypernodes(hypergraph));
    ASSERT_EQ(4, mt_kahypar_num_hyperedges(hypergraph));
    ASSERT_EQ(12, mt_kahypar_num_pins(hypergraph));
    ASSERT_EQ(28, mt_kahypar_hypergraph_weight(hypergraph));

    mt_kahypar_free_hypergraph(hypergraph);
  }

  TEST(MtKaHyPar, ConstructGraphWithNodeWeights) {
    const mt_kahypar_hypernode_id_t num_vertices = 5;
    const mt_kahypar_hyperedge_id_t num_hyperedges = 6;

    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> vertex_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(7);
    vertex_weights[0] = 1; vertex_weights[1] = 2; vertex_weights[2] = 3;
    vertex_weights[3] = 4; vertex_weights[4] = 5;

    std::unique_ptr<mt_kahypar_hypernode_id_t[]> edges =
      std::make_unique<mt_kahypar_hypernode_id_t[]>(12);
    edges[0] = 0;  edges[1] = 1;
    edges[2] = 0;  edges[3] = 2;
    edges[4] = 1;  edges[5] = 2;
    edges[6] = 1;  edges[7] = 3;
    edges[8] = 2;  edges[9] = 3;
    edges[10] = 3; edges[11] = 4;

    mt_kahypar_graph_t* graph = mt_kahypar_create_graph(
      num_vertices, num_hyperedges, edges.get(), nullptr, vertex_weights.get());

    ASSERT_EQ(5, mt_kahypar_num_nodes(graph));
    ASSERT_EQ(6, mt_kahypar_num_edges(graph));
    ASSERT_EQ(15, mt_kahypar_graph_weight(graph));

    mt_kahypar_free_graph(graph);
  }

  TEST(MtKaHyPar, CreatesPartitionedHypergraph) {
    const mt_kahypar_hypernode_id_t num_vertices = 7;
    const mt_kahypar_hyperedge_id_t num_hyperedges = 4;

    std::unique_ptr<size_t[]> hyperedge_indices = std::make_unique<size_t[]>(5);
    hyperedge_indices[0] = 0; hyperedge_indices[1] = 2; hyperedge_indices[2] = 6;
    hyperedge_indices[3] = 9; hyperedge_indices[4] = 12;

    std::unique_ptr<mt_kahypar_hyperedge_id_t[]> hyperedges = std::make_unique<mt_kahypar_hyperedge_id_t[]>(12);
    hyperedges[0] = 0;  hyperedges[1] = 2;                                        // Hyperedge 0
    hyperedges[2] = 0;  hyperedges[3] = 1; hyperedges[4] = 3;  hyperedges[5] = 4; // Hyperedge 1
    hyperedges[6] = 3;  hyperedges[7] = 4; hyperedges[8] = 6;                     // Hyperedge 2
    hyperedges[9] = 2; hyperedges[10] = 5; hyperedges[11] = 6;                    // Hyperedge 3

    mt_kahypar_hypergraph_t* hypergraph = mt_kahypar_create_hypergraph(
      num_vertices, num_hyperedges, hyperedge_indices.get(), hyperedges.get(), nullptr, nullptr);

    std::unique_ptr<mt_kahypar_partition_id_t[]> partition =
      std::make_unique<mt_kahypar_partition_id_t[]>(7);
    partition[0] = 0; partition[1] = 0; partition[2] = 0;
    partition[3] = 1; partition[4] = 1; partition[5] = 1; partition[6] = 1;

    mt_kahypar_partitioned_hypergraph_t* partitioned_hg =
      mt_kahypar_create_partitioned_hypergraph(hypergraph, 2, partition.get());

    std::unique_ptr<mt_kahypar_partition_id_t[]> actual_partition =
      std::make_unique<mt_kahypar_partition_id_t[]>(7);
    mt_kahypar_get_hypergraph_partition(partitioned_hg, actual_partition.get());

    ASSERT_EQ(2, mt_kahypar_km1(partitioned_hg));
    for ( mt_kahypar_hypernode_id_t hn = 0; hn < 7; ++hn ) {
      ASSERT_EQ(partition[hn], actual_partition[hn]);
    }

    mt_kahypar_free_hypergraph(hypergraph);
    mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
  }

  TEST(MtKaHyPar, CreatesPartitionedGraph) {
    const mt_kahypar_hypernode_id_t num_vertices = 5;
    const mt_kahypar_hyperedge_id_t num_hyperedges = 6;

    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> vertex_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(7);
    vertex_weights[0] = 1; vertex_weights[1] = 2; vertex_weights[2] = 3;
    vertex_weights[3] = 4; vertex_weights[4] = 5;

    std::unique_ptr<mt_kahypar_hypernode_id_t[]> edges =
      std::make_unique<mt_kahypar_hypernode_id_t[]>(12);
    edges[0] = 0;  edges[1] = 1;
    edges[2] = 0;  edges[3] = 2;
    edges[4] = 1;  edges[5] = 2;
    edges[6] = 1;  edges[7] = 3;
    edges[8] = 2;  edges[9] = 3;
    edges[10] = 3; edges[11] = 4;

    mt_kahypar_graph_t* graph = mt_kahypar_create_graph(
      num_vertices, num_hyperedges, edges.get(), nullptr, vertex_weights.get());

    std::unique_ptr<mt_kahypar_partition_id_t[]> partition =
      std::make_unique<mt_kahypar_partition_id_t[]>(7);
    partition[0] = 0; partition[1] = 0; partition[2] = 1;
    partition[3] = 1; partition[4] = 1;

    mt_kahypar_partitioned_graph_t* partitioned_graph =
      mt_kahypar_create_partitioned_graph(graph, 2, partition.get());

    std::unique_ptr<mt_kahypar_partition_id_t[]> actual_partition =
      std::make_unique<mt_kahypar_partition_id_t[]>(7);
    mt_kahypar_get_graph_partition(partitioned_graph, actual_partition.get());

    ASSERT_EQ(3, mt_kahypar_graph_cut(partitioned_graph));
    for ( mt_kahypar_hypernode_id_t hn = 0; hn < 7; ++hn ) {
      ASSERT_EQ(partition[hn], actual_partition[hn]);
    }

    mt_kahypar_free_graph(graph);
    mt_kahypar_free_partitioned_graph(partitioned_graph);
  }

  TEST(MtKaHyPar, WritesAndLoadsHypergraphPartitionFile) {
    const mt_kahypar_hypernode_id_t num_vertices = 7;
    const mt_kahypar_hyperedge_id_t num_hyperedges = 4;

    std::unique_ptr<size_t[]> hyperedge_indices = std::make_unique<size_t[]>(5);
    hyperedge_indices[0] = 0; hyperedge_indices[1] = 2; hyperedge_indices[2] = 6;
    hyperedge_indices[3] = 9; hyperedge_indices[4] = 12;

    std::unique_ptr<mt_kahypar_hyperedge_id_t[]> hyperedges = std::make_unique<mt_kahypar_hyperedge_id_t[]>(12);
    hyperedges[0] = 0;  hyperedges[1] = 2;                                        // Hyperedge 0
    hyperedges[2] = 0;  hyperedges[3] = 1; hyperedges[4] = 3;  hyperedges[5] = 4; // Hyperedge 1
    hyperedges[6] = 3;  hyperedges[7] = 4; hyperedges[8] = 6;                     // Hyperedge 2
    hyperedges[9] = 2; hyperedges[10] = 5; hyperedges[11] = 6;                    // Hyperedge 3

    mt_kahypar_hypergraph_t* hypergraph = mt_kahypar_create_hypergraph(
      num_vertices, num_hyperedges, hyperedge_indices.get(), hyperedges.get(), nullptr, nullptr);

    std::unique_ptr<mt_kahypar_partition_id_t[]> partition = std::make_unique<mt_kahypar_partition_id_t[]>(7);
    partition[0] = 0; partition[1] = 0; partition[2] = 0;
    partition[3] = 1; partition[4] = 1; partition[5] = 1; partition[6] = 1;

    mt_kahypar_partitioned_hypergraph_t* partitioned_hg =
      mt_kahypar_create_partitioned_hypergraph(hypergraph, 2, partition.get());

    mt_kahypar_write_hypergraph_partition_to_file(partitioned_hg, "tmp.partition");

    mt_kahypar_partitioned_hypergraph_t* partitioned_hg_2 =
      mt_kahypar_read_hypergraph_partition_from_file(hypergraph, 2, "tmp.partition");

    std::unique_ptr<mt_kahypar_partition_id_t[]> actual_partition =
      std::make_unique<mt_kahypar_partition_id_t[]>(7);
    mt_kahypar_get_hypergraph_partition(partitioned_hg_2, actual_partition.get());

    ASSERT_EQ(2, mt_kahypar_km1(partitioned_hg_2));
    for ( mt_kahypar_hypernode_id_t hn = 0; hn < 5; ++hn ) {
      ASSERT_EQ(partition[hn], actual_partition[hn]);
    }

    mt_kahypar_free_hypergraph(hypergraph);
    mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
    mt_kahypar_free_partitioned_hypergraph(partitioned_hg_2);
  }

  TEST(MtKaHyPar, WritesAndLoadsGraphPartitionFile) {
    const mt_kahypar_hypernode_id_t num_vertices = 5;
    const mt_kahypar_hyperedge_id_t num_hyperedges = 6;

    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> vertex_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(7);
    vertex_weights[0] = 1; vertex_weights[1] = 2; vertex_weights[2] = 3;
    vertex_weights[3] = 4; vertex_weights[4] = 5;

    std::unique_ptr<mt_kahypar_hypernode_id_t[]> edges =
      std::make_unique<mt_kahypar_hypernode_id_t[]>(12);
    edges[0] = 0;  edges[1] = 1;
    edges[2] = 0;  edges[3] = 2;
    edges[4] = 1;  edges[5] = 2;
    edges[6] = 1;  edges[7] = 3;
    edges[8] = 2;  edges[9] = 3;
    edges[10] = 3; edges[11] = 4;

    mt_kahypar_graph_t* graph = mt_kahypar_create_graph(
      num_vertices, num_hyperedges, edges.get(), nullptr, vertex_weights.get());

    std::unique_ptr<mt_kahypar_partition_id_t[]> partition =
      std::make_unique<mt_kahypar_partition_id_t[]>(7);
    partition[0] = 0; partition[1] = 0; partition[2] = 1;
    partition[3] = 1; partition[4] = 1;

    mt_kahypar_partitioned_graph_t* partitioned_graph =
      mt_kahypar_create_partitioned_graph(graph, 2, partition.get());

    mt_kahypar_write_graph_partition_to_file(partitioned_graph, "tmp.partition");

    mt_kahypar_partitioned_graph_t* partitioned_graph_2 =
      mt_kahypar_read_graph_partition_from_file(graph, 2, "tmp.partition");

    std::unique_ptr<mt_kahypar_partition_id_t[]> actual_partition =
      std::make_unique<mt_kahypar_partition_id_t[]>(7);
    mt_kahypar_get_graph_partition(partitioned_graph_2, actual_partition.get());

    ASSERT_EQ(3, mt_kahypar_graph_cut(partitioned_graph_2));
    for ( mt_kahypar_hypernode_id_t hn = 0; hn < 5; ++hn ) {
      ASSERT_EQ(partition[hn], actual_partition[hn]);
    }

    mt_kahypar_free_graph(graph);
    mt_kahypar_free_partitioned_graph(partitioned_graph);
    mt_kahypar_free_partitioned_graph(partitioned_graph_2);
  }

  namespace {
    mt_kahypar_partitioned_hypergraph_t* partition(const char* filename,
                                                   const mt_kahypar_file_format_type_t file_format,
                                                   const mt_kahypar_preset_type_t preset,
                                                   const mt_kahypar_partition_id_t num_blocks) {
      // Setup Partitioning Context
      mt_kahypar_context_t* context = mt_kahypar_context_new();
      mt_kahypar_load_preset(context, preset);
      mt_kahypar_set_partitioning_parameters(context, num_blocks, 0.03, KM1, 0);
      mt_kahypar_set_context_parameter(context, VERBOSE, debug ? "1" : "0");

      // Load Hypergraph
      mt_kahypar_hypergraph_t* hypergraph =
        mt_kahypar_read_hypergraph_from_file(filename, context, file_format);

      // Partition Hypergraph
      mt_kahypar_partitioned_hypergraph_t* partitioned_hg =
        mt_kahypar_partition_hypergraph(hypergraph, context);

      double imbalance = mt_kahypar_hypergraph_imbalance(partitioned_hg, context);
      mt_kahypar_hyperedge_weight_t objective = mt_kahypar_km1(partitioned_hg);
      if ( debug ) {
        LOG << " imbalance =" << imbalance << "\n"
            << "cut =" << mt_kahypar_hypergraph_cut(partitioned_hg) << "\n"
            << "km1 =" << objective << "\n"
            << "soed =" << mt_kahypar_soed(partitioned_hg);
      }
      EXPECT_LE(imbalance, 0.03);

      // Verify Partition IDs
      std::unique_ptr<mt_kahypar_partition_id_t[]> partition =
        std::make_unique<mt_kahypar_partition_id_t[]>(mt_kahypar_num_hypernodes(hypergraph));
      mt_kahypar_get_hypergraph_partition(partitioned_hg, partition.get());
      std::vector<mt_kahypar_hypernode_weight_t> expected_block_weights(num_blocks);
      for ( mt_kahypar_hypernode_id_t hn = 0; hn < mt_kahypar_num_hypernodes(hypergraph); ++hn ) {
        EXPECT_GE(partition[hn], 0);
        EXPECT_LT(partition[hn], num_blocks);
        ++expected_block_weights[partition[hn]];
      }

      // Verify Block Weights
      std::unique_ptr<mt_kahypar_hypernode_weight_t[]> block_weights =
        std::make_unique<mt_kahypar_hypernode_weight_t[]>(num_blocks);
      mt_kahypar_get_hypergraph_block_weights(partitioned_hg, block_weights.get());
      for ( mt_kahypar_partition_id_t i = 0; i < num_blocks; ++i ) {
        EXPECT_EQ(expected_block_weights[i], block_weights[i]);
      }

      mt_kahypar_free_context(context);
      mt_kahypar_free_hypergraph(hypergraph);
      return partitioned_hg;
    }
  }

  class APartitioner : public Test {
    private:
      static constexpr bool debug = false;
      static constexpr char HYPERGRAPH_FILE[] = "test_instances/ibm01.hgr";
      static constexpr char GRAPH_FILE[] = "test_instances/delaunay_n15.graph";


    public:
      APartitioner() :
        context(nullptr),
        hypergraph(nullptr),
        graph(nullptr),
        partitioned_hg(nullptr),
        partitioned_graph(nullptr) { }


    void PartitionGraph(const mt_kahypar_preset_type_t preset,
                        const mt_kahypar_partition_id_t num_blocks,
                        const double epsilon,
                        const mt_kahypar_objective_t objective,
                        const bool verbose = false) {
      SetUpContext(preset, num_blocks, epsilon, objective, verbose);
      LoadGraph();

      partitioned_graph = mt_kahypar_partition_graph(graph, context);

      double imbalance = mt_kahypar_graph_imbalance(partitioned_graph, context);
      mt_kahypar_hyperedge_weight_t cut = mt_kahypar_graph_cut(partitioned_graph);
      if ( debug ) {
        LOG << " imbalance =" << imbalance << "\n"
            << "cut =" << cut;
      }
      ASSERT_LE(imbalance, epsilon);

      // Verify Partition IDs
      std::unique_ptr<mt_kahypar_partition_id_t[]> partition =
        std::make_unique<mt_kahypar_partition_id_t[]>(mt_kahypar_num_nodes(graph));
      mt_kahypar_get_graph_partition(partitioned_graph, partition.get());
      std::vector<mt_kahypar_hypernode_weight_t> expected_block_weights(num_blocks);
      for ( mt_kahypar_hypernode_id_t hn = 0; hn < mt_kahypar_num_nodes(graph); ++hn ) {
        ASSERT_GE(partition[hn], 0);
        ASSERT_LT(partition[hn], num_blocks);
        ++expected_block_weights[partition[hn]];
      }

      // Verify Block Weights
      std::unique_ptr<mt_kahypar_hypernode_weight_t[]> block_weights =
        std::make_unique<mt_kahypar_hypernode_weight_t[]>(num_blocks);
      mt_kahypar_get_graph_block_weights(partitioned_graph, block_weights.get());
      for ( mt_kahypar_partition_id_t i = 0; i < num_blocks; ++i ) {
        EXPECT_EQ(expected_block_weights[i], block_weights[i]);
      }
    }

    void PartitionHypergraph(const mt_kahypar_preset_type_t preset,
                             const mt_kahypar_partition_id_t num_blocks,
                             const double epsilon,
                             const mt_kahypar_objective_t objective,
                             const bool verbose = false) {
      SetUpContext(preset, num_blocks, epsilon, objective, verbose);
      LoadHypergraph();

      partitioned_hg = mt_kahypar_partition_hypergraph(hypergraph, context);

      double imbalance = mt_kahypar_hypergraph_imbalance(partitioned_hg, context);
      mt_kahypar_hyperedge_weight_t km1 = mt_kahypar_km1(partitioned_hg);
      if ( debug ) {
        LOG << " imbalance =" << imbalance << "\n"
            << "cut =" << mt_kahypar_hypergraph_cut(partitioned_hg) << "\n"
            << "km1 =" << km1 << "\n"
            << "soed =" << mt_kahypar_soed(partitioned_hg);
      }
      ASSERT_LE(imbalance, epsilon);

      // Verify Partition IDs
      std::unique_ptr<mt_kahypar_partition_id_t[]> partition =
        std::make_unique<mt_kahypar_partition_id_t[]>(mt_kahypar_num_hypernodes(hypergraph));
      mt_kahypar_get_hypergraph_partition(partitioned_hg, partition.get());
      std::vector<mt_kahypar_hypernode_weight_t> expected_block_weights(num_blocks);
      for ( mt_kahypar_hypernode_id_t hn = 0; hn < mt_kahypar_num_hypernodes(hypergraph); ++hn ) {
        ASSERT_GE(partition[hn], 0);
        ASSERT_LT(partition[hn], num_blocks);
        ++expected_block_weights[partition[hn]];
      }

      // Verify Block Weights
      std::unique_ptr<mt_kahypar_hypernode_weight_t[]> block_weights =
        std::make_unique<mt_kahypar_hypernode_weight_t[]>(num_blocks);
      mt_kahypar_get_hypergraph_block_weights(partitioned_hg, block_weights.get());
      for ( mt_kahypar_partition_id_t i = 0; i < num_blocks; ++i ) {
        EXPECT_EQ(expected_block_weights[i], block_weights[i]);
      }
    }

    void ImproveGraphPartition(const mt_kahypar_preset_type_t preset,
                               const size_t num_vcycles,
                               const bool verbose = false) {
      mt_kahypar_load_preset(context, preset);
      mt_kahypar_set_context_parameter(context, VERBOSE, ( debug || verbose ) ? "1" : "0");

      mt_kahypar_hyperedge_weight_t before = mt_kahypar_graph_cut(partitioned_graph);
      mt_kahypar_improve_graph_partition(partitioned_graph, context, num_vcycles);
      mt_kahypar_hyperedge_weight_t after = mt_kahypar_graph_cut(partitioned_graph);
      ASSERT_LE(after, before);
    }

    void ImproveHypergraphPartition(const mt_kahypar_preset_type_t preset,
                                    const size_t num_vcycles,
                                    const bool verbose = false) {
      mt_kahypar_load_preset(context, preset);
      mt_kahypar_set_context_parameter(context, VERBOSE, ( debug || verbose ) ? "1" : "0");

      mt_kahypar_hyperedge_weight_t before = mt_kahypar_km1(partitioned_hg);
      mt_kahypar_improve_hypergraph_partition(partitioned_hg, context, num_vcycles);
      mt_kahypar_hyperedge_weight_t after = mt_kahypar_km1(partitioned_hg);
      ASSERT_LE(after, before);
    }

    void SetUp()  {
      mt_kahypar_initialize_thread_pool(std::thread::hardware_concurrency(), false);
      context = mt_kahypar_context_new();
    }

    void TearDown() {
      mt_kahypar_free_context(context);
      mt_kahypar_free_hypergraph(hypergraph);
      mt_kahypar_free_graph(graph);
      mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
      mt_kahypar_free_partitioned_graph(partitioned_graph);
    }

    mt_kahypar_context_t* context;
    mt_kahypar_hypergraph_t* hypergraph;
    mt_kahypar_graph_t* graph;
    mt_kahypar_partitioned_hypergraph_t* partitioned_hg;
    mt_kahypar_partitioned_graph_t* partitioned_graph;

   private:
    void LoadGraph() {
      if ( !graph ) {
        graph = mt_kahypar_read_graph_from_file(GRAPH_FILE, context, METIS);
      }
    }

    void LoadHypergraph() {
      if ( !graph ) {
        hypergraph = mt_kahypar_read_hypergraph_from_file(HYPERGRAPH_FILE, context, HMETIS);
      }
    }

    void SetUpContext(const mt_kahypar_preset_type_t preset,
                      const mt_kahypar_partition_id_t num_blocks,
                      const double epsilon,
                      const mt_kahypar_objective_t objective,
                      const bool verbose = false) {
      mt_kahypar_load_preset(context, preset);
      mt_kahypar_set_partitioning_parameters(context, num_blocks, epsilon, objective, 0);
      mt_kahypar_set_context_parameter(context, VERBOSE, ( debug || verbose ) ? "1" : "0");
    }
  };

  TEST_F(APartitioner, PartitionsAHypergraphInTwoBlocksWithSpeedPreset) {
    PartitionHypergraph(SPEED, 2, 0.03, KM1, false);
  }

  TEST_F(APartitioner, PartitionsAGraphInTwoBlocksWithSpeedPreset) {
    PartitionGraph(SPEED, 2, 0.03, CUT, false);
  }

  TEST_F(APartitioner, PartitionsAHypergraphInFourBlocksWithSpeedPreset) {
    PartitionHypergraph(SPEED, 4, 0.03, KM1, false);
  }

  TEST_F(APartitioner, PartitionsAGraphInFourBlocksWithSpeedPreset) {
    PartitionGraph(SPEED, 4, 0.03, CUT, false);
  }

  TEST_F(APartitioner, PartitionsAHypergraphInTwoBlocksWithHighQualityPreset) {
    PartitionHypergraph(HIGH_QUALITY, 2, 0.03, KM1, false);
  }

  TEST_F(APartitioner, PartitionsAGraphInTwoBlocksWithHighQualityPreset) {
    PartitionGraph(HIGH_QUALITY, 2, 0.03, CUT, false);
  }

  TEST_F(APartitioner, PartitionsAHypergraphInFourBlocksWithHighQualityPreset) {
    PartitionHypergraph(HIGH_QUALITY, 4, 0.03, KM1, false);
  }

  TEST_F(APartitioner, PartitionsAGraphInFourBlocksWithHighQualityPreset) {
    PartitionGraph(HIGH_QUALITY, 4, 0.03, CUT, false);
  }

  TEST_F(APartitioner, PartitionsAHypergraphInTwoBlocksWithDeterministicPreset) {
    PartitionHypergraph(DETERMINISTIC, 2, 0.03, KM1, false);
  }

  TEST_F(APartitioner, PartitionsAGraphInTwoBlocksWithDeterministicPreset) {
    PartitionGraph(DETERMINISTIC, 2, 0.03, CUT, false);
  }

  TEST_F(APartitioner, PartitionsAHypergraphInFourBlocksWithDeterministicPreset) {
    PartitionHypergraph(DETERMINISTIC, 4, 0.03, KM1, false);
  }

  TEST_F(APartitioner, PartitionsAGraphInFourBlocksWithDeterministicPreset) {
    PartitionGraph(DETERMINISTIC, 4, 0.03, CUT, false);
  }

  TEST_F(APartitioner, CanPartitionTwoHypergraphsSimultanously) {
    tbb::parallel_invoke([&]() {
      PartitionHypergraph(SPEED, 4, 0.03, KM1, false);
    }, [&] {
      PartitionGraph(DETERMINISTIC, 4, 0.03, CUT, false);
    });
  }

  TEST_F(APartitioner, ChecksIfDeterministicPresetProducesSameResultsForHypergraphs) {
    PartitionHypergraph(DETERMINISTIC, 8, 0.03, KM1, false);
    const double objective_1 = mt_kahypar_km1(partitioned_hg);
    PartitionHypergraph(DETERMINISTIC, 8, 0.03, KM1, false);
    const double objective_2 = mt_kahypar_km1(partitioned_hg);
    PartitionHypergraph(DETERMINISTIC, 8, 0.03, KM1, false);
    const double objective_3 = mt_kahypar_km1(partitioned_hg);
    ASSERT_EQ(objective_1, objective_2);
    ASSERT_EQ(objective_1, objective_3);
  }

  TEST_F(APartitioner, ChecksIfDeterministicPresetProducesSameResultsForGraphs) {
    PartitionGraph(DETERMINISTIC, 8, 0.03, CUT, false);
    const double objective_1 = mt_kahypar_graph_cut(partitioned_graph);
    PartitionGraph(DETERMINISTIC, 8, 0.03, CUT, false);
    const double objective_2 = mt_kahypar_graph_cut(partitioned_graph);
    PartitionGraph(DETERMINISTIC, 8, 0.03, CUT, false);
    const double objective_3 = mt_kahypar_graph_cut(partitioned_graph);
    ASSERT_EQ(objective_1, objective_2);
    ASSERT_EQ(objective_1, objective_3);
  }

  TEST_F(APartitioner, ImprovesHypergraphPartitionWithOneVCycle) {
    PartitionHypergraph(SPEED, 4, 0.03, KM1, false);
    ImproveHypergraphPartition(SPEED, 1, false);
  }

  TEST_F(APartitioner, ImprovesGraphPartitionWithOneVCycle) {
    PartitionGraph(SPEED, 4, 0.03, CUT, false);
    ImproveGraphPartition(SPEED, 1, false);
  }

  TEST_F(APartitioner, ImprovesHypergraphPartitionWithOneVCycleAndDifferentPresetType) {
    PartitionHypergraph(SPEED, 4, 0.03, KM1, false);
    ImproveHypergraphPartition(HIGH_QUALITY, 1, false);
  }

  TEST_F(APartitioner, ImprovesGraphPartitionWithOneVCycleAndDifferentPresetType) {
    PartitionGraph(SPEED, 4, 0.03, KM1, false);
    ImproveGraphPartition(HIGH_QUALITY, 1, false);
  }

  TEST_F(APartitioner, ImprovesHypergraphPartitionWithThreeVCycles) {
    PartitionHypergraph(SPEED, 4, 0.03, KM1, false);
    ImproveHypergraphPartition(HIGH_QUALITY, 3, false);
  }

  TEST_F(APartitioner, ImprovesGraphPartitionWithThreeVCycles) {
    PartitionGraph(SPEED, 4, 0.03, KM1, false);
    ImproveGraphPartition(SPEED, 3, false);
  }

  TEST_F(APartitioner, PartitionsHypergraphWithIndividualBlockWeights) {
    // Setup Individual Block Weights
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> block_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(4);
    block_weights[0] = 2131; block_weights[1] = 1213;
    block_weights[2] = 7287; block_weights[3] = 2501;
    mt_kahypar_set_individual_target_block_weights(context, 4, block_weights.get());

    PartitionHypergraph(SPEED, 4, 0.03, KM1, false);

    // Verify Block Weights
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> actual_block_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(4);
    mt_kahypar_get_hypergraph_block_weights(partitioned_hg, actual_block_weights.get());
    for ( mt_kahypar_partition_id_t i = 0; i < 4; ++i ) {
      ASSERT_LE(actual_block_weights[i], block_weights[i]);
    }
  }

  TEST_F(APartitioner, PartitionsGraphWithIndividualBlockWeights) {
    // Setup Individual Block Weights
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> block_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(4);
    block_weights[0] = 11201; block_weights[1] = 4384;
    block_weights[2] = 14174; block_weights[3] = 3989;
    mt_kahypar_set_individual_target_block_weights(context, 4, block_weights.get());

    PartitionGraph(SPEED, 4, 0.03, CUT, false);

    // Verify Block Weights
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> actual_block_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(4);
    mt_kahypar_get_graph_block_weights(partitioned_graph, actual_block_weights.get());
    for ( mt_kahypar_partition_id_t i = 0; i < 4; ++i ) {
      ASSERT_LE(actual_block_weights[i], block_weights[i]);
    }
  }

  TEST_F(APartitioner, PartitionsHypergraphWithIndividualBlockWeightsAndVCycle) {
    // Setup Individual Block Weights
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> block_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(4);
    block_weights[0] = 2131; block_weights[1] = 1213;
    block_weights[2] = 7287; block_weights[3] = 2501;
    mt_kahypar_set_individual_target_block_weights(context, 4, block_weights.get());

    PartitionHypergraph(SPEED, 4, 0.03, KM1, false);
    ImproveHypergraphPartition(SPEED, 1, false);

    // Verify Block Weights
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> actual_block_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(4);
    mt_kahypar_get_hypergraph_block_weights(partitioned_hg, actual_block_weights.get());
    for ( mt_kahypar_partition_id_t i = 0; i < 4; ++i ) {
      ASSERT_LE(actual_block_weights[i], block_weights[i]);
    }
  }

  TEST_F(APartitioner, PartitionsGraphWithIndividualBlockWeightsAndVCycle) {
    // Setup Individual Block Weights
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> block_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(4);
    block_weights[0] = 11201; block_weights[1] = 4384;
    block_weights[2] = 14174; block_weights[3] = 3989;
    mt_kahypar_set_individual_target_block_weights(context, 4, block_weights.get());

    PartitionGraph(SPEED, 4, 0.03, CUT, false);
    ImproveGraphPartition(SPEED, 1, false);

    // Verify Block Weights
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> actual_block_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(4);
    mt_kahypar_get_graph_block_weights(partitioned_graph, actual_block_weights.get());
    for ( mt_kahypar_partition_id_t i = 0; i < 4; ++i ) {
      ASSERT_LE(actual_block_weights[i], block_weights[i]);
    }
  }

  TEST(MtKaHyPar, CanSetContextParameter) {
    mt_kahypar_context_t* context = mt_kahypar_context_new();
    ASSERT_EQ(0, mt_kahypar_set_context_parameter(context, NUM_BLOCKS, "4"));
    ASSERT_EQ(0, mt_kahypar_set_context_parameter(context, EPSILON, "0.03"));
    ASSERT_EQ(0, mt_kahypar_set_context_parameter(context, OBJECTIVE, "km1"));
    ASSERT_EQ(0, mt_kahypar_set_context_parameter(context, SEED, "42"));
    ASSERT_EQ(0, mt_kahypar_set_context_parameter(context, NUM_VCYCLES, "3"));
    ASSERT_EQ(0, mt_kahypar_set_context_parameter(context, VERBOSE, "1"));


    Context& c = *reinterpret_cast<Context*>(context);
    ASSERT_EQ(4, c.partition.k);
    ASSERT_EQ(0.03, c.partition.epsilon);
    ASSERT_EQ(Objective::km1, c.partition.objective);
    ASSERT_EQ(42, c.partition.seed);
    ASSERT_EQ(3, c.partition.num_vcycles);
    ASSERT_TRUE(c.partition.verbose_output);

    mt_kahypar_free_context(context);
  }

  namespace {
    void checkIfContextAreEqual(const Context& lhs, const Context& rhs) {
      // partition
      ASSERT_EQ(lhs.partition.paradigm, rhs.partition.paradigm);
      ASSERT_EQ(lhs.partition.mode, rhs.partition.mode);
      ASSERT_EQ(lhs.partition.objective, rhs.partition.objective);
      ASSERT_EQ(lhs.partition.file_format, rhs.partition.file_format);
      ASSERT_EQ(lhs.partition.instance_type, rhs.partition.instance_type);
      ASSERT_EQ(lhs.partition.preset_type, rhs.partition.preset_type);
      ASSERT_EQ(lhs.partition.epsilon, rhs.partition.epsilon);
      ASSERT_EQ(lhs.partition.k, rhs.partition.k);
      ASSERT_EQ(lhs.partition.seed, rhs.partition.seed);
      ASSERT_EQ(lhs.partition.num_vcycles, rhs.partition.num_vcycles);
      ASSERT_EQ(lhs.partition.time_limit, rhs.partition.time_limit);
      ASSERT_EQ(lhs.partition.large_hyperedge_size_threshold_factor, rhs.partition.large_hyperedge_size_threshold_factor);
      ASSERT_EQ(lhs.partition.large_hyperedge_size_threshold, rhs.partition.large_hyperedge_size_threshold);
      ASSERT_EQ(lhs.partition.smallest_large_he_size_threshold, rhs.partition.smallest_large_he_size_threshold);
      ASSERT_EQ(lhs.partition.ignore_hyperedge_size_threshold, rhs.partition.ignore_hyperedge_size_threshold);
      ASSERT_EQ(lhs.partition.verbose_output, rhs.partition.verbose_output);
      ASSERT_EQ(lhs.partition.show_detailed_timings, rhs.partition.show_detailed_timings);
      ASSERT_EQ(lhs.partition.show_detailed_clustering_timings, rhs.partition.show_detailed_clustering_timings);
      ASSERT_EQ(lhs.partition.measure_detailed_uncontraction_timings, rhs.partition.measure_detailed_uncontraction_timings);
      ASSERT_EQ(lhs.partition.timings_output_depth, rhs.partition.timings_output_depth);
      ASSERT_EQ(lhs.partition.show_memory_consumption, rhs.partition.show_memory_consumption);
      ASSERT_EQ(lhs.partition.show_advanced_cut_analysis, rhs.partition.show_advanced_cut_analysis);
      ASSERT_EQ(lhs.partition.enable_progress_bar, rhs.partition.enable_progress_bar);
      ASSERT_EQ(lhs.partition.sp_process_output, rhs.partition.sp_process_output);
      ASSERT_EQ(lhs.partition.csv_output, rhs.partition.csv_output);
      ASSERT_EQ(lhs.partition.write_partition_file, rhs.partition.write_partition_file);
      ASSERT_EQ(lhs.partition.deterministic, rhs.partition.deterministic);

      // shared memory
      ASSERT_EQ(lhs.shared_memory.num_threads, rhs.shared_memory.num_threads);
      ASSERT_EQ(lhs.shared_memory.static_balancing_work_packages, rhs.shared_memory.static_balancing_work_packages);
      ASSERT_EQ(lhs.shared_memory.use_localized_random_shuffle, rhs.shared_memory.use_localized_random_shuffle);
      ASSERT_EQ(lhs.shared_memory.shuffle_block_size, rhs.shared_memory.shuffle_block_size);
      ASSERT_EQ(lhs.shared_memory.degree_of_parallelism, rhs.shared_memory.degree_of_parallelism);

      // preprocessing
      ASSERT_EQ(lhs.preprocessing.stable_construction_of_incident_edges,
                rhs.preprocessing.stable_construction_of_incident_edges);
      ASSERT_EQ(lhs.preprocessing.use_community_detection,
                rhs.preprocessing.use_community_detection);
      ASSERT_EQ(lhs.preprocessing.disable_community_detection_for_mesh_graphs,
                rhs.preprocessing.disable_community_detection_for_mesh_graphs);

      // community detection
      ASSERT_EQ(lhs.preprocessing.community_detection.edge_weight_function,
                rhs.preprocessing.community_detection.edge_weight_function);
      ASSERT_EQ(lhs.preprocessing.community_detection.max_pass_iterations,
                rhs.preprocessing.community_detection.max_pass_iterations);
      ASSERT_EQ(lhs.preprocessing.community_detection.low_memory_contraction,
                rhs.preprocessing.community_detection.low_memory_contraction);
      ASSERT_DOUBLE_EQ(lhs.preprocessing.community_detection.min_vertex_move_fraction,
                       rhs.preprocessing.community_detection.min_vertex_move_fraction);
      ASSERT_EQ(lhs.preprocessing.community_detection.vertex_degree_sampling_threshold,
                rhs.preprocessing.community_detection.vertex_degree_sampling_threshold);
      ASSERT_EQ(lhs.preprocessing.community_detection.num_sub_rounds_deterministic,
                rhs.preprocessing.community_detection.num_sub_rounds_deterministic);

      // coarsening
      ASSERT_EQ(lhs.coarsening.algorithm, rhs.coarsening.algorithm);
      ASSERT_EQ(lhs.coarsening.contraction_limit_multiplier, rhs.coarsening.contraction_limit_multiplier);
      ASSERT_EQ(lhs.coarsening.use_adaptive_edge_size, rhs.coarsening.use_adaptive_edge_size);
      ASSERT_EQ(lhs.coarsening.use_adaptive_max_allowed_node_weight, rhs.coarsening.use_adaptive_max_allowed_node_weight);
      ASSERT_EQ(lhs.coarsening.adaptive_node_weight_shrink_factor_threshold, rhs.coarsening.adaptive_node_weight_shrink_factor_threshold);
      ASSERT_EQ(lhs.coarsening.max_allowed_weight_multiplier, rhs.coarsening.max_allowed_weight_multiplier);
      ASSERT_EQ(lhs.coarsening.minimum_shrink_factor, rhs.coarsening.minimum_shrink_factor);
      ASSERT_EQ(lhs.coarsening.maximum_shrink_factor, rhs.coarsening.maximum_shrink_factor);
      ASSERT_EQ(lhs.coarsening.vertex_degree_sampling_threshold, rhs.coarsening.vertex_degree_sampling_threshold);
      ASSERT_EQ(lhs.coarsening.num_sub_rounds_deterministic, rhs.coarsening.num_sub_rounds_deterministic);
      ASSERT_EQ(lhs.coarsening.max_allowed_node_weight, rhs.coarsening.max_allowed_node_weight);
      ASSERT_EQ(lhs.coarsening.contraction_limit, rhs.coarsening.contraction_limit);

      // coarsening -> rating
      ASSERT_EQ(lhs.coarsening.rating.rating_function, rhs.coarsening.rating.rating_function);
      ASSERT_EQ(lhs.coarsening.rating.heavy_node_penalty_policy, rhs.coarsening.rating.heavy_node_penalty_policy);
      ASSERT_EQ(lhs.coarsening.rating.acceptance_policy, rhs.coarsening.rating.acceptance_policy);

      // initial partitioning
      ASSERT_EQ(lhs.initial_partitioning.mode, rhs.initial_partitioning.mode);
      ASSERT_EQ(lhs.initial_partitioning.runs, rhs.initial_partitioning.runs);
      ASSERT_EQ(lhs.initial_partitioning.use_adaptive_ip_runs, rhs.initial_partitioning.use_adaptive_ip_runs);
      ASSERT_EQ(lhs.initial_partitioning.min_adaptive_ip_runs, rhs.initial_partitioning.min_adaptive_ip_runs);
      ASSERT_EQ(lhs.initial_partitioning.perform_refinement_on_best_partitions, rhs.initial_partitioning.perform_refinement_on_best_partitions);
      ASSERT_EQ(lhs.initial_partitioning.fm_refinment_rounds, rhs.initial_partitioning.fm_refinment_rounds);
      ASSERT_EQ(lhs.initial_partitioning.remove_degree_zero_hns_before_ip, rhs.initial_partitioning.remove_degree_zero_hns_before_ip);
      ASSERT_EQ(lhs.initial_partitioning.lp_maximum_iterations, rhs.initial_partitioning.lp_maximum_iterations);
      ASSERT_EQ(lhs.initial_partitioning.lp_initial_block_size, rhs.initial_partitioning.lp_initial_block_size);
      ASSERT_EQ(lhs.initial_partitioning.population_size, rhs.initial_partitioning.population_size);

      // initial partitioning -> refinement
      ASSERT_EQ(lhs.initial_partitioning.refinement.refine_until_no_improvement,
                rhs.initial_partitioning.refinement.refine_until_no_improvement);
      ASSERT_EQ(lhs.initial_partitioning.refinement.relative_improvement_threshold,
                rhs.initial_partitioning.refinement.relative_improvement_threshold);
      ASSERT_EQ(lhs.initial_partitioning.refinement.max_batch_size,
                rhs.initial_partitioning.refinement.max_batch_size);
      ASSERT_EQ(lhs.initial_partitioning.refinement.min_border_vertices_per_thread,
                rhs.initial_partitioning.refinement.min_border_vertices_per_thread);


      // initial partitioning -> refinement -> label propagation
      ASSERT_EQ(lhs.initial_partitioning.refinement.label_propagation.algorithm,
                rhs.initial_partitioning.refinement.label_propagation.algorithm);
      ASSERT_EQ(lhs.initial_partitioning.refinement.label_propagation.maximum_iterations,
                rhs.initial_partitioning.refinement.label_propagation.maximum_iterations);
      ASSERT_EQ(lhs.initial_partitioning.refinement.label_propagation.rebalancing,
                rhs.initial_partitioning.refinement.label_propagation.rebalancing);
      ASSERT_EQ(lhs.initial_partitioning.refinement.label_propagation.execute_sequential,
                rhs.initial_partitioning.refinement.label_propagation.execute_sequential);
      ASSERT_EQ(lhs.initial_partitioning.refinement.label_propagation.hyperedge_size_activation_threshold,
                rhs.initial_partitioning.refinement.label_propagation.hyperedge_size_activation_threshold);


      // initial partitioning -> refinement -> fm
      ASSERT_EQ(lhs.initial_partitioning.refinement.fm.algorithm,
                rhs.initial_partitioning.refinement.fm.algorithm);
      ASSERT_EQ(lhs.initial_partitioning.refinement.fm.multitry_rounds,
                rhs.initial_partitioning.refinement.fm.multitry_rounds);
      ASSERT_EQ(lhs.initial_partitioning.refinement.fm.num_seed_nodes,
                rhs.initial_partitioning.refinement.fm.num_seed_nodes);
      ASSERT_EQ(lhs.initial_partitioning.refinement.fm.rollback_balance_violation_factor,
                rhs.initial_partitioning.refinement.fm.rollback_balance_violation_factor);
      ASSERT_EQ(lhs.initial_partitioning.refinement.fm.min_improvement,
                rhs.initial_partitioning.refinement.fm.min_improvement);
      ASSERT_EQ(lhs.initial_partitioning.refinement.fm.time_limit_factor,
                rhs.initial_partitioning.refinement.fm.time_limit_factor);
      ASSERT_EQ(lhs.initial_partitioning.refinement.fm.perform_moves_global,
                rhs.initial_partitioning.refinement.fm.perform_moves_global);
      ASSERT_EQ(lhs.initial_partitioning.refinement.fm.rollback_parallel,
                rhs.initial_partitioning.refinement.fm.rollback_parallel);
      ASSERT_EQ(lhs.initial_partitioning.refinement.fm.iter_moves_on_recalc,
                rhs.initial_partitioning.refinement.fm.iter_moves_on_recalc);
      ASSERT_EQ(lhs.initial_partitioning.refinement.fm.shuffle,
                rhs.initial_partitioning.refinement.fm.shuffle);
      ASSERT_EQ(lhs.initial_partitioning.refinement.fm.obey_minimal_parallelism,
                rhs.initial_partitioning.refinement.fm.obey_minimal_parallelism);
      ASSERT_EQ(lhs.initial_partitioning.refinement.fm.release_nodes,
                rhs.initial_partitioning.refinement.fm.release_nodes);

      // initial partitioning -> refinement -> flows
      ASSERT_EQ(lhs.initial_partitioning.refinement.flows.algorithm,
                rhs.initial_partitioning.refinement.flows.algorithm);
      ASSERT_EQ(lhs.initial_partitioning.refinement.flows.alpha,
                rhs.initial_partitioning.refinement.flows.alpha);
      ASSERT_EQ(lhs.initial_partitioning.refinement.flows.max_num_pins,
                rhs.initial_partitioning.refinement.flows.max_num_pins);
      ASSERT_EQ(lhs.initial_partitioning.refinement.flows.find_most_balanced_cut,
                rhs.initial_partitioning.refinement.flows.find_most_balanced_cut);
      ASSERT_EQ(lhs.initial_partitioning.refinement.flows.determine_distance_from_cut,
                rhs.initial_partitioning.refinement.flows.determine_distance_from_cut);
      ASSERT_EQ(lhs.initial_partitioning.refinement.flows.parallel_searches_multiplier,
                rhs.initial_partitioning.refinement.flows.parallel_searches_multiplier);
      ASSERT_EQ(lhs.initial_partitioning.refinement.flows.num_parallel_searches,
                rhs.initial_partitioning.refinement.flows.num_parallel_searches);
      ASSERT_EQ(lhs.initial_partitioning.refinement.flows.max_bfs_distance,
                rhs.initial_partitioning.refinement.flows.max_bfs_distance);
      ASSERT_EQ(lhs.initial_partitioning.refinement.flows.min_relative_improvement_per_round,
                rhs.initial_partitioning.refinement.flows.min_relative_improvement_per_round);
      ASSERT_EQ(lhs.initial_partitioning.refinement.flows.time_limit_factor,
                rhs.initial_partitioning.refinement.flows.time_limit_factor);
      ASSERT_EQ(lhs.initial_partitioning.refinement.flows.skip_small_cuts,
                rhs.initial_partitioning.refinement.flows.skip_small_cuts);
      ASSERT_EQ(lhs.initial_partitioning.refinement.flows.skip_unpromising_blocks,
                rhs.initial_partitioning.refinement.flows.skip_unpromising_blocks);
      ASSERT_EQ(lhs.initial_partitioning.refinement.flows.pierce_in_bulk,
                rhs.initial_partitioning.refinement.flows.pierce_in_bulk);

      // initial partitioning -> refinement -> deterministic
      ASSERT_EQ(lhs.initial_partitioning.refinement.deterministic_refinement.num_sub_rounds_sync_lp,
                rhs.initial_partitioning.refinement.deterministic_refinement.num_sub_rounds_sync_lp);
      ASSERT_EQ(lhs.initial_partitioning.refinement.deterministic_refinement.use_active_node_set,
                rhs.initial_partitioning.refinement.deterministic_refinement.use_active_node_set);
      ASSERT_EQ(lhs.initial_partitioning.refinement.deterministic_refinement.recalculate_gains_on_second_apply,
                rhs.initial_partitioning.refinement.deterministic_refinement.recalculate_gains_on_second_apply);

      // initial partitioning -> refinement -> global fm
      ASSERT_EQ(lhs.initial_partitioning.refinement.global_fm.use_global_fm,
                rhs.initial_partitioning.refinement.global_fm.use_global_fm);
      ASSERT_EQ(lhs.initial_partitioning.refinement.global_fm.refine_until_no_improvement,
                rhs.initial_partitioning.refinement.global_fm.refine_until_no_improvement);
      ASSERT_EQ(lhs.initial_partitioning.refinement.global_fm.num_seed_nodes,
                rhs.initial_partitioning.refinement.global_fm.num_seed_nodes);
      ASSERT_EQ(lhs.initial_partitioning.refinement.global_fm.obey_minimal_parallelism,
                rhs.initial_partitioning.refinement.global_fm.obey_minimal_parallelism);


      // refinement
      ASSERT_EQ(lhs.refinement.refine_until_no_improvement,
                rhs.refinement.refine_until_no_improvement);
      ASSERT_EQ(lhs.refinement.relative_improvement_threshold,
                rhs.refinement.relative_improvement_threshold);
      ASSERT_EQ(lhs.refinement.max_batch_size,
                rhs.refinement.max_batch_size);
      ASSERT_EQ(lhs.refinement.min_border_vertices_per_thread,
                rhs.refinement.min_border_vertices_per_thread);


      // refinement -> label propagation
      ASSERT_EQ(lhs.refinement.label_propagation.algorithm,
                rhs.refinement.label_propagation.algorithm);
      ASSERT_EQ(lhs.refinement.label_propagation.maximum_iterations,
                rhs.refinement.label_propagation.maximum_iterations);
      ASSERT_EQ(lhs.refinement.label_propagation.rebalancing,
                rhs.refinement.label_propagation.rebalancing);
      ASSERT_EQ(lhs.refinement.label_propagation.execute_sequential,
                rhs.refinement.label_propagation.execute_sequential);
      ASSERT_EQ(lhs.refinement.label_propagation.hyperedge_size_activation_threshold,
                rhs.refinement.label_propagation.hyperedge_size_activation_threshold);


      // refinement -> fm
      ASSERT_EQ(lhs.refinement.fm.algorithm,
                rhs.refinement.fm.algorithm);
      ASSERT_EQ(lhs.refinement.fm.multitry_rounds,
                rhs.refinement.fm.multitry_rounds);
      ASSERT_EQ(lhs.refinement.fm.num_seed_nodes,
                rhs.refinement.fm.num_seed_nodes);
      ASSERT_EQ(lhs.refinement.fm.rollback_balance_violation_factor,
                rhs.refinement.fm.rollback_balance_violation_factor);
      ASSERT_EQ(lhs.refinement.fm.min_improvement,
                rhs.refinement.fm.min_improvement);
      ASSERT_EQ(lhs.refinement.fm.time_limit_factor,
                rhs.refinement.fm.time_limit_factor);
      ASSERT_EQ(lhs.refinement.fm.perform_moves_global,
                rhs.refinement.fm.perform_moves_global);
      ASSERT_EQ(lhs.refinement.fm.rollback_parallel,
                rhs.refinement.fm.rollback_parallel);
      ASSERT_EQ(lhs.refinement.fm.iter_moves_on_recalc,
                rhs.refinement.fm.iter_moves_on_recalc);
      ASSERT_EQ(lhs.refinement.fm.shuffle,
                rhs.refinement.fm.shuffle);
      ASSERT_EQ(lhs.refinement.fm.obey_minimal_parallelism,
                rhs.refinement.fm.obey_minimal_parallelism);
      ASSERT_EQ(lhs.refinement.fm.release_nodes,
                rhs.refinement.fm.release_nodes);

      // refinement -> flows
      ASSERT_EQ(lhs.refinement.flows.algorithm,
                rhs.refinement.flows.algorithm);
      ASSERT_EQ(lhs.refinement.flows.alpha,
                rhs.refinement.flows.alpha);
      ASSERT_EQ(lhs.refinement.flows.max_num_pins,
                rhs.refinement.flows.max_num_pins);
      ASSERT_EQ(lhs.refinement.flows.find_most_balanced_cut,
                rhs.refinement.flows.find_most_balanced_cut);
      ASSERT_EQ(lhs.refinement.flows.determine_distance_from_cut,
                rhs.refinement.flows.determine_distance_from_cut);
      ASSERT_EQ(lhs.refinement.flows.parallel_searches_multiplier,
                rhs.refinement.flows.parallel_searches_multiplier);
      ASSERT_EQ(lhs.refinement.flows.num_parallel_searches,
                rhs.refinement.flows.num_parallel_searches);
      ASSERT_EQ(lhs.refinement.flows.max_bfs_distance,
                rhs.refinement.flows.max_bfs_distance);
      ASSERT_EQ(lhs.refinement.flows.min_relative_improvement_per_round,
                rhs.refinement.flows.min_relative_improvement_per_round);
      ASSERT_EQ(lhs.refinement.flows.time_limit_factor,
                rhs.refinement.flows.time_limit_factor);
      ASSERT_EQ(lhs.refinement.flows.skip_small_cuts,
                rhs.refinement.flows.skip_small_cuts);
      ASSERT_EQ(lhs.refinement.flows.skip_unpromising_blocks,
                rhs.refinement.flows.skip_unpromising_blocks);
      ASSERT_EQ(lhs.refinement.flows.pierce_in_bulk,
                rhs.refinement.flows.pierce_in_bulk);

      // refinement -> deterministic
      ASSERT_EQ(lhs.refinement.deterministic_refinement.num_sub_rounds_sync_lp,
                rhs.refinement.deterministic_refinement.num_sub_rounds_sync_lp);
      ASSERT_EQ(lhs.refinement.deterministic_refinement.use_active_node_set,
                rhs.refinement.deterministic_refinement.use_active_node_set);
      ASSERT_EQ(lhs.refinement.deterministic_refinement.recalculate_gains_on_second_apply,
                rhs.refinement.deterministic_refinement.recalculate_gains_on_second_apply);

      // refinement -> global fm
      ASSERT_EQ(lhs.refinement.global_fm.use_global_fm,
                rhs.refinement.global_fm.use_global_fm);
      ASSERT_EQ(lhs.refinement.global_fm.refine_until_no_improvement,
                rhs.refinement.global_fm.refine_until_no_improvement);
      ASSERT_EQ(lhs.refinement.global_fm.num_seed_nodes,
                rhs.refinement.global_fm.num_seed_nodes);
      ASSERT_EQ(lhs.refinement.global_fm.obey_minimal_parallelism,
                rhs.refinement.global_fm.obey_minimal_parallelism);
    }
  }

  TEST(MtKaHyPar, LoadDefaultPreset) {
    mt_kahypar_context_t* default_preset = mt_kahypar_context_new();
    mt_kahypar_load_preset(default_preset, SPEED);
    mt_kahypar_context_t* default_preset_ini = mt_kahypar_context_new();
    mt_kahypar_configure_context_from_file(default_preset_ini, "../../../config/default_preset.ini");

    Context& default_context = *reinterpret_cast<Context*>(default_preset);
    Context& default_context_ini = *reinterpret_cast<Context*>(default_preset_ini);

    checkIfContextAreEqual(default_context, default_context_ini);

    mt_kahypar_free_context(default_preset);
    mt_kahypar_free_context(default_preset_ini);
  }

  TEST(MtKaHyPar, LoadDefaultFlowPreset) {
    mt_kahypar_context_t* default_flow_preset = mt_kahypar_context_new();
    mt_kahypar_load_preset(default_flow_preset, HIGH_QUALITY);
    mt_kahypar_context_t* default_flow_preset_ini = mt_kahypar_context_new();
    mt_kahypar_configure_context_from_file(default_flow_preset_ini, "../../../config/default_flow_preset.ini");

    Context& default_flow_context = *reinterpret_cast<Context*>(default_flow_preset);
    Context& default_flow_context_ini = *reinterpret_cast<Context*>(default_flow_preset_ini);

    checkIfContextAreEqual(default_flow_context, default_flow_context_ini);

    mt_kahypar_free_context(default_flow_preset);
    mt_kahypar_free_context(default_flow_preset_ini);
  }

  TEST(MtKaHyPar, LoadDeterministicPreset) {
    mt_kahypar_context_t* deterministic_preset = mt_kahypar_context_new();
    mt_kahypar_load_preset(deterministic_preset, DETERMINISTIC);
    mt_kahypar_context_t* deterministic_preset_ini = mt_kahypar_context_new();
    mt_kahypar_configure_context_from_file(deterministic_preset_ini, "../../../config/deterministic_preset.ini");

    Context& deterministic_context = *reinterpret_cast<Context*>(deterministic_preset);
    Context& deterministic_context_ini = *reinterpret_cast<Context*>(deterministic_preset_ini);

    checkIfContextAreEqual(deterministic_context, deterministic_context_ini);

    mt_kahypar_free_context(deterministic_preset);
    mt_kahypar_free_context(deterministic_preset_ini);
  }

}
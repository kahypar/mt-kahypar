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

#include <tbb/parallel_invoke.h>

#include "mtkahypar.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/io/hypergraph_io.h"

using ::testing::Test;

namespace mt_kahypar {
  TEST(MtKaHyPar, LoadsContextFromFile) {
    mt_kahypar_error_t error;
    mt_kahypar_context_t* context = mt_kahypar_context_from_file("test_preset.ini", &error);
    // verbose should be false by default
    ASSERT_FALSE(reinterpret_cast<Context*>(context)->partition.verbose_output);
    ASSERT_EQ(DEFAULT, mt_kahypar_get_preset(context));
    mt_kahypar_free_context(context);

    context = mt_kahypar_context_from_file("error.ini", &error);
    ASSERT_EQ(nullptr, context);
    ASSERT_EQ(INVALID_INPUT, error.status);
    mt_kahypar_free_error_content(&error);
  }

  TEST(MtKaHyPar, CanSetAndGetPartitioningParameters) {
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
    mt_kahypar_set_partitioning_parameters(context, 8, 0.05, CUT);
    ASSERT_EQ(DEFAULT, mt_kahypar_get_preset(context));
    ASSERT_EQ(8, mt_kahypar_get_num_blocks(context));
    ASSERT_EQ(0.05, mt_kahypar_get_epsilon(context));
    ASSERT_EQ(CUT, mt_kahypar_get_objective(context));
    mt_kahypar_free_context(context);
  }

  TEST(MtKaHyPar, CanSetContextParameter) {
    mt_kahypar_error_t error;
    auto check_error_status = [&]{
      ASSERT_EQ(INVALID_PARAMETER, error.status);
      ASSERT_NE(nullptr, error.msg);
      mt_kahypar_free_error_content(&error);
      error = mt_kahypar_error_t{};
    };

    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
    ASSERT_EQ(0, mt_kahypar_set_context_parameter(context, NUM_BLOCKS, "4", &error));
    ASSERT_EQ(0, mt_kahypar_set_context_parameter(context, EPSILON, "0.03", &error));
    ASSERT_EQ(0, mt_kahypar_set_context_parameter(context, OBJECTIVE, "km1", &error));
    ASSERT_EQ(0, mt_kahypar_set_context_parameter(context, NUM_VCYCLES, "0", &error));
    ASSERT_EQ(0, mt_kahypar_set_context_parameter(context, NUM_VCYCLES, "3", &error));
    ASSERT_EQ(0, mt_kahypar_set_context_parameter(context, VERBOSE, "1", &error));

    ASSERT_EQ(INVALID_PARAMETER, mt_kahypar_set_context_parameter(context, NUM_BLOCKS, "x", &error));
    check_error_status();
    ASSERT_EQ(INVALID_PARAMETER, mt_kahypar_set_context_parameter(context, EPSILON, "", &error));
    check_error_status();
    ASSERT_EQ(INVALID_PARAMETER, mt_kahypar_set_context_parameter(context, OBJECTIVE, "wrong_objective", &error));
    check_error_status();
    ASSERT_EQ(INVALID_PARAMETER, mt_kahypar_set_context_parameter(context, NUM_VCYCLES, "a0", &error));
    check_error_status();
    ASSERT_EQ(INVALID_PARAMETER, mt_kahypar_set_context_parameter(context, VERBOSE, "2", &error));
    check_error_status();

    Context& c = *reinterpret_cast<Context*>(context);
    ASSERT_EQ(4, c.partition.k);
    ASSERT_EQ(0.03, c.partition.epsilon);
    ASSERT_EQ(Objective::km1, c.partition.objective);
    ASSERT_EQ(3, c.partition.num_vcycles);
    ASSERT_TRUE(c.partition.verbose_output);

    mt_kahypar_free_context(context);
  }

  TEST(MtKaHyPar, ReadHypergraphFile) {
    mt_kahypar_error_t error;
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
    mt_kahypar_hypergraph_t hypergraph =
      mt_kahypar_read_hypergraph_from_file("test_instances/ibm01.hgr", context, HMETIS, &error);

    ASSERT_EQ(12752, mt_kahypar_num_hypernodes(hypergraph));
    ASSERT_EQ(14111, mt_kahypar_num_hyperedges(hypergraph));
    ASSERT_EQ(50566, mt_kahypar_num_pins(hypergraph));
    ASSERT_EQ(12752, mt_kahypar_hypergraph_weight(hypergraph));

    mt_kahypar_free_hypergraph(hypergraph);
  }

  TEST(MtKaHyPar, ReadGraphFile) {
    mt_kahypar_error_t error;
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
    mt_kahypar_hypergraph_t graph =
      mt_kahypar_read_hypergraph_from_file("test_instances/delaunay_n15.graph", context, METIS, &error);

    ASSERT_EQ(32768,  mt_kahypar_num_hypernodes(graph));
    ASSERT_EQ(196548, mt_kahypar_num_hyperedges(graph));
    ASSERT_EQ(196548, mt_kahypar_num_pins(graph));
    ASSERT_EQ(32768,  mt_kahypar_hypergraph_weight(graph));

    mt_kahypar_free_hypergraph(graph);
  }

  TEST(MtKaHyPar, ConstructUnweightedStaticHypergraph) {
    mt_kahypar_error_t error;
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
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

    mt_kahypar_hypergraph_t hypergraph = mt_kahypar_create_hypergraph(
      context, num_vertices, num_hyperedges, hyperedge_indices.get(),
      hyperedges.get(), nullptr, nullptr, &error);
    ASSERT_EQ(hypergraph.type, STATIC_HYPERGRAPH);

    ASSERT_EQ(7, mt_kahypar_num_hypernodes(hypergraph));
    ASSERT_EQ(4, mt_kahypar_num_hyperedges(hypergraph));
    ASSERT_EQ(12, mt_kahypar_num_pins(hypergraph));
    ASSERT_EQ(7, mt_kahypar_hypergraph_weight(hypergraph));

    mt_kahypar_free_hypergraph(hypergraph);
  }

  TEST(MtKaHyPar, ConstructUnweightedDynamicHypergraph) {
    mt_kahypar_error_t error;
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(HIGHEST_QUALITY);
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

    mt_kahypar_hypergraph_t hypergraph = mt_kahypar_create_hypergraph(
      context, num_vertices, num_hyperedges, hyperedge_indices.get(),
      hyperedges.get(), nullptr, nullptr, &error);
    ASSERT_EQ(hypergraph.type, DYNAMIC_HYPERGRAPH);

    ASSERT_EQ(7, mt_kahypar_num_hypernodes(hypergraph));
    ASSERT_EQ(4, mt_kahypar_num_hyperedges(hypergraph));
    ASSERT_EQ(12, mt_kahypar_num_pins(hypergraph));
    ASSERT_EQ(7, mt_kahypar_hypergraph_weight(hypergraph));

    mt_kahypar_free_hypergraph(hypergraph);
  }

  TEST(MtKaHyPar, ConstructUnweightedStaticGraph) {
    mt_kahypar_error_t error;
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
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

    mt_kahypar_hypergraph_t graph = mt_kahypar_create_graph(
      context, num_vertices, num_hyperedges, edges.get(), nullptr, nullptr, &error);
    ASSERT_EQ(graph.type, STATIC_GRAPH);

    ASSERT_EQ(5, mt_kahypar_num_hypernodes(graph));
    ASSERT_EQ(12, mt_kahypar_num_hyperedges(graph));
    ASSERT_EQ(5, mt_kahypar_hypergraph_weight(graph));

    mt_kahypar_free_hypergraph(graph);
  }

  TEST(MtKaHyPar, ConstructUnweightedDynamicGraph) {
    mt_kahypar_error_t error;
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(HIGHEST_QUALITY);
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

    mt_kahypar_hypergraph_t graph = mt_kahypar_create_graph(
      context, num_vertices, num_hyperedges, edges.get(), nullptr, nullptr, &error);
    ASSERT_EQ(graph.type, DYNAMIC_GRAPH);

    ASSERT_EQ(5, mt_kahypar_num_hypernodes(graph));
    ASSERT_EQ(12, mt_kahypar_num_hyperedges(graph));
    ASSERT_EQ(5, mt_kahypar_hypergraph_weight(graph));

    mt_kahypar_free_hypergraph(graph);
  }

  TEST(MtKaHyPar, ConstructHypergraphWithNodeWeights) {
    mt_kahypar_error_t error;
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
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

    mt_kahypar_hypergraph_t hypergraph = mt_kahypar_create_hypergraph(
      context, num_vertices, num_hyperedges, hyperedge_indices.get(),
      hyperedges.get(), nullptr, vertex_weights.get(), &error);

    ASSERT_EQ(7, mt_kahypar_num_hypernodes(hypergraph));
    ASSERT_EQ(4, mt_kahypar_num_hyperedges(hypergraph));
    ASSERT_EQ(12, mt_kahypar_num_pins(hypergraph));
    ASSERT_EQ(28, mt_kahypar_hypergraph_weight(hypergraph));

    mt_kahypar_free_hypergraph(hypergraph);
  }

  TEST(MtKaHyPar, ConstructGraphWithNodeWeights) {
    mt_kahypar_error_t error;
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
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

    mt_kahypar_hypergraph_t graph = mt_kahypar_create_graph(
      context, num_vertices, num_hyperedges, edges.get(), nullptr, vertex_weights.get(), &error);

    ASSERT_EQ(5, mt_kahypar_num_hypernodes(graph));
    ASSERT_EQ(12, mt_kahypar_num_hyperedges(graph));
    ASSERT_EQ(15, mt_kahypar_hypergraph_weight(graph));

    mt_kahypar_free_hypergraph(graph);
  }

  TEST(MtKaHyPar, GetsPropertiesOfHypergraphWithNodeWeights) {
    mt_kahypar_error_t error;
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
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

    mt_kahypar_hypergraph_t hypergraph = mt_kahypar_create_hypergraph(
      context, num_vertices, num_hyperedges, hyperedge_indices.get(),
      hyperedges.get(), nullptr, vertex_weights.get(), &error);
    ASSERT_FALSE(mt_kahypar_is_graph(hypergraph));

    ASSERT_EQ(2, mt_kahypar_hypernode_degree(hypergraph, 0));
    ASSERT_EQ(1, mt_kahypar_hypernode_degree(hypergraph, 1));
    ASSERT_EQ(2, mt_kahypar_hypernode_degree(hypergraph, 2));
    ASSERT_EQ(2, mt_kahypar_hypernode_degree(hypergraph, 3));
    ASSERT_EQ(2, mt_kahypar_hypernode_degree(hypergraph, 4));
    ASSERT_EQ(1, mt_kahypar_hypernode_degree(hypergraph, 5));
    ASSERT_EQ(2, mt_kahypar_hypernode_degree(hypergraph, 6));
    ASSERT_EQ(2, mt_kahypar_hyperedge_size(hypergraph, 0));
    ASSERT_EQ(4, mt_kahypar_hyperedge_size(hypergraph, 1));
    ASSERT_EQ(3, mt_kahypar_hyperedge_size(hypergraph, 2));
    ASSERT_EQ(3, mt_kahypar_hyperedge_size(hypergraph, 3));

    for (size_t i = 0; i < 7; ++i) {
      ASSERT_EQ(vertex_weights[i], mt_kahypar_hypernode_weight(hypergraph, i));
    }
    for (size_t i = 0; i < 4; ++i) {
      ASSERT_EQ(1, mt_kahypar_hyperedge_weight(hypergraph, i));
    }

    mt_kahypar_free_hypergraph(hypergraph);
  }

  TEST(MtKaHyPar, GetsPropertiesOfGraphWithNodeWeights) {
    mt_kahypar_error_t error;
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
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

    mt_kahypar_hypergraph_t graph = mt_kahypar_create_graph(
      context, num_vertices, num_hyperedges, edges.get(), nullptr, vertex_weights.get(), &error);
    ASSERT_TRUE(mt_kahypar_is_graph(graph));

    ASSERT_EQ(2, mt_kahypar_hypernode_degree(graph, 0));
    ASSERT_EQ(3, mt_kahypar_hypernode_degree(graph, 1));
    ASSERT_EQ(3, mt_kahypar_hypernode_degree(graph, 2));
    ASSERT_EQ(3, mt_kahypar_hypernode_degree(graph, 3));
    ASSERT_EQ(1, mt_kahypar_hypernode_degree(graph, 4));

    for (size_t i = 0; i < 5; ++i) {
      ASSERT_EQ(vertex_weights[i], mt_kahypar_hypernode_weight(graph, i));
    }
    for (size_t i = 0; i < 12; ++i) {
      ASSERT_EQ(1, mt_kahypar_hyperedge_weight(graph, i));
      ASSERT_EQ(2, mt_kahypar_hyperedge_size(graph, i));
    }

    ASSERT_EQ(0, mt_kahypar_edge_source(graph, 0));  ASSERT_EQ(1, mt_kahypar_edge_target(graph, 0));
    ASSERT_EQ(0, mt_kahypar_edge_source(graph, 1));  ASSERT_EQ(2, mt_kahypar_edge_target(graph, 1));
    ASSERT_EQ(1, mt_kahypar_edge_source(graph, 2));  ASSERT_EQ(0, mt_kahypar_edge_target(graph, 2));
    ASSERT_EQ(1, mt_kahypar_edge_source(graph, 3));  ASSERT_EQ(2, mt_kahypar_edge_target(graph, 3));
    ASSERT_EQ(1, mt_kahypar_edge_source(graph, 4));  ASSERT_EQ(3, mt_kahypar_edge_target(graph, 4));
    ASSERT_EQ(2, mt_kahypar_edge_source(graph, 5));  ASSERT_EQ(0, mt_kahypar_edge_target(graph, 5));
    ASSERT_EQ(2, mt_kahypar_edge_source(graph, 6));  ASSERT_EQ(1, mt_kahypar_edge_target(graph, 6));
    ASSERT_EQ(2, mt_kahypar_edge_source(graph, 7));  ASSERT_EQ(3, mt_kahypar_edge_target(graph, 7));
    ASSERT_EQ(3, mt_kahypar_edge_source(graph, 8));  ASSERT_EQ(1, mt_kahypar_edge_target(graph, 8));
    ASSERT_EQ(3, mt_kahypar_edge_source(graph, 9));  ASSERT_EQ(2, mt_kahypar_edge_target(graph, 9));
    ASSERT_EQ(3, mt_kahypar_edge_source(graph, 10)); ASSERT_EQ(4, mt_kahypar_edge_target(graph, 10));
    ASSERT_EQ(4, mt_kahypar_edge_source(graph, 11)); ASSERT_EQ(3, mt_kahypar_edge_target(graph, 11));

    mt_kahypar_free_hypergraph(graph);
  }

  TEST(MtKaHyPar, IteratesOverUnweightedStaticHypergraph) {
    mt_kahypar_error_t error;
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
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

    mt_kahypar_hypergraph_t hypergraph = mt_kahypar_create_hypergraph(
      context, num_vertices, num_hyperedges, hyperedge_indices.get(),
      hyperedges.get(), nullptr, nullptr, &error);
    ASSERT_EQ(hypergraph.type, STATIC_HYPERGRAPH);

    std::unique_ptr<mt_kahypar_hyperedge_id_t[]> incidence_buffer = std::make_unique<mt_kahypar_hyperedge_id_t[]>(2);
    ASSERT_EQ(2, mt_kahypar_get_incident_hyperedges(hypergraph, 0, incidence_buffer.get()));
    ASSERT_EQ(0, incidence_buffer[0]);
    ASSERT_EQ(1, incidence_buffer[1]);
    ASSERT_EQ(1, mt_kahypar_get_incident_hyperedges(hypergraph, 1, incidence_buffer.get()));
    ASSERT_EQ(1, incidence_buffer[0]);
    ASSERT_EQ(2, mt_kahypar_get_incident_hyperedges(hypergraph, 2, incidence_buffer.get()));
    ASSERT_EQ(0, incidence_buffer[0]);
    ASSERT_EQ(3, incidence_buffer[1]);
    ASSERT_EQ(2, mt_kahypar_get_incident_hyperedges(hypergraph, 3, incidence_buffer.get()));
    ASSERT_EQ(1, incidence_buffer[0]);
    ASSERT_EQ(2, incidence_buffer[1]);
    ASSERT_EQ(2, mt_kahypar_get_incident_hyperedges(hypergraph, 4, incidence_buffer.get()));
    ASSERT_EQ(1, incidence_buffer[0]);
    ASSERT_EQ(2, incidence_buffer[1]);
    ASSERT_EQ(1, mt_kahypar_get_incident_hyperedges(hypergraph, 5, incidence_buffer.get()));
    ASSERT_EQ(3, incidence_buffer[0]);
    ASSERT_EQ(2, mt_kahypar_get_incident_hyperedges(hypergraph, 6, incidence_buffer.get()));
    ASSERT_EQ(2, incidence_buffer[0]);
    ASSERT_EQ(3, incidence_buffer[1]);

    std::unique_ptr<mt_kahypar_hypernode_id_t[]> pin_buffer = std::make_unique<mt_kahypar_hypernode_id_t[]>(4);
    ASSERT_EQ(2, mt_kahypar_get_hyperedge_pins(hypergraph, 0, pin_buffer.get()));
    ASSERT_EQ(0, pin_buffer[0]);
    ASSERT_EQ(2, pin_buffer[1]);
    ASSERT_EQ(4, mt_kahypar_get_hyperedge_pins(hypergraph, 1, pin_buffer.get()));
    ASSERT_EQ(0, pin_buffer[0]);
    ASSERT_EQ(1, pin_buffer[1]);
    ASSERT_EQ(3, pin_buffer[2]);
    ASSERT_EQ(4, pin_buffer[3]);
    ASSERT_EQ(3, mt_kahypar_get_hyperedge_pins(hypergraph, 2, pin_buffer.get()));
    ASSERT_EQ(3, pin_buffer[0]);
    ASSERT_EQ(4, pin_buffer[1]);
    ASSERT_EQ(6, pin_buffer[2]);
    ASSERT_EQ(3, mt_kahypar_get_hyperedge_pins(hypergraph, 3, pin_buffer.get()));
    ASSERT_EQ(2, pin_buffer[0]);
    ASSERT_EQ(5, pin_buffer[1]);
    ASSERT_EQ(6, pin_buffer[2]);

    mt_kahypar_free_hypergraph(hypergraph);
  }

  TEST(MtKaHyPar, CreatesPartitionedHypergraph) {
    mt_kahypar_error_t error;
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
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

    mt_kahypar_hypergraph_t hypergraph = mt_kahypar_create_hypergraph(
      context, num_vertices, num_hyperedges, hyperedge_indices.get(), hyperedges.get(), nullptr, nullptr, &error);

    std::unique_ptr<mt_kahypar_partition_id_t[]> partition =
      std::make_unique<mt_kahypar_partition_id_t[]>(7);
    partition[0] = 0; partition[1] = 0; partition[2] = 0;
    partition[3] = 1; partition[4] = 1; partition[5] = 1; partition[6] = 1;

    mt_kahypar_partitioned_hypergraph_t partitioned_hg =
      mt_kahypar_create_partitioned_hypergraph(hypergraph, context, 2, partition.get(), &error);
    ASSERT_EQ(partitioned_hg.type, MULTILEVEL_HYPERGRAPH_PARTITIONING);

    std::unique_ptr<mt_kahypar_partition_id_t[]> actual_partition =
      std::make_unique<mt_kahypar_partition_id_t[]>(7);
    mt_kahypar_get_partition(partitioned_hg, actual_partition.get());

    ASSERT_EQ(2, mt_kahypar_km1(partitioned_hg));
    for ( mt_kahypar_hypernode_id_t hn = 0; hn < 7; ++hn ) {
      ASSERT_EQ(partition[hn], actual_partition[hn]);
    }

    mt_kahypar_free_hypergraph(hypergraph);
    mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
  }

  TEST(MtKaHyPar, CreatesPartitionedGraph) {
    mt_kahypar_error_t error;
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
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

    mt_kahypar_hypergraph_t graph = mt_kahypar_create_graph(
      context, num_vertices, num_hyperedges, edges.get(), nullptr, vertex_weights.get(), &error);

    std::unique_ptr<mt_kahypar_partition_id_t[]> partition =
      std::make_unique<mt_kahypar_partition_id_t[]>(7);
    partition[0] = 0; partition[1] = 0; partition[2] = 1;
    partition[3] = 1; partition[4] = 1;

    mt_kahypar_partitioned_hypergraph_t partitioned_graph =
      mt_kahypar_create_partitioned_hypergraph(graph, context, 2, partition.get(), &error);
    ASSERT_EQ(partitioned_graph.type, MULTILEVEL_GRAPH_PARTITIONING);

    std::unique_ptr<mt_kahypar_partition_id_t[]> actual_partition =
      std::make_unique<mt_kahypar_partition_id_t[]>(7);
    mt_kahypar_get_partition(partitioned_graph, actual_partition.get());

    ASSERT_EQ(3, mt_kahypar_cut(partitioned_graph));
    for ( mt_kahypar_hypernode_id_t hn = 0; hn < 7; ++hn ) {
      ASSERT_EQ(partition[hn], actual_partition[hn]);
    }

    mt_kahypar_free_hypergraph(graph);
    mt_kahypar_free_partitioned_hypergraph(partitioned_graph);
  }

  TEST(MtKaHyPar, WritesAndLoadsHypergraphPartitionFile) {
    mt_kahypar_error_t error;
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
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

    mt_kahypar_hypergraph_t hypergraph = mt_kahypar_create_hypergraph(
      context, num_vertices, num_hyperedges, hyperedge_indices.get(), hyperedges.get(), nullptr, nullptr, &error);

    std::unique_ptr<mt_kahypar_partition_id_t[]> partition = std::make_unique<mt_kahypar_partition_id_t[]>(7);
    partition[0] = 0; partition[1] = 0; partition[2] = 0;
    partition[3] = 1; partition[4] = 1; partition[5] = 1; partition[6] = 1;

    mt_kahypar_partitioned_hypergraph_t partitioned_hg =
      mt_kahypar_create_partitioned_hypergraph(hypergraph, context, 2, partition.get(), &error);

    mt_kahypar_write_partition_to_file(partitioned_hg, "tmp.partition", &error);

    mt_kahypar_partitioned_hypergraph_t partitioned_hg_2 =
      mt_kahypar_read_partition_from_file(hypergraph, context, 2, "tmp.partition", &error);

    std::unique_ptr<mt_kahypar_partition_id_t[]> actual_partition =
      std::make_unique<mt_kahypar_partition_id_t[]>(7);
    mt_kahypar_get_partition(partitioned_hg_2, actual_partition.get());

    ASSERT_EQ(2, mt_kahypar_km1(partitioned_hg_2));
    for ( mt_kahypar_hypernode_id_t hn = 0; hn < 5; ++hn ) {
      ASSERT_EQ(partition[hn], actual_partition[hn]);
    }

    mt_kahypar_free_hypergraph(hypergraph);
    mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
    mt_kahypar_free_partitioned_hypergraph(partitioned_hg_2);
  }

  TEST(MtKaHyPar, ReportsPropertiesOfHypergraphPartition) {
    mt_kahypar_error_t error;
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
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

    mt_kahypar_hypergraph_t hypergraph = mt_kahypar_create_hypergraph(
      context, num_vertices, num_hyperedges, hyperedge_indices.get(), hyperedges.get(), nullptr, nullptr, &error);

    std::unique_ptr<mt_kahypar_partition_id_t[]> partition = std::make_unique<mt_kahypar_partition_id_t[]>(7);
    partition[0] = 0; partition[1] = 0; partition[2] = 0;
    partition[3] = 1; partition[4] = 1; partition[5] = 1; partition[6] = 1;

    mt_kahypar_partitioned_hypergraph_t partitioned_hg =
      mt_kahypar_create_partitioned_hypergraph(hypergraph, context, 2, partition.get(), &error);

    ASSERT_EQ(2, mt_kahypar_km1(partitioned_hg));
    ASSERT_EQ(2, mt_kahypar_num_blocks(partitioned_hg));
    ASSERT_EQ(3, mt_kahypar_block_weight(partitioned_hg, 0));
    ASSERT_EQ(4, mt_kahypar_block_weight(partitioned_hg, 1));
    for ( mt_kahypar_hypernode_id_t hn = 0; hn < 7; ++hn ) {
      ASSERT_EQ(partition[hn], mt_kahypar_block_id(partitioned_hg, hn));
    }

    ASSERT_EQ(1, mt_kahypar_num_incident_cut_hyperedges(partitioned_hg, 0));
    ASSERT_EQ(1, mt_kahypar_num_incident_cut_hyperedges(partitioned_hg, 6));

    ASSERT_EQ(1, mt_kahypar_connectivity(partitioned_hg, 0));
    ASSERT_EQ(2, mt_kahypar_connectivity(partitioned_hg, 1));
    ASSERT_EQ(1, mt_kahypar_connectivity(partitioned_hg, 2));
    ASSERT_EQ(2, mt_kahypar_connectivity(partitioned_hg, 3));

    ASSERT_EQ(2, mt_kahypar_num_pins_in_block(partitioned_hg, 0, 0));
    ASSERT_EQ(0, mt_kahypar_num_pins_in_block(partitioned_hg, 0, 1));
    ASSERT_EQ(2, mt_kahypar_num_pins_in_block(partitioned_hg, 1, 0));
    ASSERT_EQ(2, mt_kahypar_num_pins_in_block(partitioned_hg, 1, 1));
    ASSERT_EQ(0, mt_kahypar_num_pins_in_block(partitioned_hg, 2, 0));
    ASSERT_EQ(3, mt_kahypar_num_pins_in_block(partitioned_hg, 2, 1));
    ASSERT_EQ(1, mt_kahypar_num_pins_in_block(partitioned_hg, 3, 0));
    ASSERT_EQ(2, mt_kahypar_num_pins_in_block(partitioned_hg, 3, 1));

    mt_kahypar_free_hypergraph(hypergraph);
    mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
  }

  TEST(MtKaHyPar, WritesAndLoadsGraphPartitionFile) {
    mt_kahypar_error_t error;
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
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

    mt_kahypar_hypergraph_t graph = mt_kahypar_create_graph(
      context, num_vertices, num_hyperedges, edges.get(), nullptr, vertex_weights.get(), &error);

    std::unique_ptr<mt_kahypar_partition_id_t[]> partition =
      std::make_unique<mt_kahypar_partition_id_t[]>(7);
    partition[0] = 0; partition[1] = 0; partition[2] = 1;
    partition[3] = 1; partition[4] = 1;

    mt_kahypar_partitioned_hypergraph_t partitioned_graph =
      mt_kahypar_create_partitioned_hypergraph(graph, context, 2, partition.get(), &error);

    mt_kahypar_write_partition_to_file(partitioned_graph, "tmp.partition", &error);

    mt_kahypar_partitioned_hypergraph_t partitioned_graph_2 =
      mt_kahypar_read_partition_from_file(graph, context, 2, "tmp.partition", &error);

    std::unique_ptr<mt_kahypar_partition_id_t[]> actual_partition =
      std::make_unique<mt_kahypar_partition_id_t[]>(7);
    mt_kahypar_get_partition(partitioned_graph_2, actual_partition.get());

    ASSERT_EQ(3, mt_kahypar_cut(partitioned_graph_2));
    for ( mt_kahypar_hypernode_id_t hn = 0; hn < 5; ++hn ) {
      ASSERT_EQ(partition[hn], actual_partition[hn]);
    }

    mt_kahypar_free_hypergraph(graph);
    mt_kahypar_free_partitioned_hypergraph(partitioned_graph);
    mt_kahypar_free_partitioned_hypergraph(partitioned_graph_2);
  }

  TEST(MtKaHyPar, ReportErrorLoadingGraphPartitionFile) {
    mt_kahypar_error_t error{};
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
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

    mt_kahypar_hypergraph_t graph = mt_kahypar_create_graph(
      context, num_vertices, num_hyperedges, edges.get(), nullptr, nullptr, &error);

    std::unique_ptr<mt_kahypar_partition_id_t[]> partition =
      std::make_unique<mt_kahypar_partition_id_t[]>(7);
    partition[0] = 0; partition[1] = 0; partition[2] = 1;
    partition[3] = 1; partition[4] = 1;

    mt_kahypar_partitioned_hypergraph_t partitioned_graph =
      mt_kahypar_create_partitioned_hypergraph(graph, context, 2, partition.get(), &error);

    mt_kahypar_write_partition_to_file(partitioned_graph, "tmp.partition", &error);

    mt_kahypar_hypergraph_t graph_2 = mt_kahypar_create_graph(
      context, num_vertices - 1, num_hyperedges - 1, edges.get(), nullptr, nullptr, &error);
    ASSERT_EQ(error.status, 0);

    ASSERT_EQ(mt_kahypar_read_partition_from_file(graph_2, context, 2, "invalid_file", &error).partitioned_hg, nullptr);
    ASSERT_EQ(error.status, INVALID_INPUT);
    mt_kahypar_free_error_content(&error);

    ASSERT_EQ(mt_kahypar_read_partition_from_file(graph_2, context, 2, "tmp.partition", &error).partitioned_hg, nullptr);
    ASSERT_EQ(error.status, INVALID_INPUT);
    mt_kahypar_free_error_content(&error);

    mt_kahypar_free_hypergraph(graph);
    mt_kahypar_free_hypergraph(graph_2);
    mt_kahypar_free_partitioned_hypergraph(partitioned_graph);
  }

  class APartitioner : public Test {
    private:
      static constexpr bool debug = false;

    public:
      static constexpr char HYPERGRAPH_FILE[] = "test_instances/ibm01.hgr";
      static constexpr char GRAPH_FILE[] = "test_instances/delaunay_n15.graph";
      static constexpr char HYPERGRAPH_FIX_FILE[] = "test_instances/ibm01.k4.p1.fix";
      static constexpr char GRAPH_FIX_FILE[] = "test_instances/delaunay_n15.k4.p1.fix";
      static constexpr char TARGET_GRAPH_FILE[] = "test_instances/target.graph";

      APartitioner() :
        error(),
        context(nullptr),
        hypergraph(mt_kahypar_hypergraph_t { nullptr, NULLPTR_HYPERGRAPH }),
        partitioned_hg(mt_kahypar_partitioned_hypergraph_t { nullptr, NULLPTR_PARTITION }),
        target_graph(nullptr) {
      mt_kahypar_set_seed(42);
    }

    void Partition(const char* filename,
                   const mt_kahypar_file_format_type_t format,
                   const mt_kahypar_preset_type_t preset,
                   const mt_kahypar_partition_id_t num_blocks,
                   const double epsilon,
                   const mt_kahypar_objective_t objective,
                   const bool verbose = false,
                   const bool add_fixed_vertices = false) {
      SetUpContext(preset, num_blocks, epsilon, objective, verbose);
      Load(filename, format);
      if ( add_fixed_vertices ) addFixedVertices(num_blocks);
      partition(hypergraph, &partitioned_hg, context, num_blocks, epsilon, nullptr);
    }

    void PartitionNoSetup(const mt_kahypar_partition_id_t num_blocks,
                          const double epsilon,
                          const bool add_fixed_vertices = false) {
      if ( add_fixed_vertices ) addFixedVertices(num_blocks);
      partition(hypergraph, &partitioned_hg, context, num_blocks, epsilon, nullptr);
    }

    void Map(const char* filename,
            const mt_kahypar_file_format_type_t format,
            const mt_kahypar_preset_type_t preset,
            const double epsilon,
            const bool verbose = false) {
      SetUpContext(preset, 8, epsilon, KM1, verbose);
      Load(filename, format);
      partition(hypergraph, &partitioned_hg, context, 8, epsilon, target_graph);
    }

    void PartitionAnotherHypergraph(const char* filename,
                                    const mt_kahypar_file_format_type_t format,
                                    const mt_kahypar_preset_type_t preset,
                                    const mt_kahypar_partition_id_t num_blocks,
                                    const double epsilon,
                                    const mt_kahypar_objective_t objective,
                                    const bool verbose = false) {
      mt_kahypar_context_t* c = mt_kahypar_context_from_preset(preset);
      mt_kahypar_set_partitioning_parameters(c, num_blocks, epsilon, objective);
      mt_kahypar_set_context_parameter(c, VERBOSE, ( debug || verbose ) ? "1" : "0", &error);

      mt_kahypar_hypergraph_t hg = mt_kahypar_read_hypergraph_from_file(filename, c, format, &error);
      partition(hg, nullptr, c, num_blocks, epsilon, nullptr);

      mt_kahypar_free_context(c);
      mt_kahypar_free_hypergraph(hg);
    }

    void ImprovePartition(const mt_kahypar_preset_type_t preset,
                          const mt_kahypar_partition_id_t num_blocks,
                          const double epsilon,
                          const mt_kahypar_objective_t objective,
                          const size_t num_vcycles,
                          const bool verbose = false) {
      if ( preset != mt_kahypar_get_preset(context) ) {
        mt_kahypar_free_context(context);
        context = mt_kahypar_context_from_preset(preset);
      }
      mt_kahypar_set_partitioning_parameters(context, num_blocks, epsilon, objective);
      mt_kahypar_set_context_parameter(context, VERBOSE, ( debug || verbose ) ? "1" : "0", &error);

      mt_kahypar_hyperedge_weight_t before = mt_kahypar_km1(partitioned_hg);
      mt_kahypar_improve_partition(partitioned_hg, context, num_vcycles, &error);
      mt_kahypar_hyperedge_weight_t after = mt_kahypar_km1(partitioned_hg);
      ASSERT_LE(after, before);
    }

    void ImproveMapping(const mt_kahypar_preset_type_t preset,
                        const double epsilon,
                        const size_t num_vcycles,
                        const bool verbose = false) {
      if ( preset != mt_kahypar_get_preset(context) ) {
        mt_kahypar_free_context(context);
        context = mt_kahypar_context_from_preset(preset);
      }
      mt_kahypar_set_partitioning_parameters(context, 8, epsilon, KM1);
      mt_kahypar_set_context_parameter(context, VERBOSE, ( debug || verbose ) ? "1" : "0", &error);

      mt_kahypar_hyperedge_weight_t before = mt_kahypar_steiner_tree(partitioned_hg, target_graph);
      mt_kahypar_improve_mapping(partitioned_hg, target_graph, context, num_vcycles, &error);
      mt_kahypar_hyperedge_weight_t after = mt_kahypar_steiner_tree(partitioned_hg, target_graph);
      ASSERT_LE(after, before);
    }

    void Load(const char* filename,
              const mt_kahypar_file_format_type_t format) {
      if ( hypergraph.type != NULLPTR_HYPERGRAPH ) {
        mt_kahypar_free_hypergraph(hypergraph);
      }
      hypergraph = mt_kahypar_read_hypergraph_from_file(filename, context, format, &error);
    }

    void SetUpContext(const mt_kahypar_preset_type_t preset,
                      const mt_kahypar_partition_id_t num_blocks,
                      const double epsilon,
                      const mt_kahypar_objective_t objective,
                      const bool verbose = false) {
      mt_kahypar_free_context(context);
      context = mt_kahypar_context_from_preset(preset);
      mt_kahypar_set_partitioning_parameters(context, num_blocks, epsilon, objective);
      mt_kahypar_set_context_parameter(context, VERBOSE, ( debug || verbose ) ? "1" : "0", &error);
    }

    void verifyFixedVertexAssignment(const char* fixed_vertex_file) {
      HypernodeID num_nodes = mt_kahypar_num_hypernodes(hypergraph);
      std::vector<PartitionID> fixed_vertices;
      fixed_vertices.resize(num_nodes);
      mt_kahypar_read_fixed_vertices_from_file(fixed_vertex_file, num_nodes, fixed_vertices.data(), &error);
      vec<PartitionID> partition(num_nodes, kInvalidPartition);
      mt_kahypar_get_partition(partitioned_hg, partition.data());

      for ( HypernodeID hn = 0; hn < num_nodes; ++hn ) {
        if ( fixed_vertices[hn] != -1 ) {
          ASSERT_TRUE(mt_kahypar_is_fixed_vertex(hypergraph, hn));
          ASSERT_EQ(fixed_vertices[hn], mt_kahypar_fixed_vertex_block(hypergraph, hn));
          ASSERT_EQ(fixed_vertices[hn], partition[hn]);
        } else {
          ASSERT_FALSE(mt_kahypar_is_fixed_vertex(hypergraph, hn));
        }
      }
    }

    void SetUp()  {
      mt_kahypar_initialize(std::thread::hardware_concurrency(), false);
      target_graph = mt_kahypar_read_target_graph_from_file(TARGET_GRAPH_FILE, context, &error);
      ASSERT_NE(target_graph, nullptr);
    }

    void TearDown() {
      if (error.status != SUCCESS) {
        LOG << error.msg;
      }
      ASSERT_EQ(error.status, SUCCESS);
      mt_kahypar_free_context(context);
      mt_kahypar_free_hypergraph(hypergraph);
      mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
      mt_kahypar_free_target_graph(target_graph);
    }

    mt_kahypar_error_t error;
    mt_kahypar_context_t* context;
    mt_kahypar_hypergraph_t hypergraph;
    mt_kahypar_partitioned_hypergraph_t partitioned_hg;
    mt_kahypar_target_graph_t* target_graph;

   private:
    void addFixedVertices(const mt_kahypar_partition_id_t num_blocks ) {
      if ( hypergraph.type == STATIC_HYPERGRAPH ||
            hypergraph.type == DYNAMIC_HYPERGRAPH ) {
        mt_kahypar_add_fixed_vertices_from_file(hypergraph, HYPERGRAPH_FIX_FILE, num_blocks, &error);
      } else if ( hypergraph.type == STATIC_GRAPH ||
                  hypergraph.type == DYNAMIC_GRAPH ) {
        mt_kahypar_add_fixed_vertices_from_file(hypergraph, GRAPH_FIX_FILE, num_blocks, &error);
      }
    }

    void partition(mt_kahypar_hypergraph_t hg,
                   mt_kahypar_partitioned_hypergraph_t* phg,
                   mt_kahypar_context_t * c,
                   const mt_kahypar_partition_id_t num_blocks,
                   const double epsilon,
                   mt_kahypar_target_graph_t* target_graph) {
      mt_kahypar_partitioned_hypergraph_t p_hg { nullptr, NULLPTR_PARTITION };
      if ( target_graph ) {
        p_hg = mt_kahypar_map(hg, target_graph, c, &error);
      } else {
        p_hg = mt_kahypar_partition(hg, c, &error);
      }

      double imbalance = mt_kahypar_imbalance(p_hg, c);
      mt_kahypar_hyperedge_weight_t km1 = mt_kahypar_km1(p_hg);
      if ( debug ) {
        LOG << " imbalance =" << imbalance << "\n"
            << "cut =" << mt_kahypar_cut(p_hg) << "\n"
            << "km1 =" << km1 << "\n"
            << "soed =" << mt_kahypar_soed(p_hg) << "\n"
            << (target_graph ? "steiner_tree = " + std::to_string(mt_kahypar_steiner_tree(p_hg, target_graph)) : "");

      }
      ASSERT_LE(imbalance, epsilon);

      // Verify Partition IDs
      std::unique_ptr<mt_kahypar_partition_id_t[]> partition =
        std::make_unique<mt_kahypar_partition_id_t[]>(mt_kahypar_num_hypernodes(hg));
      mt_kahypar_get_partition(p_hg, partition.get());
      std::vector<mt_kahypar_hypernode_weight_t> expected_block_weights(num_blocks);
      for ( mt_kahypar_hypernode_id_t hn = 0; hn < mt_kahypar_num_hypernodes(hg); ++hn ) {
        ASSERT_GE(partition[hn], 0);
        ASSERT_LT(partition[hn], num_blocks);
        ++expected_block_weights[partition[hn]];
      }

      // Verify Block Weights
      std::unique_ptr<mt_kahypar_hypernode_weight_t[]> block_weights =
        std::make_unique<mt_kahypar_hypernode_weight_t[]>(num_blocks);
      mt_kahypar_get_block_weights(p_hg, block_weights.get());
      for ( mt_kahypar_partition_id_t i = 0; i < num_blocks; ++i ) {
        EXPECT_EQ(expected_block_weights[i], block_weights[i]);
      }

      if ( phg ) {
        *phg = p_hg;
      } else {
        mt_kahypar_free_partitioned_hypergraph(p_hg);
      }
    }
  };

  TEST_F(APartitioner, PartitionsAHypergraphInTwoBlocksWithDefaultPresetKm1) {
    Partition(HYPERGRAPH_FILE, HMETIS, DEFAULT, 2, 0.03, KM1, false);
  }

  TEST_F(APartitioner, PartitionsAHypergraphInTwoBlocksWithDefaultPresetSoed) {
    Partition(HYPERGRAPH_FILE, HMETIS, DEFAULT, 2, 0.03, SOED, false);
  }

  TEST_F(APartitioner, PartitionsAGraphInTwoBlocksWithDefaultPreset) {
    Partition(GRAPH_FILE, METIS, DEFAULT, 2, 0.03, CUT, false);
  }

  TEST_F(APartitioner, PartitionsAHypergraphInFourBlocksWithDefaultPresetKm1) {
    Partition(HYPERGRAPH_FILE, HMETIS, DEFAULT, 4, 0.03, KM1, false);
  }

  TEST_F(APartitioner, PartitionsAHypergraphInFourBlocksWithDefaultPresetSoed) {
    Partition(HYPERGRAPH_FILE, HMETIS, DEFAULT, 4, 0.03, SOED, false);
  }

  TEST_F(APartitioner, PartitionsAGraphInFourBlocksWithDefaultPreset) {
    Partition(GRAPH_FILE, METIS, DEFAULT, 4, 0.03, CUT, false);
  }

  TEST_F(APartitioner, PartitionsAHypergraphInTwoBlocksWithQualityPreset) {
    Partition(HYPERGRAPH_FILE, HMETIS, QUALITY, 2, 0.03, KM1, false);
  }

  TEST_F(APartitioner, PartitionsAGraphInTwoBlocksWithQualityPreset) {
    Partition(GRAPH_FILE, METIS, QUALITY, 2, 0.03, CUT, false);
  }

  TEST_F(APartitioner, PartitionsAHypergraphInFourBlocksWithQualityPreset) {
    Partition(HYPERGRAPH_FILE, HMETIS, QUALITY, 4, 0.03, KM1, false);
  }

  TEST_F(APartitioner, PartitionsAGraphInFourBlocksWithQualityPreset) {
    Partition(GRAPH_FILE, METIS, QUALITY, 4, 0.03, CUT, false);
  }

  TEST_F(APartitioner, PartitionsAHypergraphInTwoBlocksWithDeterministicPreset) {
    Partition(HYPERGRAPH_FILE, HMETIS, DETERMINISTIC, 2, 0.03, KM1, false);
  }

  TEST_F(APartitioner, PartitionsAGraphInTwoBlocksWithDeterministicPreset) {
    Partition(GRAPH_FILE, METIS, DETERMINISTIC, 2, 0.03, CUT, false);
  }

  TEST_F(APartitioner, PartitionsAHypergraphInFourBlocksWithDeterministicPreset) {
    Partition(HYPERGRAPH_FILE, HMETIS, DETERMINISTIC, 4, 0.03, KM1, false);
  }

  TEST_F(APartitioner, PartitionsAGraphInFourBlocksWithDeterministicPreset) {
    Partition(GRAPH_FILE, METIS, DETERMINISTIC, 4, 0.03, CUT, false);
  }

  TEST_F(APartitioner, PartitionsAHypergraphInTwoBlocksWithLargeKPreset) {
    Partition(HYPERGRAPH_FILE, HMETIS, LARGE_K, 2, 0.03, KM1, false);
  }

  TEST_F(APartitioner, PartitionsAGraphInTwoBlocksWithLargeKPreset) {
    Partition(GRAPH_FILE, METIS, LARGE_K, 2, 0.03, CUT, false);
  }

  TEST_F(APartitioner, PartitionsAHypergraphInFourBlocksWithLargeKPreset) {
    Partition(HYPERGRAPH_FILE, HMETIS, LARGE_K, 4, 0.03, KM1, false);
  }

  TEST_F(APartitioner, PartitionsAGraphInFourBlocksWithLargeKPreset) {
    Partition(GRAPH_FILE, METIS, LARGE_K, 4, 0.03, CUT, false);
  }

  TEST_F(APartitioner, PartitionsAHypergraphInTwoBlocksWithHighestQualityPreset) {
    Partition(HYPERGRAPH_FILE, HMETIS, HIGHEST_QUALITY, 2, 0.03, KM1, false);
  }

  TEST_F(APartitioner, PartitionsAGraphInTwoBlocksWithHighestQualityPreset) {
    Partition(GRAPH_FILE, METIS, HIGHEST_QUALITY, 2, 0.03, CUT, false);
  }

  TEST_F(APartitioner, PartitionsAHypergraphInFourBlocksWithHighestQualityPreset) {
    Partition(HYPERGRAPH_FILE, HMETIS, HIGHEST_QUALITY, 4, 0.03, KM1, false);
  }

  TEST_F(APartitioner, PartitionsAGraphInFourBlocksWithHighestQualityPreset) {
    Partition(GRAPH_FILE, METIS, HIGHEST_QUALITY, 4, 0.03, CUT, false);
  }

  TEST_F(APartitioner, CanPartitionTwoHypergraphsSimultanously) {
    tbb::parallel_invoke([&]() {
      PartitionAnotherHypergraph(HYPERGRAPH_FILE, HMETIS, DETERMINISTIC, 4, 0.03, KM1, false);
    }, [&] {
      PartitionAnotherHypergraph(GRAPH_FILE, METIS, DEFAULT, 8, 0.03, CUT, false);
    });
  }

  TEST_F(APartitioner, CanPartitionFourHypergraphsSimultanously) {
    tbb::parallel_invoke([&]() {
      PartitionAnotherHypergraph(HYPERGRAPH_FILE, HMETIS, DETERMINISTIC, 4, 0.03, KM1, false);
    }, [&] {
      PartitionAnotherHypergraph(GRAPH_FILE, METIS, DEFAULT, 8, 0.03, CUT, false);
    }, [&]() {
      PartitionAnotherHypergraph(HYPERGRAPH_FILE, HMETIS, DEFAULT, 4, 0.03, KM1, false);
    }, [&] {
      PartitionAnotherHypergraph(GRAPH_FILE, METIS, QUALITY, 4, 0.03, CUT, false);
    });
  }

  TEST_F(APartitioner, ChecksIfDeterministicPresetProducesSameResultsForHypergraphs) {
    Partition(HYPERGRAPH_FILE, HMETIS, DETERMINISTIC, 8, 0.03, KM1, false);
    const double objective_1 = mt_kahypar_km1(partitioned_hg);
    Partition(HYPERGRAPH_FILE, HMETIS, DETERMINISTIC, 8, 0.03, KM1, false);
    const double objective_2 = mt_kahypar_km1(partitioned_hg);
    Partition(HYPERGRAPH_FILE, HMETIS, DETERMINISTIC, 8, 0.03, KM1, false);
    const double objective_3 = mt_kahypar_km1(partitioned_hg);
    ASSERT_EQ(objective_1, objective_2);
    ASSERT_EQ(objective_1, objective_3);
  }

  TEST_F(APartitioner, ChecksIfDeterministicPresetProducesSameResultsForGraphs) {
    Partition(GRAPH_FILE, METIS, DETERMINISTIC, 8, 0.03, CUT, false);
    const double objective_1 = mt_kahypar_cut(partitioned_hg);
    Partition(GRAPH_FILE, METIS, DETERMINISTIC, 8, 0.03, CUT, false);
    const double objective_2 = mt_kahypar_cut(partitioned_hg);
    Partition(GRAPH_FILE, METIS, DETERMINISTIC, 8, 0.03, CUT, false);
    const double objective_3 = mt_kahypar_cut(partitioned_hg);
    ASSERT_EQ(objective_1, objective_2);
    ASSERT_EQ(objective_1, objective_3);
  }

  TEST_F(APartitioner, ImprovesHypergraphPartitionWithOneVCycle) {
    Partition(HYPERGRAPH_FILE, HMETIS, DEFAULT, 4, 0.03, KM1, false);
    ImprovePartition(DEFAULT, 4, 0.03, KM1, 1, false);
  }

  TEST_F(APartitioner, ImprovesGraphPartitionWithOneVCycle) {
    Partition(GRAPH_FILE, METIS, DEFAULT, 4, 0.03, CUT, false);
    ImprovePartition(DEFAULT, 4, 0.03, CUT, 1, false);
  }

  TEST_F(APartitioner, ImprovesHypergraphPartitionWithOneVCycleAndDifferentPresetType) {
    Partition(HYPERGRAPH_FILE, HMETIS, DEFAULT, 4, 0.03, KM1, false);
    ImprovePartition(QUALITY, 4, 0.03, KM1, 1, false);
  }

  TEST_F(APartitioner, ImprovesGraphPartitionWithOneVCycleAndDifferentPresetType) {
    Partition(GRAPH_FILE, METIS, DEFAULT, 4, 0.03, CUT, false);
    ImprovePartition(QUALITY, 4, 0.03, CUT, 1, false);
  }

  TEST_F(APartitioner, ImprovesHypergraphPartitionWithThreeVCycles) {
    Partition(HYPERGRAPH_FILE, HMETIS, DEFAULT, 4, 0.03, KM1, false);
    ImprovePartition(DEFAULT, 4, 0.03, KM1, 3, false);
  }

  TEST_F(APartitioner, ImprovesGraphPartitionWithThreeVCycles) {
    Partition(GRAPH_FILE, METIS, DEFAULT, 4, 0.03, CUT, false);
    ImprovePartition(DEFAULT, 4, 0.03, CUT, 3, false);
  }

  TEST_F(APartitioner, PartitionsHypergraphWithIndividualBlockWeightsAndVCycle) {
    // Setup Individual Block Weights
    SetUpContext(DEFAULT, 4, 0.03, KM1, false);
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> block_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(4);
    block_weights[0] = 2131; block_weights[1] = 1213;
    block_weights[2] = 7287; block_weights[3] = 2501;
    mt_kahypar_set_individual_target_block_weights(context, 4, block_weights.get());

    Load(HYPERGRAPH_FILE, HMETIS);
    PartitionNoSetup(4, 0.03);

    // Verify Block Weights
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> actual_block_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(4);
    mt_kahypar_get_block_weights(partitioned_hg, actual_block_weights.get());
    for ( mt_kahypar_partition_id_t i = 0; i < 4; ++i ) {
      ASSERT_LE(actual_block_weights[i], block_weights[i]);
    }

    ImprovePartition(DEFAULT, 4, 0.03, KM1, 1, false);

    // Verify Block Weights
    actual_block_weights = std::make_unique<mt_kahypar_hypernode_weight_t[]>(4);
    mt_kahypar_get_block_weights(partitioned_hg, actual_block_weights.get());
    for ( mt_kahypar_partition_id_t i = 0; i < 4; ++i ) {
      ASSERT_LE(actual_block_weights[i], block_weights[i]);
    }
  }

  TEST_F(APartitioner, PartitionsGraphWithIndividualBlockWeightsAndVCycle) {
    // Setup Individual Block Weights
    SetUpContext(DEFAULT, 4, 0.03, CUT, false);
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> block_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(4);
    block_weights[0] = 11201; block_weights[1] = 4384;
    block_weights[2] = 14174; block_weights[3] = 3989;
    mt_kahypar_set_individual_target_block_weights(context, 4, block_weights.get());

    Load(GRAPH_FILE, METIS);
    PartitionNoSetup(4, 0.03);

    // Verify Block Weights
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> actual_block_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(4);
    mt_kahypar_get_block_weights(partitioned_hg, actual_block_weights.get());
    for ( mt_kahypar_partition_id_t i = 0; i < 4; ++i ) {
      ASSERT_LE(actual_block_weights[i], block_weights[i]);
    }

    ImprovePartition(DEFAULT, 4, 0.03, CUT, 1, false);

    // Verify Block Weights
    actual_block_weights = std::make_unique<mt_kahypar_hypernode_weight_t[]>(4);
    mt_kahypar_get_block_weights(partitioned_hg, actual_block_weights.get());
    for ( mt_kahypar_partition_id_t i = 0; i < 4; ++i ) {
      ASSERT_LE(actual_block_weights[i], block_weights[i]);
    }
  }

  TEST_F(APartitioner, MapsAHypergraphOntoATargetGraphWithDefaultPreset) {
    Map(HYPERGRAPH_FILE, HMETIS, DEFAULT, 0.03, false);
  }

  TEST_F(APartitioner, MapsAHypergraphOntoATargetGraphWithQualityPreset) {
    Map(HYPERGRAPH_FILE, HMETIS, QUALITY, 0.03, false);
  }

  TEST_F(APartitioner, MapsAHypergraphOntoATargetGraphWithHighestQualityPreset) {
    Map(HYPERGRAPH_FILE, HMETIS, HIGHEST_QUALITY, 0.03, false);
  }

  TEST_F(APartitioner, MapsAHypergraphOntoATargetGraphWithDeterministicPreset) {
    Map(HYPERGRAPH_FILE, HMETIS, DETERMINISTIC, 0.03, false);
  }

  TEST_F(APartitioner, ChecksIfDeterministicMappingProducesSameResultsForHypergraphs) {
    // note: this test doesn't seem to be very successful at actually catching non-determinism
    Map(HYPERGRAPH_FILE, HMETIS, DETERMINISTIC, 0.03, false);
    const double objective_1 = mt_kahypar_steiner_tree(partitioned_hg, target_graph);
    Map(HYPERGRAPH_FILE, HMETIS, DETERMINISTIC, 0.03, false);
    const double objective_2 = mt_kahypar_steiner_tree(partitioned_hg, target_graph);
    Map(HYPERGRAPH_FILE, HMETIS, DETERMINISTIC, 0.03, false);
    const double objective_3 = mt_kahypar_steiner_tree(partitioned_hg, target_graph);
    ASSERT_EQ(objective_1, objective_2);
    ASSERT_EQ(objective_1, objective_3);
  }

  TEST_F(APartitioner, MapsAGraphOntoATargetGraphWithDefaultPreset) {
    Map(GRAPH_FILE, METIS, DEFAULT, 0.03, false);
  }

  TEST_F(APartitioner, MapsAGraphOntoATargetGraphWithQualityPreset) {
    Map(GRAPH_FILE, METIS, QUALITY, 0.03, false);
  }

  TEST_F(APartitioner, MapsAGraphOntoATargetGraphWithHighestQualityPreset) {
    Map(GRAPH_FILE, METIS, HIGHEST_QUALITY, 0.03, false);
  }

  TEST_F(APartitioner, MapsAGraphOntoATargetGraphWithDeterministicPreset) {
    Map(GRAPH_FILE, METIS, DETERMINISTIC, 0.03, false);
  }

  TEST_F(APartitioner, ChecksIfDeterministicMappingProducesSameResultsForGraphs) {
    // note: this test doesn't seem to be very successful at actually catching non-determinism
    Map(GRAPH_FILE, METIS, DETERMINISTIC, 0.03, false);
    const double objective_1 = mt_kahypar_steiner_tree(partitioned_hg, target_graph);
    Map(GRAPH_FILE, METIS, DETERMINISTIC, 0.03, false);
    const double objective_2 = mt_kahypar_steiner_tree(partitioned_hg, target_graph);
    Map(GRAPH_FILE, METIS, DETERMINISTIC, 0.03, false);
    const double objective_3 = mt_kahypar_steiner_tree(partitioned_hg, target_graph);
    ASSERT_EQ(objective_1, objective_2);
    ASSERT_EQ(objective_1, objective_3);
  }

  TEST_F(APartitioner, ImprovesHypergraphMappingWithOneVCycles) {
    Map(HYPERGRAPH_FILE, HMETIS, DEFAULT, 0.03, false);
    ImproveMapping(DEFAULT, 0.03, 1, false);
  }

  TEST_F(APartitioner, ImprovesGraphMappingWithOneVCycles) {
    Map(GRAPH_FILE, METIS, DEFAULT, 0.03, false);
    ImproveMapping(DEFAULT, 0.03, 1, false);
  }

  TEST_F(APartitioner, ImprovesHypergraphMappingWithOneVCyclesWithQualityPreset) {
    Map(HYPERGRAPH_FILE, HMETIS, DEFAULT, 0.03, false);
    ImproveMapping(QUALITY, 0.03, 1, false);
  }

  TEST_F(APartitioner, ImprovesHypergraphMappingGeneratedByOptimizingKm1Metric) {
    Partition(HYPERGRAPH_FILE, HMETIS, DEFAULT, 8, 0.03, KM1, false);
    ImproveMapping(DEFAULT, 0.03, 1, false);
  }

  TEST_F(APartitioner, PartitionsAHypergraphWithFixedVerticesAndDefaultPreset) {
    Partition(HYPERGRAPH_FILE, HMETIS, DEFAULT, 4, 0.03, KM1, false, true /* add fixed vertices */);
    verifyFixedVertexAssignment(HYPERGRAPH_FIX_FILE);
  }

  TEST_F(APartitioner, PartitionsAHypergraphWithFixedVerticesAndQualityPreset) {
    Partition(HYPERGRAPH_FILE, HMETIS, QUALITY, 4, 0.03, KM1, false, true /* add fixed vertices */);
    verifyFixedVertexAssignment(HYPERGRAPH_FIX_FILE);
  }

  TEST_F(APartitioner, PartitionsAHypergraphWithFixedVerticesAndHighestQualityPreset) {
    Partition(HYPERGRAPH_FILE, HMETIS, HIGHEST_QUALITY, 4, 0.03, KM1, false, true /* add fixed vertices */);
    verifyFixedVertexAssignment(HYPERGRAPH_FIX_FILE);
  }

  TEST_F(APartitioner, PartitionsAHypergraphWithFixedVerticesAndDeterministicPreset) {
    Partition(HYPERGRAPH_FILE, HMETIS, DETERMINISTIC, 4, 0.03, KM1, false, true /* add fixed vertices */);
    verifyFixedVertexAssignment(HYPERGRAPH_FIX_FILE);
  }

  TEST_F(APartitioner, PartitionsAGraphWithFixedVerticesAndDefaultPreset) {
    Partition(GRAPH_FILE, METIS, DEFAULT, 4, 0.03, CUT, false, true /* add fixed vertices */);
    verifyFixedVertexAssignment(GRAPH_FIX_FILE);
  }

  TEST_F(APartitioner, PartitionsGraphWithFixedVerticesAndQualityPreset) {
    Partition(GRAPH_FILE, METIS, QUALITY, 4, 0.03, CUT, false, true /* add fixed vertices */);
    verifyFixedVertexAssignment(GRAPH_FIX_FILE);
  }

  TEST_F(APartitioner, PartitionsAGraphWithFixedVerticesAndHighestQualityPreset) {
    Partition(GRAPH_FILE, METIS, HIGHEST_QUALITY, 4, 0.03, CUT, false, true /* add fixed vertices */);
    verifyFixedVertexAssignment(GRAPH_FIX_FILE);
  }

  TEST_F(APartitioner, PartitionsAGraphWithFixedVerticesAndDeterministicPreset) {
    Partition(GRAPH_FILE, METIS, DETERMINISTIC, 4, 0.03, CUT, false, true /* add fixed vertices */);
    verifyFixedVertexAssignment(GRAPH_FIX_FILE);
  }

  TEST_F(APartitioner, ImprovesPartitionWithFixedVertices) {
    Partition(HYPERGRAPH_FILE, HMETIS, DEFAULT, 4, 0.03, KM1, false, true /* add fixed vertices */);
    ImprovePartition(QUALITY, 4, 0.03, KM1, 1, false);
    verifyFixedVertexAssignment(HYPERGRAPH_FIX_FILE);
  }

  TEST_F(APartitioner, PartitionsAHypergraphAfterFixedVerticesHaveBeenRemoved) {
    Partition(HYPERGRAPH_FILE, HMETIS, DEFAULT, 4, 0.03, KM1, false, true /* add fixed vertices */);
    verifyFixedVertexAssignment(HYPERGRAPH_FIX_FILE);
    mt_kahypar_remove_fixed_vertices(hypergraph);
    PartitionNoSetup(4, 0.03, false);
  }

  TEST_F(APartitioner, PartitionsManyHypergraphsInParallel) {
    std::atomic<size_t> cnt(0);
    size_t max_runs = 100;
    tbb::parallel_for(0U, std::thread::hardware_concurrency(), [&](const int /*id*/) {
      while ( cnt.load(std::memory_order_relaxed) < max_runs ) {
        ++cnt;
        PartitionAnotherHypergraph("test_instances/test_instance.hgr", HMETIS, DEFAULT, 4, 0.03, KM1, false);
      }
    });
  }
}

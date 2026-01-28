/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2026 Nikolai Maas <nikolai.maas@kit.edu>
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

#include <thread>

#include "gmock/gmock.h"

#include "mt-kahypar/datastructures/fixed_vertex_support.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/partition/mapping/initial_mapping.h"
#include "mt-kahypar/partition/mapping/target_graph.h"
#include "mt-kahypar/partition/metrics.h"

using ::testing::Test;

namespace mt_kahypar {

template <typename TypeTraitsT>
struct TestConfig {
  using TypeTraits = TypeTraitsT;
};

template<typename Config>
class InitialMappingTest : public Test {

 public:
  using TypeTraits = typename Config::TypeTraits;
  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using HypergraphFactory = typename Hypergraph::Factory;

  InitialMappingTest() :
          hypergraph(),
          partitioned_hypergraph(),
          target_graph(),
          context() {
    context.partition.mode = Mode::direct;
    context.partition.epsilon = 0.05;
    context.partition.k = 8;

    context.mapping.max_steiner_tree_size = 3;
    context.mapping.strategy = OneToOneMappingStrategy::greedy_mapping;

    context.partition.preset_type = PresetType::default_preset;
    context.partition.instance_type = InstanceType::hypergraph;
    context.partition.objective = Objective::steiner_tree;
    context.partition.verbose_output = false;

    // Shared Memory
    context.shared_memory.original_num_threads = std::thread::hardware_concurrency();
    context.shared_memory.num_threads = std::thread::hardware_concurrency();
  }

  void constructFromValues(const HypernodeID num_hypernodes, const HyperedgeID num_hyperedges,
                           const vec<vec<HypernodeID>>& edge_vector, const vec<HypernodeWeight> hypernode_weight = {}) {
    hypergraph = HypergraphFactory::construct(num_hypernodes, num_hyperedges, edge_vector, nullptr, hypernode_weight.empty() ? nullptr : hypernode_weight.data());
  }

  void setup() {
    /**
     * Target Graph:
     *        1           2           4
     * 0  -------- 1  -------- 2  -------- 3
     * |           |           |           |
     * | 3         | 2         | 1         | 1
     * |      3    |      2    |      1    |
     * 4  -------- 5  -------- 6  -------- 7
    */
    vec<HyperedgeWeight> edge_weights =
      { 1, 2, 4,
        3, 2, 1, 1,
        3, 2, 1 };
    target_graph = std::make_unique<TargetGraph>(
      ds::StaticGraphFactory::construct(8, 10,
        { { 0, 1 }, { 1, 2 }, { 2, 3 },
          { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 },
          { 4, 5 }, { 5, 6 }, { 6, 7 } },
          edge_weights.data()));

    partitioned_hypergraph = PartitionedHypergraph(context.partition.k, hypergraph, parallel_tag_t());
    partitioned_hypergraph.setTargetGraph(target_graph.get());
    context.setupPartWeights(hypergraph.totalWeight());

    target_graph->precomputeDistances(context.mapping.max_steiner_tree_size);
  }

  void computeMapping() {
    InitialMapping<TypeTraits>::mapToTargetGraph(partitioned_hypergraph, *target_graph, context);
  }

  Hypergraph hypergraph;
  PartitionedHypergraph partitioned_hypergraph;
  std::unique_ptr<TargetGraph> target_graph;
  Context context;
};


typedef ::testing::Types<TestConfig<StaticHypergraphTypeTraits>
                         ENABLE_HIGHEST_QUALITY(COMMA TestConfig<DynamicHypergraphTypeTraits>) > TestConfigs;

TYPED_TEST_SUITE(InitialMappingTest, TestConfigs);


TYPED_TEST(InitialMappingTest, FindsSimpleImprovement) {
  this->constructFromValues(8, 1, { {0, 4, 5} });
  this->setup();

  for (HypernodeID hn = 0; hn < 8; ++hn) {
    this->partitioned_hypergraph.setOnlyNodePart(hn, hn);
  }
  this->partitioned_hypergraph.initializePartition();

  this->computeMapping();

  ASSERT_EQ(metrics::quality(this->partitioned_hypergraph, this->context), 2);
}


TYPED_TEST(InitialMappingTest, FindsImprovementsWithFixedVertices) {
  this->constructFromValues(8, 2, { {0, 4, 5}, {2, 3, 7} });
  this->setup();

  ds::FixedVertexSupport<typename TestFixture::Hypergraph> fixed_vertices(8, 8);
  fixed_vertices.setHypergraph(&this->hypergraph);
  fixed_vertices.fixToBlock(0, 0);
  fixed_vertices.fixToBlock(3, 3);
  this->hypergraph.addFixedVertexSupport(std::move(fixed_vertices));

  for (HypernodeID hn = 0; hn < 8; ++hn) {
    this->partitioned_hypergraph.setOnlyNodePart(hn, hn);
  }
  this->partitioned_hypergraph.initializePartition();

  this->computeMapping();

  ASSERT_EQ(this->partitioned_hypergraph.partID(0), 0);
  ASSERT_EQ(this->partitioned_hypergraph.partID(3), 3);
  ASSERT_EQ(metrics::quality(this->partitioned_hypergraph, this->context), 5);
}


TYPED_TEST(InitialMappingTest, WorksWithAllVerticesFixed) {
  this->constructFromValues(8, 2, { {0, 4, 5}, {2, 3, 7} });
  this->setup();

  ds::FixedVertexSupport<typename TestFixture::Hypergraph> fixed_vertices(8, 8);
  fixed_vertices.setHypergraph(&this->hypergraph);
  for (HypernodeID hn = 0; hn < 8; ++hn) {
    fixed_vertices.fixToBlock(hn, hn);
  }
  this->hypergraph.addFixedVertexSupport(std::move(fixed_vertices));

  for (HypernodeID hn = 0; hn < 8; ++hn) {
    this->partitioned_hypergraph.setOnlyNodePart(hn, hn);
  }
  this->partitioned_hypergraph.initializePartition();

  HyperedgeWeight initial_quality = metrics::quality(this->partitioned_hypergraph, this->context);

  this->computeMapping();

  for (HypernodeID hn = 0; hn < 8; ++hn) {
    ASSERT_EQ(this->partitioned_hypergraph.partID(hn), hn);
  }
  ASSERT_EQ(metrics::quality(this->partitioned_hypergraph, this->context), initial_quality);
}


TYPED_TEST(InitialMappingTest, FindsImprovementsWithIndividualBlockWeights) {
  this->constructFromValues(8, 2, { {0, 4, 5}, {2, 3, 7} }, { 1, 2, 3, 1, 2, 3, 1, 2 });
  this->context.partition.use_individual_part_weights = true;
  this->context.partition.max_part_weights = { 2, 4, 4, 1, 4, 4, 4, 2 };
  this->setup();

  for (HypernodeID hn = 0; hn < 8; ++hn) {
    this->partitioned_hypergraph.setOnlyNodePart(hn, hn);
  }
  this->partitioned_hypergraph.initializePartition();

  HyperedgeWeight initial_quality = metrics::quality(this->partitioned_hypergraph, this->context);

  this->computeMapping();

  for (PartitionID block = 0; block < 8; ++block) {
    ASSERT_LE(this->partitioned_hypergraph.partWeight(block), this->context.partition.max_part_weights[block]);
  }
  ASSERT_LT(metrics::quality(this->partitioned_hypergraph, this->context), initial_quality);
}

}

/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2025 Nikolai Maas <nikolai.maas@kit.edu>
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

#include <memory>
#include <functional>
#include <random>

#include "gmock/gmock.h"

#include "mt-kahypar/partition/mapping/target_graph.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/parallel/thread_management.h"
#include "mt-kahypar/partition/refinement/rebalancing/advanced_rebalancer.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/randomize.h"

using ::testing::Test;

namespace mt_kahypar {

template <typename TypeTraitsT, typename GainTypesT, bool is_steiner_tree_gain, bool use_gain_cache>
struct TestConfig {
  using TypeTraits = TypeTraitsT;
  using GainTypes = GainTypesT;
  static constexpr bool steiner_tree_gain = is_steiner_tree_gain;
  static constexpr bool with_gain_cache = use_gain_cache;
};

template<typename Config>
class RepairEmptyBlocksTest : public Test {

 public:
  using TypeTraits = typename Config::TypeTraits;
  using GainTypes = typename Config::GainTypes;
  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using HypergraphFactory = typename Hypergraph::Factory;
  using TargetGraphFactory = typename ds::StaticGraph::Factory;
  using GainCache = typename GainTypes::GainCache;
  using GainComputation = typename GainTypes::GainComputation;
  using MyRepairEmtpyBlocks = RepairEmtpyBlocks<GraphAndGainTypes<TypeTraits, GainTypes>>;

  RepairEmptyBlocksTest() :
          hypergraph(),
          partitioned_hypergraph(),
          context(),
          gain_cache(),
          gain_computation(),
          target_graph(),
          repair_empty_blocks(context, gain_cache) {
    parallel::initialize_tbb(std::thread::hardware_concurrency());
    context.partition.mode = Mode::direct;
    context.partition.epsilon = 1.5;
    context.partition.k = 16;

    context.partition.preset_type = PresetType::default_preset;
    context.partition.instance_type = InstanceType::hypergraph;
    context.partition.partition_type = PartitionedHypergraph::TYPE;
    context.partition.verbose_output = false;

    // Shared Memory
    context.shared_memory.original_num_threads = std::thread::hardware_concurrency();
    context.shared_memory.num_threads = std::thread::hardware_concurrency();

    if (Config::steiner_tree_gain) {
      context.partition.objective = Objective::steiner_tree;
      context.partition.gain_policy = Hypergraph::is_graph ? GainPolicy::steiner_tree_for_graphs : GainPolicy::steiner_tree;
    } else {
      context.partition.objective = Hypergraph::is_graph ? Objective::cut : Objective::km1;
      context.partition.gain_policy = Hypergraph::is_graph ? GainPolicy::cut_for_graphs : GainPolicy::km1;
    }

    gain_computation = std::make_unique<GainComputation>(context);
    context.partition.allow_empty_blocks = false;
  }

  void constructFromFile() {
    if constexpr ( Hypergraph::is_graph ) {
      hypergraph = io::readInputFile<Hypergraph>(
        "../tests/instances/delaunay_n10.graph", FileFormat::Metis, true);
    } else {
      hypergraph = io::readInputFile<Hypergraph>(
        "../tests/instances/contracted_ibm01.hgr", FileFormat::hMetis, true);
    }
  }

  void constructFromValues(const vec<HypernodeWeight>& hypernode_weight) {
    const HypernodeID num_hypernodes = 16;
    ASSERT(num_hypernodes == hypernode_weight.size());
    if constexpr ( Hypergraph::is_graph ) {
      hypergraph = HypergraphFactory::construct(num_hypernodes, 10,
        { {0, 1}, {2, 3}, {4, 5}, {1, 3}, {5, 7}, {11, 13}, {14, 15}, {12, 13}, {8, 9}, {9, 14} },
        nullptr, hypernode_weight.data(), true);
    } else {
      hypergraph = HypergraphFactory::construct(num_hypernodes, 4,
        { {0, 1, 2, 3, 4, 5}, {1, 3, 5, 7, 11, 13}, {15, 14, 13, 12, 11}, {8, 9, 12, 13} },
        nullptr, hypernode_weight.data(), true);
    }
  }

  void repairEmptyBlocks() {
    partitioned_hypergraph.initializePartition();
    if constexpr (Config::with_gain_cache) {
      gain_cache.reset(partitioned_hypergraph.initialNumNodes(), context.partition.k);
      gain_cache.initializeGainCache(partitioned_hypergraph);
    }

    mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(partitioned_hypergraph);
    repair_empty_blocks.repairEmptyBlocks(phg, *gain_computation, [&](const Move& m) {
      if constexpr (Config::with_gain_cache) {
        ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(gain_cache, m.node, m.from, m.to));
      } else {
        ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(m.node, m.from, m.to));
      }
    });
  }

  void setup() {
    partitioned_hypergraph = PartitionedHypergraph(context.partition.k, hypergraph, parallel_tag_t());
    context.setupPartWeights(hypergraph.totalWeight());

    if constexpr (Config::steiner_tree_gain) {
      vec<std::pair<HypernodeID, HypernodeID>> edges;
      vec<HyperedgeWeight> edges_weights;
      for (PartitionID i = 0; i < context.partition.k; ++i) {
        for (PartitionID j = 0; j < i; ++j) {
          edges.emplace_back(i, j);
          edges_weights.emplace_back(5 * i + j + 1);
        }
      }
      target_graph = std::make_unique<TargetGraph>(
        TargetGraphFactory::construct_from_graph_edges(context.partition.k, edges.size(), edges, edges_weights.data()));
      target_graph->precomputeDistances(2);
      partitioned_hypergraph.setTargetGraph(target_graph.get());
    }
  }

  Hypergraph hypergraph;
  PartitionedHypergraph partitioned_hypergraph;
  Context context;
  GainCache gain_cache;
  std::unique_ptr<GainComputation> gain_computation;
  std::unique_ptr<TargetGraph> target_graph;
  MyRepairEmtpyBlocks repair_empty_blocks;
};


typedef ::testing::Types<TestConfig<StaticHypergraphTypeTraits, Km1GainTypes, false, false>,
                         TestConfig<StaticHypergraphTypeTraits, Km1GainTypes, false, true>
                         ENABLE_STEINER_TREE(COMMA TestConfig<StaticHypergraphTypeTraits COMMA SteinerTreeGainTypes COMMA true COMMA false>)
                         ENABLE_STEINER_TREE(COMMA TestConfig<StaticHypergraphTypeTraits COMMA SteinerTreeGainTypes COMMA true COMMA true>)
                         ENABLE_GRAPHS(COMMA TestConfig<StaticGraphTypeTraits COMMA CutGainForGraphsTypes COMMA false COMMA false>) > TestConfigs;

TYPED_TEST_SUITE(RepairEmptyBlocksTest, TestConfigs);


TYPED_TEST(RepairEmptyBlocksTest, SuccessfullyRepairsEmptyBlocks) {
  this->constructFromValues({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
  this->setup();

  for (HyperedgeID hn = 0; hn < this->hypergraph.initialNumNodes(); ++hn) {
    this->partitioned_hypergraph.setOnlyNodePart(hn, hn % 8);
  }
  this->repairEmptyBlocks();

  for (PartitionID block = 0; block < this->context.partition.k; ++block) {
    ASSERT_GT(this->partitioned_hypergraph.partWeight(block), 0) << V(block);
  }
}


TYPED_TEST(RepairEmptyBlocksTest, SuccessfullyRepairsEmptyBlocksWithIndividualPartWeights) {
  this->constructFromValues({4, 3, 3, 2, 2, 1, 1, 1, 3, 4, 2, 3, 1, 2, 2, 2});
  this->context.partition.use_individual_part_weights = true;
  this->context.partition.max_part_weights = {10, 10, 8, 8, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1};
  this->context.partition.max_part_weights.resize(this->context.partition.k, 5);
  this->setup();

  for (HyperedgeID hn = 0; hn < this->hypergraph.initialNumNodes(); ++hn) {
    this->partitioned_hypergraph.setOnlyNodePart(hn, hn % 8);
  }
  this->repairEmptyBlocks();

  for (PartitionID block = 0; block < this->context.partition.k; ++block) {
    ASSERT_GT(this->partitioned_hypergraph.partWeight(block), 0) << V(block);
  }
}


TYPED_TEST(RepairEmptyBlocksTest, IsDeterministicWithIndividualPartWeights) {
  this->constructFromValues({4, 3, 3, 2, 2, 1, 1, 1, 3, 4, 2, 3, 1, 2, 2, 2});

  vec<PartitionID> partition(16, kInvalidPartition);

  for (size_t i = 0; i < 5; ++i) {
    this->context.partition.use_individual_part_weights = true;
    this->context.partition.max_part_weights = {10, 10, 8, 8, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1};
    this->setup();

    for (HyperedgeID hn = 0; hn < this->hypergraph.initialNumNodes(); ++hn) {
      this->partitioned_hypergraph.setOnlyNodePart(hn, hn % 8);
    }
    this->repairEmptyBlocks();

    for (HypernodeID hn: this->hypergraph.nodes()) {
      if (i == 0) {
        partition[hn] = this->partitioned_hypergraph.partID(hn);
      } else {
        ASSERT_EQ(this->partitioned_hypergraph.partID(hn), partition[hn]);
      }
    }
  }
}


TYPED_TEST(RepairEmptyBlocksTest, IsDeterministicOnRealInstance) {
  this->constructFromFile();

  vec<PartitionID> partition(this->hypergraph.initialNumNodes(), kInvalidPartition);

  for (size_t i = 0; i < 5; ++i) {
    this->setup();

    for (HyperedgeID hn = 0; hn < this->hypergraph.initialNumNodes(); ++hn) {
      this->partitioned_hypergraph.setOnlyNodePart(hn, hn % 8);
    }
    this->repairEmptyBlocks();

    for (HypernodeID hn: this->hypergraph.nodes()) {
      if (i == 0) {
        partition[hn] = this->partitioned_hypergraph.partID(hn);
      } else {
        ASSERT_EQ(this->partitioned_hypergraph.partID(hn), partition[hn]);
      }
    }
  }
}


TYPED_TEST(RepairEmptyBlocksTest, IsDeterministicOnRealInstanceWithIndividualPartWeights) {
  this->constructFromFile();

  vec<PartitionID> partition(this->hypergraph.initialNumNodes(), kInvalidPartition);

  for (size_t i = 0; i < 5; ++i) {
    const HypernodeWeight avg_weight = this->hypergraph.totalWeight() / 16 + 1;
    auto map_weight = [&](const double factor) { return static_cast<HypernodeWeight>(factor * avg_weight); };

    this->context.partition.use_individual_part_weights = true;
    this->context.partition.max_part_weights = {
      map_weight(2.0), map_weight(3.5), map_weight(3.0), map_weight(2.5),
      map_weight(2.0), map_weight(3.5), map_weight(3.0), map_weight(2.5),
      map_weight(0.8), map_weight(0.3), map_weight(0.1), map_weight(0.01),
      map_weight(0.8), map_weight(0.3), map_weight(0.1), map_weight(0.01)
    };
    this->setup();

    for (HyperedgeID hn = 0; hn < this->hypergraph.initialNumNodes(); ++hn) {
      this->partitioned_hypergraph.setOnlyNodePart(hn, hn % 8);
    }
    this->repairEmptyBlocks();

    for (HypernodeID hn: this->hypergraph.nodes()) {
      if (i == 0) {
        partition[hn] = this->partitioned_hypergraph.partID(hn);
      } else {
        ASSERT_EQ(this->partitioned_hypergraph.partID(hn), partition[hn]);
      }
    }
  }
}

}

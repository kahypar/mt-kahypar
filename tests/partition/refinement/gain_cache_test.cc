/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <atomic>
#include "gmock/gmock.h"

#include "tbb/parallel_for.h"
#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/static_graph_factory.h"
#include "mt-kahypar/datastructures/thread_safe_fast_reset_flag_array.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/process_mapping/process_graph.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/randomize.h"

using ::testing::Test;

namespace mt_kahypar {

template<typename TypeTraitsT, typename GainTypesT>
struct TestConfig {
  using TypeTraits = TypeTraitsT;
  using GainTypes = GainTypesT;
};

template<typename Config>
class AGainCache : public Test {
  using TypeTraits = typename Config::TypeTraits;
  using Hypergraph = typename TypeTraits::Hypergraph;
  using ParallelHyperedge = typename Hypergraph::ParallelHyperedge;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using GainTypes = typename Config::GainTypes;
  using GainCache = typename GainTypes::GainCache;
  using AttributedGains = typename GainTypes::AttributedGains;

 public:
  static constexpr PartitionID k = 8;
  static constexpr HypernodeID contraction_limit = 160;
  static constexpr size_t max_batch_size = 25;

  AGainCache() :
    hypergraph(),
    partitioned_hg(),
    process_graph(nullptr),
    gain_cache(),
    was_moved() {

    if constexpr ( Hypergraph::is_graph ) {
      hypergraph = io::readInputFile<Hypergraph>(
        "../tests/instances/delaunay_n10.graph", FileFormat::Metis, true);
    } else {
      hypergraph = io::readInputFile<Hypergraph>(
        "../tests/instances/contracted_unweighted_ibm01.hgr", FileFormat::hMetis, true);
    }
    partitioned_hg = PartitionedHypergraph(k, hypergraph, parallel_tag_t { });

    if ( GainCache::TYPE == GainPolicy::process_mapping ) {
      /**
       * Process Graph:
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
      process_graph = std::make_unique<ProcessGraph>(
        ds::StaticGraphFactory::construct(8, 10,
          { { 0, 1 }, { 1, 2 }, { 2, 3 },
            { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 },
            { 4, 5 }, { 5, 6 }, { 6, 7 } },
            edge_weights.data()));
      process_graph->precomputeDistances(3);
      partitioned_hg.setProcessGraph(process_graph.get());
    }

    was_moved.setSize(hypergraph.initialNumNodes());
  }

  void initializePartition() {
    std::vector<PartitionID> partition;
    if constexpr ( Hypergraph::is_graph ) {
      io::readPartitionFile("../tests/instances/delaunay_n10.graph.part8", partition);
    } else {
      io::readPartitionFile("../tests/instances/contracted_unweighted_ibm01.hgr.part8", partition);
    }
    partitioned_hg.doParallelForAllNodes([&](const HypernodeID& hn) {
      partitioned_hg.setOnlyNodePart(hn, partition[hn]);
    });
    partitioned_hg.initializePartition();
  }

  void moveAllNodesAtRandom() {
    utils::Randomize& rand = utils::Randomize::instance();
    was_moved.reset();
    partitioned_hg.doParallelForAllNodes([&](const HypernodeID& hn) {
      if ( rand.flipCoin(SCHED_GETCPU) ) {
        const PartitionID from = partitioned_hg.partID(hn);
        const PartitionID to = rand.getRandomInt(0, k - 1, SCHED_GETCPU);
        if ( from != to && was_moved.compare_and_set_to_true(hn) ) {
          partitioned_hg.changeNodePart(gain_cache, hn, from, to);
        }
      }
    });

    partitioned_hg.doParallelForAllNodes([&](const HypernodeID& hn) {
      if ( was_moved[hn] ) {
        gain_cache.recomputeInvalidTerms(partitioned_hg, hn);
      }
    });
  }

  void moveAllNodesOfBatchAtRandom(const Batch& batch) {
    utils::Randomize& rand = utils::Randomize::instance();
    was_moved.reset();
    auto move_node = [&](const HypernodeID hn) {
      if ( rand.flipCoin(SCHED_GETCPU) )  {
        const PartitionID from = partitioned_hg.partID(hn);
        const PartitionID to = rand.getRandomInt(0, k - 1, SCHED_GETCPU);
        if ( from != to && was_moved.compare_and_set_to_true(hn) ) {
          partitioned_hg.changeNodePart(gain_cache, hn, from, to);
        }
      }
    };

    tbb::parallel_for(UL(0), batch.size(),
      [&](const size_t i) {
        move_node(batch[i].u);
        move_node(batch[i].v);
      });

    tbb::parallel_for(UL(0), batch.size(),
      [&](const size_t i) {
        if ( was_moved[batch[i].u] ) gain_cache.recomputeInvalidTerms(partitioned_hg, batch[i].u);
        if ( was_moved[batch[i].v] ) gain_cache.recomputeInvalidTerms(partitioned_hg, batch[i].v);
      });
  }

  void simulateNLevelWithGainCacheUpdates(const bool simulate_localized_refinement) {
    if constexpr ( !Hypergraph::is_static_hypergraph ) {
      // Coarsening
      utils::Randomize& rand = utils::Randomize::instance();
      std::atomic<HypernodeID> current_num_nodes(hypergraph.initialNumNodes());
      tbb::enumerable_thread_specific<vec<HypernodeID>> local_rep;
      vec<vec<ParallelHyperedge>> parallel_hes;
      while ( current_num_nodes > contraction_limit ) {
        const HypernodeID num_nodes_before_pass = current_num_nodes.load();
        tbb::parallel_for(ID(0), current_num_nodes.load(),
          [&](const HypernodeID& hn) {
            if ( hypergraph.nodeIsEnabled(hn) && current_num_nodes.load() > contraction_limit ) {
              vec<HypernodeID> representatives = local_rep.local();
              for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
                for ( const HypernodeID& pin : hypergraph.pins(he) ) {
                  if ( hn != pin ) {
                    representatives.push_back(pin);
                  }
                }
              }
              // Choose contraction partner at random
              if ( representatives.size() > 0 ) {
                const HypernodeID rep = representatives[rand.getRandomInt(
                  0, static_cast<int>(representatives.size() - 1), SCHED_GETCPU)];
                if ( hypergraph.registerContraction(hn, rep) ) {
                  current_num_nodes -= hypergraph.contract(
                    rep, std::numeric_limits<HypernodeWeight>::max());
                }
                representatives.clear();
              }
            }
          });
        if ( current_num_nodes.load() == num_nodes_before_pass ) break;
        parallel_hes.emplace_back(
          hypergraph.removeSinglePinAndParallelHyperedges());
      }

      // Initial Partitioning
      initializePartition();
      gain_cache.initializeGainCache(partitioned_hg);

      // Uncoarsening
      VersionedBatchVector hierarchy =
        hypergraph.createBatchUncontractionHierarchy(max_batch_size);
      while ( !hierarchy.empty() ) {
        BatchVector& batches = hierarchy.back();
        while ( !batches.empty() ) {
          Batch& batch = batches.back();
          if ( batch.size() > 0 ) {
            // Uncontract batch with gain cache update
            partitioned_hg.uncontract(batch, gain_cache);
            if ( simulate_localized_refinement ) {
              moveAllNodesOfBatchAtRandom(batch);
            }
          }
          batches.pop_back();
        }

        if ( !parallel_hes.empty() ) {
          // Restore single-pin and parallel nets
          partitioned_hg.restoreSinglePinAndParallelNets(
            parallel_hes.back(), gain_cache);
          parallel_hes.pop_back();
        }
        hierarchy.pop_back();
      }
    } else {
      initializePartition();
      gain_cache.initializeGainCache(partitioned_hg);
    }
  }

  void verifyGainCacheEntries() {
    partitioned_hg.doParallelForAllNodes([&](const HypernodeID& hn) {
      const PartitionID from = partitioned_hg.partID(hn);
      ASSERT_EQ(gain_cache.penaltyTerm(hn, partitioned_hg.partID(hn)),
        gain_cache.recomputePenaltyTerm(partitioned_hg, hn));
      for ( const PartitionID to : gain_cache.adjacentBlocks(hn) ) {
        if ( from != to ) {
          EXPECT_EQ(gain_cache.benefitTerm(hn, to),
            gain_cache.recomputeBenefitTerm(partitioned_hg, hn, to))
              << V(hn) << " " << V(from) << " " << V(to);
        }
      }
    });
  }

  Gain attributedGain(const SyncronizedEdgeUpdate& sync_update) {
    return -AttributedGains::gain(sync_update);
  }

  Hypergraph hypergraph;
  PartitionedHypergraph partitioned_hg;
  std::unique_ptr<ProcessGraph> process_graph;
  GainCache gain_cache;
  ds::ThreadSafeFastResetFlagArray<> was_moved;
};

typedef ::testing::Types<TestConfig<StaticHypergraphTypeTraits, Km1GainTypes>,
                         TestConfig<StaticHypergraphTypeTraits, CutGainTypes>,
                         TestConfig<StaticHypergraphTypeTraits, SoedGainTypes>,
                         TestConfig<StaticHypergraphTypeTraits, ProcessMappingGainTypes>,
                         TestConfig<StaticGraphTypeTraits, CutGainForGraphsTypes>
                         ENABLE_N_LEVEL(COMMA TestConfig<DynamicHypergraphTypeTraits COMMA Km1GainTypes>)
                         ENABLE_N_LEVEL(COMMA TestConfig<DynamicHypergraphTypeTraits COMMA CutGainTypes>)
                         ENABLE_N_LEVEL(COMMA TestConfig<DynamicHypergraphTypeTraits COMMA SoedGainTypes>)
                         ENABLE_N_LEVEL(COMMA TestConfig<DynamicHypergraphTypeTraits COMMA ProcessMappingGainTypes>)
                         ENABLE_N_LEVEL(COMMA TestConfig<DynamicGraphTypeTraits COMMA CutGainForGraphsTypes>)
                         ENABLE_LARGE_K(COMMA TestConfig<LargeKHypergraphTypeTraits COMMA Km1GainTypes>)
                         ENABLE_LARGE_K(COMMA TestConfig<LargeKHypergraphTypeTraits COMMA CutGainTypes>)
                         ENABLE_LARGE_K(COMMA TestConfig<LargeKHypergraphTypeTraits COMMA SoedGainTypes>)
                         ENABLE_LARGE_K(COMMA TestConfig<LargeKHypergraphTypeTraits COMMA ProcessMappingGainTypes>)> TestConfigs;

TYPED_TEST_CASE(AGainCache, TestConfigs);

TYPED_TEST(AGainCache, HasCorrectInitialGains) {
  this->initializePartition();
  this->gain_cache.initializeGainCache(this->partitioned_hg);
  this->verifyGainCacheEntries();
}

TYPED_TEST(AGainCache, HasCorrectGainsAfterMovingAllNodesAtRandom) {
  this->initializePartition();
  this->gain_cache.initializeGainCache(this->partitioned_hg);
  this->moveAllNodesAtRandom();
  this->verifyGainCacheEntries();
}

TYPED_TEST(AGainCache, ComparesGainsWithAttributedGains) {
  this->initializePartition();
  this->gain_cache.initializeGainCache(this->partitioned_hg);

  utils::Randomize& rand = utils::Randomize::instance();
  Gain attributed_gain = 0;
  auto delta = [&](const SyncronizedEdgeUpdate& sync_update) {
    attributed_gain += this->attributedGain(sync_update);
  };
  vec<PartitionID> adjacent_blocks;
  for ( const HypernodeID& hn : this->partitioned_hg.nodes() ) {
    adjacent_blocks.clear();
    for ( const PartitionID to : this->gain_cache.adjacentBlocks(hn) ) {
      adjacent_blocks.push_back(to);
    }
    const PartitionID from = this->partitioned_hg.partID(hn);
    const PartitionID to = adjacent_blocks[rand.getRandomInt(0,
      static_cast<int>(adjacent_blocks.size()) - 1, SCHED_GETCPU)];
    if ( from != to ) {
      const Gain expected_gain = this->gain_cache.gain(hn, from, to);
      this->partitioned_hg.changeNodePart(this->gain_cache, hn, from, to,
        std::numeric_limits<HyperedgeWeight>::max(), []{}, delta);
      ASSERT_EQ(expected_gain, attributed_gain);
    }
    attributed_gain = 0;
  }
}

#ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES

TYPED_TEST(AGainCache, HasCorrectGainsAfterNLevelUncontraction) {
  this->simulateNLevelWithGainCacheUpdates(false);
  this->verifyGainCacheEntries();
}

TYPED_TEST(AGainCache, HasCorrectGainsAfterNLevelUncontractionWithLocalizedRefinement) {
  this->simulateNLevelWithGainCacheUpdates(true);
  this->verifyGainCacheEntries();
}

#endif

}
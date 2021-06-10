/*******************************************************************************
 * This file is part of KaHyPar.
 *
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

#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/refinement/advanced/problem_construction.h"
#include "tests/partition/refinement/advanced_refiner_mock.h"

using ::testing::Test;

namespace mt_kahypar {

class AProblemConstruction : public Test {
 public:
  AProblemConstruction() :
    hg(),
    phg(),
    context(),
    max_part_weights(8, std::numeric_limits<HypernodeWeight>::max()) {

    context.partition.graph_filename = "../tests/instances/ibm01.hgr";
    context.partition.k = 8;
    context.partition.epsilon = 0.03;
    context.partition.mode = kahypar::Mode::direct_kway;
    context.partition.objective = kahypar::Objective::km1;
    context.shared_memory.num_threads = std::thread::hardware_concurrency();
    context.refinement.advanced.algorithm = AdvancedRefinementAlgorithm::mock;
    context.refinement.advanced.num_threads_per_search = 1;
    context.refinement.advanced.num_cut_edges_per_block_pair = 50;
    context.refinement.advanced.max_bfs_distance = 2;

    // Read hypergraph
    hg = io::readHypergraphFile(context.partition.graph_filename);
    phg = PartitionedHypergraph(
      context.partition.k, hg, parallel_tag_t());
    context.setupPartWeights(hg.totalWeight());

    // Read Partition
    std::vector<PartitionID> partition;
    io::readPartitionFile("../tests/instances/ibm01.hgr.part8", partition);
    phg.doParallelForAllNodes([&](const HypernodeID& hn) {
      phg.setOnlyNodePart(hn, partition[hn]);
    });
    phg.initializePartition();

    AdvancedRefinerMockControl::instance().reset();
    AdvancedRefinerMockControl::instance().max_prob_size_func = [&](ProblemStats& stats) {
      bool limit_reached = true;
      for ( const PartitionID block : stats.containedBlocks() ) {
        bool block_limit_reached = stats.nodeWeightOfBlock(block) >= max_part_weights[block];
        if ( block_limit_reached ) stats.lockBlock(block);
        limit_reached &= block_limit_reached;
      }
      return limit_reached;
    };
  }

  void verifyThatPartWeightsAreLessEqualToMaxPartWeight(const vec<HypernodeID> nodes,
                                                        const SearchID search_id,
                                                        const QuotientGraph& qg) {
    vec<HypernodeWeight> part_weights(context.partition.k, 0);
    for ( const HypernodeID& hn : nodes ) {
      part_weights[phg.partID(hn)] += phg.nodeWeight(hn);
    }

    vec<bool> used_blocks(context.partition.k, false);
      for ( const BlockPair& blocks : qg.getBlockPairs(search_id) ) {
      used_blocks[blocks.i] = true;
      used_blocks[blocks.j] = true;
    }

    for ( PartitionID i = 0; i < context.partition.k; ++i ) {
      if ( used_blocks[i] ) {
        ASSERT_LE(part_weights[i], max_part_weights[i]);
      } else {
        ASSERT_EQ(0, part_weights[i]);
      }
    }
  }

  Hypergraph hg;
  PartitionedHypergraph phg;
  Context context;
  vec<HypernodeWeight> max_part_weights;
};

template <class F, class K>
void executeConcurrent(F f1, K f2) {
  std::atomic<int> cnt(0);

  tbb::parallel_invoke([&] {
    cnt++;
    while (cnt < 2) { }
    f1();
  }, [&] {
    cnt++;
    while (cnt < 2) { }
    f2();
  });
}

void verifyThatVertexSetAreDisjoint(const vec<HypernodeID>& nodes_1, const vec<HypernodeID>& nodes_2) {
  std::set<HypernodeID> nodes;
  for ( const HypernodeID& hn : nodes_1 ) {
    nodes.insert(hn);
  }
  for ( const HypernodeID& hn : nodes_2 ) {
    ASSERT_TRUE(nodes.find(hn) == nodes.end());
  }
}

TEST_F(AProblemConstruction, GrowAnAdvancedRefinementProblemAroundTwoBlocks1) {
  ProblemConstruction constructor(hg, context);
  AdvancedRefinerAdapter refiner(hg, context);
  QuotientGraph qg(hg, context);
  qg.initialize(phg);

  max_part_weights.assign(context.partition.k, 400);
  max_part_weights[2] = 300;
  SearchID search_id = qg.requestNewSearch(refiner);
  vec<HypernodeID> nodes = constructor.construct(
    search_id, qg, refiner, phg);

  verifyThatPartWeightsAreLessEqualToMaxPartWeight(nodes, search_id, qg);
}

TEST_F(AProblemConstruction, GrowAnAdvancedRefinementProblemAroundTwoBlocks2) {
  ProblemConstruction constructor(hg, context);
  AdvancedRefinerAdapter refiner(hg, context);
  QuotientGraph qg(hg, context);
  qg.initialize(phg);

  max_part_weights.assign(context.partition.k, 800);
  max_part_weights[2] = 500;
  SearchID search_id = qg.requestNewSearch(refiner);
  vec<HypernodeID> nodes = constructor.construct(
    search_id, qg, refiner, phg);

  verifyThatPartWeightsAreLessEqualToMaxPartWeight(nodes, search_id, qg);
}

TEST_F(AProblemConstruction, GrowTwoAdvancedRefinementProblemAroundTwoBlocksSimultanously) {
  ProblemConstruction constructor(hg, context);
  AdvancedRefinerAdapter refiner(hg, context);
  QuotientGraph qg(hg, context);
  qg.initialize(phg);

  max_part_weights.assign(context.partition.k, 400);

  vec<HypernodeID> nodes_1;
  vec<HypernodeID> nodes_2;
  executeConcurrent([&] {
    SearchID search_id = qg.requestNewSearch(refiner);
     nodes_1 = constructor.construct(
      search_id, qg, refiner, phg);
    verifyThatPartWeightsAreLessEqualToMaxPartWeight(nodes_1, search_id, qg);
  }, [&] {
    SearchID search_id = qg.requestNewSearch(refiner);
    nodes_2 = constructor.construct(
      search_id, qg, refiner, phg);
    verifyThatPartWeightsAreLessEqualToMaxPartWeight(nodes_2, search_id, qg);
  });
  verifyThatVertexSetAreDisjoint(nodes_1, nodes_2);
}

TEST_F(AProblemConstruction, GrowAnAdvancedRefinementProblemAroundFourBlocks1) {
  AdvancedRefinerMockControl::instance().max_num_blocks = 4;
  ProblemConstruction constructor(hg, context);
  AdvancedRefinerAdapter refiner(hg, context);
  QuotientGraph qg(hg, context);
  qg.initialize(phg);

  max_part_weights.assign(context.partition.k, 800);
  max_part_weights[2] = 500;
  SearchID search_id = qg.requestNewSearch(refiner);
  vec<HypernodeID> nodes = constructor.construct(
    search_id, qg, refiner, phg);

  verifyThatPartWeightsAreLessEqualToMaxPartWeight(nodes, search_id, qg);
}

TEST_F(AProblemConstruction, GrowAnAdvancedRefinementProblemAroundFourBlocks2) {
  AdvancedRefinerMockControl::instance().max_num_blocks = 4;
  ProblemConstruction constructor(hg, context);
  AdvancedRefinerAdapter refiner(hg, context);
  QuotientGraph qg(hg, context);
  qg.initialize(phg);

  max_part_weights.assign(context.partition.k, 800);
  max_part_weights[2] = 500;
  max_part_weights[6] = 300;
  SearchID search_id = qg.requestNewSearch(refiner);
  vec<HypernodeID> nodes = constructor.construct(
    search_id, qg, refiner, phg);

  verifyThatPartWeightsAreLessEqualToMaxPartWeight(nodes, search_id, qg);
}

TEST_F(AProblemConstruction, GrowTwoAdvancedRefinementProblemAroundFourBlocksSimultanously) {
  AdvancedRefinerMockControl::instance().max_num_blocks = 4;
  ProblemConstruction constructor(hg, context);
  AdvancedRefinerAdapter refiner(hg, context);
  QuotientGraph qg(hg, context);
  qg.initialize(phg);

  max_part_weights.assign(context.partition.k, 500);
  vec<HypernodeID> nodes_1;
  vec<HypernodeID> nodes_2;
  executeConcurrent([&] {
    SearchID search_id = qg.requestNewSearch(refiner);
    nodes_1 = constructor.construct(
      search_id, qg, refiner, phg);
    verifyThatPartWeightsAreLessEqualToMaxPartWeight(nodes_1, search_id, qg);
  }, [&] {
    SearchID search_id = qg.requestNewSearch(refiner);
    nodes_2 = constructor.construct(
      search_id, qg, refiner, phg);
    verifyThatPartWeightsAreLessEqualToMaxPartWeight(nodes_2, search_id, qg);
  });
  verifyThatVertexSetAreDisjoint(nodes_1, nodes_2);
}



}
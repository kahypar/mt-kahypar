/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "gmock/gmock.h"

#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/refinement/flows/problem_construction.h"
#include "tests/partition/refinement/flow_refiner_mock.h"

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
    context.partition.mode = Mode::direct;
    context.partition.objective = kahypar::Objective::km1;
    context.shared_memory.num_threads = std::thread::hardware_concurrency();
    context.refinement.flows.algorithm = FlowAlgorithm::mock;
    context.refinement.flows.max_bfs_distance = 2;

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

    FlowRefinerMockControl::instance().reset();
  }

  void verifyThatPartWeightsAreLessEqualToMaxPartWeight(const Subhypergraph& sub_hg,
                                                        const SearchID search_id,
                                                        const QuotientGraph& qg) {
    vec<HypernodeWeight> part_weights(context.partition.k, 0);
    for ( const HypernodeID& hn : sub_hg.nodes_of_block_0 ) {
      part_weights[phg.partID(hn)] += phg.nodeWeight(hn);
    }
    for ( const HypernodeID& hn : sub_hg.nodes_of_block_1 ) {
      part_weights[phg.partID(hn)] += phg.nodeWeight(hn);
    }

    vec<bool> used_blocks(context.partition.k, false);
    const BlockPair blocks = qg.getBlockPair(search_id);
    used_blocks[blocks.i] = true;
    used_blocks[blocks.j] = true;
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

void verifyThatVertexSetAreDisjoint(const Subhypergraph& sub_hg_1, const Subhypergraph& sub_hg_2) {
  std::set<HypernodeID> nodes;
  for ( const HypernodeID& hn : sub_hg_1.nodes_of_block_0 ) {
    nodes.insert(hn);
  }
  for ( const HypernodeID& hn : sub_hg_1.nodes_of_block_1 ) {
    nodes.insert(hn);
  }
  for ( const HypernodeID& hn : sub_hg_2.nodes_of_block_0 ) {
    ASSERT_TRUE(nodes.find(hn) == nodes.end());
  }
  for ( const HypernodeID& hn : sub_hg_2.nodes_of_block_1 ) {
    ASSERT_TRUE(nodes.find(hn) == nodes.end());
  }
}

TEST_F(AProblemConstruction, GrowAnFlowProblemAroundTwoBlocks1) {
  ProblemConstruction constructor(hg, context);
  FlowRefinerAdapter refiner(hg, context);
  QuotientGraph qg(hg, context);
  refiner.initialize(context.shared_memory.num_threads);
  qg.initialize(phg);

  max_part_weights.assign(context.partition.k, 400);
  max_part_weights[2] = 300;
  SearchID search_id = qg.requestNewSearch(refiner);
  Subhypergraph sub_hg = constructor.construct(search_id, qg, phg);

  verifyThatPartWeightsAreLessEqualToMaxPartWeight(sub_hg, search_id, qg);
}

TEST_F(AProblemConstruction, GrowAnFlowProblemAroundTwoBlocks2) {
  ProblemConstruction constructor(hg, context);
  FlowRefinerAdapter refiner(hg, context);
  QuotientGraph qg(hg, context);
  refiner.initialize(context.shared_memory.num_threads);
  qg.initialize(phg);

  max_part_weights.assign(context.partition.k, 800);
  max_part_weights[2] = 500;
  SearchID search_id = qg.requestNewSearch(refiner);
  Subhypergraph sub_hg = constructor.construct(search_id, qg, phg);

  verifyThatPartWeightsAreLessEqualToMaxPartWeight(sub_hg, search_id, qg);
}

TEST_F(AProblemConstruction, GrowTwoFlowProblemAroundTwoBlocksSimultanously) {
  ProblemConstruction constructor(hg, context);
  FlowRefinerAdapter refiner(hg, context);
  QuotientGraph qg(hg, context);
  refiner.initialize(context.shared_memory.num_threads);
  qg.initialize(phg);

  max_part_weights.assign(context.partition.k, 400);

  Subhypergraph sub_hg_1;
  Subhypergraph sub_hg_2;
  executeConcurrent([&] {
    SearchID search_id = qg.requestNewSearch(refiner);
     sub_hg_1 = constructor.construct(search_id, qg, phg);
    verifyThatPartWeightsAreLessEqualToMaxPartWeight(sub_hg_1, search_id, qg);
  }, [&] {
    SearchID search_id = qg.requestNewSearch(refiner);
    sub_hg_2 = constructor.construct(search_id, qg, phg);
    verifyThatPartWeightsAreLessEqualToMaxPartWeight(sub_hg_2, search_id, qg);
  });
  verifyThatVertexSetAreDisjoint(sub_hg_1, sub_hg_2);
}

}
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
#include "mt-kahypar/partition/refinement/advanced/quotient_graph.h"
#include "tests/partition/refinement/advanced_refiner_mock.h"

using ::testing::Test;

#define MOVE(HN, FROM, TO) Move { FROM, TO, HN, 0 }

namespace mt_kahypar {

class AQuotientGraph : public Test {
 public:
  AQuotientGraph() :
    hg(),
    phg(),
    context() {

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
    hg = io::readHypergraphFile(
      context.partition.graph_filename, TBBNumaArena::GLOBAL_TASK_GROUP);
    phg = PartitionedHypergraph(
      context.partition.k, TBBNumaArena::GLOBAL_TASK_GROUP, hg);
    context.setupPartWeights(hg.totalWeight());

    // Read Partition
    std::vector<PartitionID> partition;
    io::readPartitionFile("../tests/instances/ibm01.hgr.part8", partition);
    phg.doParallelForAllNodes([&](const HypernodeID& hn) {
      phg.setOnlyNodePart(hn, partition[hn]);
    });
    phg.initializePartition(TBBNumaArena::GLOBAL_TASK_GROUP);
  }

  Hypergraph hg;
  PartitionedHypergraph phg;
  Context context;
};

TEST_F(AQuotientGraph, SimulatesBlockScheduling) {
  AdvancedRefinerAdapter refiner(hg, context, TBBNumaArena::GLOBAL_TASK_GROUP);
  QuotientGraph qg(context);
  qg.initialize(phg);
  const bool debug = false;

  vec<vec<CAtomic<HyperedgeWeight>>> cut_he_weights(
    context.partition.k, vec<CAtomic<HyperedgeWeight>>(
      context.partition.k, CAtomic<HyperedgeWeight>(0)));
  for ( PartitionID i = 0; i < context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < context.partition.k; ++j ) {
      cut_he_weights[i][j] += qg.getCutHyperedgeWeightOfBlockPair(i, j);
    }
  }

  tbb::parallel_for(0U, std::thread::hardware_concurrency(), [&](const unsigned int) {
    while ( !qg.terminate() ) {
      SearchID search_id = qg.requestNewSearch(refiner);
      if ( search_id != QuotientGraph::INVALID_SEARCH_ID ) {
        vec<BlockPairCutHyperedges> block_pair_cut_hes = qg.requestCutHyperedges(search_id, 10);
        ASSERT(qg.getBlockPairs(search_id).size() == 1);
        const PartitionID i = qg.getBlockPairs(search_id)[0].i;
        const PartitionID j = qg.getBlockPairs(search_id)[0].j;
        size_t num_edges = 0;
        for ( const HyperedgeID& he : block_pair_cut_hes[0].cut_hes ) {
          cut_he_weights[i][j] -= phg.edgeWeight(he);
          ++num_edges;
        }
        ASSERT_LE(num_edges, 10);
        qg.finalizeConstruction(search_id);
        qg.finalizeSearch(search_id, false);
        refiner.finalizeSearch(search_id);

        if ( debug ) {
          LOG << "Thread" << sched_getcpu() << "executes search on block pair (" << i << "," << j << ")"
              << "with" << block_pair_cut_hes[0].cut_hes.size() << "cut hyperedges ( Search ID:" << search_id << ")";
        }
      }
    }
  });

  // Each edge should be scheduled once
  for ( PartitionID i = 0; i < context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < context.partition.k; ++j ) {
      ASSERT_EQ(0, cut_he_weights[i][j]);
    }
  }
}

TEST_F(AQuotientGraph, SimulatesBlockSchedulingWithSuccessfulSearches) {
  AdvancedRefinerAdapter refiner(hg, context, TBBNumaArena::GLOBAL_TASK_GROUP);
  QuotientGraph qg(context);
  qg.initialize(phg);
  const bool debug = false;

  vec<vec<CAtomic<HyperedgeWeight>>> cut_he_weights(
    context.partition.k, vec<CAtomic<HyperedgeWeight>>(
      context.partition.k, CAtomic<HyperedgeWeight>(0)));
  for ( PartitionID i = 0; i < context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < context.partition.k; ++j ) {
      cut_he_weights[i][j] += qg.getCutHyperedgeWeightOfBlockPair(i, j);
    }
  }

  tbb::parallel_for(0U, std::thread::hardware_concurrency(), [&](const unsigned int cpu_id) {
    while ( !qg.terminate() ) {
      SearchID search_id = qg.requestNewSearch(refiner);
      if ( search_id != QuotientGraph::INVALID_SEARCH_ID ) {
        vec<BlockPairCutHyperedges> block_pair_cut_hes = qg.requestCutHyperedges(search_id, 10);
        ASSERT(qg.getBlockPairs(search_id).size() == 1);
        const PartitionID i = qg.getBlockPairs(search_id)[0].i;
        const PartitionID j = qg.getBlockPairs(search_id)[0].j;
        HyperedgeWeight cut_he_weight = 0;
        size_t num_edges = 0;
        for ( const HyperedgeID& he : block_pair_cut_hes[0].cut_hes ) {
          cut_he_weight += phg.edgeWeight(he);
          ++num_edges;
        }
        ASSERT_LE(num_edges, 10);
        cut_he_weights[i][j] -= cut_he_weight;
        qg.finalizeConstruction(search_id);

        bool success = utils::Randomize::instance().flipCoin(cpu_id);
        if ( success ) {
          cut_he_weights[i][j] += cut_he_weight;
        }
        qg.finalizeSearch(search_id, success);
        refiner.finalizeSearch(search_id);

        if ( debug ) {
          LOG << "Thread" << sched_getcpu() << "executes search on block pair (" << i << "," << j << ")"
              << "with" << block_pair_cut_hes[0].cut_hes.size() << "cut hyperedges ( Search ID:" << search_id << ", Success:"
              << std::boolalpha << success << ")";
        }
      }
    }
  });

  // Each edge should be scheduled once
  for ( PartitionID i = 0; i < context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < context.partition.k; ++j ) {
      ASSERT_EQ(0, cut_he_weights[i][j]);
    }
  }
}

TEST_F(AQuotientGraph, SimulatesBlockSchedulingWithSearchesThatRequestFourBlocks) {
  AdvancedRefinerMockControl::instance().max_num_blocks = 4;
  AdvancedRefinerAdapter refiner(hg, context, TBBNumaArena::GLOBAL_TASK_GROUP);
  QuotientGraph qg(context);
  qg.initialize(phg);

  vec<vec<CAtomic<HyperedgeWeight>>> cut_he_weights(
    context.partition.k, vec<CAtomic<HyperedgeWeight>>(
      context.partition.k, CAtomic<HyperedgeWeight>(0)));
  for ( PartitionID i = 0; i < context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < context.partition.k; ++j ) {
      cut_he_weights[i][j] += qg.getCutHyperedgeWeightOfBlockPair(i, j);
    }
  }

  tbb::parallel_for(0U, std::thread::hardware_concurrency(), [&](const unsigned int) {
    while ( !qg.terminate() ) {
      SearchID search_id = qg.requestNewSearch(refiner);
      if ( search_id != QuotientGraph::INVALID_SEARCH_ID ) {
        vec<BlockPairCutHyperedges> block_pair_cut_hes = qg.requestCutHyperedges(search_id, 10);
        size_t num_edges = 0;
        for ( const BlockPairCutHyperedges& bpch : block_pair_cut_hes ) {
          const BlockPair& blocks = bpch.blocks;
          for ( const HyperedgeID& he : bpch.cut_hes ) {
            cut_he_weights[blocks.i][blocks.j] -= phg.edgeWeight(he);
            ++num_edges;
          }
        }
        ASSERT_LE(num_edges, 10);
        qg.finalizeConstruction(search_id);
        qg.finalizeSearch(search_id, false);
        refiner.finalizeSearch(search_id);
      }
    }
  });

  // Each edge should be scheduled once
  for ( PartitionID i = 0; i < context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < context.partition.k; ++j ) {
      ASSERT_EQ(0, cut_he_weights[i][j]);
    }
  }
}

TEST_F(AQuotientGraph, SimulatesBlockSchedulingWithSearchesThatRequestFourBlocksWithSuccessfullSearches) {
  AdvancedRefinerMockControl::instance().max_num_blocks = 4;
  AdvancedRefinerAdapter refiner(hg, context, TBBNumaArena::GLOBAL_TASK_GROUP);
  QuotientGraph qg(context);
  qg.initialize(phg);

  vec<vec<CAtomic<HyperedgeWeight>>> cut_he_weights(
    context.partition.k, vec<CAtomic<HyperedgeWeight>>(
      context.partition.k, CAtomic<HyperedgeWeight>(0)));
  for ( PartitionID i = 0; i < context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < context.partition.k; ++j ) {
      cut_he_weights[i][j] += qg.getCutHyperedgeWeightOfBlockPair(i, j);
    }
  }

  tbb::parallel_for(0U, std::thread::hardware_concurrency(), [&](const unsigned int cpu_id) {
    while ( !qg.terminate() ) {
      SearchID search_id = qg.requestNewSearch(refiner);
      if ( search_id != QuotientGraph::INVALID_SEARCH_ID ) {
        vec<BlockPairCutHyperedges> block_pair_cut_hes = qg.requestCutHyperedges(search_id, 10);
        vec<HyperedgeWeight> cut_he_weight(block_pair_cut_hes.size(), 0);
        size_t num_edges = 0;
        for ( size_t i = 0; i < block_pair_cut_hes.size(); ++i ) {
          for ( const HyperedgeID& he : block_pair_cut_hes[i].cut_hes ) {
            cut_he_weight[i] += phg.edgeWeight(he);
            ++num_edges;
          }
        }
        ASSERT_LE(num_edges, 10);
        qg.finalizeConstruction(search_id);

        bool success = utils::Randomize::instance().flipCoin(cpu_id);
        if ( !success ) {
          for ( size_t i = 0; i < block_pair_cut_hes.size(); ++i ) {
            const BlockPair& blocks = block_pair_cut_hes[i].blocks;
            cut_he_weights[blocks.i][blocks.j] -= cut_he_weight[i];
          }
        }
        qg.finalizeSearch(search_id, success);
        refiner.finalizeSearch(search_id);
      }
    }
  });

  // Each edge should be scheduled once
  for ( PartitionID i = 0; i < context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < context.partition.k; ++j ) {
      ASSERT_EQ(0, cut_he_weights[i][j]);
    }
  }
}

}
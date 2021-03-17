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

TEST_F(AQuotientGraph, TEST) {
  QuotientGraph qg(context);
  qg.initialize(phg);

  SearchID search_1 = qg.requestNewSearch();
  SearchID search_2 = qg.requestNewSearch();
  SearchID search_3 = qg.requestNewSearch();
  SearchID search_4 = qg.requestNewSearch();
  LOG << V(qg.getBlockPair(search_1).i) << V(qg.getBlockPair(search_1).j);
  LOG << V(qg.getBlockPair(search_2).i) << V(qg.getBlockPair(search_2).j);
  LOG << V(qg.getBlockPair(search_3).i) << V(qg.getBlockPair(search_3).j);
  LOG << V(qg.getBlockPair(search_4).i) << V(qg.getBlockPair(search_4).j);

  vec<HyperedgeID> cut_hes = qg.requestCutHyperedges(search_1, 50);
  LOG << V(cut_hes.size());
  for ( const HyperedgeID he : cut_hes ) {
    std::cout << he << ", ";
  }
  std::cout << std::endl;

  qg.finalizeConstruction(search_1);
  qg.finalizeSearch(search_1, true);
}

}
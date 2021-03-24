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
#include "mt-kahypar/partition/refinement/advanced/advanced_refinement_problem_construction.h"
#include "tests/partition/refinement/advanced_refiner_mock.h"

using ::testing::Test;

namespace mt_kahypar {

class AAdvancedRefinementProblemConstruction : public Test {
 public:
  AAdvancedRefinementProblemConstruction() :
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

TEST_F(AAdvancedRefinementProblemConstruction, TEST) {
  // AdvancedRefinerMockControl::instance().max_num_nodes = 200;
  AdvancedRefinementProblemConstruction constructor(hg, context);
  AdvancedRefinerAdapter refiner(hg, context, TBBNumaArena::GLOBAL_TASK_GROUP);
  QuotientGraph qg(context);
  qg.initialize(phg);

  SearchID search_id = qg.requestNewSearch(refiner);

  vec<HypernodeID> nodes = constructor.construct(
    search_id, qg, refiner, phg);
  qg.finalizeConstruction(search_id);

  // LOG << V(hg.initialNumNodes()) << V(phg.partWeight(2)) << V(phg.partWeight(4)) << V(nodes.size());

  qg.finalizeSearch(search_id, false);
  refiner.finalizeSearch(search_id);
}


}
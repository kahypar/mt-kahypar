/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "tests/partition/coarsening/coarsener_fixtures.h"

using ::testing::Test;

namespace mt_kahypar {
using TestTypeTraits = ds::TestTypeTraits<2>;
using PartitionedHyperGraph = typename TestTypeTraits::PartitionedHyperGraph;
using Coarsener = MultilevelCoarsenerT<TestTypeTraits, HeavyEdgeScore,
                                       NoWeightPenalty, BestRatingWithoutTieBreaking>;

TEST_F(ACoarsener, DecreasesNumberOfPins) {
  context.coarsening.contraction_limit = 4;
  Coarsener coarsener(hypergraph, context, TBB::GLOBAL_TASK_GROUP);
  decreasesNumberOfPins(coarsener, 6);
}

TEST_F(ACoarsener, DecreasesNumberOfHyperedges) {
  context.coarsening.contraction_limit = 4;
  Coarsener coarsener(hypergraph, context, TBB::GLOBAL_TASK_GROUP);
  decreasesNumberOfHyperedges(coarsener, 3);
}

TEST_F(ACoarsener, RemovesHyperedgesOfSizeOneDuringCoarsening) {
  context.coarsening.contraction_limit = 4;
  Coarsener coarsener(hypergraph, context, TBB::GLOBAL_TASK_GROUP);
  doCoarsening(coarsener);
  auto& hypergraph = coarsener.coarsestHypergraph();
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_GE(hypergraph.edgeSize(he), 2);
  }
}

TEST_F(ACoarsener, RemovesParallelHyperedgesDuringCoarsening) {
  context.coarsening.contraction_limit = 4;
  Coarsener coarsener(hypergraph, context, TBB::GLOBAL_TASK_GROUP);
  doCoarsening(coarsener);
  auto& hypergraph = coarsener.coarsestHypergraph();
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(hypergraph.edgeWeight(he), 2);
  }
}

TEST_F(ACoarsener, ProjectsPartitionBackToOriginalHypergraph) {
  context.coarsening.contraction_limit = 4;
  Coarsener coarsener(hypergraph, context, TBB::GLOBAL_TASK_GROUP);
  doCoarsening(coarsener);
  PartitionedHyperGraph& coarsest_partitioned_hypergraph =
    coarsener.coarsestPartitionedHypergraph();
  assignPartitionIDs(coarsest_partitioned_hypergraph);
  PartitionedHyperGraph partitioned_hypergraph = coarsener.uncoarsen(nullptr_refiner);
  for ( const HypernodeID& hn : partitioned_hypergraph.nodes() ) {
    PartitionID part_id = common::get_numa_node_of_vertex(hn);
    ASSERT_EQ(part_id, partitioned_hypergraph.partID(hn));
  }
}

}  // namespace mt_kahypar

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

#include "mt-kahypar/partition/coarsening/multilevel_coarsener.h"
#include "tests/partition/coarsening/coarsener_fixtures.h"

using ::testing::Test;

namespace mt_kahypar {
using Coarsener = MultilevelCoarsenerT<ds::TestTypeTraits<2>, HeavyEdgeScore,
                                       NoWeightPenalty, BestRatingWithoutTieBreaking>;

TEST_F(ACoarsener, DecreasesNumberOfPins) {
  context.coarsening.contraction_limit = 7;
  Coarsener coarsener(hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  decreasesNumberOfPins(coarsener, 12);
}

TEST_F(ACoarsener, DecreasesNumberOfHyperedges) {
  context.coarsening.contraction_limit = 7;
  Coarsener coarsener(hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  decreasesNumberOfHyperedges(coarsener, 6);
}

TEST_F(ACoarsener, RemovesHyperedgesOfSizeOneDuringCoarsening) {
  context.coarsening.contraction_limit = 7;
  Coarsener coarsener(hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  doCoarsening(coarsener);
  auto& hypergraph = coarsener.coarsestHypergraph();
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_GE(hypergraph.edgeSize(he), 2);
  }
}

TEST_F(ACoarsener, RemovesParallelHyperedgesDuringCoarsening) {
  context.coarsening.contraction_limit = 7;
  Coarsener coarsener(hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  doCoarsening(coarsener);
  auto& hypergraph = coarsener.coarsestHypergraph();
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(hypergraph.edgeWeight(he), 2);
  }
}

TEST_F(ACoarsener, ProjectsPartitionBackToOriginalHypergraph) {
  context.coarsening.contraction_limit = 7;
  Coarsener coarsener(hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  doCoarsening(coarsener);
  assignPartitionIDs(coarsener.coarsestHypergraph());
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_EQ(-1, hypergraph.partID(hn));
  }
  coarsener.uncoarsen(nullptr_refiner);
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    PartitionID part_id = StreamingHypergraph::get_numa_node_of_vertex(hn);
    ASSERT_EQ(part_id, hypergraph.partID(hn));
  }
}

}  // namespace mt_kahypar

/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "tests/partition/coarsening/coarsener_fixtures.h"
#include "mt-kahypar/partition/coarsening/multilevel_uncoarsener.h"
#include "mt-kahypar/partition/coarsening/nlevel_uncoarsener.h"
#include "mt-kahypar/partition/coarsening/coarsening_commons.h"

using ::testing::Test;

namespace mt_kahypar {
#ifdef USE_STRONG_PARTITIONER
using Coarsener = NLevelCoarsener<HeavyEdgeScore, NoWeightPenalty, BestRatingWithoutTieBreaking>;
using Uncoarsener = NLevelUncoarsener;
bool nlevel = true;
#else
using Coarsener = MultilevelCoarsener<HeavyEdgeScore, NoWeightPenalty, BestRatingWithoutTieBreaking>;
using Uncoarsener = MultilevelUncoarsener;
bool nlevel = false;
#endif

TEST_F(ACoarsener, DecreasesNumberOfPins) {
  context.coarsening.contraction_limit = 4;
  UncoarseningData uncoarseningData(nlevel, hypergraph, context);
  Coarsener coarsener(hypergraph, context, uncoarseningData);
  decreasesNumberOfPins(coarsener, 6);
}

TEST_F(ACoarsener, DecreasesNumberOfHyperedges) {
  context.coarsening.contraction_limit = 4;
  UncoarseningData uncoarseningData(nlevel, hypergraph, context);
  Coarsener coarsener(hypergraph, context, uncoarseningData);
  decreasesNumberOfHyperedges(coarsener, 3);
}

TEST_F(ACoarsener, RemovesHyperedgesOfSizeOneDuringCoarsening) {
  context.coarsening.contraction_limit = 4;
  UncoarseningData uncoarseningData(nlevel, hypergraph, context);
  Coarsener coarsener(hypergraph, context, uncoarseningData);
  doCoarsening(coarsener);
  auto& hypergraph = coarsener.coarsestHypergraph();
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_GE(hypergraph.edgeSize(he), 2);
  }
}

TEST_F(ACoarsener, RemovesParallelHyperedgesDuringCoarsening) {
  context.coarsening.contraction_limit = 4;
  UncoarseningData uncoarseningData(nlevel, hypergraph, context);
  Coarsener coarsener(hypergraph, context, uncoarseningData);
  doCoarsening(coarsener);
  auto& hypergraph = coarsener.coarsestHypergraph();
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(hypergraph.edgeWeight(he), 2);
  }
}

TEST_F(ACoarsener, ProjectsPartitionBackToOriginalHypergraph) {
  context.coarsening.contraction_limit = 4;
  UncoarseningData uncoarseningData(nlevel, hypergraph, context);
  Coarsener coarsener(hypergraph, context, uncoarseningData);
  Uncoarsener uncoarsener(hypergraph, context, uncoarseningData);
  context.type = kahypar::ContextType::initial_partitioning;
  doCoarsening(coarsener);
  PartitionedHyperGraph& coarsest_partitioned_hypergraph =
    coarsener.coarsestPartitionedHypergraph();
  assignPartitionIDs(coarsest_partitioned_hypergraph);
  PartitionedHyperGraph partitioned_hypergraph = uncoarsener.uncoarsen(nullptr_refiner, nullptr_refiner);
  for ( const HypernodeID& hn : partitioned_hypergraph.nodes() ) {
    PartitionID part_id = 0;
    ASSERT_EQ(part_id, partitioned_hypergraph.partID(hn));
  }
}

}  // namespace mt_kahypar

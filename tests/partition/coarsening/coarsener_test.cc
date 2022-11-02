/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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
  context.type = ContextType::initial_partitioning;
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

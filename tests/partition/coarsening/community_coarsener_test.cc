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

#include "mt-kahypar/partition/coarsening/community_coarsener.h"
#include "tests/partition/coarsening/community_coarsener_fixtures.h"

using ::testing::Test;

namespace mt_kahypar {
using Coarsener = CommunityCoarsenerT<ds::TestTypeTraits<2>, HeavyEdgeScore,
                                      NoWeightPenalty, BestRatingPreferringUnmatched>;

TEST_F(ACommunityCoarsener, DecreasesNumberOfPins) {
  Coarsener coarsener(hypergraph, context);
  decreasesNumberOfPins(coarsener, hypergraph, 14);
}

TEST_F(ACommunityCoarsener, DecreasesNumberOfHyperedges) {
  Coarsener coarsener(hypergraph, context);
  decreasesNumberOfHyperedges(coarsener, hypergraph, 7);
}

TEST_F(ACommunityCoarsener, RemovesHyperedgesOfSizeOneDuringCoarsening) {
  Coarsener coarsener(hypergraph, context);
  removesHyperedgesOfSizeOneDuringCoarsening(
    coarsener, hypergraph, { 0, 5, 281474976710656, 281474976710661 });
}

TEST_F(ACommunityCoarsener, ReAddsHyperedgesOfSizeOneDuringUncoarsening) {
  Coarsener coarsener(hypergraph, context);
  reAddsHyperedgesOfSizeOneDuringUncoarsening(
    coarsener, hypergraph, { 0, 5, 281474976710656, 281474976710661 });
}

TEST_F(ACommunityCoarsener, RemovesParallelHyperedgesDuringCoarsening) {
  Coarsener coarsener(hypergraph, context);
  removesParallelHyperedgesDuringCoarsening(
    coarsener, hypergraph, { 2, 4, 7, 9, 281474976710658, 281474976710660, 281474976710663 });
}

TEST_F(ACommunityCoarsener, UpdatesEdgeWeightOfRepresentativeHyperedgeOnParallelHyperedgeRemoval) {
  Coarsener coarsener(hypergraph, context);
  updatesEdgeWeightOfRepresentativeHyperedgeOnParallelHyperedgeRemoval(
    coarsener, hypergraph, { { 1, 2 }, { 3, 2 }, { 6, 2 }, { 8, 2 },
      { 281474976710657, 2 }, { 281474976710659, 2 }, { 281474976710662, 2 } });
}

TEST_F(ACommunityCoarsener, RestoresParallelHyperedgesDuringUncoarsening) {
  Coarsener coarsener(hypergraph, context);
  restoresParallelHyperedgesDuringUncoarsening(
    coarsener, hypergraph, { 2, 4, 7, 9, 281474976710658, 281474976710660, 281474976710663 });
}

TEST_F(ACommunityCoarsener, DoesNotCoarsenUntilCoarseningLimit) {
  Coarsener coarsener(hypergraph, context);
  doesNotCoarsenUntilCoarseningLimit(coarsener, hypergraph, context, 4, 3, 8);
}
}  // namespace mt_kahypar

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

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/zero_gain_cache.h"
#include "tests/datastructures/hypergraph_fixtures.h"

using ::testing::Test;

namespace mt_kahypar {

class AZeroGainCache : public ds::AHypergraph<1> {
 private:
  using Base = AHypergraph<1>;

 public:
  using Base::TestHypergraph;
  using Cache = ZeroGainCache<TestHypergraph>;

  AZeroGainCache() :
    Base(),
    hypergraph(construct_hypergraph(16,
                                    { { 0, 1,  2,  3,  4,  5,  6,  7 },
                                      { 8, 9, 10, 11, 12, 13, 14, 15 },
                                      { 1, 4 }, { 3, 6 }, { 9, 12 }, {11, 14} },
                                    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                                    { 0, 0, 0, 0, 0, 0 },
                                    { }, 4)),
    context(),
    zero_gain_cache(nullptr),
    id(16),
    delta(0) {

    context.partition.k = 4;
    context.partition.epsilon = 0.0;
    context.setupPartWeights(16);

    // Assign part ids
    for ( HypernodeID original_hn = 0; original_hn < 16; ++original_hn ) {
      const HypernodeID hn = hypergraph.globalNodeID(original_hn);
      hypergraph.setNodePart(hn, original_hn / 4);
      id[original_hn] = hn;
    }
    hypergraph.updateGlobalPartInfos();
    hypergraph.initializeNumCutHyperedges();

    zero_gain_cache = std::make_unique<Cache>(hypergraph.initialNumNodes(), context);
  }

  void insertAllZeroGainMoves() {
    zero_gain_cache->insert(hypergraph, id[0],  0, 1);
    zero_gain_cache->insert(hypergraph, id[2],  0, 1);
    zero_gain_cache->insert(hypergraph, id[5],  1, 0);
    zero_gain_cache->insert(hypergraph, id[7],  1, 0);
    zero_gain_cache->insert(hypergraph, id[8],  2, 3);
    zero_gain_cache->insert(hypergraph, id[10], 2, 3);
    zero_gain_cache->insert(hypergraph, id[13], 3, 2);
    zero_gain_cache->insert(hypergraph, id[15], 3, 2);
  }

  TestHypergraph hypergraph;
  Context context;
  std::unique_ptr<Cache> zero_gain_cache;
  std::vector<HypernodeID> id;
  Gain delta;
};

TEST_F(AZeroGainCache, InsertsAZeroGainMove) {
  Cache& cache = *zero_gain_cache;
  cache.insert(hypergraph, id[0], 0, 1);

  ASSERT_EQ(0, cache._cache_entry[0].from);
  ASSERT_TRUE(cache._cache_entry[0].valid_to[1]);
  ASSERT_EQ(id[0], cache._cache[0][1][0]);
}

TEST_F(AZeroGainCache, InsertsTwoZeroGainMoves) {
  Cache& cache = *zero_gain_cache;
  cache.insert(hypergraph, id[0], 0, 1);
  cache.insert(hypergraph, id[3], 0, 1);

  ASSERT_EQ(0, cache._cache_entry[3].from);
  ASSERT_TRUE(cache._cache_entry[3].valid_to[1]);
  ASSERT_EQ(id[3], cache._cache[0][1][1]);
}

TEST_F(AZeroGainCache, InsertsTwoZeroGainMovesToDifferentBlocks) {
  Cache& cache = *zero_gain_cache;
  cache.insert(hypergraph, id[0], 0, 1);
  cache.insert(hypergraph, id[5], 1, 0);

  ASSERT_EQ(1, cache._cache_entry[5].from);
  ASSERT_TRUE(cache._cache_entry[5].valid_to[0]);
  ASSERT_EQ(id[5], cache._cache[1][0][0]);
}

TEST_F(AZeroGainCache, InsertsTwoZeroGainMovesForSameVertex1) {
  Cache& cache = *zero_gain_cache;
  cache.insert(hypergraph, id[0], 0, 1);
  cache.insert(hypergraph, id[0], 0, 2);

  ASSERT_EQ(0, cache._cache_entry[0].from);
  ASSERT_TRUE(cache._cache_entry[0].valid_to[1]);
  ASSERT_TRUE(cache._cache_entry[0].valid_to[2]);
  ASSERT_EQ(id[0], cache._cache[0][1][0]);
  ASSERT_EQ(id[0], cache._cache[0][2][0]);
}

TEST_F(AZeroGainCache, InsertsTwoZeroGainMovesForSameVertex2) {
  Cache& cache = *zero_gain_cache;
  cache.insert(hypergraph, id[0], 0, 1);
  cache.insert(hypergraph, id[3], 0, 1);
  cache.insert(hypergraph, id[0], 0, 2);
  ASSERT_EQ(0, cache._cache_entry[0].from);
  ASSERT_TRUE(cache._cache_entry[0].valid_to[1]);
  ASSERT_TRUE(cache._cache_entry[0].valid_to[2]);
  ASSERT_EQ(id[0], cache._cache[0][1][0]);
  ASSERT_EQ(id[0], cache._cache[0][2][0]);
}

TEST_F(AZeroGainCache, ReinsertZeroGainMoveAfterChangeNodePart1) {
  Cache& cache = *zero_gain_cache;
  cache.insert(hypergraph, id[0], 0, 1);
  hypergraph.changeNodePart(id[0], 0, 1);
  cache.insert(hypergraph, id[0], 1, 0);

  ASSERT_EQ(1, cache._cache_entry[0].from);
  ASSERT_TRUE(cache._cache_entry[0].valid_to[0]);
  ASSERT_FALSE(cache._cache_entry[0].valid_to[1]);
  ASSERT_EQ(id[0], cache._cache[0][1][0]); // Will be removed lazily
  ASSERT_EQ(id[0], cache._cache[1][0][0]);
}

TEST_F(AZeroGainCache, ReinsertZeroGainMoveAfterChangeNodePart2) {
  Cache& cache = *zero_gain_cache;
  cache.insert(hypergraph, id[0], 0, 1);
  cache.insert(hypergraph, id[3], 0, 1);
  hypergraph.changeNodePart(id[0], 0, 1);
  cache.insert(hypergraph, id[0], 1, 0);

  ASSERT_EQ(1, cache._cache_entry[0].from);
  ASSERT_EQ(0, cache._cache_entry[3].from);
  ASSERT_TRUE(cache._cache_entry[0].valid_to[0]);
  ASSERT_FALSE(cache._cache_entry[0].valid_to[1]);
  ASSERT_TRUE(cache._cache_entry[3].valid_to[1]);
  ASSERT_FALSE(cache._cache_entry[3].valid_to[0]);
  ASSERT_EQ(id[0], cache._cache[0][1][0]); // Will be removed lazily
  ASSERT_EQ(id[3], cache._cache[0][1][1]);
  ASSERT_EQ(id[0], cache._cache[1][0][0]);
}


TEST_F(AZeroGainCache, CheckIfMoveIsPossible1) {
  insertAllZeroGainMoves();
  ASSERT_TRUE(zero_gain_cache->isMovePossible(hypergraph, id[1], 0, 1));
}

TEST_F(AZeroGainCache, CheckIfMoveIsPossible2) {
  insertAllZeroGainMoves();
  hypergraph.changeNodePart(id[5], 1, 0);
  hypergraph.changeNodePart(id[7], 1, 0);
  // There is no zero gain move left to move hn 1 from block 0 to 1
  ASSERT_FALSE(zero_gain_cache->isMovePossible(hypergraph, id[1], 0, 1));
}

TEST_F(AZeroGainCache, PerformRebalancingMoves1) {
  insertAllZeroGainMoves();
  ASSERT_TRUE(zero_gain_cache->performMove(hypergraph, id[1], 0, 1, delta));
  ASSERT_EQ(1, hypergraph.partID(id[1]));
  ASSERT_EQ(0, hypergraph.partID(id[5]));
  ASSERT_EQ(1, hypergraph.partID(id[7]));
}

TEST_F(AZeroGainCache, PerformRebalancingMoves2) {
  insertAllZeroGainMoves();
  hypergraph.changeNodePart(id[5], 1, 0);
  hypergraph.changeNodePart(id[2], 0, 1);
  ASSERT_TRUE(zero_gain_cache->performMove(hypergraph, id[1], 0, 1, delta));
  ASSERT_EQ(1, hypergraph.partID(id[1]));
  ASSERT_EQ(0, hypergraph.partID(id[5]));
  ASSERT_EQ(0, hypergraph.partID(id[7]));
}

TEST_F(AZeroGainCache, PerformRebalancingMoves3) {
  insertAllZeroGainMoves();
  ASSERT_TRUE(zero_gain_cache->performMove(hypergraph, id[1], 0, 1, delta));
  ASSERT_TRUE(zero_gain_cache->performMove(hypergraph, id[12], 3, 2, delta));
  ASSERT_EQ(1, hypergraph.partID(id[1]));
  ASSERT_EQ(0, hypergraph.partID(id[5]));
  ASSERT_EQ(1, hypergraph.partID(id[7]));
  ASSERT_EQ(2, hypergraph.partID(id[12]));
  ASSERT_EQ(3, hypergraph.partID(id[8]));
  ASSERT_EQ(2, hypergraph.partID(id[10]));
}


} // namespace mt_kahyper
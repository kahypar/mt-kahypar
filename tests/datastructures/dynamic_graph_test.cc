/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2022 Nikolai Maas <nikolai.maas@student.kit.edu>
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

#include <atomic>

#include "mt-kahypar/definitions.h"
#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/datastructures/dynamic_graph.h"
#include "mt-kahypar/datastructures/dynamic_graph_factory.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
namespace ds {

using ADynamicGraph = HypergraphFixture<DynamicGraph, DynamicGraphFactory, true>;

template<typename F, typename K>
void executeParallel(const F& f1, const K& f2) {
  std::atomic<size_t> cnt(0);
  tbb::parallel_invoke([&] {
    ++cnt;
    while ( cnt < 2 ) { }
    f1();
  }, [&] {
    ++cnt;
    while ( cnt < 2 ) { }
    f2();
  });
}

TEST_F(ADynamicGraph, HasCorrectStats) {
  ASSERT_EQ(7,  hypergraph.initialNumNodes());
  ASSERT_EQ(12,  hypergraph.initialNumEdges());
  ASSERT_EQ(12, hypergraph.initialNumPins());
  ASSERT_EQ(12, hypergraph.initialTotalVertexDegree());
  ASSERT_EQ(7,  hypergraph.totalWeight());
  ASSERT_EQ(2,  hypergraph.maxEdgeSize());
}

TEST_F(ADynamicGraph, HasCorrectInitialNodeIterator) {
  HypernodeID expected_hn = 0;
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_EQ(expected_hn++, hn);
  }
  ASSERT_EQ(7, expected_hn);
}

TEST_F(ADynamicGraph, HasCorrectNodeIteratorIfVerticesAreDisabled) {
  hypergraph.removeDegreeZeroHypernode(0);
  hypergraph.disableHypernode(5);
  const std::vector<HypernodeID> expected_iter =
    { 1, 2, 3, 4, 6 };
  HypernodeID pos = 0;
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_EQ(expected_iter[pos++], hn);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(ADynamicGraph, HasCorrectInitialEdgeIterator) {
  const HyperedgeID offset = DynamicAdjacencyArray::index_offset_per_node;
  std::vector<HyperedgeID> expected_iter;
  HyperedgeID current_offset = 0;
  // node 0
  current_offset += offset;
  // node 1
  current_offset += offset;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  // node 2
  current_offset += offset;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  // node 3
  current_offset += offset;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  // node 4
  current_offset += offset;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  // node 5
  current_offset += offset;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  // node 6
  current_offset += offset;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  expected_iter.push_back(current_offset);
  current_offset += 1;

  HypernodeID pos = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(expected_iter[pos++], he);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

} // namespace ds
} // namespace mt_kahypar

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

#include "tests/datastructures/hypergraph_fixtures.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

using AHypergraphWithTwoStreamingHypergraphs = AHypergraph<2>;
using TestHypergraph = typename AHypergraphWithTwoStreamingHypergraphs::TestHypergraph;

TestHypergraph construct_test_hypergraph(const AHypergraphWithTwoStreamingHypergraphs& test) {
  return test.construct_hypergraph(7, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} },
                                      { 0, 0, 0, 1, 1, 1, 1 },
                                      { 0, 0, 1, 1 } );
}

template< typename IDType >
auto identity = [](const IDType& id) { return id; };

template< typename IDType, typename F, typename K = decltype(identity<IDType>) >
void verifyIterator(const std::set<IDType>& reference, F&& it_func, K map_func = identity<IDType>) {
  size_t count = 0;
  for ( const IDType& id : it_func() ) {
    ASSERT_TRUE(reference.find(map_func(id)) != reference.end()) << V(map_func(id));
    count++;
  }
  ASSERT_EQ(count, reference.size());
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ContainsCorrectNumberofNodesEdgesAndPins) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  ASSERT_EQ(7, hypergraph.initialNumNodes());
  ASSERT_EQ(4, hypergraph.initialNumEdges());
  ASSERT_EQ(12, hypergraph.initialNumPins());
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalNodeIterators1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  verifyIterator<HypernodeID>({0, 1, 2}, [&] {
    return hypergraph.nodes(0);
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalNodeIterators2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  verifyIterator<HypernodeID>({281474976710656, 281474976710657, 281474976710658, 281474976710659}, [&] {
    return hypergraph.nodes(1);
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalNodeIteratorsWithDisabledHypernodes1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHypernode(1);

  verifyIterator<HypernodeID>({0, 2}, [&] {
    return hypergraph.nodes(0);
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalNodeIteratorsWithDisabledHypernodes2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHypernode(281474976710657);
  hypergraph.disableHypernode(281474976710658);
  
  verifyIterator<HypernodeID>({281474976710656, 281474976710659}, [&] {
    return hypergraph.nodes(1);
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksOriginalNodeIds1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  verifyIterator<HypernodeID>({0, 1, 2}, [&] {
    return hypergraph.nodes(0);
  }, [&](const HypernodeID& hn) {
    return hypergraph.originalNodeID(hn);
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksOriginalNodeIds2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  verifyIterator<HypernodeID>({3, 4, 5, 6}, [&] {
    return hypergraph.nodes(1);
  }, [&](const HypernodeID& hn) {
    return hypergraph.originalNodeID(hn);
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalNodeIterator) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  verifyIterator<HypernodeID>({0, 1, 2, 281474976710656, 281474976710657, 
    281474976710658, 281474976710659}, [&] {
    return hypergraph.nodes();
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalNodeIteratorWithDisabledHypernodes1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHypernode(2);
  hypergraph.disableHypernode(281474976710657);
  hypergraph.disableHypernode(281474976710659);
  
  verifyIterator<HypernodeID>({0, 1, 281474976710656, 281474976710658}, [&] {
    return hypergraph.nodes();
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalNodeIteratorWithDisabledHypernodes2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHypernode(0);
  hypergraph.disableHypernode(2);
  hypergraph.disableHypernode(281474976710656);
  hypergraph.disableHypernode(281474976710659);
  
  verifyIterator<HypernodeID>({1, 281474976710657, 281474976710658}, [&] {
    return hypergraph.nodes();
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalEdgeIterators1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  verifyIterator<HyperedgeID>({0, 1}, [&] {
    return hypergraph.edges(0);
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalEdgeIterators2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HyperedgeID>({281474976710656, 281474976710657}, [&] {
    return hypergraph.edges(1);
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalEdgeIteratorsWithDisabledHyperedges1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHyperedge(1);
  
  verifyIterator<HyperedgeID>({0}, [&] {
    return hypergraph.edges(0);
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalEdgeIteratorsWithDisabledHyperedges2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHyperedge(281474976710656);
  
  verifyIterator<HyperedgeID>({281474976710657}, [&] {
    return hypergraph.edges(1);
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalEdgeIterators) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  verifyIterator<HyperedgeID>({0, 1, 281474976710656, 281474976710657}, [&] {
    return hypergraph.edges();
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalEdgeIteratorsWithDisabledHyperedges1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHyperedge(1);
  hypergraph.disableHyperedge(281474976710656);
  
  verifyIterator<HyperedgeID>({0, 281474976710657}, [&] {
    return hypergraph.edges();
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalEdgeIteratorsWithDisabledHyperedges2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHyperedge(0);
  hypergraph.disableHyperedge(281474976710657);
  
  verifyIterator<HyperedgeID>({1, 281474976710656}, [&] {
    return hypergraph.edges();
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksIncidentEdgesOfHypernode1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  verifyIterator<HyperedgeID>({0, 1}, [&] {
    return hypergraph.incidentEdges(GLOBAL_ID(hypergraph, 0));
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksIncidentEdgesOfHypernode2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  verifyIterator<HyperedgeID>({0, 281474976710657}, [&] {
    return hypergraph.incidentEdges(GLOBAL_ID(hypergraph, 2));
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksIncidentEdgesOfHypernode3) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  verifyIterator<HyperedgeID>({1, 281474976710656}, [&] {
    return hypergraph.incidentEdges(GLOBAL_ID(hypergraph, 3));
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksIncidentEdgesOfHypernode4) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  verifyIterator<HyperedgeID>({281474976710657}, [&] {
    return hypergraph.incidentEdges(GLOBAL_ID(hypergraph, 5));
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksIncidentEdgesOfHypernode5) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  verifyIterator<HyperedgeID>({281474976710656, 281474976710657}, [&] {
    return hypergraph.incidentEdges(GLOBAL_ID(hypergraph, 6));
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksPinsOfHyperedge1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  verifyIterator<HypernodeID>({GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 2)}, [&] {
    return hypergraph.pins(0);
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksPinsOfHyperedge2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  verifyIterator<HypernodeID>({GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4)}, [&] {
    return hypergraph.pins(1);
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksPinsOfHyperedge3) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  verifyIterator<HypernodeID>({GLOBAL_ID(hypergraph, 3), 
    GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 6)}, [&] {
    return hypergraph.pins(281474976710656);
  });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksPinsOfHyperedge4) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  verifyIterator<HypernodeID>({GLOBAL_ID(hypergraph, 2), 
    GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)}, [&] {
    return hypergraph.pins(281474976710657);
  });
}

} // namespace ds
} // namespace mt_kahypar

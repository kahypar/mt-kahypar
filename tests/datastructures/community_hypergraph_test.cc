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

#include <atomic>

#include "gmock/gmock.h"

#include "tests/datastructures/hypergraph_fixtures.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

using ACommunityHypergraph = AHypergraph<2>;
using TestHypergraph = typename ACommunityHypergraph::TestHypergraph;
using TestStreamingHypergraph = typename ACommunityHypergraph::TestStreamingHypergraph;
using Memento = typename TestStreamingHypergraph::Memento;

TestHypergraph construct_test_hypergraph(const ACommunityHypergraph& test) {
  return test.construct_hypergraph(7, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} },
                                      { 0, 0, 0, 1, 1, 1, 1 },
                                      { 0, 0, 1, 1 },
                                      { 0, 0, 0, 1, 1, 2, 1 } );
}

void assignPartitionIDs(TestHypergraph& hypergraph) {
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    PartitionID part_id = TestStreamingHypergraph::get_numa_node_of_vertex(hn);
    hypergraph.setNodePart(hn, part_id);
  }
}

template< typename IDType, typename F >
void verifyIterator(const std::set<IDType>& reference, F&& it_func, bool log = false) {
  size_t count = 0;
  for ( const IDType& id : it_func() ) {
    if ( log ) LOG << V(id);
    ASSERT_TRUE(reference.find(id) != reference.end()) << V(id);
    count++;
  }
  ASSERT_EQ(count, reference.size());
}

void verifyPinIterators(const TestHypergraph& hypergraph,
                        const std::vector<HyperedgeID> hyperedges,
                        const std::vector<std::set<HypernodeID>>& references,
                        bool log = false) {
  ASSERT(hyperedges.size() == references.size());
  for ( size_t i = 0; i < hyperedges.size(); ++i ) {
    const HyperedgeID he = hyperedges[i];
    const std::set<HypernodeID>& reference = references[i];
    size_t count = 0;
    for ( const HypernodeID& pin : hypergraph.pins(he) ) {
      if ( log ) LOG << V(he) << V(pin);
      ASSERT_TRUE(reference.find(pin) != reference.end()) << V(he) << V(pin);
      count++;
    }
    ASSERT_EQ(count, reference.size());
  }
}

TEST_F(ACommunityHypergraph, InitializesCommunityHyperedges) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.initializeCommunityHyperedges();

  ASSERT_EQ(1, hypergraph.numCommunitiesOfHyperedge(0));
  ASSERT_EQ(2, hypergraph.numCommunitiesOfHyperedge(1));
  ASSERT_EQ(1, hypergraph.numCommunitiesOfHyperedge(281474976710656));
  ASSERT_EQ(3, hypergraph.numCommunitiesOfHyperedge(281474976710657));
}

TEST_F(ACommunityHypergraph, VerifiesCommunityHyperedgeSizes) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.initializeCommunityHyperedges();

  ASSERT_EQ(2, hypergraph.edgeSize(0, 0));
  ASSERT_EQ(2, hypergraph.edgeSize(1, 0));
  ASSERT_EQ(2, hypergraph.edgeSize(1, 1));
  ASSERT_EQ(3, hypergraph.edgeSize(281474976710656, 1));
  ASSERT_EQ(1, hypergraph.edgeSize(281474976710657, 0));
  ASSERT_EQ(1, hypergraph.edgeSize(281474976710657, 1));
  ASSERT_EQ(1, hypergraph.edgeSize(281474976710657, 2));
}

TEST_F(ACommunityHypergraph, VerifiesPinsOfCommunityHyperedges) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  hypergraph.initializeCommunityHyperedges();

  verifyIterator<HypernodeID>({id[0], id[2]}, [&] {
    return hypergraph.pins(0, 0);
  });

  verifyIterator<HypernodeID>({id[0], id[1]}, [&] {
    return hypergraph.pins(1, 0);
  });

  verifyIterator<HypernodeID>({id[3], id[4]}, [&] {
    return hypergraph.pins(1, 1);
  });

  verifyIterator<HypernodeID>({id[3], id[4], id[6]}, [&] {
    return hypergraph.pins(281474976710656, 1);
  });

  verifyIterator<HypernodeID>({id[2]}, [&] {
    return hypergraph.pins(281474976710657, 0);
  });

  verifyIterator<HypernodeID>({id[6]}, [&] {
    return hypergraph.pins(281474976710657, 1);
  });

  verifyIterator<HypernodeID>({id[5]}, [&] {
    return hypergraph.pins(281474976710657, 2);
  });
}

TEST_F(ACommunityHypergraph, VerifiesIncidentNetsOfHypernodes) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  hypergraph.initializeCommunityHyperedges();

  // Note, in case we iterate over incident edges of the community hyperedges
  // we skip single-pin nets

  verifyIterator<HyperedgeID>({0, 1}, [&] {
    return hypergraph.incidentEdges(id[0], 0);
  });

  verifyIterator<HyperedgeID>({1}, [&] {
    return hypergraph.incidentEdges(id[1], 0);
  });

  verifyIterator<HyperedgeID>({0, 281474976710657}, [&] {
    return hypergraph.incidentEdges(id[2], 0);
  });

  verifyIterator<HyperedgeID>({1, 281474976710656}, [&] {
    return hypergraph.incidentEdges(id[3], 1);
  });

  verifyIterator<HyperedgeID>({1, 281474976710656}, [&] {
    return hypergraph.incidentEdges(id[4], 1);
  });

  verifyIterator<HyperedgeID>({281474976710657}, [&] {
    return hypergraph.incidentEdges(id[5], 0);
  });

  verifyIterator<HyperedgeID>({281474976710656, 281474976710657}, [&] {
    return hypergraph.incidentEdges(id[6], 0);
  });
}

TEST_F(ACommunityHypergraph, ContractsTwoHypernodes1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  hypergraph.initializeCommunityHyperedges();

  HypernodeID u = id[0];
  HypernodeID v = id[2];
  hypergraph.contract(u, v, 0);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(2, hypergraph.nodeWeight(u));

  verifyIterator<HyperedgeID>({0, 1, 281474976710657}, [&] {
    return hypergraph.incidentEdges(u, 0);
  });

  verifyIterator<HyperedgeID>({id[0]}, [&] {
    return hypergraph.pins(0, 0);
  });

  verifyIterator<HyperedgeID>({id[0], id[1]}, [&] {
    return hypergraph.pins(1, 0);
  });

  verifyIterator<HyperedgeID>({id[0]}, [&] {
    return hypergraph.pins(281474976710657, 0);
  });
}

TEST_F(ACommunityHypergraph, ContractsTwoHypernodes2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  hypergraph.initializeCommunityHyperedges();

  HypernodeID u = id[0];
  HypernodeID v = id[1];
  hypergraph.contract(u, v, 0);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(2, hypergraph.nodeWeight(u));

  verifyIterator<HyperedgeID>({0, 1}, [&] {
    return hypergraph.incidentEdges(u, 0);
  });

  verifyIterator<HyperedgeID>({id[0], id[2]}, [&] {
    return hypergraph.pins(0, 0);
  });

  verifyIterator<HyperedgeID>({id[0]}, [&] {
    return hypergraph.pins(1, 0);
  });
}

TEST_F(ACommunityHypergraph, ContractsTwoHypernodes3) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  hypergraph.initializeCommunityHyperedges();

  HypernodeID u = id[1];
  HypernodeID v = id[2];
  hypergraph.contract(u, v, 0);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(2, hypergraph.nodeWeight(u));

  verifyIterator<HyperedgeID>({0, 1, 281474976710657}, [&] {
    return hypergraph.incidentEdges(u, 0);
  });

  verifyIterator<HyperedgeID>({id[0], id[1]}, [&] {
    return hypergraph.pins(0, 0);
  });

  verifyIterator<HyperedgeID>({id[0], id[1]}, [&] {
    return hypergraph.pins(1, 0);
  });

  verifyIterator<HyperedgeID>({id[1]}, [&] {
    return hypergraph.pins(281474976710657, 0);
  });
}

TEST_F(ACommunityHypergraph, ContractsTwoHypernodes4) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  hypergraph.initializeCommunityHyperedges();

  HypernodeID u = id[3];
  HypernodeID v = id[4];
  hypergraph.contract(u, v, 1);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(2, hypergraph.nodeWeight(u));

  verifyIterator<HyperedgeID>({1, 281474976710656}, [&] {
    return hypergraph.incidentEdges(u, 1);
  });

  verifyIterator<HyperedgeID>({id[3]}, [&] {
    return hypergraph.pins(1, 1);
  });

  verifyIterator<HyperedgeID>({id[3], id[6]}, [&] {
    return hypergraph.pins(281474976710656, 1);
  });
}

TEST_F(ACommunityHypergraph, ContractsTwoHypernodes5) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  hypergraph.initializeCommunityHyperedges();

  HypernodeID u = id[3];
  HypernodeID v = id[6];
  hypergraph.contract(u, v, 1);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(2, hypergraph.nodeWeight(u));

  verifyIterator<HyperedgeID>({1, 281474976710656, 281474976710657}, [&] {
    return hypergraph.incidentEdges(u, 1);
  });

  verifyIterator<HyperedgeID>({id[3], id[4]}, [&] {
    return hypergraph.pins(1, 1);
  });

  verifyIterator<HyperedgeID>({id[3], id[4]}, [&] {
    return hypergraph.pins(281474976710656, 1);
  });

  verifyIterator<HyperedgeID>({id[3]}, [&] {
    return hypergraph.pins(281474976710657, 1);
  });
}


TEST_F(ACommunityHypergraph, ContractsTwoHypernodes6) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  hypergraph.initializeCommunityHyperedges();

  HypernodeID u = id[4];
  HypernodeID v = id[6];
  hypergraph.contract(u, v, 1);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(2, hypergraph.nodeWeight(u));

  verifyIterator<HyperedgeID>({1, 281474976710656, 281474976710657}, [&] {
    return hypergraph.incidentEdges(u, 1);
  });

  verifyIterator<HyperedgeID>({id[3], id[4]}, [&] {
    return hypergraph.pins(1, 1);
  });

  verifyIterator<HyperedgeID>({id[3], id[4]}, [&] {
    return hypergraph.pins(281474976710656, 1);
  });

  verifyIterator<HyperedgeID>({id[4]}, [&] {
    return hypergraph.pins(281474976710657, 1);
  });
}

TEST_F(ACommunityHypergraph, ResetCommunityHyperedgesAfterContraction1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());
  hypergraph.initializeCommunityHyperedges();

  HypernodeID u = id[0];
  HypernodeID v = id[2];
  std::vector<Memento> mementos = { hypergraph.contract(u, v, 0) };

  hypergraph.resetCommunityHyperedges(mementos);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0] }, {id[0], id[1], id[3], id[4]}, {id[3], id[4], id[6]}, {id[0], id[5], id[6]} });

  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(mementos[0], parallel_he_representative);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0], id[2] }, {id[0], id[1], id[3], id[4]}, {id[3], id[4], id[6]}, {id[2], id[5], id[6]} });
}

TEST_F(ACommunityHypergraph, ResetCommunityHyperedgesAfterContraction2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());
  hypergraph.initializeCommunityHyperedges();

  HypernodeID u = id[0];
  HypernodeID v = id[1];
  std::vector<Memento> mementos = { hypergraph.contract(u, v, 0) };

  hypergraph.resetCommunityHyperedges(mementos);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0], id[2] }, {id[0], id[3], id[4]}, {id[3], id[4], id[6]}, {id[2], id[5], id[6]} });

  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(mementos[0], parallel_he_representative);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0], id[2] }, {id[0], id[1], id[3], id[4]}, {id[3], id[4], id[6]}, {id[2], id[5], id[6]} });
}

TEST_F(ACommunityHypergraph, ResetCommunityHyperedgesAfterContraction3) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());
  hypergraph.initializeCommunityHyperedges();

  HypernodeID u = id[3];
  HypernodeID v = id[4];
  std::vector<Memento> mementos = { hypergraph.contract(u, v, 1) };

  hypergraph.resetCommunityHyperedges(mementos);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0], id[2] }, {id[0], id[1], id[3]}, {id[3], id[6]}, {id[2], id[5], id[6]} });

  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(mementos[0], parallel_he_representative);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0], id[2] }, {id[0], id[1], id[3], id[4]}, {id[3], id[4], id[6]}, {id[2], id[5], id[6]} });
}

TEST_F(ACommunityHypergraph, ResetCommunityHyperedgesAfterContraction4) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());
  hypergraph.initializeCommunityHyperedges();

  HypernodeID u = id[4];
  HypernodeID v = id[6];
  std::vector<Memento> mementos = { hypergraph.contract(u, v, 1) };

  hypergraph.resetCommunityHyperedges(mementos);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0], id[2] }, {id[0], id[1], id[3], id[4]}, {id[3], id[4]}, {id[2], id[4], id[5]} });

  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(mementos[0], parallel_he_representative);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0], id[2] }, {id[0], id[1], id[3], id[4]}, {id[3], id[4], id[6]}, {id[2], id[5], id[6]} });
}

TEST_F(ACommunityHypergraph, ResetCommunityHyperedgesAfterSeveralContractions1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());
  hypergraph.initializeCommunityHyperedges();

  std::vector<Memento> mementos = { hypergraph.contract(id[0], id[2], 0),
                                    hypergraph.contract(id[3], id[4], 1) };

  hypergraph.resetCommunityHyperedges(mementos);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0] }, {id[0], id[1], id[3]}, {id[3], id[6]}, {id[0], id[5], id[6]} });

  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(mementos[1], parallel_he_representative);
  hypergraph.uncontract(mementos[0], parallel_he_representative);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0], id[2] }, {id[0], id[1], id[3], id[4]}, {id[3], id[4], id[6]}, {id[2], id[5], id[6]} });
}

TEST_F(ACommunityHypergraph, ResetCommunityHyperedgesAfterSeveralContractions2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());
  hypergraph.initializeCommunityHyperedges();

  std::vector<Memento> mementos = { hypergraph.contract(id[0], id[1], 0),
                                    hypergraph.contract(id[3], id[6], 1) };

  hypergraph.resetCommunityHyperedges(mementos);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0], id[2] }, {id[0], id[3], id[4]}, {id[3], id[4]}, {id[2], id[3], id[5]} });

  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(mementos[1], parallel_he_representative);
  hypergraph.uncontract(mementos[0], parallel_he_representative);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0], id[2] }, {id[0], id[1], id[3], id[4]}, {id[3], id[4], id[6]}, {id[2], id[5], id[6]} });
}

TEST_F(ACommunityHypergraph, ResetCommunityHyperedgesAfterSeveralContractions3) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());
  hypergraph.initializeCommunityHyperedges();

  std::vector<Memento> mementos = { hypergraph.contract(id[0], id[1], 0),
                                    hypergraph.contract(id[3], id[6], 1),
                                    hypergraph.contract(id[0], id[2], 0) };

  hypergraph.resetCommunityHyperedges(mementos);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0] }, {id[0], id[3], id[4]}, {id[3], id[4]}, {id[0], id[3], id[5]} });

  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(mementos[2], parallel_he_representative);
  hypergraph.uncontract(mementos[1], parallel_he_representative);
  hypergraph.uncontract(mementos[0], parallel_he_representative);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0], id[2] }, {id[0], id[1], id[3], id[4]}, {id[3], id[4], id[6]}, {id[2], id[5], id[6]} });
}

TEST_F(ACommunityHypergraph, ResetCommunityHyperedgesAfterSeveralContractions4) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());
  hypergraph.initializeCommunityHyperedges();

  std::vector<Memento> mementos = { hypergraph.contract(id[1], id[0], 0),
                                    hypergraph.contract(id[3], id[6], 1),
                                    hypergraph.contract(id[2], id[1], 0),
                                    hypergraph.contract(id[4], id[3], 1) };

  hypergraph.resetCommunityHyperedges(mementos);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { { id[2] }, {id[2], id[4]}, {id[4]}, {id[2], id[4], id[5]} });

  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(mementos[3], parallel_he_representative);
  hypergraph.uncontract(mementos[2], parallel_he_representative);
  hypergraph.uncontract(mementos[1], parallel_he_representative);
  hypergraph.uncontract(mementos[0], parallel_he_representative);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0], id[2] }, {id[0], id[1], id[3], id[4]}, {id[3], id[4], id[6]}, {id[2], id[5], id[6]} });
}

void doParallelContractions(TestHypergraph& hypergraph,
                            const std::vector<std::vector<Memento>>& mementos,
                            const std::vector<std::vector<HyperedgeID>>& remove_edge = {}) {
  ASSERT((PartitionID) mementos.size() == hypergraph.numCommunities());
  ASSERT(mementos.size() < std::thread::hardware_concurrency());
  ASSERT(remove_edge.size() == 0 || (PartitionID) remove_edge.size() == hypergraph.numCommunities());

  std::atomic<size_t> cnt(0);
  size_t num_threads = mementos.size();
  tbb::task_group group;
  tbb::task_arena arena(num_threads, 0);
  for ( PartitionID community_id = 0; community_id < (PartitionID) mementos.size(); ++community_id ) {
    arena.execute([&, community_id] {
      group.run([&, community_id] {
        cnt++;
        while( cnt < num_threads ) { }

        for ( size_t i = 0; i < mementos[community_id].size(); ++i ) {
          HypernodeID u = mementos[community_id][i].u;
          HypernodeID v = mementos[community_id][i].v;
          hypergraph.contract(u, v, community_id);
        }

        if ( remove_edge.size() > 0 ) {
          for ( size_t i = 0; i < remove_edge[community_id].size(); ++i ) {
            hypergraph.removeSinglePinEdge(remove_edge[community_id][i], community_id);
          }
        }
      });
    });
  }
  arena.execute([&] {
    group.wait();
  });
}

TEST_F(ACommunityHypergraph, DoesParallelContractionsOnHypergraph1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());
  hypergraph.initializeCommunityHyperedges();

  std::vector<Memento> mementos = { Memento { id[1], id[0] },
                                    Memento { id[3], id[6] },
                                    Memento { id[2], id[1] },
                                    Memento { id[4], id[3] } };

  std::vector<std::vector<Memento>> m = { {mementos[0], mementos[2]}, {mementos[1], mementos[3]}, {} };
  doParallelContractions(hypergraph, m);

  hypergraph.resetCommunityHyperedges(mementos);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { { id[2] }, {id[2], id[4]}, {id[4]}, {id[2], id[4], id[5]} });

  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(mementos[3], parallel_he_representative);
  hypergraph.uncontract(mementos[2], parallel_he_representative);
  hypergraph.uncontract(mementos[1], parallel_he_representative);
  hypergraph.uncontract(mementos[0], parallel_he_representative);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0], id[2] }, {id[0], id[1], id[3], id[4]}, {id[3], id[4], id[6]}, {id[2], id[5], id[6]} });
}

TEST_F(ACommunityHypergraph, DoesParallelContractionsOnHypergraph2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());
  hypergraph.initializeCommunityHyperedges();

  std::vector<Memento> mementos = { Memento { id[0], id[2] },
                                    Memento { id[3], id[4] },
                                    Memento { id[6], id[3] },
                                    Memento { id[0], id[1] } };

  std::vector<std::vector<Memento>> m = { {mementos[0], mementos[3]}, {mementos[1], mementos[2]}, {} };
  doParallelContractions(hypergraph, m);

  hypergraph.resetCommunityHyperedges(mementos);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0] }, {id[0], id[6]}, {id[6]}, {id[0], id[5], id[6]} });

  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(mementos[3], parallel_he_representative);
  hypergraph.uncontract(mementos[2], parallel_he_representative);
  hypergraph.uncontract(mementos[1], parallel_he_representative);
  hypergraph.uncontract(mementos[0], parallel_he_representative);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0], id[2] }, {id[0], id[1], id[3], id[4]}, {id[3], id[4], id[6]}, {id[2], id[5], id[6]} });
}

TEST_F(ACommunityHypergraph, DoesParallelContractionsOnHypergraph3) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());
  hypergraph.initializeCommunityHyperedges();

  std::vector<Memento> mementos = { Memento { id[2], id[1] },
                                    Memento { id[2], id[0] },
                                    Memento { id[6], id[4] },
                                    Memento { id[6], id[3] } };

  std::vector<std::vector<Memento>> m = { {mementos[0], mementos[1]}, {mementos[2], mementos[3]}, {} };
  doParallelContractions(hypergraph, m);

  hypergraph.resetCommunityHyperedges(mementos);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { { id[2] }, {id[2], id[6]}, {id[6]}, {id[2], id[5], id[6]} });

  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(mementos[3], parallel_he_representative);
  hypergraph.uncontract(mementos[2], parallel_he_representative);
  hypergraph.uncontract(mementos[1], parallel_he_representative);
  hypergraph.uncontract(mementos[0], parallel_he_representative);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0], id[2] }, {id[0], id[1], id[3], id[4]}, {id[3], id[4], id[6]}, {id[2], id[5], id[6]} });
}

TEST_F(ACommunityHypergraph, DoesParallelContractionsOnHypergraph4) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());
  hypergraph.initializeCommunityHyperedges();

  std::vector<Memento> mementos = { Memento { id[4], id[3] },
                                    Memento { id[1], id[2] },
                                    Memento { id[1], id[0] },
                                    Memento { id[4], id[6] } };

  std::vector<std::vector<Memento>> m = { {mementos[1], mementos[2]}, {mementos[0], mementos[3]}, {} };
  doParallelContractions(hypergraph, m);

  hypergraph.resetCommunityHyperedges(mementos);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[1] }, {id[1], id[4]}, {id[4]}, {id[1], id[4], id[5]} });

  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(mementos[3], parallel_he_representative);
  hypergraph.uncontract(mementos[2], parallel_he_representative);
  hypergraph.uncontract(mementos[1], parallel_he_representative);
  hypergraph.uncontract(mementos[0], parallel_he_representative);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0], id[2] }, {id[0], id[1], id[3], id[4]}, {id[3], id[4], id[6]}, {id[2], id[5], id[6]} });
}

TEST_F(ACommunityHypergraph, DoesParallelContractionsOnHypergraphWithEdgeRemoval) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = {GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
    GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6)};
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());
  hypergraph.initializeCommunityHyperedges();

  std::vector<Memento> mementos = { Memento { id[4], id[3] },
                                    Memento { id[1], id[2] },
                                    Memento { id[1], id[0] },
                                    Memento { id[4], id[6] } };

  std::vector<std::vector<Memento>> m = { {mementos[1], mementos[2]}, {mementos[0], mementos[3]}, {} };
  doParallelContractions(hypergraph, m, { {0}, {}, {} });

  hypergraph.resetCommunityHyperedges(mementos);

  verifyPinIterators(hypergraph, {1, 281474976710656, 281474976710657},
   { {id[1], id[4]}, {id[4]}, {id[1], id[4], id[5]} });

  assignPartitionIDs(hypergraph);
  hypergraph.restoreSinglePinHyperedge(0);
  hypergraph.uncontract(mementos[3], parallel_he_representative);
  hypergraph.uncontract(mementos[2], parallel_he_representative);
  hypergraph.uncontract(mementos[1], parallel_he_representative);
  hypergraph.uncontract(mementos[0], parallel_he_representative);

  verifyPinIterators(hypergraph, {0, 1, 281474976710656, 281474976710657},
   { {id[0], id[2] }, {id[0], id[1], id[3], id[4]}, {id[3], id[4], id[6]}, {id[2], id[5], id[6]} });
}

} // namespace ds
} // namespace mt_kahypar
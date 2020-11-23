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
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/sparsifier_hypergraph.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

using ASparsifierHypergraph = HypergraphFixture<Hypergraph, HypergraphFactory>;
using SparsifierHyperGraph = SparsifierHypergraph<Hypergraph, HypergraphFactory>;

void verifyPinsOfSparsifiedHypergraph(SparsifierHyperGraph& s_hypergraph,
                                      const HyperedgeID he,
                                      const parallel::scalable_vector<HypernodeID>& pins) {
  size_t pos = 0;
  for ( const HypernodeID& pin : s_hypergraph.pins(he) ) {
    ASSERT_EQ(pins[pos++], pin);
  }
  ASSERT_EQ(pins.size(), pos);
}

template <class F1, class F2>
void executeConcurrent(const F1& f1, const F2& f2) {
  std::atomic<int> cnt(0);
  int num_threads = std::min(2, TBBNumaArena::instance().total_number_of_threads());
  tbb::parallel_invoke([&] {
    cnt++;
    while (cnt < num_threads) { }
    f1();
  }, [&] {
    cnt++;
    while (cnt < num_threads) { }
    f2();
  });
}

template <class F1, class F2, class F3>
void executeConcurrent(const F1& f1, const F2& f2, const F3& f3) {
  std::atomic<int> cnt(0);
  int num_threads = std::min(3, TBBNumaArena::instance().total_number_of_threads());
  tbb::parallel_invoke([&] {
    cnt++;
    while (cnt < num_threads) { }
    f1();
  }, [&] {
    cnt++;
    while (cnt < num_threads) { }
    f2();
  }, [&] {
    cnt++;
    while (cnt < num_threads) { }
    f3();
  });
}

TEST_F(ASparsifierHypergraph, HasCorrectNumNodesAndEdges) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(7, s_hg.numNodes());
  ASSERT_EQ(4, s_hg.numEdges());
}

TEST_F(ASparsifierHypergraph, HasCorrectNodeWeights) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, s_hg.nodeWeight(0));
  ASSERT_EQ(1, s_hg.nodeWeight(1));
  ASSERT_EQ(1, s_hg.nodeWeight(2));
  ASSERT_EQ(1, s_hg.nodeWeight(3));
  ASSERT_EQ(1, s_hg.nodeWeight(4));
  ASSERT_EQ(1, s_hg.nodeWeight(5));
  ASSERT_EQ(1, s_hg.nodeWeight(6));
}

TEST_F(ASparsifierHypergraph, HasCorrectNodesDegrees) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(2, s_hg.nodeDegree(0));
  ASSERT_EQ(1, s_hg.nodeDegree(1));
  ASSERT_EQ(2, s_hg.nodeDegree(2));
  ASSERT_EQ(2, s_hg.nodeDegree(3));
  ASSERT_EQ(2, s_hg.nodeDegree(4));
  ASSERT_EQ(1, s_hg.nodeDegree(5));
  ASSERT_EQ(2, s_hg.nodeDegree(6));
}

TEST_F(ASparsifierHypergraph, HasCorrectEdgeWeight) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, s_hg.edgeWeight(0));
  ASSERT_EQ(1, s_hg.edgeWeight(1));
  ASSERT_EQ(1, s_hg.edgeWeight(2));
  ASSERT_EQ(1, s_hg.edgeWeight(3));
}

TEST_F(ASparsifierHypergraph, HasCorrectPins1) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyPinsOfSparsifiedHypergraph(s_hg, ID(0), {0, 2});
}

TEST_F(ASparsifierHypergraph, HasCorrectPins2) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyPinsOfSparsifiedHypergraph(s_hg, ID(1), {0, 1, 3, 4});
}

TEST_F(ASparsifierHypergraph, HasCorrectPins3) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyPinsOfSparsifiedHypergraph(s_hg, ID(2), {3, 4, 6});
}

TEST_F(ASparsifierHypergraph, HasCorrectPins4) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyPinsOfSparsifiedHypergraph(s_hg, ID(3), {2, 5, 6});
}

TEST_F(ASparsifierHypergraph, ContractsTwoVertices1) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  s_hg.contract(0, 2);
  ASSERT_EQ(1, s_hg.numRemovedNodes());
  ASSERT_EQ(2, s_hg.nodeWeight(0));
  ASSERT_TRUE(s_hg.nodeIsEnabled(0));
  ASSERT_FALSE(s_hg.nodeIsEnabled(2));

  Hypergraph hg = s_hg.sparsify();
  ASSERT_EQ(6, hg.initialNumNodes());
  ASSERT_EQ(2, hg.nodeWeight(0));
  verifyPins(hg, {0, 1, 2, 3},
    { { 0 }, { 0, 1, 2, 3 }, {2, 3, 5}, { 0, 4, 5 } });
}

TEST_F(ASparsifierHypergraph, ContractsTwoVertices2) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  s_hg.contract(3, 4);
  ASSERT_EQ(1, s_hg.numRemovedNodes());
  ASSERT_EQ(2, s_hg.nodeWeight(3));
  ASSERT_TRUE(s_hg.nodeIsEnabled(3));
  ASSERT_FALSE(s_hg.nodeIsEnabled(4));

  Hypergraph hg = s_hg.sparsify();
  ASSERT_EQ(6, hg.initialNumNodes());
  ASSERT_EQ(2, hg.nodeWeight(3));
  verifyPins(hg, {0, 1, 2, 3},
    { { 0, 2 }, { 0, 1, 3 }, {3, 5}, { 2, 4, 5 } });
}

TEST_F(ASparsifierHypergraph, ContractsSeveralVertices1) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  s_hg.contract(0, 2);
  s_hg.contract(3, 4);
  ASSERT_EQ(2, s_hg.numRemovedNodes());

  Hypergraph hg = s_hg.sparsify();
  ASSERT_EQ(5, hg.initialNumNodes());
  ASSERT_EQ(2, hg.nodeWeight(0));
  ASSERT_EQ(2, hg.nodeWeight(2));
  verifyPins(hg, {0, 1, 2, 3},
    { { 0 }, { 0, 1, 2 }, {2, 4}, { 0, 3, 4 } });
}

TEST_F(ASparsifierHypergraph, ContractsSeveralVertices2) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  s_hg.contract(0, 2);
  s_hg.contract(0, 1);
  ASSERT_EQ(2, s_hg.numRemovedNodes());

  Hypergraph hg = s_hg.sparsify();
  ASSERT_EQ(5, hg.initialNumNodes());
  ASSERT_EQ(3, hg.nodeWeight(0));
  verifyPins(hg, {0, 1, 2, 3},
    { { 0 }, { 0, 1, 2 }, {1, 2, 4}, { 0, 3, 4 } });
}

TEST_F(ASparsifierHypergraph, ContractsSeveralVerticesConcurrently1) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  executeConcurrent([&] {
    s_hg.contract(0, 2);
  }, [&] {
    s_hg.contract(3, 4);
  });
  ASSERT_EQ(2, s_hg.numRemovedNodes());

  Hypergraph hg = s_hg.sparsify();
  ASSERT_EQ(5, hg.initialNumNodes());
  ASSERT_EQ(2, hg.nodeWeight(0));
  ASSERT_EQ(2, hg.nodeWeight(2));
  verifyPins(hg, {0, 1, 2, 3},
    { { 0 }, { 0, 1, 2 }, {2, 4}, { 0, 3, 4 } });
}

TEST_F(ASparsifierHypergraph, ContractsSeveralVerticesConcurrently2) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  executeConcurrent([&] {
    s_hg.contract(0, 2);
  }, [&] {
    s_hg.contract(3, 4);
  }, [&] {
    s_hg.contract(5, 6);
  });
  ASSERT_EQ(3, s_hg.numRemovedNodes());

  Hypergraph hg = s_hg.sparsify();
  ASSERT_EQ(4, hg.initialNumNodes());
  ASSERT_EQ(2, hg.nodeWeight(0));
  ASSERT_EQ(2, hg.nodeWeight(2));
  ASSERT_EQ(2, hg.nodeWeight(3));
  verifyPins(hg, {0, 1, 2, 3},
    { { 0 }, { 0, 1, 2 }, {2, 3}, { 0, 3 } });
}

TEST_F(ASparsifierHypergraph, ContractsSeveralVerticesConcurrently3) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  executeConcurrent([&] {
    s_hg.contract(0, 2);
    s_hg.contract(0, 1);
  }, [&] {
    s_hg.contract(3, 4);
  }, [&] {
    s_hg.contract(5, 6);
  });
  ASSERT_EQ(4, s_hg.numRemovedNodes());

  Hypergraph hg = s_hg.sparsify();
  ASSERT_EQ(3, hg.initialNumNodes());
  ASSERT_EQ(3, hg.nodeWeight(0));
  ASSERT_EQ(2, hg.nodeWeight(1));
  ASSERT_EQ(2, hg.nodeWeight(2));
  verifyPins(hg, {0, 1, 2, 3},
    { { 0 }, { 0, 1 }, {1, 2}, { 0, 2 } });
}

TEST_F(ASparsifierHypergraph, ModifiesEdgeWeight1) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  s_hg.setEdgeWeight(0, 3);
  ASSERT_EQ(3, s_hg.edgeWeight(0));

  Hypergraph hg = s_hg.sparsify();
  ASSERT_EQ(3, hg.edgeWeight(0));
}

TEST_F(ASparsifierHypergraph, ModifiesEdgeWeight2) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  s_hg.setEdgeWeight(2, 4);
  ASSERT_EQ(4, s_hg.edgeWeight(2));

  Hypergraph hg = s_hg.sparsify();
  ASSERT_EQ(4, hg.edgeWeight(2));
}

TEST_F(ASparsifierHypergraph, ModifiesEdgeWeightConcurrently) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  s_hg.setEdgeWeight(2, 4);
  executeConcurrent([&] {
    s_hg.setEdgeWeight(0, 2);
  }, [&] {
    s_hg.setEdgeWeight(1, 5);
  }, [&] {
    s_hg.setEdgeWeight(3, 4);
  });
  ASSERT_EQ(2, s_hg.edgeWeight(0));
  ASSERT_EQ(5, s_hg.edgeWeight(1));
  ASSERT_EQ(4, s_hg.edgeWeight(3));

  Hypergraph hg = s_hg.sparsify();
  ASSERT_EQ(2, hg.edgeWeight(0));
  ASSERT_EQ(5, hg.edgeWeight(1));
  ASSERT_EQ(4, hg.edgeWeight(3));
}

TEST_F(ASparsifierHypergraph, RemovesAHyperedge1) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  s_hg.remove(0);
  ASSERT_EQ(1, s_hg.nodeDegree(0));
  ASSERT_EQ(1, s_hg.nodeDegree(2));

  Hypergraph hg = s_hg.sparsify();
  ASSERT_EQ(3, hg.initialNumEdges());
  verifyPins(hg, {0, 1, 2},
    { { 0, 1, 3, 4 }, {3, 4, 6}, { 2, 5, 6 } });
}

TEST_F(ASparsifierHypergraph, RemovesAHyperedge2) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  s_hg.remove(1);
  ASSERT_EQ(1, s_hg.nodeDegree(0));
  ASSERT_EQ(0, s_hg.nodeDegree(1));
  ASSERT_EQ(1, s_hg.nodeDegree(3));
  ASSERT_EQ(1, s_hg.nodeDegree(4));

  Hypergraph hg = s_hg.sparsify();
  ASSERT_EQ(3, hg.initialNumEdges());
  verifyPins(hg, {0, 1, 2},
    { { 0, 2 }, {3, 4, 6}, { 2, 5, 6 } });
}

TEST_F(ASparsifierHypergraph, RemovesHyperedgesConcurrently) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  executeConcurrent([&] {
    s_hg.remove(2);
  }, [&] {
    s_hg.remove(3);
  });
  ASSERT_EQ(2, s_hg.nodeDegree(0));
  ASSERT_EQ(1, s_hg.nodeDegree(1));
  ASSERT_EQ(1, s_hg.nodeDegree(2));
  ASSERT_EQ(1, s_hg.nodeDegree(3));
  ASSERT_EQ(1, s_hg.nodeDegree(4));
  ASSERT_EQ(0, s_hg.nodeDegree(5));
  ASSERT_EQ(0, s_hg.nodeDegree(6));

  Hypergraph hg = s_hg.sparsify();
  ASSERT_EQ(2, hg.initialNumEdges());
  verifyPins(hg, {0, 1},
    { { 0, 2 }, { 0, 1, 3, 4 } });
}

TEST_F(ASparsifierHypergraph, ReplacesAHyperedge1) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  parallel::scalable_vector<HyperedgeID> edge = { 4, 5, 6 };
  s_hg.replace(1, std::move(edge));
  ASSERT_EQ(1, s_hg.nodeDegree(0));
  ASSERT_EQ(0, s_hg.nodeDegree(1));
  ASSERT_EQ(2, s_hg.nodeDegree(2));
  ASSERT_EQ(1, s_hg.nodeDegree(3));
  ASSERT_EQ(2, s_hg.nodeDegree(4));
  ASSERT_EQ(2, s_hg.nodeDegree(5));
  ASSERT_EQ(3, s_hg.nodeDegree(6));

  Hypergraph hg = s_hg.sparsify();
  verifyPins(hg, {0, 1, 2, 3},
    { { 0, 2 }, {4, 5, 6}, {3, 4, 6}, { 2, 5, 6 } });
}

TEST_F(ASparsifierHypergraph, ReplacesAHyperedge2) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  parallel::scalable_vector<HyperedgeID> edge = { 0, 1, 2, 3, 4, 5, 6 };
  s_hg.replace(3, std::move(edge));
  ASSERT_EQ(3, s_hg.nodeDegree(0));
  ASSERT_EQ(2, s_hg.nodeDegree(1));
  ASSERT_EQ(2, s_hg.nodeDegree(2));
  ASSERT_EQ(3, s_hg.nodeDegree(3));
  ASSERT_EQ(3, s_hg.nodeDegree(4));
  ASSERT_EQ(1, s_hg.nodeDegree(5));
  ASSERT_EQ(2, s_hg.nodeDegree(6));

  Hypergraph hg = s_hg.sparsify();
  verifyPins(hg, {0, 1, 2, 3},
    { { 0, 2 }, {0, 1, 3, 4}, {3, 4, 6}, { 0, 1, 2, 3, 4, 5, 6 } });
}

TEST_F(ASparsifierHypergraph, ReplacesAHyperedgesConcurrently) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  executeConcurrent([&] {
    parallel::scalable_vector<HyperedgeID> edge = { 0, 4, 6 };
    s_hg.replace(1, std::move(edge));
  }, [&] {
    parallel::scalable_vector<HyperedgeID> edge = { 2, 3, 5 };
    s_hg.replace(2, std::move(edge));
  });
  ASSERT_EQ(2, s_hg.nodeDegree(0));
  ASSERT_EQ(0, s_hg.nodeDegree(1));
  ASSERT_EQ(3, s_hg.nodeDegree(2));
  ASSERT_EQ(1, s_hg.nodeDegree(3));
  ASSERT_EQ(1, s_hg.nodeDegree(4));
  ASSERT_EQ(2, s_hg.nodeDegree(5));
  ASSERT_EQ(2, s_hg.nodeDegree(6));

  Hypergraph hg = s_hg.sparsify();
  verifyPins(hg, {0, 1, 2, 3},
    { { 0, 2 }, {0, 4, 6}, {2, 3, 5}, { 2, 5, 6} });
}

TEST_F(ASparsifierHypergraph, ExecutesMixedModificationsConcurrently1) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  executeConcurrent([&] {
    parallel::scalable_vector<HyperedgeID> edge = { 0, 4, 6 };
    s_hg.replace(1, std::move(edge));
  }, [&] {
    s_hg.contract(3, 4);
  });
  ASSERT_EQ(1, s_hg.numRemovedNodes());
  ASSERT_EQ(2, s_hg.nodeDegree(0));
  ASSERT_EQ(0, s_hg.nodeDegree(1));
  ASSERT_EQ(2, s_hg.nodeDegree(2));
  ASSERT_EQ(1, s_hg.nodeDegree(3));
  ASSERT_EQ(1, s_hg.nodeDegree(5));
  ASSERT_EQ(3, s_hg.nodeDegree(6));

  Hypergraph hg = s_hg.sparsify();
  verifyPins(hg, {0, 1, 2, 3},
    { { 0, 2 }, {0, 3, 5}, {3, 5}, { 2, 4, 5} });
}

TEST_F(ASparsifierHypergraph, ExecutesMixedModificationsConcurrently2) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  executeConcurrent([&] {
    s_hg.remove(3);
  }, [&] {
    s_hg.contract(0, 1);
  });
  ASSERT_EQ(1, s_hg.numRemovedNodes());
  ASSERT_EQ(2, s_hg.nodeDegree(0));
  ASSERT_EQ(1, s_hg.nodeDegree(2));
  ASSERT_EQ(2, s_hg.nodeDegree(3));
  ASSERT_EQ(2, s_hg.nodeDegree(4));
  ASSERT_EQ(0, s_hg.nodeDegree(5));
  ASSERT_EQ(1, s_hg.nodeDegree(6));

  Hypergraph hg = s_hg.sparsify();
  verifyPins(hg, {0, 1, 2},
    { { 0, 1 }, {0, 2, 3}, {2, 3, 5} });
}

TEST_F(ASparsifierHypergraph, ExecutesMixedModificationsConcurrently3) {
  SparsifierHyperGraph s_hg(hypergraph, TBBNumaArena::GLOBAL_TASK_GROUP);
  executeConcurrent([&] {
    s_hg.remove(3);
  }, [&] {
    parallel::scalable_vector<HyperedgeID> edge = { 1, 2, 3 };
    s_hg.replace(2, std::move(edge));
  });
  ASSERT_EQ(2, s_hg.nodeDegree(0));
  ASSERT_EQ(2, s_hg.nodeDegree(1));
  ASSERT_EQ(2, s_hg.nodeDegree(2));
  ASSERT_EQ(2, s_hg.nodeDegree(3));
  ASSERT_EQ(1, s_hg.nodeDegree(4));
  ASSERT_EQ(0, s_hg.nodeDegree(5));
  ASSERT_EQ(0, s_hg.nodeDegree(6));

  Hypergraph hg = s_hg.sparsify();
  verifyPins(hg, {0, 1, 2},
    { { 0, 2 }, {0, 1, 3, 4}, {1, 2, 3} });
}

TEST_F(ASparsifierHypergraph, SparsifiesHypergraphWithDisabledHypernodes) {
  Hypergraph degree_zero_hg = HypergraphFactory::construct(
    TBBNumaArena::GLOBAL_TASK_GROUP, 7, 2,
    { { 3, 4, 5 }, {4, 5, 6} } );
  degree_zero_hg.removeHypernode(1);
  degree_zero_hg.removeHypernode(2);

  SparsifierHyperGraph s_hg(degree_zero_hg, TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(2, s_hg.numRemovedNodes());
  Hypergraph hg = s_hg.sparsify();
  ASSERT_EQ(5, hg.initialNumNodes());
  verifyPins(hg, {0, 1},
    { {1, 2, 3}, {2, 3, 4} });
}

TEST_F(ASparsifierHypergraph, SparsifiesHypergraphWithDisabledHyperedges) {
  Hypergraph single_pin_hg = HypergraphFactory::construct(
    TBBNumaArena::GLOBAL_TASK_GROUP, 7, 3,
    { {0}, { 3, 4, 5 }, {4, 5, 6} } );
  single_pin_hg.removeEdge(0);

  SparsifierHyperGraph s_hg(single_pin_hg, TBBNumaArena::GLOBAL_TASK_GROUP);
  Hypergraph hg = s_hg.sparsify();
  ASSERT_EQ(2, hg.initialNumEdges());
  verifyPins(hg, {0, 1},
    { {3, 4, 5}, {4, 5, 6} });
}

TEST_F(ASparsifierHypergraph, SparsifiesHypergraphWithDisabledHypernodesAndHyperedges) {
  Hypergraph single_pin_hg = HypergraphFactory::construct(
    TBBNumaArena::GLOBAL_TASK_GROUP, 7, 3,
    { {0}, { 3, 4, 5 }, {4, 5, 6} } );
  single_pin_hg.removeEdge(0);
  single_pin_hg.removeHypernode(1);
  single_pin_hg.removeHypernode(2);

  SparsifierHyperGraph s_hg(single_pin_hg, TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(2, s_hg.numRemovedNodes());
  Hypergraph hg = s_hg.sparsify();
  ASSERT_EQ(5, hg.initialNumNodes());
  ASSERT_EQ(2, hg.initialNumEdges());
  verifyPins(hg, {0, 1},
    { {1, 2, 3}, {2, 3, 4} });
}

} // namespace ds
} // namespace mt_kahypar
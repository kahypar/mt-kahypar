/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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
#include <atomic>

#include "mt-kahypar/datastructures/dynamic_adjacency_array.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

void verifyNeighbors(const HypernodeID u,
                     const HypernodeID num_nodes,
                     const DynamicAdjacencyArray& adjacency_array,
                     const std::set<HypernodeID>& _expected_neighbors,
                     bool strict = false) {
  size_t num_neighbors = 0;
  size_t degree = 0;
  std::vector<bool> actual_neighbors(num_nodes, false);
  for ( const HyperedgeID& he : adjacency_array.incidentEdges(u) ) {
    const HypernodeID neighbor = adjacency_array.edge(he).target;
    ASSERT_NE(_expected_neighbors.find(neighbor), _expected_neighbors.end())
      << "Vertex " << neighbor << " should not be neighbor of vertex " << u;
    ASSERT_EQ(u, adjacency_array.edge(he).source)
      << "Source of " << he << " (target: " << adjacency_array.edge(he).target << ") should be "
      << u << " but is " << adjacency_array.edge(he).source;
    ASSERT_TRUE(!strict || !actual_neighbors[neighbor])
      << "Vertex " << u << " contain duplicate edge with target " << neighbor;
    if (!actual_neighbors[neighbor]) {
      ++num_neighbors;
    }
    ++degree;
    actual_neighbors[neighbor] = true;
  }
  ASSERT_EQ(num_neighbors, _expected_neighbors.size());
  ASSERT_EQ(degree, adjacency_array.nodeDegree(u));
  ASSERT_TRUE(!strict || num_neighbors == degree);
}

void verifyEdges(HyperedgeID expected_num_edges, const DynamicAdjacencyArray& adjacency_array,
                 const std::set<std::pair<HypernodeID, HypernodeID>>& _expected_edges) {
  size_t num_edges = 0;
  for ( const HyperedgeID& he : adjacency_array.edges() ) {
    const HypernodeID source = adjacency_array.edge(he).source;
    const HypernodeID target = adjacency_array.edge(he).target;
    ASSERT_TRUE(_expected_edges.find({source, target}) != _expected_edges.end() ||
                _expected_edges.find({target, source}) != _expected_edges.end())
      << "Edge (" << source << ", " << target << ") is invalid.";
    ++num_edges;
  }
  ASSERT_EQ(num_edges, 2 * expected_num_edges);
}

kahypar::ds::FastResetFlagArray<> createFlagArray(const HypernodeID num_nodes,
                                                  const std::vector<HypernodeID>& contained_nodes) {
  kahypar::ds::FastResetFlagArray<> flag_array(num_nodes);
  for ( const HypernodeID& node : contained_nodes ) {
    flag_array.set(node, true);
  }
  return flag_array;
}

TEST(ADynamicAdjacencyArray, VerifyInitialEdges) {
  DynamicAdjacencyArray adjacency_array(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  verifyEdges(6, adjacency_array, { {1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6} });
}

TEST(ADynamicAdjacencyArray, VerifyEdgesAfterContractions1) {
  DynamicAdjacencyArray adjacency_array(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  adjacency_array.contract(2, 3);
  adjacency_array.contract(4, 6);
  verifyEdges(4, adjacency_array, { {1, 2}, {1, 4}, {4, 5} });
}

TEST(ADynamicAdjacencyArray, VerifyEdgesAfterContractions2) {
  DynamicAdjacencyArray adjacency_array(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  adjacency_array.contract(1, 4);
  adjacency_array.contract(1, 5);
  verifyEdges(4, adjacency_array, {{1, 2}, {2, 3}, {1, 6} });
}

TEST(ADynamicAdjacencyArray, VerifyEdgesAfterContractions3) {
  DynamicAdjacencyArray adjacency_array(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  adjacency_array.contract(2, 1);
  adjacency_array.contract(4, 6);
  adjacency_array.contract(2, 3);
  adjacency_array.contract(4, 0);
  adjacency_array.contract(2, 4);
  verifyEdges(2, adjacency_array, { {2, 5} });
  adjacency_array.contract(2, 5);
  verifyEdges(0, adjacency_array, { });
}

TEST(ADynamicAdjacencyArray, VerifyInitialNeighborsOfEachVertex) {
  DynamicAdjacencyArray adjacency_array(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  verifyNeighbors(0, 7, adjacency_array, { });
  verifyNeighbors(1, 7, adjacency_array, { 2, 4 });
  verifyNeighbors(2, 7, adjacency_array, { 1, 3 });
  verifyNeighbors(3, 7, adjacency_array, { 2 });
  verifyNeighbors(4, 7, adjacency_array, { 1, 5, 6 });
  verifyNeighbors(5, 7, adjacency_array, { 4, 6 });
  verifyNeighbors(6, 7, adjacency_array, { 4, 5 });
}

TEST(ADynamicAdjacencyArray, ContractTwoVertices1) {
  DynamicAdjacencyArray adjacency_array(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  adjacency_array.contract(0, 1);
  verifyNeighbors(0, 7, adjacency_array, { 2, 4 });
}

TEST(ADynamicAdjacencyArray, ContractTwoVertices2) {
  DynamicAdjacencyArray adjacency_array(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  adjacency_array.contract(1, 2);
  verifyNeighbors(1, 7, adjacency_array, { 3, 4 });
}

TEST(ADynamicAdjacencyArray, ContractTwoVertices3) {
  DynamicAdjacencyArray adjacency_array(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  adjacency_array.contract(1, 5);
  verifyNeighbors(1, 7, adjacency_array, { 2, 4, 6 });
  verifyNeighbors(2, 7, adjacency_array, { 1, 3 });
  verifyNeighbors(4, 7, adjacency_array, { 1, 6 });
  verifyNeighbors(6, 7, adjacency_array, { 1, 4 });
}

TEST(ADynamicAdjacencyArray, ContractSeveralVertices1) {
  DynamicAdjacencyArray adjacency_array(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  adjacency_array.contract(0, 1);
  adjacency_array.contract(0, 2);
  verifyNeighbors(0, 7, adjacency_array, { 3, 4 });
}

TEST(ADynamicAdjacencyArray, ContractSeveralVertices2) {
  DynamicAdjacencyArray adjacency_array(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  adjacency_array.contract(1, 0);
  adjacency_array.contract(1, 2);
  verifyNeighbors(1, 7, adjacency_array, { 3, 4 });
  adjacency_array.contract(4, 5);
  adjacency_array.contract(4, 6);
  verifyNeighbors(4, 7, adjacency_array, { 1 });
  adjacency_array.contract(1, 3);
  adjacency_array.contract(1, 4);
  verifyNeighbors(1, 7, adjacency_array, { });
}

TEST(ADynamicAdjacencyArray, UncontractTwoVertices1) {
  DynamicAdjacencyArray adjacency_array(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  adjacency_array.contract(1, 2);
  adjacency_array.uncontract(1, 2);
  verifyNeighbors(1, 7, adjacency_array, { 2, 4 });
  verifyNeighbors(2, 7, adjacency_array, { 1, 3 });
}

TEST(ADynamicAdjacencyArray, UncontractTwoVertices2) {
  DynamicAdjacencyArray adjacency_array(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  adjacency_array.contract(1, 5);
  adjacency_array.uncontract(1, 5);
  verifyNeighbors(1, 7, adjacency_array, { 2, 4 });
  verifyNeighbors(5, 7, adjacency_array, { 4, 6 });
  verifyNeighbors(2, 7, adjacency_array, { 1, 3 });
  verifyNeighbors(4, 7, adjacency_array, { 1, 5, 6 });
  verifyNeighbors(6, 7, adjacency_array, { 4, 5 });
}

TEST(ADynamicAdjacencyArray, UncontractSeveralVertices1) {
  DynamicAdjacencyArray adjacency_array(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  adjacency_array.contract(1, 2);
  adjacency_array.contract(1, 0);
  adjacency_array.contract(4, 5);
  adjacency_array.contract(4, 6);
  adjacency_array.contract(4, 3);
  adjacency_array.contract(4, 1);
  verifyNeighbors(4, 7, adjacency_array, { });
  adjacency_array.uncontract(4, 1);
  verifyNeighbors(1, 7, adjacency_array, { 4 });
  verifyNeighbors(4, 7, adjacency_array, { 1 });
  adjacency_array.uncontract(4, 3);
  verifyNeighbors(4, 7, adjacency_array, { 1 });
  verifyNeighbors(3, 7, adjacency_array, { 1 });
  adjacency_array.uncontract(4, 6);
  verifyNeighbors(4, 7, adjacency_array, { 1, 6 });
  verifyNeighbors(6, 7, adjacency_array, { 4 });
  adjacency_array.uncontract(4, 5);
  verifyNeighbors(4, 7, adjacency_array, { 1, 5, 6 });
  verifyNeighbors(5, 7, adjacency_array, { 4, 6 });
  adjacency_array.uncontract(1, 0);
  verifyNeighbors(0, 7, adjacency_array, { });
  verifyNeighbors(1, 7, adjacency_array, { 3, 4 });
  adjacency_array.uncontract(1, 2);
  verifyNeighbors(0, 7, adjacency_array, { });
  verifyNeighbors(1, 7, adjacency_array, { 2, 4 });
  verifyNeighbors(2, 7, adjacency_array, { 1, 3 });
  verifyNeighbors(3, 7, adjacency_array, { 2 });
  verifyNeighbors(4, 7, adjacency_array, { 1, 5, 6 });
  verifyNeighbors(5, 7, adjacency_array, { 4, 6 });
  verifyNeighbors(6, 7, adjacency_array, { 4, 5 });
}

TEST(ADynamicAdjacencyArray, UncontractSeveralVertices2) {
  DynamicAdjacencyArray adjacency_array(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  adjacency_array.contract(3, 1);
  adjacency_array.contract(3, 4);
  adjacency_array.contract(5, 6);
  adjacency_array.contract(3, 5);
  adjacency_array.contract(0, 2);
  adjacency_array.contract(0, 3);
  verifyNeighbors(0, 7, adjacency_array, { });
  adjacency_array.uncontract(0, 3);
  verifyNeighbors(0, 7, adjacency_array, { 3 });
  verifyNeighbors(3, 7, adjacency_array, { 0 });
  adjacency_array.uncontract(0, 2);
  verifyNeighbors(0, 7, adjacency_array, { });
  verifyNeighbors(2, 7, adjacency_array, { 3 });
  adjacency_array.uncontract(3, 5);
  verifyNeighbors(3, 7, adjacency_array, { 2, 5 });
  verifyNeighbors(5, 7, adjacency_array, { 3 });
  adjacency_array.uncontract(5, 6);
  verifyNeighbors(5, 7, adjacency_array, { 3, 6 });
  verifyNeighbors(6, 7, adjacency_array, { 3, 5 });
  adjacency_array.uncontract(3, 4);
  verifyNeighbors(3, 7, adjacency_array, { 2, 4 });
  verifyNeighbors(4, 7, adjacency_array, { 3, 5, 6 });
  adjacency_array.uncontract(3, 1);
  verifyNeighbors(0, 7, adjacency_array, { });
  verifyNeighbors(1, 7, adjacency_array, { 2, 4 });
  verifyNeighbors(2, 7, adjacency_array, { 1, 3 });
  verifyNeighbors(3, 7, adjacency_array, { 2 });
  verifyNeighbors(4, 7, adjacency_array, { 1, 5, 6 });
  verifyNeighbors(5, 7, adjacency_array, { 4, 6 });
  verifyNeighbors(6, 7, adjacency_array, { 4, 5 });
}

TEST(ADynamicAdjacencyArray, RemovesParrallelEdges1) {
  DynamicAdjacencyArray adjacency_array(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  adjacency_array.contract(2, 4);
  adjacency_array.removeParallelEdges();
  verifyNeighbors(1, 7, adjacency_array, { 2 }, true);
  verifyNeighbors(2, 7, adjacency_array, { 1, 3, 5, 6 }, true);
}

TEST(ADynamicAdjacencyArray, RemovesParrallelEdges2) {
  DynamicAdjacencyArray adjacency_array(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  adjacency_array.contract(5, 6);
  adjacency_array.removeParallelEdges();
  verifyNeighbors(5, 7, adjacency_array, { 4 }, true);
  verifyNeighbors(4, 7, adjacency_array, { 1, 5 }, true);
}

TEST(ADynamicAdjacencyArray, RemovesParrallelEdges3) {
  DynamicAdjacencyArray adjacency_array(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  adjacency_array.contract(4, 2);
  adjacency_array.contract(4, 5);
  adjacency_array.removeParallelEdges();
  verifyNeighbors(1, 7, adjacency_array, { 4 }, true);
  verifyNeighbors(3, 7, adjacency_array, { 4 }, true);
  verifyNeighbors(4, 7, adjacency_array, { 1, 3, 6 }, true);
  verifyNeighbors(6, 7, adjacency_array, { 4 }, true);
}


// using OwnershipVector = parallel::scalable_vector<SpinLock>;

// template<typename F, typename K>
// void executeParallel(const F& f1, const K& f2) {
//   std::atomic<size_t> cnt(0);
//   tbb::parallel_invoke([&] {
//     ++cnt;
//     while ( cnt < 2 ) { }
//     f1();
//   }, [&] {
//     ++cnt;
//     while ( cnt < 2 ) { }
//     f2();
//   });
// }


}  // namespace ds
}  // namespace mt_kahypar

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

#include "mt-kahypar/datastructures/incident_edge_array.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

void verifyNeighbors(const HypernodeID u,
                     const HypernodeID num_nodes,
                     const IncidentEdgeArray& incident_edges,
                     const std::set<HypernodeID>& _expected_neighbors) {
  size_t num_neighbors = 0;
  std::vector<bool> actual_neighbors(num_nodes, false);
  for ( const HyperedgeID& he : incident_edges.incidentEdges(u) ) {
    const HypernodeID source = incident_edges.edge(he).source;
    const HypernodeID neighbor = incident_edges.edge(he).target;
    ASSERT_TRUE(_expected_neighbors.find(neighbor) != _expected_neighbors.end())
      << "Vertex " << neighbor << " should not be neighbor of vertex " << u;
    ASSERT_FALSE(actual_neighbors[he])
      << "Vertex " << neighbor << " occurs more than once neighbors of vertex " << u;
    actual_neighbors[he] = true;
    ++num_neighbors;
  }
  ASSERT_EQ(num_neighbors, _expected_neighbors.size());
  ASSERT_EQ(num_neighbors, incident_edges.nodeDegree(u));
}

kahypar::ds::FastResetFlagArray<> createFlagArray(const HyperedgeID num_hyperedges,
                                                  const std::vector<HyperedgeID>& contained_hes) {
  kahypar::ds::FastResetFlagArray<> flag_array(num_hyperedges);
  for ( const HyperedgeID& he : contained_hes ) {
    flag_array.set(he, true);
  }
  return flag_array;
}

TEST(AIncidentEdgeArray, VerifyInitialNeighborsOfEachVertex) {
  IncidentEdgeArray incident_edges(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  verifyNeighbors(0, 7, incident_edges, { });
  verifyNeighbors(1, 7, incident_edges, { 2, 4 });
  verifyNeighbors(2, 7, incident_edges, { 1, 3 });
  verifyNeighbors(3, 7, incident_edges, { 2 });
  verifyNeighbors(4, 7, incident_edges, { 1, 5, 6 });
  verifyNeighbors(5, 7, incident_edges, { 4, 6 });
  verifyNeighbors(6, 7, incident_edges, { 4, 5 });
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

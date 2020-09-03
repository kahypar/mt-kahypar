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

#include "mt-kahypar/datastructures/incident_net_array.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

void verifyIncidentNets(const HypernodeID u,
                        const HyperedgeID num_hyperedges,
                        const IncidentNetArray& incident_nets,
                        const std::set<HyperedgeID>& _expected_incident_nets) {
  size_t num_incident_edges = 0;
  std::vector<bool> actual_incident_edges(num_hyperedges, false);
  for ( const HyperedgeID& he : incident_nets.incidentEdges(u) ) {
    ASSERT_TRUE(_expected_incident_nets.find(he) != _expected_incident_nets.end())
      << "Hyperedge " << he << " should be not part of incident nets of vertex " << u;
    ASSERT_FALSE(actual_incident_edges[he])
      << "Hyperedge " << he << " occurs more than once in incident nets of vertex " << u;
    actual_incident_edges[he] = true;
    ++num_incident_edges;
  }
  ASSERT_EQ(num_incident_edges, _expected_incident_nets.size());
}

TEST(AIncidentNetArray, VerifyInitialIncidentNetsOfEachVertex) {
  IncidentNetArray incident_nets(
    7, {{0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6}});
  verifyIncidentNets(0, 4, incident_nets, { 0, 1 });
  verifyIncidentNets(1, 4, incident_nets, { 1 });
  verifyIncidentNets(2, 4, incident_nets, { 0, 3 });
  verifyIncidentNets(3, 4, incident_nets, { 1, 2 });
  verifyIncidentNets(4, 4, incident_nets, { 1, 2 });
  verifyIncidentNets(5, 4, incident_nets, { 3 });
  verifyIncidentNets(6, 4, incident_nets, { 2, 3 });
}

}  // namespace ds
}  // namespace mt_kahypar

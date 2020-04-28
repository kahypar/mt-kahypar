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

#define HYPERGRAPH_UNIT_TEST true

#include "kahypar/definitions.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/hardware_topology.h"
#include "mt-kahypar/parallel/tbb_numa_arena.h"
#include "tests/parallel/topology_mock.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {
#define GLOBAL_ID(hypergraph, id) hypergraph.globalNodeID(id)
#define GLOBAL_NODE_ID(hypergraph, id) hypergraph.globalNodeID(id)
#define GLOBAL_EDGE_ID(hypergraph, id) hypergraph.globalEdgeID(id)

auto identity = [](const HypernodeID& id) { return id; };

class HypergraphFixture : public Test {
 public:
  HypergraphFixture() :
    hypergraph(HypergraphFactory::construct(TBBNumaArena::GLOBAL_TASK_GROUP,
      7 , 4, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} })) {
    id.resize(7);
    for ( const HypernodeID& hn : hypergraph.nodes() ) {
      id[hypergraph.originalNodeID(hn)] = hn;
    }
  }

  static void SetUpTestSuite() {
    TBBNumaArena::instance(HardwareTopology::instance().num_cpus());
  }

  template <typename K = decltype(identity)>
  void verifyIncidentNets(const Hypergraph& hg,
                          const HypernodeID hn,
                          const std::set<HypernodeID>& reference,
                          K map_func = identity,
                          bool log = false) {
    size_t count = 0;
    for (const HyperedgeID& he : hg.incidentEdges(hn)) {
      if (log) LOG << V(he) << V(map_func(he));
      ASSERT_TRUE(reference.find(map_func(he)) != reference.end()) << V(map_func(he));
      count++;
    }
    ASSERT_EQ(count, reference.size());
  }

  template <typename K = decltype(identity)>
  void verifyIncidentNets(const HypernodeID hn,
                          const std::set<HypernodeID>& reference,
                          K map_func = identity,
                          bool log = false) {
    verifyIncidentNets(hypergraph, hn, reference, map_func, log);
  }

  void verifyPins(const Hypergraph& hg,
                  const std::vector<HyperedgeID> hyperedges,
                  const std::vector< std::set<HypernodeID> >& references,
                  bool log = false) {
    ASSERT(hyperedges.size() == references.size());
    for (size_t i = 0; i < hyperedges.size(); ++i) {
      const HyperedgeID he = hyperedges[i];
      const std::set<HypernodeID>& reference = references[i];
      size_t count = 0;
      for (const HypernodeID& pin : hg.pins(he)) {
        if (log) LOG << V(he) << V(pin);
        ASSERT_TRUE(reference.find(pin) != reference.end()) << V(he) << V(pin);
        count++;
      }
      ASSERT_EQ(count, reference.size());
    }
  }

  void verifyPins(const std::vector<HyperedgeID> hyperedges,
                  const std::vector< std::set<HypernodeID> >& references,
                  bool log = false) {
    verifyPins(hypergraph, hyperedges, references, log);
  }

  void verifyCommunityPins(const Hypergraph& hg,
                           const PartitionID community_id,
                           const std::vector<HyperedgeID> hyperedges,
                           const std::vector< std::set<HypernodeID> >& references,
                           bool log = false) {
    ASSERT(hyperedges.size() == references.size());
    for (size_t i = 0; i < hyperedges.size(); ++i) {
      const HyperedgeID he = hyperedges[i];
      const std::set<HypernodeID>& reference = references[i];
      size_t count = 0;
      for (const HypernodeID& pin : hg.pins(he, community_id)) {
        if (log) LOG << V(he) << V(pin);
        ASSERT_TRUE(reference.find(pin) != reference.end()) << V(he) << V(pin);
        count++;
      }
      ASSERT_EQ(count, reference.size());
    }
  }

  void verifyCommunityPins(const PartitionID community_id,
                           const std::vector<HyperedgeID> hyperedges,
                           const std::vector< std::set<HypernodeID> >& references,
                           bool log = false) {
    verifyCommunityPins(hypergraph, community_id, hyperedges, references, log);
  }

  void assignCommunityIds() {
    hypergraph.setCommunityID(id[0], 0);
    hypergraph.setCommunityID(id[1], 0);
    hypergraph.setCommunityID(id[2], 0);
    hypergraph.setCommunityID(id[3], 1);
    hypergraph.setCommunityID(id[4], 1);
    hypergraph.setCommunityID(id[5], 2);
    hypergraph.setCommunityID(id[6], 2);
    hypergraph.initializeCommunities(TBBNumaArena::GLOBAL_TASK_GROUP);
  }

  Hypergraph hypergraph;
  std::vector<HypernodeID> id;
};

}  // namespace ds
}  // namespace mt_kahypar

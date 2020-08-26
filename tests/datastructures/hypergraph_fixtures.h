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

#include "kahypar/definitions.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/hardware_topology.h"
#include "mt-kahypar/parallel/tbb_numa_arena.h"
#include "tests/parallel/topology_mock.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

static auto identity = [](const HypernodeID& id) { return id; };

template<typename HyperGraph, typename HyperGraphFactory>
class HypergraphFixture : public Test {
 public:
  HypergraphFixture() :
    hypergraph(HyperGraphFactory::construct(TBBNumaArena::GLOBAL_TASK_GROUP,
      7 , 4, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} }, nullptr, nullptr, true)) {
  }

  template <typename K = decltype(identity)>
  void verifyIncidentNets(const HyperGraph& hg,
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

  void verifyPins(const HyperGraph& hg,
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

  void assignCommunityIds() {
    hypergraph.setCommunityID(0, 0);
    hypergraph.setCommunityID(1, 0);
    hypergraph.setCommunityID(2, 0);
    hypergraph.setCommunityID(3, 1);
    hypergraph.setCommunityID(4, 1);
    hypergraph.setCommunityID(5, 2);
    hypergraph.setCommunityID(6, 2);
  }

  HyperGraph hypergraph;
};

}  // namespace ds
}  // namespace mt_kahypar

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

#include "tests/parallel/topology_mock.h"
#include "mt-kahypar/parallel/hardware_topology.h"
#include "mt-kahypar/parallel/tbb_numa_arena.h"
#include "mt-kahypar/datastructures/hypergraph.h"
#include "mt-kahypar/definitions.h"
#include "kahypar/definitions.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

#define GLOBAL_ID(hypergraph, id) hypergraph.globalNodeID(id)

template < int NUM_NUMA_NODES >
struct TestTypeTraits {
  using TopoMock = mt_kahypar::parallel::TopologyMock<NUM_NUMA_NODES>;
  using HwTopology = mt_kahypar::parallel::HardwareTopology<TopoMock, parallel::topology_t, parallel::node_t>;
  using TBB = mt_kahypar::parallel::TBBNumaArena<HwTopology>;
  using HyperGraph = mt_kahypar::ds::Hypergraph<HypernodeID, HyperedgeID, 
    HypernodeWeight, HyperedgeWeight, PartitionID, HwTopology, TBB>;
  using StreamingHyperGraph = mt_kahypar::ds::StreamingHypergraph<HypernodeID, HyperedgeID, 
    HypernodeWeight, HyperedgeWeight, PartitionID, HwTopology, TBB>;
};

template < int NUM_NUMA_NODES >
class AHypergraph : public Test {
 private:
  using HyperedgeVector = parallel::scalable_vector<HyperedgeID>;
  using TBBArena = typename TestTypeTraits<NUM_NUMA_NODES>::TBB;

 public:
  using TestStreamingHypergraph = typename TestTypeTraits<NUM_NUMA_NODES>::StreamingHyperGraph;
  using TestHypergraph = typename TestTypeTraits<NUM_NUMA_NODES>::HyperGraph;

  AHypergraph() { 
    TBBArena::instance(std::thread::hardware_concurrency());
  }

  TestHypergraph construct_hypergraph(const HypernodeID num_hypernodes,
                                      const std::vector<HyperedgeVector>& hyperedges,
                                      const std::vector<HypernodeID>& node_mapping,
                                      const std::vector<HyperedgeID>& edge_mapping,
                                      const std::vector<PartitionID>& communities = {}) const {
    ASSERT(num_hypernodes == node_mapping.size());
    ASSERT(hyperedges.size() == edge_mapping.size());
    
    // Create hypergraphs
    std::vector<TestStreamingHypergraph> numa_hypergraphs;
    for ( int node = 0; node < NUM_NUMA_NODES; ++node ) {
      TBBArena::instance().numa_task_arena(node).execute([&] {
        numa_hypergraphs.emplace_back(node);
      });
    }             

    // Stream hyperedges
    for ( HyperedgeID node = 0; node < NUM_NUMA_NODES; ++node ) {
      TBBArena::instance().numa_task_arena(node).execute([&] {
        for ( size_t i = 0; i < hyperedges.size(); ++i ) {
          ASSERT(edge_mapping[i] < NUM_NUMA_NODES);
          if ( edge_mapping[i] == node ) {
            numa_hypergraphs[node].streamHyperedge(hyperedges[i], 1);
          }
        }
        numa_hypergraphs[node].initializeHyperedges(num_hypernodes);
      });
    }     

    // Create hypergraph (that also initialize hypernodes)
    TestHypergraph hypergraph(num_hypernodes, std::move(numa_hypergraphs), node_mapping);    

    if ( communities.size() > 0 ) {
      ASSERT(num_hypernodes == communities.size());
      for ( HypernodeID hn = 0; hn < num_hypernodes; ++hn ) {
        hypergraph.streamCommunityID(hypergraph.globalNodeID(hn), communities[hn]);
      }
      hypergraph.initializeCommunities();
    }
    
    return hypergraph; 
  }

};


} // namespace ds
} // namespace mt_kahypar

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

#pragma once

#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_invoke.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace ds {

template<class Hypergraph>
class CommunitySupport {

 static_assert(!Hypergraph::is_numa_aware,  "Only non-numa-aware hypergraphs are allowed");
 static_assert(!Hypergraph::is_partitioned,  "Only unpartitioned hypergraphs are allowed");

 using Counter = parallel::scalable_vector<HypernodeID>;
 using ThreadLocalCounter = tbb::enumerable_thread_specific<Counter>;

 public:
  explicit CommunitySupport() :
    _is_initialized(false),
    _num_communities(0),
    _communities_num_hypernodes(),
    _communities_num_pins(),
    _community_degree() { }

  CommunitySupport(const CommunitySupport&) = delete;
  CommunitySupport & operator= (const CommunitySupport &) = delete;

  CommunitySupport(CommunitySupport&& other) :
    _is_initialized(other._is_initialized),
    _num_communities(other._num_communities),
    _communities_num_hypernodes(std::move(other._communities_num_hypernodes)),
    _communities_num_pins(std::move(other._communities_num_pins)),
    _community_degree(std::move(other._community_degree)) { }

  CommunitySupport & operator= (CommunitySupport&& other) {
    _is_initialized = other._is_initialized;
    _num_communities = other._num_communities;
    _communities_num_hypernodes = std::move(other._communities_num_hypernodes);
    _communities_num_pins = std::move(other._communities_num_pins);
    _community_degree = std::move(other._community_degree);
    return *this;
  }

  /*!
   * Initializes community-related information after all vertices are assigned to a community.
   * This includes:
   *  1.) Number of Communities
   *  2.) Number of Vertices per Community
   *  3.) Number of Pins per Community
   *  4.) For each hypernode v of community C, we compute a unique id within
   *      that community in the range [0, |C|)
   * Note, in case 'hypergraph' is part of numa hypergraph, than 'hypergraphs' is
   * not empty and used to gather some information about communities, which 'hypergraph'
   * is not able to.
   */
  void initialize(const Hypergraph& hypergraph,
                  const int node,
                  const parallel::scalable_vector<Hypergraph>& hypergraphs = {}) {
    // Compute number of communities
    if ( hypergraphs.empty() ) {
      computeNumberOfCommunities(hypergraph);
    } else {
      computeNumberOfCommunities(hypergraphs);
    }

    if ( _num_communities > 1 ) {
      ThreadLocalCounter local_communities_num_hypernodes(_num_communities, 0);
      ThreadLocalCounter local_communities_num_pins(_num_communities, 0);
      ThreadLocalCounter local_community_degree(_num_communities, 0);
      // Iterate in parallel over all edges and vertices and gather stats about
      // communities. Stats are first aggregated in thread locals and afterwards
      // sequential in the member vector.
      tbb::parallel_invoke([&] {
        hypergraph.doParallelForAllNodes(TBBNumaArena::GLOBAL_TASK_GROUP, [&](const HypernodeID hn) {
          Counter& communities_num_hypernodes = local_communities_num_hypernodes.local();
          Counter& community_degree = local_community_degree.local();
          const PartitionID community_id = hypergraph.communityID(hn);
          ASSERT(community_id < _num_communities);
          ++communities_num_hypernodes[community_id];
          community_degree[community_id] += hypergraph.nodeDegree(hn);
        });
      }, [&] {
        hypergraph.doParallelForAllEdges(TBBNumaArena::GLOBAL_TASK_GROUP, [&](const HyperedgeID he) {
          Counter& communities_num_pins = local_communities_num_pins.local();
          for ( const HypernodeID& pin : hypergraph.pins(he) ) {
            const int hn_node = common::get_numa_node_of_vertex(pin);
            PartitionID community_id = kInvalidPartition;
            if ( hn_node != node ) {
              // In case the current pin is on an other NUMA node, we look up
              // its community id on the corresponding hypergraph, otherwise ...
              ASSERT(hn_node < static_cast<int>(hypergraphs.size()));
              community_id = hypergraphs[hn_node].communityID(pin);
            } else {
              // we look it up on the hypergraph this class is responsible for
              community_id = hypergraph.communityID(pin);
            }
            ASSERT(community_id != kInvalidPartition && community_id < _num_communities);
            ++communities_num_pins[community_id];
          }
        });
      }, [&] {
        _communities_num_hypernodes.assign(_num_communities, 0);
        _communities_num_pins.assign(_num_communities, 0);
        _community_degree.assign(_num_communities, 0);
      });

      // Aggregate thread locals in member vectors
      tbb::parallel_invoke([&] {
        for ( const Counter& communities_num_hypernodes : local_communities_num_hypernodes ) {
          for ( PartitionID community_id = 0; community_id < _num_communities; ++community_id ) {
            _communities_num_hypernodes[community_id] += communities_num_hypernodes[community_id];
          }
        }
      }, [&] {
        for ( const Counter& communities_num_pins : local_communities_num_pins ) {
          for ( PartitionID community_id = 0; community_id < _num_communities; ++community_id ) {
            _communities_num_pins[community_id] += communities_num_pins[community_id];
          }
        }
      }, [&] {
        for ( const Counter& community_degree : local_community_degree ) {
          for ( PartitionID community_id = 0; community_id < _num_communities; ++community_id ) {
            _community_degree[community_id] += community_degree[community_id];
          }
        }
      });
    } else {
      // Special case, if community information are not available
      // => all are assigned to the same community
      _communities_num_hypernodes.assign(1, hypergraph.initialNumNodes());
      _communities_num_pins.assign(1, hypergraph.initialNumPins());
      _community_degree.assign(1, hypergraph.initialTotalVertexDegree());
    }

    _is_initialized = true;
  }

  PartitionID numCommunities() const {
    ASSERT(_is_initialized);
    return _num_communities;
  }

  // ! Number of hypernodes in community
  HypernodeID numCommunityHypernodes(const PartitionID community) const {
    ASSERT(_is_initialized);
    ASSERT(community < _num_communities);
    return _communities_num_hypernodes[community];
  }

  // ! Number of pins in community
  HypernodeID numCommunityPins(const PartitionID community) const {
    ASSERT(_is_initialized);
    ASSERT(community < _num_communities);
    return _communities_num_pins[community];
  }

  HyperedgeID communityDegree(const PartitionID community) const {
    ASSERT(_is_initialized);
    ASSERT(community < _num_communities);
    return _community_degree[community];
  }

 private:
  void computeNumberOfCommunities(const Hypergraph& hypergraph,
                                  const int node = 0) {
    // The number of communities is the maximum community id plus 1
    _num_communities = tbb::parallel_reduce(tbb::blocked_range<HypernodeID>(0UL, hypergraph.initialNumNodes()),
      _num_communities, [&](const tbb::blocked_range<HypernodeID>& range, PartitionID init) {
        PartitionID num_communities = init;
        for (HypernodeID id = range.begin(); id < range.end(); ++id) {
          const HypernodeID hn = common::get_global_vertex_id(node, id);
          if ( hypergraph.nodeIsEnabled(hn) ) {
            num_communities = std::max(num_communities, hypergraph.communityID(hn) + 1);
          }
        }
        return num_communities;
      },
      [](const PartitionID lhs, const PartitionID rhs) {
        return std::max(lhs, rhs);
      });
    _num_communities = std::max(_num_communities, 1);
  }

  void computeNumberOfCommunities(const parallel::scalable_vector<Hypergraph>& hypergraphs) {
    ASSERT(hypergraphs.size() > 0);
    int node = 0;
    for ( const Hypergraph& hypergraph : hypergraphs ) {
      computeNumberOfCommunities(hypergraph, node);
      ++node;
    }
  }

  // ! Indicates, if community information are initialized
  bool _is_initialized;
  // ! Number of communities
  PartitionID _num_communities;
  // ! Number of hypernodes in a community
  parallel::scalable_vector<HypernodeID> _communities_num_hypernodes;
  // ! Number of pins in a community
  parallel::scalable_vector<HypernodeID> _communities_num_pins;
  // ! Total degree of a community
  parallel::scalable_vector<HyperedgeID> _community_degree;
};

} // namespace ds
} // namespace mt_kahypar
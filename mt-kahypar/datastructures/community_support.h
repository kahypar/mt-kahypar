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

#include <type_traits>

#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_invoke.h"
#include "tbb/parallel_reduce.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/memory_tree.h"
#include "mt-kahypar/utils/range.h"
#include <external_tools/kahypar/kahypar/utils/math.h>

namespace mt_kahypar {
namespace ds {

template<class Hypergraph>
class CommunitySupport {

 static constexpr bool enable_heavy_assert = false;

 static_assert(!Hypergraph::is_numa_aware,  "Only non-numa-aware hypergraphs are allowed");
 static_assert(!Hypergraph::is_partitioned, "Only unpartitioned hypergraphs are allowed");

 using Counter = parallel::scalable_vector<HypernodeID>;
 using AtomicCounter = parallel::scalable_vector<parallel::IntegralAtomicWrapper<HypernodeID>>;
 using ThreadLocalCounter = tbb::enumerable_thread_specific<Counter>;

  /**
   * Represents a community hyperedge of the hypergraph and contains all information
   * associated with it.
   */
  class CommunityHyperedge {
   public:
    CommunityHyperedge() :
      _community_id(kInvalidPartition),
      _begin(0),
      _size(0),
      _weight(0),
      _hash(kEdgeHashSeed),
      _valid(false) { }

    CommunityHyperedge(const PartitionID community_id,
                       const size_t begin,
                       const size_t size,
                       const HypernodeWeight weight) :
      _community_id(community_id),
      _begin(begin),
      _size(size),
      _weight(weight),
      _hash(kEdgeHashSeed),
      _valid(true) { }

    void disable() {
      ASSERT(!isDisabled());
      _valid = false;
    }

    bool isDisabled() const {
      return _valid == false;
    }

    PartitionID communityID() const {
      return _community_id;
    }

    // ! Returns the index of the first element in _incidence_array
    size_t firstEntry() const {
      return _begin;
    }

    // ! Sets the index of the first element in _incidence_array to begin
    void setFirstEntry(size_t begin) {
      _begin = begin;
    }

    // ! Returns the index of the first element in _incidence_array
    size_t firstInvalidEntry() const {
      return _begin + _size;
    }

    size_t size() const {
      return _size;
    }

    void setSize(size_t size) {
      _size = size;
    }

    void incrementSize() {
      ++_size;
    }

    void decrementSize() {
      ASSERT(_size > 0);
      --_size;
    }

    HyperedgeWeight weight() const {
      return _weight;
    }

    void setWeight(HyperedgeWeight weight) {
      _weight = weight;
    }

    size_t& hash() {
      return _hash;
    }

    size_t hash() const {
      return _hash;
    }

    bool operator== (const CommunityHyperedge& rhs) const {
      return _begin == rhs._begin && _size == rhs._size;
    }

    bool operator!= (const CommunityHyperedge& rhs) const {
      return !operator== (this, rhs);
    }

   private:
    // ! Community id of hyperedge
    PartitionID _community_id;
    // ! Index of the first element in _incidence_array
    size_t _begin;
    // ! Number of _incidence_array elements
    size_t _size;
    // ! hyperedge weight
    HyperedgeWeight _weight;
    // ! Hash of pins
    size_t _hash;
    // ! Flag indicating whether or not the element is active.
    bool _valid;
  };

  static_assert(std::is_trivially_copyable<CommunityHyperedge>::value,
    "Community Hyperedge is not trivally copyable");

  using CommunitiesOfHyperedges = parallel::scalable_vector<parallel::scalable_vector<PartitionID> >;
  using CommunityHyperedges = parallel::scalable_vector<parallel::scalable_vector<CommunityHyperedge> >;
  // ! Iterator to iterate over the pins of a hyperedge
  using IncidenceIterator = typename Hypergraph::IncidenceIterator;

 public:
  using CommunityIterator = parallel::scalable_vector<PartitionID>::const_iterator;

  explicit CommunitySupport() :
    _node(0),
    _is_initialized(false),
    _num_communities(0),
    _communities_num_hypernodes(),
    _communities_num_pins(),
    _community_degree(),
    _are_community_hyperedges_initialized(false),
    _community_hyperedge_ids(),
    _community_hyperedges(),
    _vertex_to_community_node_id() { }

  CommunitySupport(const CommunitySupport&) = delete;
  CommunitySupport & operator= (const CommunitySupport &) = delete;

  CommunitySupport(CommunitySupport&& other) :
    _node(other._node),
    _is_initialized(other._is_initialized),
    _num_communities(other._num_communities),
    _communities_num_hypernodes(std::move(other._communities_num_hypernodes)),
    _communities_num_pins(std::move(other._communities_num_pins)),
    _community_degree(std::move(other._community_degree)),
    _are_community_hyperedges_initialized(other._are_community_hyperedges_initialized),
    _community_hyperedge_ids(std::move(other._community_hyperedge_ids)),
    _community_hyperedges(std::move(other._community_hyperedges)),
    _vertex_to_community_node_id(std::move(other._vertex_to_community_node_id)) { }

  CommunitySupport & operator= (CommunitySupport&& other) {
    _node = other._node;
    _is_initialized = other._is_initialized;
    _num_communities = other._num_communities;
    _communities_num_hypernodes = std::move(other._communities_num_hypernodes);
    _communities_num_pins = std::move(other._communities_num_pins);
    _community_degree = std::move(other._community_degree);
    _are_community_hyperedges_initialized = other._are_community_hyperedges_initialized;
    _community_hyperedge_ids = std::move(other._community_hyperedge_ids);
    _community_hyperedges = std::move(other._community_hyperedges);
    _vertex_to_community_node_id = std::move(other._vertex_to_community_node_id);
    return *this;
  }

  bool isInitialized() const {
    return _is_initialized;
  }

  bool areCommunityHyperedgesInitialized() const {
    return _are_community_hyperedges_initialized;
  }

  // ! Number of communities
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

  // ! Returns a range to loop over the pins of hyperedge e that belong to a certain community.
  // ! Note, this function fails if community hyperedges are not initialized.
  IteratorRange<IncidenceIterator> pins(const Hypergraph& hypergraph,
                                        const HyperedgeID e,
                                        const PartitionID community_id) const {
    ASSERT(_are_community_hyperedges_initialized);
    const CommunityHyperedge& community_he = community_hyperedge(e, community_id);
    return IteratorRange<IncidenceIterator>(
      hypergraph._incidence_array.cbegin() + community_he.firstEntry(),
      hypergraph._incidence_array.cbegin() + community_he.firstInvalidEntry());
  }

  // ! Returns a range to loop over the set of communities contained in hyperedge e.
  IteratorRange<CommunityIterator> communities(const HyperedgeID e) const {
    ASSERT(_are_community_hyperedges_initialized);
    HyperedgeID local_id = common::get_local_position_of_edge(e);
    ASSERT(local_id < _community_hyperedge_ids.size(), "Hyperedge" << e << "does not exist");
    ASSERT(_node == common::get_numa_node_of_edge(e), "Hyperedge" << e << "is not part of NUMA node" << _node);
    return IteratorRange<CommunityIterator>(
      _community_hyperedge_ids[local_id].cbegin(),
      _community_hyperedge_ids[local_id].cend());
  }

  // ! Consider hypernode u is part of community C = {v_1, ..., v_n},
  // ! than this function returns a unique id for hypernode u in the
  // ! range [0,n).
  HypernodeID communityNodeId(const HypernodeID u) const {
    const HypernodeID local_id = common::get_local_position_of_vertex(u);
    ASSERT(local_id < _vertex_to_community_node_id.size(), "Hypernode" << u << "does not exist");
    ASSERT(_node == common::get_numa_node_of_vertex(u), "Hypernode" << u << "is not part of NUMA node" << _node);
    return _vertex_to_community_node_id[local_id];
  }

  // ! Weight of a community hyperedge
  HypernodeWeight edgeWeight(const HyperedgeID e, const PartitionID community_id) const {
    ASSERT(_are_community_hyperedges_initialized);
    return community_hyperedge(e, community_id).weight();
  }

  // ! Sets the weight of a community hyperedge
  void setEdgeWeight(const HyperedgeID e, const PartitionID community_id, const HyperedgeWeight weight) {
    ASSERT(_are_community_hyperedges_initialized);
    community_hyperedge(e, community_id).setWeight(weight);
  }

  // ! Number of pins of a hyperedge that are assigned to a community
  HypernodeID edgeSize(const HyperedgeID e, const PartitionID community_id) const {
    ASSERT(_are_community_hyperedges_initialized);
    return community_hyperedge(e, community_id).size();
  }

  // ! Hash value defined over the pins of a hyperedge that belongs to a community
  size_t edgeHash(const HyperedgeID e, const PartitionID community_id) const {
    ASSERT(_are_community_hyperedges_initialized);
    return community_hyperedge(e, community_id).hash();
  }

  // ! Number of communities which pins of hyperedge belongs to
  size_t numCommunitiesInHyperedge(const HyperedgeID e) const {
    ASSERT(_are_community_hyperedges_initialized);
    const HyperedgeID local_pos = common::get_local_position_of_edge(e);
    ASSERT(local_pos < _community_hyperedges.size());
    return _community_hyperedges[local_pos].size();
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
                  const parallel::scalable_vector<Hypergraph>& hypergraphs,
                  const TaskGroupID task_group_id) {
    _node = hypergraph.numaNode();
    // Compute number of communities
    if ( hypergraphs.empty() ) {
      computeNumberOfCommunities(hypergraph);
    } else {
      computeNumberOfCommunities(hypergraphs);
    }

    if ( _num_communities > 1 ) {
      AtomicCounter tmp_communities_num_hypernodes(_num_communities,
        parallel::IntegralAtomicWrapper<HypernodeID>(0));
      ThreadLocalCounter local_communities_num_pins(_num_communities, 0);
      ThreadLocalCounter local_community_degree(_num_communities, 0);
      // Iterate in parallel over all edges and vertices and gather stats about
      // communities. Stats are first aggregated in thread locals and afterwards
      // sequential in the member vector.
      tbb::parallel_invoke([&] {
        _vertex_to_community_node_id.resize(hypergraph.initialNumNodes());
        hypergraph.doParallelForAllNodes(task_group_id, [&](const HypernodeID hn) {
          Counter& community_degree = local_community_degree.local();
          const PartitionID community_id = hypergraph.communityID(hn);
          ASSERT(community_id < _num_communities);
          const HypernodeID local_id = common::get_local_position_of_vertex(hn);
          ASSERT(local_id < _vertex_to_community_node_id.size());
          _vertex_to_community_node_id[local_id] = tmp_communities_num_hypernodes[community_id]++;
          community_degree[community_id] += hypergraph.nodeDegree(hn);
        });
      }, [&] {
        hypergraph.doParallelForAllEdges(task_group_id, [&](const HyperedgeID he) {
          Counter& communities_num_pins = local_communities_num_pins.local();
          for ( const HypernodeID& pin : hypergraph.pins(he) ) {
            const int hn_node = common::get_numa_node_of_vertex(pin);
            PartitionID community_id = kInvalidPartition;
            if ( hn_node != _node ) {
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
        for ( PartitionID community_id = 0; community_id < _num_communities; ++community_id ) {
          _communities_num_hypernodes[community_id] = tmp_communities_num_hypernodes[community_id].load();
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
        _vertex_to_community_node_id.resize(hypergraph.initialNumNodes());
      _communities_num_hypernodes.assign(1, hypergraph.initialNumNodes());
      _communities_num_pins.assign(1, hypergraph.initialNumPins());
      _community_degree.assign(1, hypergraph.initialTotalVertexDegree());
    }

    // In case 'hypergraph' is part of a numa-aware hypergraph, finalizeCommunityNodeIds
    // have to be called in order to initialize the community node ids.
    _is_initialized = true;
  }

  // ! In order to get unique community node ids in case 'hypergraph' is part of a numa-aware
  // ! hypergraph, we have to add the prefix sum over the number of nodes in each community
  // ! of hypergraphs on a numa node with an id smaller than the current hypergraph to all
  // ! local community node ids.
  void finalizeCommunityNodeIds(const Hypergraph& hypergraph,
                                const parallel::scalable_vector<Hypergraph>& hypergraphs,
                                const TaskGroupID task_group_id) {
    ASSERT(_is_initialized);
    ASSERT(!hypergraphs.empty());
    if ( hypergraph.numaNode() == 0 ) {
      return;
    }

    parallel::scalable_vector<HypernodeID> num_hypernodes_prefix_sum(_num_communities, 0);
    tbb::parallel_for(0, _num_communities, [&](const PartitionID community_id) {
      for ( const Hypergraph& hg : hypergraphs ) {
        if ( hg.numaNode() < hypergraph.numaNode() ) {
          num_hypernodes_prefix_sum[community_id] += hg.numCommunityHypernodes(community_id);
        }
      }
    });

    hypergraph.doParallelForAllNodes(task_group_id, [&](const HypernodeID hn) {
      const PartitionID community_id = hypergraph.communityID(hn);
      const HypernodeID local_id = common::get_local_position_of_vertex(hn);
      ASSERT(local_id < _vertex_to_community_node_id.size());
      _vertex_to_community_node_id[local_id] += num_hypernodes_prefix_sum[community_id];
    });

    _is_initialized = true;
  }

  /*!
  * Initializes community hyperedges.
  * This includes:
  *   1.) Sort the pins of each hyperedge in increasing order of their community id
  *   2.) Introduce for each community id contained in a hyperedge a seperate
  *       community hyperedge pointing to a range of consecutive pins with
  *       same community in that hyperedge
  */
  void initializeCommunityHyperedges(Hypergraph& hypergraph,
                                     const parallel::scalable_vector<Hypergraph>& hypergraphs) {
    auto get_community_id =
      [&](const HypernodeID& hn) {
      if ( hypergraphs.empty() ) {
        return hypergraph.communityID(hn);
      } else {
        return common::hypergraph_of_vertex(hn, hypergraphs).communityID(hn);
      }
    };

    auto add_community_hyperedge =
      [&](const HyperedgeID he,
          const PartitionID community_id,
          const size_t start,
          const size_t end,
          const HyperedgeWeight weight) {
        ASSERT(community_id != kInvalidPartition);
        ASSERT(he < _community_hyperedges.size());
        ASSERT(start < end);
        _community_hyperedge_ids[he].push_back(community_id);
        _community_hyperedges[he].emplace_back(community_id, start, end - start, weight);

        // Compute community hyperedge hash
        for (size_t pos = start; pos < end; ++pos) {
          _community_hyperedges[he].back().hash() +=
            kahypar::math::hash(hypergraph._incidence_array[pos]);
        }
      };

    _community_hyperedge_ids.resize(hypergraph.initialNumEdges());
    _community_hyperedges.resize(hypergraph.initialNumEdges());
    tbb::parallel_for(ID(0), hypergraph.initialNumEdges(), [&](const HyperedgeID& he) {
      auto& e = hypergraph._hyperedges[he];
      if ( !e.isDisabled() ) {
        // Sort pins of hyperedge in increasing order of their community ids
        size_t incidence_array_start = e.firstEntry();
        size_t incidence_array_end = e.firstInvalidEntry();
        std::sort(hypergraph._incidence_array.begin() + incidence_array_start,
                  hypergraph._incidence_array.begin() + incidence_array_end,
                  [&](const HypernodeID& lhs, const HypernodeID& rhs) {
          return get_community_id(lhs) < get_community_id(rhs);
        });

        // Add community hyperedges for each consecutive range of pins with
        // the same community id
        size_t last_community_start = incidence_array_start;
        PartitionID last_community_id = get_community_id(
          hypergraph._incidence_array[last_community_start]);
        for (size_t incidence_array_pos = incidence_array_start + 1;
            incidence_array_pos < incidence_array_end;
            ++incidence_array_pos) {
          const HypernodeID pin = hypergraph._incidence_array[incidence_array_pos];
          const PartitionID community_id = get_community_id(pin);
          if (community_id != last_community_id) {
            add_community_hyperedge(he, last_community_id, last_community_start, incidence_array_pos, e.weight());
            last_community_start = incidence_array_pos;
            last_community_id = community_id;
          }
        }
        add_community_hyperedge(he, last_community_id, last_community_start, incidence_array_end, e.weight());

        HEAVY_COARSENING_ASSERT([&] {
          if (e.firstEntry() != _community_hyperedges[he][0].firstEntry()) {
            return false;
          }
          for (size_t i = 1; i < _community_hyperedges[he].size(); ++i) {
            if (_community_hyperedges[he][i - 1].firstInvalidEntry() !=
                _community_hyperedges[he][i].firstEntry()) {
              return false;
            }
          }
          if (e.firstInvalidEntry() != _community_hyperedges[he].back().firstInvalidEntry()) {
            return false;
          }
          return true;
        } (), "Initialization of community hyperedges failed!");
      }
    });
    _are_community_hyperedges_initialized = true;
  }

  void removeCommunityHyperedges(const parallel::scalable_vector<HypernodeID>& contraction_index,
                                 const parallel::scalable_vector<Hypergraph>& hypergraphs) {
    ASSERT(_are_community_hyperedges_initialized);
    unused(contraction_index);
    unused(hypergraphs);
    if (!Hypergraph::is_static_hypergraph) {
      // TODO(heuer): implement removal of community hyperedges for dynamic hypergraph here
    }

    CommunitiesOfHyperedges tmp_community_hyperedge_ids;
    CommunityHyperedges tmp_community_hyperedges;
    _community_hyperedge_ids = std::move(tmp_community_hyperedge_ids);
    _community_hyperedges = std::move(tmp_community_hyperedges);
    _are_community_hyperedges_initialized = false;
  }

  // ! Copy community support in parallel
  CommunitySupport copy(const TaskGroupID) {
    CommunitySupport community_support;
    community_support._node = _node;
    community_support._is_initialized = _is_initialized;
    community_support._num_communities = _num_communities;
    community_support._are_community_hyperedges_initialized =
      _are_community_hyperedges_initialized;

    // Copy all members in parallel
    tbb::parallel_invoke([&] {
      community_support._communities_num_hypernodes.resize(
        _communities_num_hypernodes.size());
      memcpy(community_support._communities_num_hypernodes.data(),
        _communities_num_hypernodes.data(), sizeof(HypernodeID) * _communities_num_hypernodes.size());
    }, [&] {
      community_support._communities_num_pins.resize(
        _communities_num_pins.size());
      memcpy(community_support._communities_num_pins.data(),
        _communities_num_pins.data(), sizeof(HypernodeID) * _communities_num_pins.size());
    }, [&] {
      community_support._community_degree.resize(
        _community_degree.size());
      memcpy(community_support._community_degree.data(),
        _community_degree.data(), sizeof(HyperedgeID) * _community_degree.size());
    }, [&] {
      const size_t size = _community_hyperedge_ids.size();
      community_support._community_hyperedge_ids.resize(size);
      tbb::parallel_for(0UL, size, [&](const size_t i) {
        const size_t he_size = _community_hyperedge_ids[i].size();
        community_support._community_hyperedge_ids[i].resize(he_size);
        memcpy(community_support._community_hyperedge_ids[i].data(),
          _community_hyperedge_ids[i].data(), sizeof(PartitionID) * he_size);
      });
    }, [&] {
      const size_t size = _community_hyperedges.size();
      community_support._community_hyperedges.resize(size);
      tbb::parallel_for(0UL, size, [&](const size_t i) {
        const size_t he_size = _community_hyperedges[i].size();
        community_support._community_hyperedges[i].resize(he_size);
        memcpy(community_support._community_hyperedges[i].data(),
          _community_hyperedges[i].data(), sizeof(CommunityHyperedge) * he_size);
      });
    }, [&] {
      community_support._vertex_to_community_node_id.resize(
        _vertex_to_community_node_id.size());
      memcpy(community_support._vertex_to_community_node_id.data(),
        _vertex_to_community_node_id.data(), sizeof(HypernodeID) *
        _vertex_to_community_node_id.size());
    });

    return community_support;
  }

  // Copy community support sequential
  CommunitySupport copy() {
    CommunitySupport community_support;
    community_support._node = _node;
    community_support._is_initialized = _is_initialized;
    community_support._num_communities = _num_communities;

    community_support._communities_num_hypernodes.resize(
      _communities_num_hypernodes.size());
    memcpy(community_support._communities_num_hypernodes.data(),
      _communities_num_hypernodes.data(), sizeof(HypernodeID) * _num_communities);
    community_support._communities_num_pins.resize(
      _communities_num_pins.size());
    memcpy(community_support._communities_num_pins.data(),
      _communities_num_pins.data(), sizeof(HypernodeID) * _num_communities);
    community_support._community_degree.resize(
      _community_degree.size());
    memcpy(community_support._community_degree.data(),
      _community_degree.data(), sizeof(HyperedgeID) * _num_communities);

    community_support._are_community_hyperedges_initialized =
      _are_community_hyperedges_initialized;
    community_support._community_hyperedge_ids.resize(_community_hyperedge_ids.size());
    std::copy(_community_hyperedge_ids.begin(), _community_hyperedge_ids.end(),
      community_support._community_hyperedge_ids.begin());
    community_support._community_hyperedges.resize(_community_hyperedges.size());
    std::copy(_community_hyperedges.begin(), _community_hyperedges.end(),
      community_support._community_hyperedges.begin());

    community_support._vertex_to_community_node_id.resize(
      _vertex_to_community_node_id.size());
    memcpy(community_support._vertex_to_community_node_id.data(),
      _vertex_to_community_node_id.data(), sizeof(HypernodeID) *
      _vertex_to_community_node_id.size());

    return community_support;
  }

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    parent->addChild("Num Hypernodes per Community",
      sizeof(HypernodeID) * _communities_num_hypernodes.size());
    parent->addChild("Num Pins per Community",
      sizeof(HypernodeID) * _communities_num_pins.size());
    parent->addChild("Community Degree",
      sizeof(HyperedgeID) * _community_degree.size());

    size_t size_community_hyperedge_ids = 0;
    for ( const auto& community_he_ids : _community_hyperedge_ids ) {
      size_community_hyperedge_ids += sizeof(PartitionID) * community_he_ids.size();
    }
    parent->addChild("Community IDs of HEs", size_community_hyperedge_ids);

    size_t size_community_hyperedges = 0;
    for ( const auto& community_hes : _community_hyperedges ) {
      size_community_hyperedges += sizeof(CommunityHyperedge) * community_hes.size();
    }
    parent->addChild("Community Hyperedges of HEs", size_community_hyperedges);

    parent->addChild("Vertex to Community Node ID",
      sizeof(HypernodeID) * _vertex_to_community_node_id.size());
  }

 private:
  void computeNumberOfCommunities(const Hypergraph& hypergraph,
                                  const int node = 0) {
    // The number of communities is the maximum community id plus 1
    _num_communities = tbb::parallel_reduce(tbb::blocked_range<HypernodeID>(ID(0), hypergraph.initialNumNodes()),
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

  // ! Accessor for community hyperedge-related information
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const CommunityHyperedge& community_hyperedge(const HyperedgeID e, const PartitionID community_id) const {
    const HypernodeID local_id = common::get_local_position_of_edge(e);
    ASSERT(local_id < _community_hyperedges.size(), "Hyperedge" << e << "does not exist");
    ASSERT(_node == common::get_numa_node_of_edge(e), "Hyperedge" << e << "is not part of NUMA node" << _node);

    size_t community_hyperedge_pos = 0;
    for ( ; community_hyperedge_pos < _community_hyperedges[local_id].size(); ++community_hyperedge_pos ) {
      if ( _community_hyperedges[local_id][community_hyperedge_pos].communityID() == community_id ) {
        break;
      }
    }

    ASSERT(community_hyperedge_pos < _community_hyperedges[local_id].size(),
           "Community hyperedge" << e << "with community id" << community_id << "not found");
    return _community_hyperedges[local_id][community_hyperedge_pos];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE CommunityHyperedge& community_hyperedge(const HyperedgeID e, const PartitionID community_id) {
    return const_cast<CommunityHyperedge&>(static_cast<const CommunitySupport&>(*this).community_hyperedge(e, community_id));
  }

  // ! NUMA node over which this class is constructed
  int _node;
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

  // ! Indicates, if community hyperedges are initialized
  bool _are_community_hyperedges_initialized;
  // ! Community Ids contained in a hyperedge
  CommunitiesOfHyperedges _community_hyperedge_ids;
  // ! For each hyperedge this structure contains all community hyperedges
  CommunityHyperedges _community_hyperedges;

  // ! Mapping that maps a vertex u to a unique id in the range [0,|C_u|)
  // ! where C_v is the community id of vertex u
  parallel::scalable_vector<HypernodeID> _vertex_to_community_node_id;

};

} // namespace ds
} // namespace mt_kahypar
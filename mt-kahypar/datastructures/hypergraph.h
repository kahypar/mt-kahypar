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

#include <algorithm>
#include <type_traits>
#include <chrono>
#include <functional>
#include <set>

#include "tbb/enumerable_thread_specific.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/streaming_hypergraph.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
namespace ds {

template <typename HypernodeType_ = Mandatory,
          typename HyperedgeType_ = Mandatory,
          typename HypernodeWeightType_ = Mandatory,
          typename HyperedgeWeightType_ = Mandatory,
          typename PartitionIDType_ = Mandatory,
          typename HardwareTopology = Mandatory,
          typename TBBNumaArena = Mandatory>
class Hypergraph {

 private:
  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  using HypernodeID = HypernodeType_;
  using HyperedgeID = HyperedgeType_;
  using HypernodeWeight = HypernodeWeightType_;
  using HyperedgeWeight = HyperedgeWeightType_;
  using PartitionID = PartitionIDType_;

  static constexpr PartitionID kInvalidPartition = -1;
  static constexpr HyperedgeID kInvalidHyperedge = std::numeric_limits<HyperedgeID>::max();


 public:
  using StreamingHypergraph = mt_kahypar::ds::StreamingHypergraph<HypernodeID,
                                                                  HyperedgeID,
                                                                  HypernodeWeight,
                                                                  HyperedgeWeight,
                                                                  PartitionID,
                                                                  HardwareTopology,
                                                                  TBBNumaArena>;

 private:
  using HypernodeIterator = typename StreamingHypergraph::HypernodeIterator;
  using HyperedgeIterator = typename StreamingHypergraph::HyperedgeIterator;
  using IncidenceIterator = typename StreamingHypergraph::IncidenceIterator;
  using ConnectivitySetIterator = typename StreamingHypergraph::ConnectivitySetIterator;
  using CommunityIterator = typename StreamingHypergraph::CommunityIterator;
  using Memento = typename StreamingHypergraph::Memento;

  // ! Generic function that will be called if a hypernode v moves from a block from to a block to for
  // ! each incident net of the moved vertex v.
  // ! It will be called with the following arguments
  // !  1.) Hyperedge Weight
  // !  2.) Hyperedge Size
  // !  3.) Pin count in block from after move
  // !  4.) Pin count in block to after move
  // ! This function can be used to compute e.g. the delta in cut or km1 metric after a move
  using DeltaFunction = std::function<void(const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID)>;
  #define NOOP_FUNC [](const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }

  using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

  /*!
   * Iterator for HypergraphElements (Hypernodes/Hyperedges)
   *
   * The iterator is used in for-each loops over all hypernodes/hyperedges.
   * In order to support iteration over coarsened hypergraphs, this iterator
   * skips over HypergraphElements marked as invalid.
   * Iterating over the set of vertices \f$V\f$ therefore is linear in the
   * size \f$|V|\f$ of the original hypergraph - even if it has been coarsened
   * to much smaller size. The same also holds true for for-each loops over
   * the set of hyperedges.
   *
   * In order to be as generic as possible, the iterator does not expose the
   * internal Hypernode/Hyperedge representations. Instead only handles to
   * the respective elements are returned, i.e. the IDs of the corresponding
   * hypernodes/hyperedges.
   *
   */
  template <typename Iterator>
  class GlobalHypergraphElementIterator :
    public std::iterator<std::forward_iterator_tag,    // iterator_category
                         typename Iterator::IDType,   // value_type
                         std::ptrdiff_t,   // difference_type
                         const typename Iterator::IDType*,   // pointer
                         typename Iterator::IDType>{   // reference
    using IDType = typename Iterator::IDType;
    using Iterators = std::vector<std::pair<Iterator, Iterator>>;

   public:
    /*!
     * Construct a GlobalHypergraphElementIterator
     * See Hypergraph::nodes() or Hypergraph::edges() for usage.
     *
     * \param iterators container of iterators (see StreamingHypergraph::HypergraphElementIterator)
     */
    GlobalHypergraphElementIterator(Iterators&& iterators) :
      _iterators(std::move(iterators)),
      _idx(0) {
      ASSERT(_iterators.size() > 0);
      while ( _idx < _iterators.size() - 1 &&
              *_iterators[_idx].first == *_iterators[_idx].second  ) {
        ++_idx;
      }
    }

    // ! Returns the id of the element the iterator currently points to.
    IDType operator* () const {
      return *_iterators[_idx].first;
    }

    // ! Prefix increment. The iterator advances to the next valid element.
    GlobalHypergraphElementIterator& operator++ () {
      ++_iterators[_idx].first;
      while ( *_iterators[_idx].first == *_iterators[_idx].second &&
              _idx < _iterators.size() - 1 ) {
        ++_idx;
      }
      return *this;
    }

    // ! Postfix increment. The iterator advances to the next valid element.
    GlobalHypergraphElementIterator operator++ (int) {
      GlobalHypergraphElementIterator copy = *this;
      operator++ ();
      return copy;
    }

    // ! Convenience function for range-based for-loops
    friend GlobalHypergraphElementIterator end<>(const std::pair<GlobalHypergraphElementIterator,
                                                           GlobalHypergraphElementIterator>& iter_pair);
    // ! Convenience function for range-based for-loops
    friend GlobalHypergraphElementIterator begin<>(const std::pair<GlobalHypergraphElementIterator,
                                                             GlobalHypergraphElementIterator>& iter_pair);

    bool operator!= (const GlobalHypergraphElementIterator& rhs) {
      return *_iterators[_idx].first != *rhs._iterators[rhs._idx].first;
    }

   private:
    Iterators _iterators;
    size_t _idx;
  };

  /*!
   * For each block \f$V_i\f$ of the \f$k\f$-way partition \f$\mathrm{\Pi} = \{V_1, \dots, V_k\}\f$,
   * a PartInfo object stores the number of hypernodes currently assigned to block \f$V_i\f$
   * as well as the sum of their weights.
   */
  class PartInfo {
    public:
      bool operator== (const PartInfo& other) const {
        return weight == other.weight && size == other.size;
      }

      HypernodeWeight weight;
      int64_t size;
  };


  /**
   * Each thread contains its local part weight and size information. If a hypernode changes
   * its block, the modification to part weights and sizes are only applied to
   * the block weights of the calling thread. The part weights are stored
   * relative to the initial global block weights as a delta. The current block weights
   * can be calculated by summing up the deltas of all threads (relative to the initial
   * global block weight).
   */
  class ThreadPartInfos {

    public:
      static ThreadPartInfos construct( const PartitionID k, const std::vector<PartInfo>& global ) {
        return ThreadPartInfos(k, global);
      }

      void apply( const PartitionID id, const PartInfo& delta ) {
        ASSERT(id >= 0 && id < _k, V(id) << V(_k));
        _delta[id].weight += delta.weight;
        _delta[id].size += delta.size;
        _current[id].weight += delta.weight;
        _current[id].size += delta.size;
      }

      /**
       * Appliying the deltas of all threads to local block weights and sizes.
       */
      void snapshot( const tbb::enumerable_thread_specific<ThreadPartInfos>& infos ) {
        // Reset current block weights
        for ( PartitionID k = 0; k < _k; ++k ) {
          _current[k] = _global[k];
        }

        // Applying deltas of all threads
        for ( const ThreadPartInfos& thread_info : infos ) {
          // It can happen that in some situations (when frequently updating local part
          // weights) that the current thread info is initialized and therefore iterating
          // over all k's would fail
          if (thread_info._delta.size() == (size_t) _k ) {
            for ( PartitionID k = 0; k < _k; ++k ) {
              _current[k].weight += thread_info._delta[k].weight;
              _current[k].size += thread_info._delta[k].size;
            }
          }
        }
      }

      const std::vector<PartInfo>& delta() const {
        return _delta;
      }

      HyperedgeWeight weight(const PartitionID id) const {
        ASSERT(id >= 0 && id < _k);
        return _current[id].weight;
      }

      size_t size(const PartitionID id) const {
        ASSERT(id >= 0 && id < _k);
        return _current[id].size;
      }

      void reset() {
        for ( PartitionID k = 0; k < _k; ++k ) {
          _current[k] = _global[k];
          _delta[k] = PartInfo {0, 0};
        }
      }

    private:
      ThreadPartInfos(const PartitionID k, const std::vector<PartInfo>& global) :
        _k(k),
        _global(global),
        _delta(global.size()),
        _current(global) { }

      // ! Number of blocks
      const PartitionID _k;
      // ! Initial global part infos
      const std::vector<PartInfo>& _global;
      // ! Delta relative to initial global part infos (changes performed by local thread)
      std::vector<PartInfo> _delta;
      // ! Represents current part infos for local thread
      std::vector<PartInfo> _current;
  };

  // ! TBB Thread Local Storage
  using ThreadLocalPartInfos = tbb::enumerable_thread_specific<ThreadPartInfos>;

  // ! Iterator to iterator over the hypernodes
  using GlobalHypernodeIterator = GlobalHypergraphElementIterator<HypernodeIterator>;
  // ! Iterator to iterator over the hyperedges
  using GlobalHyperedgeIterator = GlobalHypergraphElementIterator<HyperedgeIterator>;

 public:

  constexpr static size_t kEdgeHashSeed = StreamingHypergraph::kEdgeHashSeed;

  // ! Empty Hypergraph
  explicit Hypergraph() :
    _num_hypernodes(0),
    _num_hyperedges(0),
    _num_pins(0),
    _num_communities(0),
    _k(0),
    _communities_num_hypernodes(),
    _communities_num_pins(),
    _part_info(),
    _local_part_info([&] { return ThreadPartInfos::construct(_k, _part_info); }),
    _hypergraphs(),
    _node_mapping(),
    _edge_mapping(),
    _community_node_mapping() { }

  // ! Constructs a hypergraph based on the given numa hypergraphs
  // ! and additionaly computes a vertex to numa node mapping
  // ! based on the current hyperedge distribution
  Hypergraph(const HypernodeID num_hypernodes,
             std::vector<StreamingHypergraph>&& hypergraphs,
             PartitionID k) :
    _num_hypernodes(num_hypernodes),
    _num_hyperedges(0),
    _num_pins(0),
    _num_communities(0),
    _k(k),
    _communities_num_hypernodes(),
    _communities_num_pins(),
    _part_info(k),
    _local_part_info([&] { return ThreadPartInfos::construct(_k, _part_info); }),
    _hypergraphs(std::move(hypergraphs)),
    _node_mapping(num_hypernodes, 0),
    _edge_mapping(),
    _community_node_mapping() {
    computeNodeMapping();
    initializeHypernodes();
  }
  // ! Constructs a hypergraph based on the given numa hypergraphs
  // ! and node mapping
  Hypergraph(const HypernodeID num_hypernodes,
             std::vector<StreamingHypergraph>&& hypergraphs,
             std::vector<HypernodeID>&& node_mapping,
             PartitionID k) :
    _num_hypernodes(num_hypernodes),
    _num_hyperedges(0),
    _num_pins(0),
    _num_communities(0),
    _k(k),
    _communities_num_hypernodes(),
    _communities_num_pins(),
    _part_info(k),
    _local_part_info([&] { return ThreadPartInfos::construct(_k, _part_info); }),
    _hypergraphs(std::move(hypergraphs)),
    _node_mapping(std::move(node_mapping)),
    _edge_mapping(),
    _community_node_mapping() {
    initializeHypernodes();
  }

  Hypergraph(const Hypergraph&) = delete;
  Hypergraph& operator= (const Hypergraph&) = delete;

  Hypergraph(Hypergraph&&) = default;

  Hypergraph& operator= (Hypergraph&& other) {
    _num_hypernodes = other._num_hypernodes;
    _num_hyperedges = other._num_hyperedges;
    _num_pins = other._num_pins;
    _num_communities = other._num_communities;
    _k = other._k;
    _communities_num_hypernodes = std::move(other._communities_num_hypernodes);
    _communities_num_pins = std::move(other._communities_num_hypernodes);
    _part_info = std::move(other._part_info);
    _local_part_info = ThreadLocalPartInfos([&] { return ThreadPartInfos::construct(_k, _part_info); });
    _hypergraphs = std::move(other._hypergraphs);
    _node_mapping = std::move(other._node_mapping);
    _edge_mapping = std::move(other._edge_mapping);
    _community_node_mapping = std::move(other._community_node_mapping);
    return *this;
  }

  // ####################### General Hypergraph Stats #######################

  // ! Initial number of hypernodes
  HypernodeID initialNumNodes() const {
    return _num_hypernodes;
  }

  // ! Initial number of hypernodes on numa node
  HypernodeID initialNumNodes(const int node) const {
    ASSERT(node < (int)_hypergraphs.size());
    return _hypergraphs[node].initialNumNodes();
  }

  // ! Initial number of hyperedges
  HyperedgeID initialNumEdges() const {
    return _num_hyperedges;
  }

  // ! Initial number of hyperedges on numa node
  HyperedgeID initialNumEdges(const int node) const {
    ASSERT(node < (int)_hypergraphs.size());
    return _hypergraphs[node].initialNumEdges();
  }

  // ! Initial number of pins
  HypernodeID initialNumPins() const {
    return _num_pins;
  }

  // ! Initial number of pins on numa node
  HypernodeID initialNumPins(const int node) const {
    ASSERT(node < (int)_hypergraphs.size());
    return _hypergraphs[node].initialNumPins();
  }

  // ! Total weight of hypergraph
  HypernodeWeight totalWeight() const {
    HypernodeWeight weight = 0;
    for ( const StreamingHypergraph& hypergraph : _hypergraphs ) {
      weight += hypergraph.totalWeight();
    }
    return weight;
  }

  // ! Recomputes the total weight of the hypergraph (in parallel)
  void updateTotalWeight() {
    TBBNumaArena::instance().execute_on_all_numa_nodes([&](const int node) {
      _hypergraphs[node].updateTotalWeight();
    });
  }

  // ! Number of blocks this hypergraph is partitioned into
  PartitionID k() const {
    return _k;
  }


  // ####################### Iterators #######################

  // ! Returns an iterator over the set of active nodes of the hypergraph
  std::pair<GlobalHypernodeIterator, GlobalHypernodeIterator> nodes() const {
    ASSERT(_hypergraphs.size() > 0);
    std::vector<std::pair<HypernodeIterator, HypernodeIterator>> iterators;
    std::vector<std::pair<HypernodeIterator, HypernodeIterator>> end;
    for ( size_t node = 0; node < _hypergraphs.size(); ++node ) {
      iterators.emplace_back(_hypergraphs[node].nodes());
    }
    size_t last = iterators.size() - 1;
    end.emplace_back(std::make_pair(iterators[last].second, iterators[last].second));
    return std::make_pair(GlobalHypernodeIterator(std::move(iterators)),
                          GlobalHypernodeIterator(std::move(end)));
  }

  // ! Returns an iterator over the set of active nodes of the hypergraph on a numa node
  std::pair<HypernodeIterator, HypernodeIterator>  nodes(const int node) const {
    ASSERT(node < (int) _hypergraphs.size());
    return _hypergraphs[node].nodes();
  }

  // ! Returns an iterator over the set of active edges of the hypergraph
  std::pair<GlobalHyperedgeIterator, GlobalHyperedgeIterator> edges() const {
    ASSERT(_hypergraphs.size() > 0);
    std::vector<std::pair<HyperedgeIterator, HyperedgeIterator>> iterators;
    std::vector<std::pair<HyperedgeIterator, HyperedgeIterator>> end;
    for ( size_t node = 0; node < _hypergraphs.size(); ++node ) {
      iterators.emplace_back(_hypergraphs[node].edges());
    }
    size_t last = iterators.size() - 1;
    end.emplace_back(std::make_pair(iterators[last].second, iterators[last].second));
    return std::make_pair(GlobalHyperedgeIterator(std::move(iterators)),
                          GlobalHyperedgeIterator(std::move(end)));
  }

  // ! Returns an iterator over the set of active nodes of the hypergraph on a numa node
  std::pair<HyperedgeIterator, HyperedgeIterator> edges(const int node) const {
    ASSERT(node < (int) _hypergraphs.size());
    return _hypergraphs[node].edges();
  }


  /*!
   * Illustration for different incidentEdges iterators:
   *
   * Structure of incident nets for a vertex u during community coarsening:
   *
   * | <-- single-pin community hyperedges --> | <--   valid hyperedges     --> | <-- invalid hyperedges --> |
   *
   *                                            <-- validIncidentEdges(u,c) -->
   *
   *  <-------------------------- incidentEdges(u,c) ------------------------->
   *
   *  <----------------------------------------- incidentEdges(u) ------------------------------------------>
   */


  // ! Returns a for-each iterator-pair to loop over the set of incident hyperedges of hypernode u.
  // ! During parallel community coarsening we do not remove parallel hyperedges (which spans more than
  // ! one community) incident to hypernode u. Instead, we just invalidate those hyperedges
  // ! (which means swapping them to the end of incident nets and storing a pointer to all
  // ! invalidated hyperedges). This function returns an iterator over all VALID and INVALID
  // ! hyperedges incident of vertex u.
  std::pair<IncidenceIterator, IncidenceIterator> incidentEdges(const HypernodeID u) const {
    return hypergraph_of_vertex(u).incidentEdges(u);
  }

  // ! Returns a for-each iterator-pair to loop over the set of incident hyperedges of hypernode u.
  // ! During parallel community coarsening we do not remove parallel hyperedges (which spans more than
  // ! one community) incident to hypernode u. Instead, we just invalidate those hyperedges
  // ! (which means swapping them to the end of incident nets and storing a pointer to all
  // ! invalidated hyperedges). This function returns an iterator over all VALID hyperedges incident
  // ! of vertex u.
  std::pair<IncidenceIterator, IncidenceIterator> validIncidentEdges(const HypernodeID u, const PartitionID community_id) const {
    return hypergraph_of_vertex(u).validIncidentEdges(u, community_id);
  }

  // ! Returns a for-each iterator-pair to loop over the set of incident hyperedges of hypernode u.
  // ! Note, this function returns the same incident edges than validIncidentEdges(u, community_id), but
  // ! additionally it skips all hyperedges that only contains a single vertex of a community (single-pin
  // ! community hyperedges)
  std::pair<IncidenceIterator, IncidenceIterator> incidentEdges(const HypernodeID u, const PartitionID community_id) const {
    return hypergraph_of_vertex(u).incidentEdges(u, community_id);
  }

  // ! Returns a for-each iterator-pair to loop over the set pins of hyperedge e.
  // ! Note, this function fails if community hyperedges are initialized.
  std::pair<IncidenceIterator, IncidenceIterator> pins(const HyperedgeID e) const {
    return hypergraph_of_edge(e).pins(e);
  }

  // ! Returns a for-each iterator-pair to loop over the set pins of hyperedge e that belongs to a certain community.
  // ! Note, this function fails if community hyperedges are not initialized.
  std::pair<IncidenceIterator, IncidenceIterator> pins(const HyperedgeID e, const PartitionID community_id) const {
    return hypergraph_of_edge(e).pins(e, community_id);
  }

  // ! Returns a for-each iterator-pair to loop over the set of communities contained in hyperedge e.
  std::pair<CommunityIterator, CommunityIterator> communities(const HyperedgeID e) const {
    return hypergraph_of_edge(e).communities(e);
  }

  // ! Returns a for-each iterator-pair to loop over the set of block ids contained in hyperedge e.
  std::pair<ConnectivitySetIterator, ConnectivitySetIterator> connectivitySet(const HyperedgeID e) const {
    return hypergraph_of_edge(e).connectivitySet(e);
  }

  // ####################### Hypernode Information #######################

  // ! Returns for a vertex of the hypergraph its original vertex id
  // ! defined in the input file
  HypernodeID originalNodeID(const HypernodeID u) const {
    return hypergraph_of_vertex(u).originalNodeId(u);
  }

  // ! Reverse operation of originalNodeID(u)
  HypernodeID globalNodeID(const HypernodeID u) const {
    ASSERT(u < _node_mapping.size());
    return _node_mapping[u];
  }

  // ! Weight of a vertex
  HypernodeWeight nodeWeight(const HypernodeID u) const {
    return hypergraph_of_vertex(u).nodeWeight(u);
  }

  // ! Sets the weight of a vertex
  void setNodeWeight(const HypernodeID u, const HypernodeWeight weight) {
    hypergraph_of_vertex(u).setNodeWeight(u, weight);
  }

  // ! Degree of a hypernode
  HyperedgeID nodeDegree(const HypernodeID u) const {
    return hypergraph_of_vertex(u).nodeDegree(u);
  }

  // ! Returns, whether a hypernode is enabled or not
  bool nodeIsEnabled(const HypernodeID u) const {
    return hypergraph_of_vertex(u).nodeIsEnabled(u);
  }

  // ! Enables a hypernode (must be disabled before)
  void enableHypernode(const HypernodeID u) {
    hypergraph_of_vertex(u).enableHypernode(u);
  }

  // ! Disable a hypernode (must be enabled before)
  void disableHypernode(const HypernodeID u) {
    hypergraph_of_vertex(u).disableHypernode(u);
  }

  // ####################### Hyperedge Information #######################

  // ! Returns for a hyperedge of the hypergraph its original hyperedge id
  // ! defined in the input file
  HypernodeID originalEdgeID(const HyperedgeID e) const {
    return hypergraph_of_edge(e).originalEdgeId(e);
  }

  // ! Reverse operation of originalEdgeID(e)
  HypernodeID globalEdgeID(const HyperedgeID e) const {
    ASSERT(e < _edge_mapping.size());
    return _edge_mapping[e];
  }

  // ! Weight of a hyperedge
  HypernodeWeight edgeWeight(const HyperedgeID e) const {
    return hypergraph_of_edge(e).edgeWeight(e);
  }

  // ! Sets the weight of a hyperedge
  void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
    hypergraph_of_edge(e).setEdgeWeight(e, weight);
  }

  // ! Number of pins of a hyperedge
  HypernodeID edgeSize(const HyperedgeID e) const {
    return hypergraph_of_edge(e).edgeSize(e);
  }

  // ! Hash value defined over the pins of a hyperedge
  size_t edgeHash(const HyperedgeID e) const {
    return hypergraph_of_edge(e).edgeHash(e);
  }

  // ! Returns, whether a hyperedge is enabled or not
  bool edgeIsEnabled(const HyperedgeID e) const {
    return hypergraph_of_edge(e).edgeIsEnabled(e);
  }

  // ! Enables a hyperedge (must be disabled before)
  void enableHyperedge(const HyperedgeID e) {
    hypergraph_of_edge(e).enableHyperedge(e);
  }

  // ! Disabled a hyperedge (must be enabled before)
  void disableHyperedge(const HyperedgeID e) {
    hypergraph_of_edge(e).disableHyperedge(e);
  }

  // ####################### Community Hyperedge Information #######################

  // ! Weight of a community hyperedge
  HypernodeWeight edgeWeight(const HyperedgeID e, const PartitionID community_id) const {
    return hypergraph_of_edge(e).edgeWeight(e, community_id);
  }

  // ! Sets the weight of a community hyperedge
  void setEdgeWeight(const HyperedgeID e, const PartitionID community_id, const HyperedgeWeight weight) {
    hypergraph_of_edge(e).setEdgeWeight(e, community_id, weight);
  }

  // ! Number of pins of a hyperedge that are assigned to a community
  HypernodeID edgeSize(const HyperedgeID e, const PartitionID community_id) const {
    return hypergraph_of_edge(e).edgeSize(e, community_id);
  }

  // ! Hash value defined over the pins of a hyperedge that belongs to a community
  size_t edgeHash(const HyperedgeID e, const PartitionID community_id) const {
    return hypergraph_of_edge(e).edgeHash(e, community_id);
  }

  // ####################### Community Information #######################

  // ! Number of communities
  PartitionID numCommunities() const {
    return _num_communities;
  }

  // ! Community id which hypernode u is assigned to
  PartitionID communityID(const HypernodeID u) const {
    return hypergraph_of_vertex(u).communityID(u);
  }

  // ! Assign a community to a hypernode
  // ! Note, in order to use all community-related functions, initializeCommunities()
  // ! have to be called after assigning to each vertex a community id
  void setCommunityID(const HypernodeID hn, const PartitionID community_id) {
    hypergraph_of_vertex(hn).setCommunityID(hn, community_id);
  }

  // ! Consider hypernode u is part of community C = {v_1, ..., v_n},
  // ! than this function returns a unique id for hypernode u in the
  // ! range [0,n).
  HypernodeID communityNodeId(const HypernodeID u) const {
    return hypergraph_of_vertex(u).communityNodeId(u);
  }

  // ! Number of hypernodes in community
  HypernodeID numCommunityHypernodes(const PartitionID community) const {
    ASSERT(community < (PartitionID) _communities_num_hypernodes.size());
    return _communities_num_hypernodes[community];
  }

  // ! Number of pins in community
  HypernodeID numCommunityPins(const PartitionID community) const {
    ASSERT(community < (PartitionID) _communities_num_pins.size());
    return _communities_num_pins[community];
  }

  // ! Number of communities which pins of hyperedge belongs to
  size_t numCommunitiesInHyperedge(const HyperedgeID e) const {
    return hypergraph_of_edge(e).numCommunitiesInHyperedge(e);
  }

  // ! Numa node to which community is assigned to
  PartitionID communityNumaNode(const PartitionID community_id) const {
    ASSERT(community_id < (PartitionID) _community_node_mapping.size());
    return _community_node_mapping[community_id];
  }

  // ! Sets the community to numa node mapping
  void setCommunityNodeMapping(std::vector<PartitionID>&& community_node_mapping) {
    _community_node_mapping = std::move(community_node_mapping);
  }

  // ####################### Partition Information #######################

  /**
   * About moving a vertex:
   * The functions setNodePart(...) and changeNodePart(...) are implemented thread-safe
   * and lock-free. All involved data structures are based on atomics. We maintain
   * the following information when moving a node:
   *  1.) Block id of a vertex
   *  2.) Pin count in a block of a hyperedge
   *  3.) Connectivity of a hyperedge (number of blocks to which the pins of a hyperedge
   *      belongs to)
   *  4.) Connectivity set of a hyperedge (blocks to which the pins of a hyperedge belongs
   *      to)
   *  5.) Number of incident cut hyperedges of a vertex (number with hyperedges with
   *      connectivity greater than 1)
   *  6.) Block weights and sizes.
   *
   * In order, to change the block id of a vertex, we start by performing a CAS operation
   * on the block id of a vertex. This prevents that two threads are moving one vertex to
   * different blocks concurrently. If the operation succeeds we start updating all involved
   * data structures. Otherwise, we abort the operation and return false (indicating that
   * the move of the vertex failed). On succees, we update the local block weights of the
   * calling thread (see description below) and update 2.) - 5.) by iterating over each
   * hyperedge and increment or decrement the pin count of a block of a hyperedge. Both
   * operations return the pin count after the atomic update and based on that we can decide
   * if the hyperedge is still cut or not.
   *
   * It is guaranteed, that if all parallel vertex move operations are finished, the values
   * for 1.) - 5.) reflects the current state of the hypergraph. However, this is not the case
   * in a parallel setting. Since the update of all incident edges of the moved vertex is not
   * performed in a transactional/exclusive fashion, it can happen that threads that read one
   * of the values 2.) to 5.) have an intermediate view on those. Algorithms that
   * build on top of that have to consider that behavior and resolve or deal with that issue on
   * a higher layer.
   *
   * Another feature of changeNodePart is that someone can pass a generic delta function to compute
   * e.g. the delta in the cut- or km1-metric when moving vertices in parallel. E.g. one
   * can pass the following function to changeNodePart (as delta_func) to compute the delta
   * for the cut-metric of the move:
   *
   * HyperedgeWeight delta = 0;
   * auto cut_delta = [&](const HyperedgeWeight edge_weight,
   *                             const HypernodeID edge_size,
   *                             const HypernodeID pin_count_in_from_part_after,
   *                             const HypernodeID pin_count_in_to_part_after) {
   *   delta += mt_kahypar::Hypergraph::cutDelta(
   *     edge_weight, edge_size, pin_count_in_from_part_after, pin_count_in_to_part_after);
   * };
   * hypergraph.changeNodePart(hn, from, to, cut_delta); // 'delta' will contain afterwards the delta
   *                                                     // for the cut metric
   *
   * Summing up all deltas of all parallel moves will be delta after the parallel execution phase.
   *
   */

  // ! Sets the block id an unassigned vertex u.
  // ! Returns true, if assignment of unassigned vertex u to block id succeeds.
  bool setNodePart(const HypernodeID u, const PartitionID id) {
    StreamingHypergraph& hypergraph_of_u = hypergraph_of_vertex(u);
    ASSERT(id < _k && id != kInvalidPartition, "Part ID" << id << "is invalid");

    // Sets the node part of vertex u to id. The operation succeeds
    // if CAS operation on part_id of vertex u succeeds
    if ( hypergraph_of_u.setNodePart(u, id) ) {
      // Update local block weights of calling thread
      _local_part_info.local().apply(id, PartInfo{ nodeWeight(u), 1 });

      for ( const HyperedgeID& he : incidentEdges(u) ) {
        hypergraph_of_edge(he).incrementPinCountInPart(he, id);
      }
      return true;
    }

    return false;
  }

  // ! Changes the block id of vertex u from block 'from' to block 'to'
  // ! Returns true, if move of vertex u to corresponding block succeeds.
  bool changeNodePart(const HypernodeID u,
                      const PartitionID from,
                      const PartitionID to,
                      const DeltaFunction& delta_func = NOOP_FUNC) {
    StreamingHypergraph& hypergraph_of_u = hypergraph_of_vertex(u);
    ASSERT(to < _k && to != kInvalidPartition, "Part ID" << to << "is invalid");

    // Changes the node part of vertex u to id. The operation succeeds
    // if CAS operation on part_id of vertex u succeeds
    if ( hypergraph_of_u.changeNodePart(u, from, to) ) {
      // Update local block weights of calling thread
      _local_part_info.local().apply(from, PartInfo{ -nodeWeight(u), -1 });
      _local_part_info.local().apply(to, PartInfo{ nodeWeight(u), 1 });

      for ( const HyperedgeID& he : incidentEdges(u) ) {
        HyperedgeID pin_count_in_from_part_after = hypergraph_of_edge(he).decrementPinCountInPart(he, from);
        HyperedgeID pin_count_in_to_part_after = hypergraph_of_edge(he).incrementPinCountInPart(he, to);
        bool no_pins_left_in_source_part = pin_count_in_from_part_after == 0;
        bool only_one_pin_in_to_part = pin_count_in_to_part_after == 1;
        HypernodeID edge_size = edgeSize(he);

        delta_func(edgeWeight(he), edge_size, pin_count_in_from_part_after, pin_count_in_to_part_after);

        if ( no_pins_left_in_source_part && !only_one_pin_in_to_part &&
             pin_count_in_to_part_after == edge_size ) {
          // In that case, hyperedge he becomes an internal hyperedge
          for ( const HypernodeID& pin : pins(he) ) {
            hypergraph_of_vertex(pin).decrementIncidentNumCutHyperedges(pin);
          }
        } else if ( !no_pins_left_in_source_part && only_one_pin_in_to_part &&
                    pin_count_in_from_part_after == edge_size - 1 ) {
          // In that case, hyperedge he becomes an cut hyperede
          for ( const HypernodeID& pin : pins(he) ) {
            hypergraph_of_vertex(pin).incrementIncidentNumCutHyperedges(pin);
          }
        }
      }

      return true;
    }

    return false;
  }

  // ! Helper function to compute delta for cut-metric after changeNodePart
  static HyperedgeWeight cutDelta(const HyperedgeWeight edge_weight,
                                  const HypernodeID edge_size,
                                  const HypernodeID pin_count_in_from_part_after,
                                  const HypernodeID pin_count_in_to_part_after ) {
    if ( pin_count_in_to_part_after == edge_size ) {
      return -edge_weight;
    } else if ( pin_count_in_from_part_after == edge_size - 1 &&
                pin_count_in_to_part_after == 1 ) {
      return edge_weight;
    } else {
      return 0;
    }
  }

  // ! Helper function to compute delta for km1-metric after changeNodePart
  static HyperedgeWeight km1Delta(const HyperedgeWeight edge_weight,
                                  const HypernodeID,
                                  const HypernodeID pin_count_in_from_part_after,
                                  const HypernodeID pin_count_in_to_part_after ) {
    return ( pin_count_in_to_part_after == 1 ? edge_weight : 0 ) +
           ( pin_count_in_from_part_after == 0 ? -edge_weight : 0 );
  }

  // ! Block which vertex u belongs to
  PartitionID partID(const HypernodeID u) const {
    return hypergraph_of_vertex(u).partID(u);
  }

  // ! Returns, whether hypernode u is adjacent to a least one cut hyperedge.
  bool isBorderNode(const HypernodeID u) const {
    return hypergraph_of_vertex(u).isBorderNode(u);
  }

  // ! Number of incident cut hyperedges of vertex u
  HyperedgeID numIncidentCutHyperedges(const HypernodeID u) const {
    return hypergraph_of_vertex(u).numIncidentCutHyperedges(u);
  }

  // ! Initializes the number of cut hyperedges for each vertex
  // ! NOTE, this function have to be called after initial partitioning
  // ! and before local search.
  void initializeNumCutHyperedges() {
    for ( StreamingHypergraph& hypergraph : _hypergraphs ) {
      hypergraph.initializeNumCutHyperedges(_hypergraphs);
    }
  }

  // ! Number of blocks which pins of hyperedge e belongs to
  PartitionID connectivity(const HyperedgeID e) const {
    return hypergraph_of_edge(e).connectivity(e);
  }

  // ! Returns the number pins of hyperedge e that are part of block id
  HypernodeID pinCountInPart(const HyperedgeID e, const PartitionID id) const {
    return hypergraph_of_edge(e).pinCountInPart(e, id);
  }

  /**
   * About Block Weights:
   * We implement the concept of global and thread local block weights. If a vertex changes
   * its block, changes are only reflected on the local block weights of the calling thread. The local and
   * global block weights can be updated with updateLocalPartInfos() resp. updateGlobalPartInfos().
   * The rationale behind this is, that we want to prevent congestion on the block weights (when implemented
   * as atomic update), when multiple threads are moving vertices around. We leave it to the algorithm
   * implemented on top of these primitives to update local and global block weights frequently.
   */

  // ! Global weight of a block
  HypernodeWeight partWeight(const PartitionID id) const {
    ASSERT(id < _k && id != kInvalidPartition, "Part ID" << id << "is invalid");
    return _part_info[id].weight;
  }

  // ! Global size of a block
  size_t partSize(const PartitionID id) const {
    ASSERT(id < _k && id != kInvalidPartition, "Part ID" << id << "is invalid");
    return _part_info[id].size;
  }

  // ! Local weight of a block (of calling thread)
  HypernodeWeight localPartWeight(const PartitionID id) {
    ASSERT(id < _k && id != kInvalidPartition, "Part ID" << id << "is invalid");
    return _local_part_info.local().weight(id);
  }

  // ! Local size of a block (of calling thread)
  size_t localPartSize(const PartitionID id) {
    ASSERT(id < _k && id != kInvalidPartition, "Part ID" << id << "is invalid");
    return _local_part_info.local().size(id);
  }

  // ! Updates the local block weights of the calling thread
  void updateLocalPartInfos() {
    _local_part_info.local().snapshot(_local_part_info);
  }

  // ! Updates the global block weights
  // ! Note, this function is not thread safe and should be only called in
  // ! a single-threaded setting.
  void updateGlobalPartInfos() {
    for ( const ThreadPartInfos& thread_part_info : _local_part_info  ) {
      // Applying deltas of each local part information to global part information
      const std::vector<PartInfo>& delta = thread_part_info.delta();
      ASSERT(delta.size() == (size_t) _k);
      for ( PartitionID k = 0; k < _k; ++k ) {
        _part_info[k].weight += delta[k].weight;
        _part_info[k].size += delta[k].size;
      }
    }

    // Reset local block weights
    for ( ThreadPartInfos& thread_part_info : _local_part_info  ) {
      thread_part_info.reset();
    }
  }

  // ####################### Contract / Uncontract #######################

  /*!
   * Contracts the vertex pair (u,v). The representative u remains in the hypergraph.
   * The contraction partner v is removed from the hypergraph.
   *
   * For each hyperedge e incident to v, a contraction lead to one of two operations:
   * 1.) If e contained both u and v, then v is removed from e.
   * 2.) If e only contained v, than the slot of v in the incidence structure of e
   *     is reused to store u.
   *
   * The returned Memento can be used to undo the contraction via an uncontract operation.
   *
   * NOTE, this function is not thread-safe and should be only called in a single-threaded
   * setting.
   *
   * \param u Representative hypernode that will remain in the hypergraph
   * \param v Contraction partner that will be removed from the hypergraph
   */
  Memento contract(const HypernodeID u,
                   const HypernodeID v) {
    ASSERT(nodeIsEnabled(u), "Hypernode" << u << "is disabled");
    ASSERT(nodeIsEnabled(v), "Hypernode" << v << "is disabled");
    // TODO(heuer): Assertions verifies that both node have same part id

    DBG << "Contracting (" << u << "," << v << ")";
    setNodeWeight(u, nodeWeight(u) + nodeWeight(v));

    StreamingHypergraph& hypergraph_of_u = hypergraph_of_vertex(u);
    for ( const HyperedgeID& he : incidentEdges(v) ) {
      hypergraph_of_edge(he).contract(u, v, he, hypergraph_of_u);
    }

    disableHypernode(v);
    return Memento { u, v };
  }

  /*!
   * Contracts the vertex pair (u,v) belonging to the same community.
   * The representative u remains in the hypergraph. The contraction partner
   * v is removed from the hypergraph.
   *
   * For each hyperedge e incident to v, a contraction lead to one of two operations:
   * 1.) If e contained both u and v, then v is removed from e.
   * 2.) If e only contained v, than the slot of v in the incidence structure of e
   *     is reused to store u.
   *
   * The returned Memento can be used to undo the contraction via an uncontract operation.
   *
   * NOTE, in order that this function works correct, community hyperedges have to be
   * initialized beforehand. This function is thread-safe as long as only one thread
   * performs contractions in one community.
   *
   * \param u Representative hypernode that will remain in the hypergraph
   * \param v Contraction partner that will be removed from the hypergraph
   * \param community_id Community to which u and v belongs to
   */
  Memento contract(const HypernodeID u,
                   const HypernodeID v,
                   const PartitionID community_id) {
    ASSERT(nodeIsEnabled(u), "Hypernode" << u << "is disabled");
    ASSERT(nodeIsEnabled(v), "Hypernode" << v << "is disabled");
    ASSERT(communityID(u) == community_id);
    ASSERT(communityID(v) == community_id);
    // TODO(heuer): Assertions verifies that both node have same part id

    DBG << "Contracting (" << u << "," << v << ")";
    setNodeWeight(u, nodeWeight(u) + nodeWeight(v));

    StreamingHypergraph& hypergraph_of_u = hypergraph_of_vertex(u);
    for ( const HyperedgeID& he : incidentEdges(v) ) {
      hypergraph_of_edge(he).contract(u, v, he, community_id, hypergraph_of_u);
    }

    disableHypernode(v);
    return Memento { u, v, community_id };
  }


  /*!
  * Undoes a contraction operation that was remembered by the memento.
  *
  * During parallel community coarsening, we detect parallel hyperedges. If a contraction
  * is performed, we check if some of the incident hyperedges become parallel. Note, that
  * we perform contractions in parallel (one thread per community). To not interfer with
  * other threads (community coarsener), we only invalidate (not remove) a parallel hyperedge
  * for all pins that belongs to the community which detects it. Thus, it can happen that
  * in some communities a parallel hyperedge is detected and in some not. Parallel hyperedges
  * that are detected during contraction of two vertices are restored before the call to
  * the uncontraction function. However, there are still cases that due to an other
  * uncontraction, that did not detect that a hyperedge is parallel, a parallel hyperedge
  * becomes non-parallel to its representative. In such cases, the restore operation have to
  * be handled by the uncontraction function.
  *
  * After coarsening we resort the incident nets array such that the first part contains
  * all invalidated and second all enabled hyperedges. To detect whether a invalidated
  * hyperedge e' becomes non-parallel to its enabled representative e'', we check the following
  * cases (u and v are the vertices which are uncontracted):
  *   1.) v was part of e', but not of e'' before contraction (and vice versa)
  *   2.) e' was part of I(v), but e'' not before contraction (and vice versa)
  *
  * NOTE, this function is not thread-safe and should only be called in a single-threaded
  * setting.
  *
  * \param memento Memento remembering the contraction operation that should be reverted
  * \param parallel_he_representative for each disabled hyperedge it contains it
  *                                   parallel representative.
  */
  void uncontract(const Memento& memento,
                  parallel::scalable_vector<HyperedgeID>& parallel_he_representative) {
    ASSERT(nodeIsEnabled(memento.u), "Hypernode" << memento.u << "is disabled");
    ASSERT(!nodeIsEnabled(memento.v), "Hypernode" << memento.v << "is not invalid");

    DBG << "uncontracting (" << memento.u << "," << memento.v << ")";
    reverseContraction(memento);
    markAllIncidentNetsOf(memento.v);

    auto find_representative = [&](const HyperedgeID e) {
      HyperedgeID representative = originalEdgeID(e);
      ASSERT(representative < parallel_he_representative.size());
      while( parallel_he_representative[representative] != kInvalidHyperedge ) {
        representative = parallel_he_representative[representative];
      }
      return globalEdgeID(representative);
    };

    // Uncontraction starts by checking if a disabled parallel hyperedge becomes non-parallel
    // to one of its representatives. Usually all parallel hyperedges are enabled before
    // uncontraction (see hypergraph pruner). However, in some cases this is not the case.
    // Therefore, we perform an explicit check here if two hyperedges become non-parallel
    // after uncontraction.
    const auto& incident_hes_of_u = hypergraph_of_vertex(memento.u).incident_nets(memento.u);
    size_t incident_hes_start = hypergraph_of_vertex(memento.u).hypernode(memento.u).invalidIncidentNets();
    std::vector<HyperedgeID> disabled_hyperedges;
    for ( size_t incident_hes_it = 0; incident_hes_it != incident_hes_start; ++incident_hes_it ) {
      disabled_hyperedges.push_back(incident_hes_of_u[incident_hes_it]);
    }
    // All disabled hyperedges have to be traversed in decreasing order of their edge id
    // when checking if they become non-parallel to one of its representatives.
    std::sort(disabled_hyperedges.begin(), disabled_hyperedges.begin() + incident_hes_start,
              [&](const HyperedgeID& lhs, const HyperedgeID& rhs) {
                return lhs > rhs;
              });

    // Check if a disabled parallel hyperedges will become non-parallel to
    // its representative.
    for ( const HyperedgeID& he : disabled_hyperedges ) {
      if (!edgeIsEnabled(he)) {
        const StreamingHypergraph& hypergraph_of_he = hypergraph_of_edge(he);
        HyperedgeID representative = find_representative(he);
        const size_t edge_size = edgeSize(find_representative(he));
        bool becomes_non_parallel = false;

        HyperedgeID last_representative = he;
        HyperedgeID current_representative = originalEdgeID(he);
        representative = globalEdgeID(current_representative);
        // Verify if hyperedge becomes non-parallel to one of its representatives (all hyperedges
        // on the path to the root in the hyperedge representative tree)
        while ( !becomes_non_parallel && parallel_he_representative[current_representative] != kInvalidHyperedge ) {
          last_representative = representative;
          current_representative = parallel_he_representative[current_representative];
          representative = globalEdgeID(current_representative);

          const StreamingHypergraph& hypergraph_of_rep = hypergraph_of_edge(representative);

          // In case, both hyperedges fall into different uncontraction cases, than both become
          // non parallel afterwards.
          bool is_case_1_rep = hypergraph_of_rep.get_uncontraction_case(representative, edge_size, memento.v) ==
                              StreamingHypergraph::UncontractionCase::CASE_1;
          bool is_case_1_he = hypergraph_of_he.get_uncontraction_case(he, edge_size, memento.v) ==
                              StreamingHypergraph::UncontractionCase::CASE_1;

          // In case, the contraction partner v contains either the disabled hyperedge or
          // the representative (but not both), than both become non parallel afterwards.
          becomes_non_parallel = becomes_non_parallel ||
                                 ( hypergraph_of_he.containsIncidentNet(he) &&
                                   !hypergraph_of_rep.containsIncidentNet(representative) ) ||
                                 ( !hypergraph_of_he.containsIncidentNet(he) &&
                                   hypergraph_of_rep.containsIncidentNet(representative) ) ||
                                 ( !is_case_1_rep && is_case_1_he ) ||
                                 ( is_case_1_rep && !is_case_1_he );
        }

        if ( becomes_non_parallel ) {
          restoreParallelHyperedge(last_representative, parallel_he_representative);
          incident_hes_start = hypergraph_of_vertex(memento.u).hypernode(memento.u).invalidIncidentNets();
        }
      }
    }

    // Uncontract all disabled parallel hyperedges
    #ifndef NDEBUG
    // In case of debug mode, we uncontract each disabled parallel hyperedge also in order
    // to check that each hyperedge is parallel to its representative
    for (size_t incident_hes_it = 0; incident_hes_it != incident_hes_start; ++incident_hes_it) {
      const HyperedgeID he = incident_hes_of_u[incident_hes_it];
      const HyperedgeID representative = find_representative(he);

      if ( hypergraph_of_edge(he).uncontract(memento.u, memento.v, he, representative,
            incident_hes_it, _hypergraphs) ) {
        --incident_hes_it;
        --incident_hes_start;
      }
    }
    #endif

    // Uncontract all active hyperedges
    size_t incident_hes_end = incident_hes_of_u.size();
    for (size_t incident_hes_it = incident_hes_start; incident_hes_it != incident_hes_end; ++incident_hes_it) {
      const HyperedgeID he = incident_hes_of_u[incident_hes_it];
      if ( hypergraph_of_edge(he).uncontract(memento.u, memento.v, he, incident_hes_it, _hypergraphs) ) {
        --incident_hes_it;
        --incident_hes_end;
      }
    }

    // Remove all previously enabled parallel hyperedges from invalid part of incident nets of v
    StreamingHypergraph& hypergraph_of_v = hypergraph_of_vertex(memento.v);
    auto& incident_hes_of_v = hypergraph_of_v.incident_nets(memento.v);
    incident_hes_start = hypergraph_of_v.hypernode(memento.v).invalidIncidentNets();
    for (size_t incident_hes_it = 0; incident_hes_it != incident_hes_start; ++incident_hes_it) {
      const HyperedgeID he = incident_hes_of_v[incident_hes_it];
      if ( edgeIsEnabled(he) ) {
        std::swap(incident_hes_of_v[incident_hes_it--], incident_hes_of_v[--incident_hes_start]);
        hypergraph_of_v.hypernode(memento.v).decrementInvalidIncidentNets();
      }
    }

    setNodeWeight(memento.u, nodeWeight(memento.u) - nodeWeight(memento.v));

    HEAVY_REFINEMENT_ASSERT( numIncidentCutHyperedges(memento.u) == numIncidentCutHEs(memento.u),
      V(memento.u) << V(numIncidentCutHyperedges(memento.u)) << V(numIncidentCutHEs(memento.u)) );
    HEAVY_REFINEMENT_ASSERT( numIncidentCutHyperedges(memento.v) == numIncidentCutHEs(memento.v),
      V(memento.v) << V(numIncidentCutHyperedges(memento.v)) << V(numIncidentCutHEs(memento.v)) );
  }

  // ####################### Remove / Restore Hyperedges #######################

  /*!
  * Removes a hyperedge from the hypergraph. This includes the removal of he from all
  * of its pins and to disable the hyperedge.
  *
  * NOTE, this function is not thread-safe and should only be called in a single-threaded
  * setting.
  */
  void removeEdge(const HyperedgeID he) {
    ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");
    for ( const HypernodeID& pin : pins(he) ) {
      hypergraph_of_vertex(pin).removeIncidentEdgeFromHypernode(he, pin);
    }
    hypergraph_of_edge(he).disableHyperedge(he);
  }

  /*!
  * Removes a single-pin community hyperedge from the hypergraph.
  *
  * NOTE, in order that this function works correct, community hyperedges have to be
  * initialized beforehand. This function is thread-safe as long as only one thread
  * performs contractions in one community.
  */
  void removeSinglePinCommunityEdge(const HyperedgeID he, const PartitionID community_id) {
    ASSERT(hypergraph_of_edge(he).numCommunitiesInHyperedge(he) == 1,
           "Only allowed to remove hyperedges that contains pins from one community");
    ASSERT(hypergraph_of_edge(he).edgeSize(he) == 1, "Hyperedge is not a single-pin hyperedge");
    removeEdge(he, community_id);
    hypergraph_of_edge(he).disableHyperedge(he);
  }

  /*!
  * Removes a parallel hyperedge from the hypergraph.
  * Note, that in case the pins of the hyperedge belongs only to one community,
  * the hyperedge is directly removed from the hypergraph. Otherwise, it is only
  * invalidated (for an explanation of the term "invalidated", see illustration over
  * function incidentEdges and documentation of function uncontraction(...)).
  *
  * NOTE, in order that this function works correct, community hyperedges have to be
  * initialized beforehand. This function is thread-safe as long as only one thread
  * performs contractions in one community.
  */
  void removeParallelEdge(const HyperedgeID he, const PartitionID community_id) {
    if ( numCommunitiesInHyperedge(he) == 1 ) {
      removeEdge(he, community_id);
    } else {
      invalidateEdge(he, community_id);
    }
  }

  // ! Restores an hyperedge of a certain size.
  void restoreEdge(const HyperedgeID he, const size_t size,
                   const HyperedgeID representative = kInvalidHyperedge) {
    ASSERT(!edgeIsEnabled(he), "Hyperedge" << he << "already enabled");
    ASSERT(representative == kInvalidHyperedge || edgeIsEnabled(representative),
           "Hyperedge" << representative << "is disabled");
    enableHyperedge(he);
    hypergraph_of_edge(he).hyperedge(he).setSize(size);
    StreamingHypergraph& hypergraph_of_he = hypergraph_of_edge(he);
    bool representative_is_cut = representative != kInvalidHyperedge && connectivity(representative) > 1;
    for ( const HypernodeID& pin : pins(he) ) {
      hypergraph_of_vertex(pin).insertIncidentEdgeToHypernode(he, pin);
      ASSERT(partID(pin) != kInvalidPartition, V(pin) << V(partID(pin)));
      hypergraph_of_he.incrementPinCountInPart(he, partID(pin));
      if ( representative_is_cut ) {
        hypergraph_of_vertex(pin).incrementIncidentNumCutHyperedges(pin);
      }
    }
  }

  // ! Restores a single-pin hyperedge
  void restoreSinglePinHyperedge(const HyperedgeID he) {
    ASSERT(!edgeIsEnabled(he), "Hyperedge" << he << "already enabled");
    restoreEdge(he, 1);
  }

  /*!
   * Restores a parallel hyperedge.
   *
   * The representative of a parallel hyperedge he is stored in _parallel_he_representative.
   * However, the representative can also be parallel to an other hyperedge. The array
   * _parallel_he_representative forms a tree of parallel hyperedges. In order to restore
   * a parallel hyperedge, we have to restore all parallel hyperedges from he to the root of
   * the tree recursively.
   *
   * @param he parallel hyperedge to restore
   * @param _parallel_he_representative parallel hyperedge representative tree
   */
  void restoreParallelHyperedge(const HyperedgeID he,
                                parallel::scalable_vector<HyperedgeID>& _parallel_he_representative) {
    if ( !edgeIsEnabled(he) ) {
      DBG << "restore parallel HE" << he << "in hypergraph";
      const HyperedgeID representative = globalEdgeID(_parallel_he_representative[originalEdgeID(he)]);
      restoreParallelHyperedge(representative, _parallel_he_representative);

      #ifdef NDEBUG
      // If we are not in debug mode, we copy the content/pins of the representative to
      // to its parallel hyperedge we want to restore.
      StreamingHypergraph& hypergraph_of_he = hypergraph_of_edge(he);
      StreamingHypergraph& hypergraph_of_rep = hypergraph_of_edge(representative);
      memcpy(hypergraph_of_he._incidence_array.data() + hypergraph_of_he.hyperedge(he).firstEntry(),
              hypergraph_of_rep._incidence_array.data() + hypergraph_of_rep.hyperedge(representative).firstEntry(),
              edgeSize(representative) * sizeof(HypernodeID));
      #endif

      HEAVY_REFINEMENT_ASSERT(verifyThatHyperedgesAreParallel(representative, he),
                              "HE" << he << "is not parallel to" << representative);
      restoreEdge(he, edgeSize(representative), representative);
      setEdgeWeight(representative, edgeWeight(representative) - edgeWeight(he));
      _parallel_he_representative[originalEdgeID(he)] = kInvalidHyperedge;
    }
  }

  // ####################### Initialization / Reset Functions #######################

  /*!
   * Adds community hyperedges to the hypergraphs. Community hyperedges
   * point to a consecutive range of pins of the original hyperedge that
   * belongs to the same community. It allows to perform parallel contractions
   * on the hypergraph within each community.
   * Additionally the incident nets array of each hypernode is reordered such
   * that the first part contains all single-pin community hyperedges and second
   * all others.
   *
   * Note, this function have to be called before parallel community coarsening.
   */
  void initializeCommunityHyperedges() {
    TBBNumaArena::instance().execute_on_all_numa_nodes([&](const int node) {
      _hypergraphs[node].initializeCommunityHyperedges(_hypergraphs);
    });

    TBBNumaArena::instance().execute_on_all_numa_nodes([&](const int node) {
      _hypergraphs[node].initializeCommunityHypernodes(_hypergraphs);
    });
  }

  /*!
   * Removes all community hyperedges from the hypergraph after parallel community
   * coarsening terminates.
   *
   * The pins of the original hyperedge are sorted in decreasing order of their
   * contraction index. The contraction index of a vertex v is defined as the index
   * of the contraction (u,v) in the contraction history, where v occurs as contraction
   * partner. This is done to fullfil the invariants required by the uncontraction method.
   *
   * Note this function have to be called after parallel community coarsening such
   * that uncontractions can be performed correctly.
   *
   * @params history contraction history
   */
  void removeCommunityHyperedges(const std::vector<Memento>& history) {
    TBBNumaArena::instance().execute_on_all_numa_nodes([&](const int node) {
      _hypergraphs[node].removeCommunityHyperedges(history, _num_hypernodes, _hypergraphs);
    });
  }

  /*!
   * Initializes community-related information after all vertices are assigned to a community.
   * This includes:
   *  1.) Number of Communities
   *  2.) Number of Vertices per Community
   *  3.) Number of Pins per Community
   *  4.) For each hypernode v of community C, we compute a unique id within
   *      that community in the range [0, |C|)
   */
  void initializeCommunities() {
    // Compute number of communities
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    _num_communities = tbb::parallel_reduce(tbb::blocked_range<HypernodeID>(0UL, _num_hypernodes), 0,
      [this](const tbb::blocked_range<HypernodeID>& range, PartitionID init) {
        PartitionID num_communities = init;
        for ( HypernodeID hn = range.begin(); hn < range.end(); ++hn ) {
          num_communities = std::max(num_communities, communityID(globalNodeID(hn)) + 1);
        }
        return num_communities;
      },
      [](const PartitionID lhs, const PartitionID rhs) {
        return std::max(lhs, rhs);
      });
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("compute_number_of_communities", "Compute Num of Communities",
      "initialize_communities", mt_kahypar::utils::Timer::Type::PREPROCESSING, 0, std::chrono::duration<double>(end - start).count());

    // Compute number of hypernodes per community and also for each node
    // a unique node id within each community
    start = std::chrono::high_resolution_clock::now();
    _communities_num_hypernodes.assign(_num_communities, 0);
    for ( const HypernodeID& hn : nodes() ) {
      PartitionID community_id = communityID(hn);
      ASSERT(community_id < _num_communities);
      hypergraph_of_vertex(hn).hypernode(hn).setCommunityNodeId(_communities_num_hypernodes[community_id]);
      ++_communities_num_hypernodes[community_id];
    }
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("compute_num_community_hns", "Compute Num Community HNs",
      "initialize_communities", mt_kahypar::utils::Timer::Type::PREPROCESSING, 1, std::chrono::duration<double>(end - start).count());

    // Compute number of pins per community
    start = std::chrono::high_resolution_clock::now();
    _communities_num_pins.assign(_num_communities, 0);
    for ( const HyperedgeID& he : edges() ) {
      for ( const HypernodeID& pin : pins(he) ) {
        ASSERT(communityID(pin) < _num_communities);
        ++_communities_num_pins[communityID(pin)];
      }
    }
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("compute_num_community_pins", "Compute Num Community Pins",
      "initialize_communities", mt_kahypar::utils::Timer::Type::PREPROCESSING, 2, std::chrono::duration<double>(end - start).count());
  }

  // ! Resets the ids of all pins in the incidence array to its original node id
  void resetPinsToOriginalNodeIds() {
    TBBNumaArena::instance().execute_on_all_numa_nodes([&](const int node) {
      _hypergraphs[node].resetPinsToOriginalNodeIds(_hypergraphs);
    });
  }

  // ! Invalidates all disabled hyperedges from the incident nets array of each node
  // ! For further details please take a look at the documentation of uncontraction(...)
  void invalidateDisabledHyperedgesFromIncidentNets() {
    TBBNumaArena::instance().execute_on_all_numa_nodes([&](const int node) {
      _hypergraphs[node].invalidateDisabledHyperedgesFromIncidentNets(_hypergraphs);
    });
  }

 private:

  // ####################### Uncontraction Functions #######################

  // ! Restores the contraction partner of the contraction
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void reverseContraction(const Memento& memento) {
    enableHypernode(memento.v);
    PartitionID part_id = partID(memento.u);
    ASSERT(part_id != kInvalidPartition);
    hypergraph_of_vertex(memento.v).setNodePart(memento.v, part_id);
    _local_part_info.local().apply(part_id, PartInfo { 0, 1 });
  }

  // ! Mark all incident nets of hypernode v in bit set.
  // ! Note, during a uncontraction of two vertices u and v only nets
  // ! that are part of I(u) and I(v) are uncontracted.
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void markAllIncidentNetsOf(const HypernodeID v) {
    hypergraph_of_vertex(v).markAllIncidentNetsOf(v, _hypergraphs);
  }

  // ####################### Remove Hyperedges #######################

  // ! Removes a community hyperedge from incident net arrays of all its pins
  void removeEdge(const HyperedgeID he, const PartitionID community_id) {
    ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");
    for ( const HypernodeID& pin : pins(he, community_id) ) {
      hypergraph_of_vertex(pin).removeIncidentEdgeFromHypernode(
        he, pin, community_id, _hypergraphs);
    }
  }

  // ! Invalidates a community hyperedge in all incident net arrays of all its pins
  void invalidateEdge(const HyperedgeID he, const PartitionID community_id) {
    ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");
    hypergraph_of_edge(he).disableHyperedge(he, community_id);
    for ( const HypernodeID& pin : pins(he, community_id) ) {
      hypergraph_of_vertex(pin).removeIncidentEdgeFromHypernode(
        he, pin, community_id, _hypergraphs, true /* invalidate only */);
    }
  }

  // ####################### Helper Functions #######################

  // ! Computes a mapping from vertex to numa hypernode.
  // ! A vertex is assigned to the numa node where it occurs most as pin.
  void computeNodeMapping() {
    size_t num_streaming_hypergraphs = _hypergraphs.size();
    // Computes mapping for each node to a streaming hypergraph
    // A node is assigned to the streaming hypergraph where it occurs
    // most as pin.
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(0UL, _num_hypernodes),
      [&](const tbb::blocked_range<HypernodeID>& range) {
      for ( HypernodeID hn = range.begin(); hn < range.end(); ++hn ) {
        size_t max_pins = 0;
        HypernodeID max_node_id = 0;
        for ( HypernodeID node = 1; node < num_streaming_hypergraphs; ++node ) {
          size_t num_pins = _hypergraphs[node].vertexPinCount(hn);
          if ( num_pins > max_pins ) {
            max_pins = num_pins;
            max_node_id = node;
          }
        }
        ASSERT(max_node_id < _hypergraphs.size());
        _node_mapping[hn] = max_node_id;
      }
    });
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("compute_node_mapping", "Compute Node Mapping",
      "initialize_hypernodes", mt_kahypar::utils::Timer::Type::IMPORT, 0, std::chrono::duration<double>(end - start).count());
  }

  /*!
   * Initializes the hypernodes of the hypergraph.
   * Hypernodes are streamed into its corresponding numa hypergraphs (defined
   * by the vertex to numa node mapping) and afterwards the incident nets data
   * structure is initialized.
   */
  void initializeHypernodes() {
    // Verify that node mapping is valid
    ASSERT([&]() {
      for ( HypernodeID hn = 0; hn < _num_hypernodes; ++hn ) {
        if ( _node_mapping[hn] >= _hypergraphs.size() ) {
          LOG << "Hypernode" << hn << "should be mapped to hypergraph on node"
              << _node_mapping[hn] << ", but there are only" << _hypergraphs.size()
              << "nodes";
          return false;
        }
      }
      return true;
    }(), "Invalid node mapping");

    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    size_t num_streaming_hypergraphs = _hypergraphs.size();
    // Stream hypernodes into corresponding streaming hypergraph, where it
    // is assigned to
    std::vector<HypernodeID> tmp_node_mapping(_num_hypernodes);
    for ( HypernodeID node = 0; node < num_streaming_hypergraphs; ++node ) {
      TBBNumaArena::instance().numa_task_arena(node).execute([&] {
        TBBNumaArena::instance().numa_task_group(node).run([&, node] {
          tbb::parallel_for(tbb::blocked_range<HypernodeID>(0UL, _num_hypernodes),
            [&](const tbb::blocked_range<HypernodeID>& range) {
            for ( HypernodeID hn = range.begin(); hn < range.end(); ++hn ) {
              if ( _node_mapping[hn] == node ) {
                tmp_node_mapping[hn] = _hypergraphs[node].streamHypernode(hn, 1);
              }
            }
          });
        });
      });
    }
    TBBNumaArena::instance().wait();
    _node_mapping = std::move(tmp_node_mapping);
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("stream_hypernodes", "Stream Hypernodes",
      "initialize_hypernodes", mt_kahypar::utils::Timer::Type::IMPORT, 1, std::chrono::duration<double>(end - start).count());

    // Initialize hypernodes on each streaming hypergraph
    // NOTE, that also involves streaming local incident nets to other
    // streaming hypergraphs
    start = std::chrono::high_resolution_clock::now();
    for ( size_t node = 0; node < num_streaming_hypergraphs; ++node ) {
      TBBNumaArena::instance().numa_task_arena(node).execute([&] {
        TBBNumaArena::instance().numa_task_group(node).run([&, node] {
          _hypergraphs[node].initializeHypernodes(_hypergraphs, _node_mapping);
        });
      });
    }
    TBBNumaArena::instance().wait();
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("initialize_numa_hypernodes", "Initialize Numa Hypernodes",
      "initialize_hypernodes", mt_kahypar::utils::Timer::Type::IMPORT, 2, std::chrono::duration<double>(end - start).count());

    // Verify that number of hypernodes is equal to number of hypernodes
    // in streaming hypergraphs
    ASSERT([&]{
      HypernodeID actual_number_of_nodes = 0;
      for ( size_t node = 0; node < num_streaming_hypergraphs; ++node ) {
        actual_number_of_nodes += _hypergraphs[node].initialNumNodes();
      }
      if ( actual_number_of_nodes == _num_hypernodes ) {
        return true;
      } else {
        LOG << V(actual_number_of_nodes) << V(_num_hypernodes);
        for ( size_t node = 0; node < num_streaming_hypergraphs; ++node ) {
          LOG << V(node) << V(_hypergraphs[node].initialNumNodes());
        }
        return false;
      }
    }(), "Invalid number hypernodes in streaming hypergraph");

    // Initialize incident nets of hypernodes
    start = std::chrono::high_resolution_clock::now();
    for ( size_t node = 0; node < num_streaming_hypergraphs; ++node ) {
    TBBNumaArena::instance().numa_task_arena(node).execute([&] {
      TBBNumaArena::instance().numa_task_group(node).run([&, node] {
          _hypergraphs[node].initializeIncidentNets();
        });
      });
    }
    TBBNumaArena::instance().wait();
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("initialize_incident_nets", "Initialize Incident Nets",
      "initialize_hypernodes", mt_kahypar::utils::Timer::Type::IMPORT, 3, std::chrono::duration<double>(end - start).count());


    ASSERT([&] {
      // Internally verify that incident nets are constructed correctly
      for ( size_t node = 0; node < num_streaming_hypergraphs; ++node ) {
        if ( !_hypergraphs[node].verify_incident_nets_of_hypergraph(_hypergraphs) ) {
          return false;
        }
      }
      return true;
    }(), "Initialization of incident nets failed");

    for ( size_t node = 0; node < num_streaming_hypergraphs; ++node ) {
      _num_hyperedges += _hypergraphs[node].initialNumEdges();
      _num_pins += _hypergraphs[node].initialNumPins();
    }

    start = std::chrono::high_resolution_clock::now();
    _edge_mapping.resize(_num_hyperedges);
    for ( const HyperedgeID& he : edges() ) {
      HyperedgeID original_id = hypergraph_of_edge(he).originalEdgeId(he);
      ASSERT(original_id < _edge_mapping.size());
      _edge_mapping[original_id] = he;
    }
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("initialize_he_mapping", "Initialize HE Mapping",
      "initialize_hypernodes", mt_kahypar::utils::Timer::Type::IMPORT, 4, std::chrono::duration<double>(end - start).count());
  }

  // ####################### Helper Functions #######################

  const StreamingHypergraph& hypergraph_of_vertex(const HypernodeID u) const {
    int node = StreamingHypergraph::get_numa_node_of_vertex(u);
    ASSERT(node < (int) _hypergraphs.size());
    return _hypergraphs[node];
  }

  const StreamingHypergraph& hypergraph_of_edge(const HyperedgeID e) const {
    int node = StreamingHypergraph::get_numa_node_of_hyperedge(e);
    ASSERT(node < (int) _hypergraphs.size());
    return _hypergraphs[node];
  }

  StreamingHypergraph& hypergraph_of_vertex(const HypernodeID u) {
    return const_cast<StreamingHypergraph&>(static_cast<const Hypergraph&>(*this).hypergraph_of_vertex(u));
  }

  StreamingHypergraph& hypergraph_of_edge(const HyperedgeID e) {
    return const_cast<StreamingHypergraph&>(static_cast<const Hypergraph&>(*this).hypergraph_of_edge(e));
  }

  // ! Returns the number of incident hyperedges of a hypernode that connect more than on block.
  // ! This method is used internally to verify the consistency of the public isBorderNode method.
  HyperedgeID numIncidentCutHEs(const HypernodeID hn) const {
    HyperedgeID num_cut_hes = 0;
    for (const HyperedgeID& he : incidentEdges(hn)) {
      if (connectivity(he) > 1) {
        ++num_cut_hes;
      }
    }
    return num_cut_hes;
  }

  // ! Only for debugging
  bool verifyThatHyperedgesAreParallel(const HyperedgeID representative, const HyperedgeID parallel_he) {
    if ( !edgeIsEnabled(representative) || edgeIsEnabled(parallel_he) ) {
      LOG << "HE" << representative << "must be enabled and HE" << parallel_he << "disabled";
      return false;
    }

    std::set<HypernodeID> contained_pins;
    for ( const HypernodeID& pin : pins(representative) ) {
      contained_pins.insert(pin);
    }

    size_t edge_size = contained_pins.size();
    StreamingHypergraph& hypergraph_of_he = hypergraph_of_edge(parallel_he);
    size_t incidence_array_start = hypergraph_of_he.hyperedge(parallel_he).firstEntry();
    for ( size_t pos = incidence_array_start; pos < incidence_array_start + edge_size; ++pos ) {
      const HypernodeID pin = hypergraph_of_he._incidence_array[pos];
      if ( contained_pins.find(pin) == contained_pins.end() ) {
        LOG << "Pin" << pin << "of HE" << parallel_he << "is not contained in HE" << representative;
        hypergraph_of_he.printHyperedgeInfo(parallel_he);
        hypergraph_of_edge(representative).printHyperedgeInfo(representative);
        return false;
      }
    }

    return true;
  }

  // ! Number of hypernodes
  HypernodeID _num_hypernodes;
  // ! Number of hyperedges
  HyperedgeID _num_hyperedges;
  // ! Number of pins
  HypernodeID _num_pins;
  // ! Number of communities
  PartitionID _num_communities;
  // ! Number of blocks
  PartitionID _k;

  // ! Number of hypernodes in a community
  parallel::scalable_vector<HypernodeID> _communities_num_hypernodes;
  // ! Number of pins in a community
  parallel::scalable_vector<HypernodeID> _communities_num_pins;
  // ! Global weight and size information for all blocks.
  std::vector<PartInfo> _part_info;
  // ! Thread local weight and size information for all blocks.
  ThreadLocalPartInfos _local_part_info;

  // ! Numa Hypergraphs
  // ! Hypergraph on index i is assigned to numa node i
  std::vector<StreamingHypergraph> _hypergraphs;
  // ! Mapping from original node id to its hypergraph node id
  std::vector<HypernodeID> _node_mapping;
  // ! Mapping from original edge id to its hypergraph edge id
  std::vector<HyperedgeID> _edge_mapping;
  // ! Mapping from community id to numa node
  std::vector<PartitionID> _community_node_mapping;
};

} // namespace ds
} // namespace mt_kahypar
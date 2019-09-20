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

#include "kahypar/macros.h"
#include "kahypar/meta/mandatory.h"

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

  static constexpr bool debug = false;

  using HypernodeID = HypernodeType_;
  using HyperedgeID = HyperedgeType_;
  using HypernodeWeight = HypernodeWeightType_;
  using HyperedgeWeight = HyperedgeWeightType_;
  using PartitionID = PartitionIDType_;

  using StreamingHypergraph = mt_kahypar::ds::StreamingHypergraph<HypernodeID, 
                                                                  HyperedgeID, 
                                                                  HypernodeWeight, 
                                                                  HyperedgeWeight, 
                                                                  PartitionID,
                                                                  HardwareTopology,
                                                                  TBBNumaArena>;

  using HypernodeIterator = typename StreamingHypergraph::HypernodeIterator;
  using HyperedgeIterator = typename StreamingHypergraph::HyperedgeIterator;
  using IncidenceIterator = typename StreamingHypergraph::IncidenceIterator;
  
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
    GlobalHypergraphElementIterator() = default;

    GlobalHypergraphElementIterator(const GlobalHypergraphElementIterator& other) = default;
    GlobalHypergraphElementIterator& operator= (const GlobalHypergraphElementIterator& other) = default;

    GlobalHypergraphElementIterator(GlobalHypergraphElementIterator&& other) = default;
    GlobalHypergraphElementIterator& operator= (GlobalHypergraphElementIterator&& other) = default;

    ~GlobalHypergraphElementIterator() = default;

    /*!
     * Construct a GlobalHypergraphElementIterator
     * See GenericHypergraph::nodes() or GenericHypergraph::edges() for usage.
     *
     * If start_element is invalid, the iterator advances to the first valid
     * element.
     *
     * \param start_element A pointer to the starting position
     * \param id The index of the element the pointer points to
     * \param max_id The maximum index allowed
     */
    GlobalHypergraphElementIterator(Iterators&& iterators) :
      _iterators(std::move(iterators)),
      _idx(0) { 
      ASSERT(_iterators.size() > 0);
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

  // ! Iterator to iterator over the hypernodes
  using GlobalHypernodeIterator = GlobalHypergraphElementIterator<HypernodeIterator>;
  // ! Iterator to iterator over the hyperedges
  using GlobalHyperedgeIterator = GlobalHypergraphElementIterator<HyperedgeIterator>;

 public:
  /*!
  * A memento stores all information necessary to undo the contraction operation
  * of a vertex pair \f$(u,v)\f$.
  *
  * A contraction operations can increase the set \f$I(u)\f$ of nets incident
  * to \f$u\f$. This in turn leads to changes in _incidence_array. Therefore
  * the memento stores the initial starting index of \f$u\f$'s incident nets
  * as well as the old size, i.e. \f$|I(u)|\f$.
  *
  */
  struct Memento {
    // ! The representative hypernode that remains in the hypergraph
    HypernodeID u;
    // ! The contraction partner of u that is removed from the hypergraph after the contraction.
    HypernodeID v;
  };

  explicit Hypergraph() :
    _num_hypernodes(0),
    _num_hyperedges(0),
    _num_pins(0),
    _num_communities(0),
    _communities_num_hypernodes(),
    _communities_num_pins(),
    _hypergraphs(),
    _node_mapping() { }

  Hypergraph(const HypernodeID num_hypernodes,
             std::vector<StreamingHypergraph>&& hypergraphs) :
    _num_hypernodes(num_hypernodes),
    _num_hyperedges(0),
    _num_pins(0),
    _num_communities(0),
    _communities_num_hypernodes(),
    _communities_num_pins(),
    _hypergraphs(std::move(hypergraphs)),
    _node_mapping(num_hypernodes, 0) { 
    computeNodeMapping();
    initializeHypernodes();
  }

  Hypergraph(const HypernodeID num_hypernodes,
             std::vector<StreamingHypergraph>&& hypergraphs,
             const std::vector<HypernodeID>& node_mapping) :
    _num_hypernodes(num_hypernodes),    
    _num_hyperedges(0),
    _num_pins(0),
    _num_communities(0),
    _communities_num_hypernodes(),
    _communities_num_pins(),
    _hypergraphs(std::move(hypergraphs)),
    _node_mapping(node_mapping) { 
    initializeHypernodes();
  }

  Hypergraph(const HypernodeID num_hypernodes,
             std::vector<StreamingHypergraph>&& hypergraphs,
             std::vector<HypernodeID>&& node_mapping) :
    _num_hypernodes(num_hypernodes),
    _num_hyperedges(0),
    _num_pins(0),
    _num_communities(0),
    _communities_num_hypernodes(),
    _communities_num_pins(),
    _hypergraphs(std::move(hypergraphs)),
    _node_mapping(std::move(node_mapping)) { 
    initializeHypernodes();
  }

  Hypergraph(const Hypergraph&) = delete;
  Hypergraph& operator= (const Hypergraph&) = delete;

  Hypergraph(Hypergraph&& other) = default;
  Hypergraph& operator= (Hypergraph&&) = default;

  ~Hypergraph() = default;


  HypernodeID initialNumNodes() const {
    return _num_hypernodes;
  }

  HyperedgeID initialNumEdges() const {
    return _num_hyperedges;
  }

  HypernodeID initialNumPins() const {
    return _num_pins;
  }

  HypernodeID initialNumNodes(const int node) const {
    ASSERT(node < (int)_hypergraphs.size());
    return _hypergraphs[node].initialNumNodes();
  }

  HyperedgeID initialNumEdges(const int node) const {
    ASSERT(node < (int)_hypergraphs.size());
    return _hypergraphs[node].initialNumEdges();
  }

  HypernodeID initialNumPins(const int node) const {
    ASSERT(node < (int)_hypergraphs.size());
    return _hypergraphs[node].initialNumPins();
  }

  PartitionID numCommunities() const {
    return _num_communities;
  }

  HypernodeID initialNumCommunityHypernodes(const PartitionID community) const {
    ASSERT(community < (PartitionID) _communities_num_hypernodes.size());
    return _communities_num_hypernodes[community];
  }

  HypernodeID initialNumCommunityPins(const PartitionID community) const {
    ASSERT(community < (PartitionID) _communities_num_pins.size());
    return _communities_num_pins[community];
  }

  // TODO(heuer): Replace with correct counter
  HypernodeID currentNumNodes() const {
    return _num_hypernodes;
  }

  // TODO(heuer): Replace with correct counter
  HyperedgeID currentNumEdges() const {
    return _num_hyperedges;
  }

  // TODO(heuer): Replace with correct counter
  HypernodeID currentNumPins() const {
    return _num_pins;
  }

  HypernodeWeight totalWeight() const {
    HypernodeWeight weight = 0;
    for ( const StreamingHypergraph& hypergraph : _hypergraphs ) {
      weight += hypergraph.totalWeight();
    }
    return weight;
  }

  // ! Returns a for-each iterator-pair to loop over the set of incident hyperedges of hypernode u.
  std::pair<IncidenceIterator, IncidenceIterator> incidentEdges(const HypernodeID u) const {
    return hypergraph_of_vertex(u).incidentEdges(u);
  }

  // ! Returns a for-each iterator-pair to loop over the set pins of hyperedge e.
  std::pair<IncidenceIterator, IncidenceIterator> pins(const HyperedgeID e) const {
    return hypergraph_of_edge(e).pins(e);
  }

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

  std::pair<HypernodeIterator, HypernodeIterator>  nodes(const int node) const {
    ASSERT(node < (int) _hypergraphs.size());
    return _hypergraphs[node].nodes();
  }

  std::pair<HyperedgeIterator, HyperedgeIterator> edges(const int node) const {
    ASSERT(node < (int) _hypergraphs.size());
    return _hypergraphs[node].edges();
  }

  HypernodeID originalNodeID(const HypernodeID u) const {
    return hypergraph_of_vertex(u).originalNodeId(u);
  }

  HypernodeID globalNodeID(const HypernodeID u) const {
    ASSERT(u < _node_mapping.size());
    return _node_mapping[u];
  }

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
   * \param u Representative hypernode that will remain in the hypergraph
   * \param v Contraction partner that will be removed from the hypergraph
   *
   */
  Memento contract(const HypernodeID u, const HypernodeID v) {
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
  * Undoes a contraction operation that was remembered by the memento.
  * This is the default uncontract method.
  *
  * \param memento Memento remembering the contraction operation that should be reverted
  */
  void uncontract(const Memento& memento) {
    ASSERT(nodeIsEnabled(memento.u), "Hypernode" << memento.u << "is disabled");
    ASSERT(!nodeIsEnabled(memento.v), "Hypernode" << memento.v << "is not invalid");

    DBG << "uncontracting (" << memento.u << "," << memento.v << ")";
    reverseContraction(memento);
    markAllIncidentNetsOf(memento.v);

    const auto& incident_hes_of_u = hypergraph_of_vertex(memento.u).incidentNets(memento.u);
    size_t incident_hes_end = incident_hes_of_u.size();

    for (size_t incident_hes_it = 0; incident_hes_it != incident_hes_end; ++incident_hes_it) {
      const HyperedgeID he = incident_hes_of_u[incident_hes_it];
      if ( hypergraph_of_edge(he).uncontract(memento.u, memento.v, he, incident_hes_it, _hypergraphs) ) {
        --incident_hes_it;
        --incident_hes_end;
      }
    }
    setNodeWeight(memento.u, nodeWeight(memento.u) - nodeWeight(memento.v));
  }

  void removeEdge(const HyperedgeID he) {
    ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");
    for ( const HypernodeID& pin : pins(he) ) {
      hypergraph_of_vertex(pin).removeIncidentEdgeFromHypernode(he, pin);
    }
    hypergraph_of_edge(he).disableHyperedge(he);
    // TODO(heuer): invalidate pin counts of he
  }

  void restoreEdge(const HyperedgeID he) {
    ASSERT(!edgeIsEnabled(he), "Hyperedge" << he << "already enabled");
    enableHyperedge(he);
    // TODO(heuer): reset partition pin counts of he
    for ( const HypernodeID& pin : pins(he) ) {
      hypergraph_of_vertex(pin).insertIncidentEdgeToHypernode(he, pin);
    }
  }

  // TODO(heuer): restore operation for parallel hyperedges
  // if part ids are integrated

  HypernodeWeight nodeWeight(const HypernodeID u) const {
    return hypergraph_of_vertex(u).nodeWeight(u);
  }

  void setNodeWeight(const HypernodeID u, const HypernodeWeight weight) {
    hypergraph_of_vertex(u).setNodeWeight(u, weight);
  }

  HypernodeWeight edgeWeight(const HyperedgeID e) const {
    return hypergraph_of_edge(e).edgeWeight(e);
  }

  void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
    hypergraph_of_edge(e).setEdgeWeight(e, weight);
  }

  HyperedgeID nodeDegree(const HypernodeID u) const {
    return hypergraph_of_vertex(u).nodeDegree(u);
  }

  HypernodeID edgeSize(const HyperedgeID e) const {
    return hypergraph_of_edge(e).edgeSize(e);
  }

  size_t edgeHash(const HyperedgeID e) const {
    return hypergraph_of_edge(e).edgeHash(e);
  }

  PartitionID communityID(const HypernodeID u) const {
    return hypergraph_of_vertex(u).communityID(u);
  }

  bool nodeIsEnabled(const HypernodeID u) const {
    return hypergraph_of_vertex(u).nodeIsEnabled(u);
  }

  bool edgeIsEnabled(const HyperedgeID e) const {
    return hypergraph_of_edge(e).edgeIsEnabled(e);
  }

  void enableHypernode(const HypernodeID u) {
    hypergraph_of_vertex(u).enableHypernode(u);
  }

  void disableHypernode(const HypernodeID u) {
    hypergraph_of_vertex(u).disableHypernode(u);
  }

  void enableHyperedge(const HyperedgeID e) {
    hypergraph_of_edge(e).enableHyperedge(e);
  }

  void streamCommunityID(const HypernodeID hn, const PartitionID community_id) {
    hypergraph_of_vertex(hn).streamCommunityID(hn, community_id);
  }

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

    start = std::chrono::high_resolution_clock::now();
    _communities_num_hypernodes.assign(_num_communities, 0);
    for ( const HypernodeID& hn : nodes() ) {
      ASSERT(communityID(hn) < _num_communities);
      ++_communities_num_hypernodes[communityID(hn)];
    }
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("compute_num_community_hns", "Compute Num Community HNs",
      "initialize_communities", mt_kahypar::utils::Timer::Type::PREPROCESSING, 1, std::chrono::duration<double>(end - start).count());

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

  void resetPinsToOriginalNodeIds() {
    tbb::task_group group;
    for ( int node = 0; node < (int)_hypergraphs.size(); ++node ) {
      TBBNumaArena::instance().numa_task_arena(node).execute([&] {
        group.run([&, node] {
          _hypergraphs[node].resetPinsToOriginalNodeIds(_hypergraphs);
        });
      });
    }
    group.wait();
  }

  // ! Only for testing
  void disableHyperedge(const HyperedgeID e) {
    hypergraph_of_edge(e).disableHyperedge(e);
  }

 private:

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void reverseContraction(const Memento& memento) {
    enableHypernode(memento.v);
    // TODO(heuer): update part id of v
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void markAllIncidentNetsOf(const HypernodeID v) {
    hypergraph_of_vertex(v).markAllIncidentNetsOf(v, _hypergraphs);
  }

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
    tbb::task_group group;
    std::vector<HypernodeID> tmp_node_mapping(_num_hypernodes);
    for ( HypernodeID node = 0; node < num_streaming_hypergraphs; ++node ) {
      TBBNumaArena::instance().numa_task_arena(node).execute([&] {
        group.run([&, node] {
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
    group.wait();
    _node_mapping = std::move(tmp_node_mapping);
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("stream_hypernodes", "Stream Hypernodes",
      "initialize_hypernodes", mt_kahypar::utils::Timer::Type::IMPORT, 1, std::chrono::duration<double>(end - start).count());

    // Initialize hypernodes on each streaming hypergraph
    // NOTE, that also involves streaming local incident nets to other
    // streaming hypergraphs
    start = std::chrono::high_resolution_clock::now();
    for ( size_t node = 0; node < num_streaming_hypergraphs; ++node ) {
      group.run([&, node] {
        TBBNumaArena::instance().numa_task_arena(node).execute([&] {
          _hypergraphs[node].initializeHypernodes(_hypergraphs, _node_mapping);
        });
      });
    }
    group.wait();
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
      group.run([&, node] {
        TBBNumaArena::instance().numa_task_arena(node).execute([&] {
          _hypergraphs[node].initializeIncidentNets();
        });
      });
    }
    group.wait();
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
  }

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

  // ! Number of hypernodes
  HypernodeID _num_hypernodes;
  // ! Number of hyperedges
  HyperedgeID _num_hyperedges;
  // ! Number of pins
  HypernodeID _num_pins;
  // ! Number of communities
  PartitionID _num_communities;

  // ! Number of hypernodes in a community
  parallel::scalable_vector<HypernodeID> _communities_num_hypernodes;
  // ! Number of pins in a community
  parallel::scalable_vector<HypernodeID> _communities_num_pins;

  std::vector<StreamingHypergraph> _hypergraphs;
  std::vector<HypernodeID> _node_mapping;

 
};

} // namespace ds
} // namespace mt_kahypar
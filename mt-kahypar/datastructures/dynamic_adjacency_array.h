/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
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

#include <cstddef>

#include "tbb/enumerable_thread_specific.h"

#include "kahypar/datastructure/fast_reset_flag_array.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/stl/scalable_unique_ptr.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/range.h"

namespace mt_kahypar {
namespace ds {

// forward declaration
class DynamicAdjacencyArray;

// Iterator over the incident edges of a vertex u
class IncidentEdgeIterator :
  public std::iterator<std::forward_iterator_tag,    // iterator_category
                        HyperedgeID,   // value_type
                        std::ptrdiff_t,   // difference_type
                        const HyperedgeID*,   // pointer
                        HyperedgeID> {   // reference
  public:
  IncidentEdgeIterator(const HypernodeID u,
                      const DynamicAdjacencyArray* dynamic_adjacency_array,
                      const size_t pos,
                      const bool end);

  HyperedgeID operator* () const;

  IncidentEdgeIterator & operator++ ();

  IncidentEdgeIterator operator++ (int) {
    IncidentEdgeIterator copy = *this;
    operator++ ();
    return copy;
  }

  bool operator!= (const IncidentEdgeIterator& rhs);

  bool operator== (const IncidentEdgeIterator& rhs);

  private:
  void traverse_headers();

  void skip_invalid();

  HypernodeID _u;
  HypernodeID _current_u;
  HypernodeID _current_size;
  HyperedgeID _current_pos;
  const DynamicAdjacencyArray* _dynamic_adjacency_array;
  bool _end;
};

// Iterator over all edges
class EdgeIterator :
  public std::iterator<std::forward_iterator_tag,    // iterator_category
                        HyperedgeID,   // value_type
                        std::ptrdiff_t,   // difference_type
                        const HyperedgeID*,   // pointer
                        HyperedgeID> {   // reference
  public:
  EdgeIterator(const HypernodeID u,
               const DynamicAdjacencyArray* dynamic_adjacency_array,
               std::function<bool (const HypernodeID)> filter);

  HyperedgeID operator* () const;

  EdgeIterator & operator++ ();

  EdgeIterator operator++ (int) {
    EdgeIterator copy = *this;
    operator++ ();
    return copy;
  }

  bool operator!= (const EdgeIterator& rhs);

  bool operator== (const EdgeIterator& rhs);

  private:
  void update_ids();

  void traverse_headers();

  HypernodeID _current_head;
  HypernodeID _current_u;
  HyperedgeID _current_id;
  HyperedgeID _current_last_id;
  const DynamicAdjacencyArray* _dynamic_adjacency_array;
  std::function<bool (const HypernodeID)> _filter;
};

class DynamicAdjacencyArray {
  using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
  using EdgeVector = parallel::scalable_vector<std::pair<HypernodeID, HypernodeID>>;
  using ThreadLocalCounter = tbb::enumerable_thread_specific<parallel::scalable_vector<size_t>>;
  using AtomicCounter = parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>>;

  using NodeFilterFunc = std::function<bool (const HypernodeID)>;
  using AcquireLockFunc = std::function<void (const HypernodeID)>;
  using ReleaseLockFunc = std::function<void (const HypernodeID)>;
  using CaseOneFunc = std::function<void (const HyperedgeID)>;
  using CaseTwoFunc = std::function<void (const HyperedgeID)>;
  #define FILTER_NONE_FUNC [] (const HypernodeID) { return true; }
  #define NOOP_LOCK_FUNC [] (const HypernodeID) { }

  static_assert(sizeof(char) == 1);

 public:
  // Represents one edge of a vertex.
  // An edge is associated with a version number. Edges
  // with a version number greater or equal than the version number in
  // header (see Header -> current_version) are active.
  struct Edge {

    bool operator== (const Edge& rhs) const {
      return target == rhs.target && source == rhs.source &&
             weight == rhs.weight && version == rhs.version;
    }

    bool operator!= (const Edge& rhs) const {
      return !operator==(rhs);
    }

    bool isSinglePin() const {
      return source == target;
    }

    bool isValid() const {
      return target != kInvalidHypernode;
    }

    void setValid(bool valid) {
      if (valid) {
        target = source;
      } else {
        ASSERT(isSinglePin());
        target = kInvalidHypernode;
      }
    }

    // ! Index of target node
    HypernodeID target;
    // ! Index of source node
    HypernodeID source;
    // ! edge weight
    HyperedgeWeight weight;
    // ! version for undoing contractions
    HypernodeID version; // TODO(maas): is it better to keep all edges?
    // ! the header of the original target
    HypernodeID original_target;
    HypernodeID original_source; // TODO(maas): we wouldn't need that without `edgeIsEnabled`
    HyperedgeID unique_id; // TODO(maas): any other solution?
  };

  struct RemovedEdgesOrWeight {
    bool is_weight;
    HypernodeID header;

    RemovedEdgesOrWeight() { };

    explicit RemovedEdgesOrWeight(bool is_weight, HypernodeID header, HyperedgeWeight weight, HypernodeID target):
      is_weight(is_weight),
      header(header),
      _num_removed_or_weight(static_cast<HyperedgeID>(weight)),
      _degree_diff_or_target(target) {
      ASSERT(weight >= 0);
    }

    static RemovedEdgesOrWeight asEdges(HypernodeID header, HyperedgeID num_removed, HyperedgeID degree_diff) {
      return RemovedEdgesOrWeight(false, header, static_cast<HyperedgeWeight>(num_removed), static_cast<HypernodeID>(degree_diff));
    }

    static RemovedEdgesOrWeight asWeight(HypernodeID header, HyperedgeWeight weight, HypernodeID target) {
      return RemovedEdgesOrWeight(true, header, weight, target);
    }

    HyperedgeID numRemoved() const {
      ASSERT(!is_weight);
      return _num_removed_or_weight;
    }

    HyperedgeID degreeDiff() const {
      ASSERT(!is_weight);
      return static_cast<HyperedgeID>(_degree_diff_or_target);
    }

    HyperedgeWeight weight() const {
      ASSERT(is_weight);
      return static_cast<HyperedgeWeight>(_num_removed_or_weight);
    }

    HypernodeID target() const {
      ASSERT(is_weight);
      return _degree_diff_or_target;
    }

   private:
    // some hand-made "union" fields
    HyperedgeID _num_removed_or_weight;
    HypernodeID _degree_diff_or_target;
  };

 private:
  // Header of the incident edge list of a vertex. The incident edge lists
  // contracted into one vertex are concatenated in a double linked list.
  struct Header {
    Header() :
      prev(0),
      next(0),
      it_prev(0),
      it_next(0),
      tail(0),
      first_active(0),
      first_inactive(0),
      degree(0),
      current_version(0),
      is_head(false) { }

    explicit Header(const HypernodeID u) :
      prev(u),
      next(u),
      it_prev(u),
      it_next(u),
      tail(u),
      first_active(0),
      first_inactive(0),
      degree(0),
      current_version(0),
      is_head(true) { }
    
    HyperedgeID size() const {
      return first_inactive - first_active;
    }

    // TODO: remove
    void print() const {
      LOG << V(prev) << ", " << V(next) << ", " << V(it_prev) << ", " << V(it_next) << ", " << V(is_head);
      // LOG << V(header.tail) << ", " << V(header.first_active) << ", " << V(header.first_inactive) << V(header.degree);
      // LOG << V(header.current_version) << V(header.is_head);
    }

    // ! Previous incident edge list
    HypernodeID prev;
    // ! Next incident edge list
    HypernodeID next;
    // ! Previous non-empty incident edge list
    HypernodeID it_prev;
    // ! Next non-empty incident edge list
    HypernodeID it_next;
    // ! If we append a vertex v to the incident edge list of a vertex u, we store
    // ! the previous tail of vertex v, such that we can restore the list of v
    // ! during uncontraction
    HypernodeID tail;
    // ! Id of the first edge
    HyperedgeID first;
    // ! All incident edges between [first_active, first_inactive) are active
    HyperedgeID first_active;
    // ! All incident edges between [first_active, first_inactive) are active
    HyperedgeID first_inactive;
    // ! Degree of the vertex
    HyperedgeID degree;
    // ! Current version of the incident edge list
    HypernodeID current_version;
    // ! True, if the vertex is the head of a incident edge list
    bool is_head;
  };

  // Used for detecting parallel edges.
  // Represents one edge with the required information
  // for detecting duplicates and removing the represented edge.
  struct ParallelEdgeInformation {
    ParallelEdgeInformation() = default;

    ParallelEdgeInformation(HypernodeID target, HypernodeID original_target,
                            HyperedgeID edge_id, HypernodeID header_id):
        target(target), original_target(original_target), edge_id(edge_id), header_id(header_id) { }

    // ! Index of target node
    HypernodeID target;
    // ! heade of target node
    HypernodeID original_target;
    // ! Index of corresponding edge
    HyperedgeID edge_id;
    // ! header
    HypernodeID header_id;
  };

  using ThreadLocalParallelEdgeVector = tbb::enumerable_thread_specific<vec<ParallelEdgeInformation>>;

 public:
  using const_iterator = IncidentEdgeIterator;

  DynamicAdjacencyArray() :
    _num_nodes(0),
    _header_array(),
    _edges(),
    _degree_diffs() { }

  DynamicAdjacencyArray(const HypernodeID num_nodes,
                        const EdgeVector& edge_vector,
                        const HyperedgeWeight* edge_weight = nullptr) :
    _num_nodes(num_nodes),
    _header_array(),
    _edges(),
    _thread_local_vec(),
    _degree_diffs() {
    tbb::parallel_invoke([&] {
        construct(edge_vector, edge_weight);
      }, [&] {
      _degree_diffs.resize(_num_nodes);
      }
    );
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Edge& edge(const HyperedgeID e) const {
    ASSERT(e < _edges.size(), "Edge" << e << "does not exist");
    return _edges[e];
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Edge& edge(const HyperedgeID e) {
    ASSERT(e <= _edges.size(), "Edge" << e << "does not exist");
    return _edges[e];
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeID numNodes() const {
    return _num_nodes;
  }

  // ! Degree of the vertex
  HypernodeID nodeDegree(const HypernodeID u) const {
    ASSERT(u < _num_nodes, "Hypernode" << u << "does not exist");
    return header(u).degree;
  }

  bool edgeIsEnabled(const HyperedgeID e) const {
    const Header& head = header(edge(e).original_source);
    return head.first_active <= e && e < head.first_inactive;
  }

  // ! Returns a range to loop over the incident edges of hypernode u.
  IteratorRange<IncidentEdgeIterator> incidentEdges(const HypernodeID u) const {
    ASSERT(u < _num_nodes, "Hypernode" << u << "does not exist");
    return IteratorRange<IncidentEdgeIterator>(
      IncidentEdgeIterator(u, this, 0UL, false),
      IncidentEdgeIterator(u, this, 0UL, true));
  }

  // ! Returns a range to loop over the incident edges of hypernode u.
  IteratorRange<IncidentEdgeIterator> incidentEdges(const HypernodeID u,
                                                   const size_t pos) const {
    ASSERT(u < _num_nodes, "Hypernode" << u << "does not exist");
    return IteratorRange<IncidentEdgeIterator>(
      IncidentEdgeIterator(u, this, pos, false),
      IncidentEdgeIterator(u, this, 0UL, true));
  }

  // ! Returns a range to loop over all edges.
  IteratorRange<EdgeIterator> edges(NodeFilterFunc filter = FILTER_NONE_FUNC) const {
    return IteratorRange<EdgeIterator>(
      EdgeIterator(0, this, filter),
      EdgeIterator(_num_nodes, this, filter));
  }


  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdgesOfNode(HypernodeID hn, const F& f) const {
    // TODO(maas): might not be optimal
    for (HypernodeID current_u: headers(hn)) {
      tbb::parallel_for(firstActiveEdge(current_u), firstInactiveEdge(current_u), [&](const HyperedgeID& e) {
        f(e);
      });
    }
  }

  // ! Contracts two incident list of u and v, whereby u is the representative and
  // ! v the contraction partner of the contraction. The contraction involves to remove
  // ! all incident edges shared between u and v from the incident edge list of v and append
  // ! the list of v to u.
  void contract(const HypernodeID u,
                const HypernodeID v,
                const AcquireLockFunc& acquire_lock = NOOP_LOCK_FUNC,
                const ReleaseLockFunc& release_lock = NOOP_LOCK_FUNC);

  // ! Uncontract two previously contracted vertices u and v.
  // ! Uncontraction involves to decrement the version number of all incident lists contained
  // ! in v and restore all incident edges with a version number equal to the new version.
  // ! Note, uncontraction must be done in relative contraction order
  void uncontract(const HypernodeID u,
                  const HypernodeID v,
                  const AcquireLockFunc& acquire_lock = NOOP_LOCK_FUNC,
                  const ReleaseLockFunc& release_lock = NOOP_LOCK_FUNC);

  // ! Uncontract two previously contracted vertices u and v.
  // ! Uncontraction involves to decrement the version number of all incident lists contained
  // ! in v and restore all incident edges with a version number equal to the new version.
  // ! Additionally it calls case_one_func for a hyperedge he, if u and v were previously both
  // ! adjacent to he and case_two_func if only v was previously adjacent to he.
  // ! Note, uncontraction must be done in relative contraction order
  void uncontract(const HypernodeID u,
                  const HypernodeID v,
                  const CaseOneFunc& case_one_func,
                  const CaseTwoFunc& case_two_func,
                  const AcquireLockFunc& acquire_lock,
                  const ReleaseLockFunc& release_lock);

  parallel::scalable_vector<RemovedEdgesOrWeight> removeSinglePinAndParallelEdges();

  void restoreSinglePinAndParallelEdges(const parallel::scalable_vector<RemovedEdgesOrWeight>& edges_to_restore);

  DynamicAdjacencyArray copy(parallel_tag_t) const;

  DynamicAdjacencyArray copy() const;

  void reset();

  void sortIncidentEdges();

  size_t size_in_bytes() const {
    return _edges.size() * sizeof(Edge)
      + _degree_diffs.size() * sizeof(int32_t)
      + _header_array.size() * sizeof(Header);
  }

 private:
  friend class IncidentEdgeIterator;
  friend class EdgeIterator;

  class HeaderIterator :
    // TODO
    public std::iterator<std::forward_iterator_tag,    // iterator_category
                          HypernodeID,   // value_type
                          std::ptrdiff_t,   // difference_type
                          const HypernodeID*,   // pointer
                          HypernodeID> {   // reference
    public:
    HeaderIterator(const HypernodeID u,
                   const DynamicAdjacencyArray* dynamic_adjacency_array,
                   const bool end):
      _u(u),
      _current_u(u),
      _dynamic_adjacency_array(dynamic_adjacency_array),
      _end(end) { }

    HypernodeID operator* () const {
      return _current_u;
    }

    HeaderIterator & operator++ () {
      _current_u = _dynamic_adjacency_array->header(_current_u).next;
      if (_current_u == _u) {
        _end = true;
      }
      return *this;
    }

    HeaderIterator operator++ (int) {
      HeaderIterator copy = *this;
      operator++ ();
      return copy;
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool operator== (const HeaderIterator& rhs) {
      return _u == rhs._u && _end == rhs._end;
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool operator!= (const HeaderIterator& rhs) {
      return !(*this == rhs);
    }

    private:
    HypernodeID _u;
    HypernodeID _current_u;
    const DynamicAdjacencyArray* _dynamic_adjacency_array;
    bool _end;
  };

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Header& header(const HypernodeID u) const {
    ASSERT(u <= _num_nodes, "Hypernode" << u << "does not exist");
    return _header_array[u];
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Header& header(const HypernodeID u) {
    ASSERT(u <= _num_nodes, "Hypernode" << u << "does not exist");
    return _header_array[u];
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HyperedgeID firstEdge(const HypernodeID u) const {
    ASSERT(u <= _num_nodes, "Hypernode" << u << "does not exist");
    return header(u).first;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HyperedgeID firstActiveEdge(const HypernodeID u) const {
    ASSERT(u <= _num_nodes, "Hypernode" << u << "does not exist");
    return header(u).first_active;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HyperedgeID firstInactiveEdge(const HypernodeID u) const {
    ASSERT(u <= _num_nodes, "Hypernode" << u << "does not exist");
    return header(u).first_inactive;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HyperedgeID lastEdge(const HypernodeID u) const {
    ASSERT(u <= _num_nodes, "Hypernode" << u << "does not exist");
    return header(u + 1).first;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void swap(Edge& lhs, Edge& rhs) {
    Edge tmp_lhs = lhs;
    lhs = rhs;
    rhs = tmp_lhs;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void swap_to_front(const HypernodeID u, const HyperedgeID e) {
    ASSERT(u <= _num_nodes, "Hypernode" << u << "does not exist");
    swap(edge(e), edge(firstActiveEdge(u)));
    ++header(u).first_active;
    ASSERT(header(u).first_active <= header(u).first_inactive);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void swap_to_back(const HypernodeID u, const HyperedgeID e) {
    ASSERT(u <= _num_nodes, "Hypernode" << u << "does not exist");
    swap(edge(e), edge(firstInactiveEdge(u) - 1));
    --header(u).first_inactive;
    ASSERT(header(u).first_active <= header(u).first_inactive);
  }

  // ! Returns a range to loop over the headers of node u.
  IteratorRange<HeaderIterator> headers(const HypernodeID u) const {
    ASSERT(u < _num_nodes, "Hypernode" << u << "does not exist");
    return IteratorRange<HeaderIterator>(
      HeaderIterator(u, this, false),
      HeaderIterator(u, this, true));
  }

  // ! Restores all previously removed incident edges
  // ! Note, function must be called in reverse order of calls to
  // ! removeIncidentEdges(...) and all uncontraction that happens
  // ! between two consecutive calls to removeIncidentEdges(...) must
  // ! be processed.
  void restoreIncidentEdges(const HypernodeID u,
                           const CaseOneFunc& case_one_func,
                           const CaseTwoFunc& case_two_func);

  void append(const HypernodeID u, const HypernodeID v);

  void splice(const HypernodeID u, const HypernodeID v);

  void removeEmptyIncidentEdgeList(const HypernodeID u);

  void restoreIteratorPointers(const HypernodeID u);

  void restoreItLink(const HypernodeID u, const HypernodeID prev, const HypernodeID current);

  HyperedgeID findBackwardsEdge(const Edge& forward, HypernodeID source) const;

  void streamWeight(StreamingVector<RemovedEdgesOrWeight>& tmp_removed_edges,
                    const ParallelEdgeInformation& e, HyperedgeWeight w);

  void construct(const EdgeVector& edge_vector, const HyperedgeWeight* edge_weight = nullptr);

  bool verifyIteratorPointers(const HypernodeID u) const;

  HypernodeID _num_nodes;
  Array<Header> _header_array;
  Array<Edge> _edges;
  // data used during parallel edge removal
  ThreadLocalParallelEdgeVector _thread_local_vec;
  Array<int32_t> _degree_diffs;
};

}  // namespace ds
}  // namespace mt_kahypar

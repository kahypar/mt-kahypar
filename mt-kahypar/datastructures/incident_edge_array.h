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
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/stl/scalable_unique_ptr.h"
#include "mt-kahypar/utils/range.h"

namespace mt_kahypar {
namespace ds {

// forward declaration
class IncidentEdgeArray;

// Iterator over the incident edges of a vertex u
class IncidentEdgeIterator :
  // TODO
  public std::iterator<std::forward_iterator_tag,    // iterator_category
                        HyperedgeID,   // value_type
                        std::ptrdiff_t,   // difference_type
                        const HyperedgeID*,   // pointer
                        HyperedgeID> {   // reference
  public:
  IncidentEdgeIterator(const HypernodeID u,
                      const IncidentEdgeArray* incident_edge_array,
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

  HypernodeID _u;
  HypernodeID _current_u;
  HypernodeID _current_size;
  HyperedgeID _current_pos;
  const IncidentEdgeArray* _incident_edge_array;
  bool _end;
};

class IncidentEdgeArray {
  using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
  using EdgeVector = parallel::scalable_vector<std::pair<HypernodeID, HypernodeID>>;
  using ThreadLocalCounter = tbb::enumerable_thread_specific<parallel::scalable_vector<size_t>>;
  using AtomicCounter = parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>>;

  using AcquireLockFunc = std::function<void (const HypernodeID)>;
  using ReleaseLockFunc = std::function<void (const HypernodeID)>;
  using CaseOneFunc = std::function<void (const HyperedgeID)>;
  using CaseTwoFunc = std::function<void (const HyperedgeID)>;
  #define NOOP_LOCK_FUNC [] (const HypernodeID) { }

  static_assert(sizeof(char) == 1);

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

    // ! Index of target node
    HypernodeID target;
    // ! Index of source node
    HypernodeID source;
    // ! edge weight
    HyperedgeWeight weight;
    // ! version for undoing contractions
    HypernodeID version;
  };
 
 public:
  // Header of the incident edge list of a vertex. The incident edge lists
  // contracted into one vertex are concatenated in a double linked list.
  struct Header {
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

   private:
    // ensure that sizeof(Header) is a multiple of sizeof(Edge)
    uint32_t __padding_0;
    uint32_t __padding_1;
  };

  static_assert(alignof(Header) == alignof(Edge));
  static_assert(sizeof(Header) % sizeof(Edge) == 0);

 public:
  using const_iterator = IncidentEdgeIterator;

  IncidentEdgeArray() :
    _num_nodes(0),
    _size_in_bytes(0),
    _index_array(),
    _incident_edge_array(nullptr) { }

  IncidentEdgeArray(const HypernodeID num_hypernodes,
                   const EdgeVector& edge_vector) :
    _num_nodes(num_hypernodes),
    _size_in_bytes(0),
    _index_array(),
    _incident_edge_array(nullptr)  {
    construct(edge_vector);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Edge& edge(const HyperedgeID e) const {
    ASSERT(e <= _size_in_bytes / sizeof(Edge), "Edge" << e << "does not exist");
    return _incident_edge_array.get()[e];
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Edge& edge(const HyperedgeID e) {
    ASSERT(e <= _size_in_bytes / sizeof(Edge), "Edge" << e << "does not exist");
    return _incident_edge_array.get()[e];
  }

  // ! Degree of the vertex
  HypernodeID nodeDegree(const HypernodeID u) const {
    ASSERT(u < _num_nodes, "Hypernode" << u << "does not exist");
    return header(u)->degree;
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

  // ! Contracts two incident list of u and v, whereby u is the representative and
  // ! v the contraction partner of the contraction. The contraction involves to remove
  // ! all incident edges shared between u and v from the incident edge list of v and append
  // ! the list of v to u.
  void contract(const HypernodeID u,
                const HypernodeID v,
                kahypar::ds::FastResetFlagArray<>& edge_bitset,
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

  // ! Removes all incidents nets of u flagged in hes_to_remove.
  void removeIncidentEdges(const HypernodeID u,
                          const kahypar::ds::FastResetFlagArray<>& hes_to_remove);

  // ! Restores all previously removed incident edges
  // ! Note, function must be called in reverse order of calls to
  // ! removeIncidentEdges(...) and all uncontraction that happens
  // ! between two consecutive calls to removeIncidentEdges(...) must
  // ! be processed.
  void restoreIncidentEdges(const HypernodeID u);

  IncidentEdgeArray copy(parallel_tag_t);

  IncidentEdgeArray copy();

  void reset();

  size_t size_in_bytes() const {
    return _size_in_bytes + sizeof(size_t) * _index_array.size();
  }

 private:
  friend class IncidentEdgeIterator;

  class HeaderIterator :
    // TODO
    public std::iterator<std::forward_iterator_tag,    // iterator_category
                          HypernodeID,   // value_type
                          std::ptrdiff_t,   // difference_type
                          const HypernodeID*,   // pointer
                          HypernodeID> {   // reference
    public:
    HeaderIterator(const HypernodeID u,
                   const IncidentEdgeArray* incident_edge_array,
                   const bool end):
      _u(u),
      _current_u(u),
      _incident_edge_array(incident_edge_array),
      _end(end) { }

    HypernodeID operator* () const {
      return _current_u;
    }

    HeaderIterator & operator++ () {
      const Header* header = _incident_edge_array->header(_current_u);
      _current_u = header->next;
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
    const IncidentEdgeArray* _incident_edge_array;
    bool _end;
  };

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Header* header(const HypernodeID u) const {
    ASSERT(u <= _num_nodes, "Hypernode" << u << "does not exist");
    return reinterpret_cast<const Header*>(_incident_edge_array.get() + _index_array[u]);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Header* header(const HypernodeID u) {
    return const_cast<Header*>(static_cast<const IncidentEdgeArray&>(*this).header(u));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HyperedgeID firstEdge(const HypernodeID u) const {
    ASSERT(u <= _num_nodes, "Hypernode" << u << "does not exist");
    return _index_array[u] + sizeof(Header) / sizeof(Edge);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HyperedgeID firstActiveEdge(const HypernodeID u) const {
    ASSERT(u <= _num_nodes, "Hypernode" << u << "does not exist");
    return firstEdge(u) + header(u)->first_active;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HyperedgeID firstInactiveEdge(const HypernodeID u) const {
    ASSERT(u <= _num_nodes, "Hypernode" << u << "does not exist");
    return firstEdge(u) + header(u)->first_inactive;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HyperedgeID lastEdge(const HypernodeID u) const {
    ASSERT(u <= _num_nodes, "Hypernode" << u << "does not exist");
    return _index_array[u + 1];
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void swap(Edge& lhs, Edge& rhs) {
    Edge tmp_lhs = lhs;
    lhs = rhs;
    rhs = tmp_lhs;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void swap_to_front(const HypernodeID u, const HyperedgeID e) {
    ASSERT(u <= _num_nodes, "Hypernode" << u << "does not exist");
    swap(edge(e), edge(firstActiveEdge(u)));
    ++header(u)->first_active;
    ASSERT(header(u)->first_active <= header(u)->first_inactive);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void swap_to_back(const HypernodeID u, const HyperedgeID e) {
    ASSERT(u <= _num_nodes, "Hypernode" << u << "does not exist");
    swap(edge(e), edge(firstInactiveEdge(u) - 1));
    --header(u)->first_inactive;
    ASSERT(header(u)->first_active <= header(u)->first_inactive);
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

  void restoreItLink(const HypernodeID u, const HypernodeID prev, const HypernodeID current);

  std::pair<HypernodeID, HyperedgeID> findOutgoingEdge(const HypernodeID u, const HypernodeID target) const;

  std::pair<HypernodeID, HyperedgeID> findRemovedOutgoingEdge(const HypernodeID u, const HypernodeID target) const;

  void construct(const EdgeVector& edge_vector, const HyperedgeWeight* edge_weight = nullptr);

  bool verifyIteratorPointers(const HypernodeID u) const;

  HypernodeID _num_nodes;
  size_t _size_in_bytes;
  Array<HyperedgeID> _index_array;
  parallel::tbb_unique_ptr<Edge> _incident_edge_array;
};

}  // namespace ds
}  // namespace mt_kahypar

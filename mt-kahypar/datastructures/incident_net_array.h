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
#pragma once

#include <stddef.h>

#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_invoke.h"
#include "tbb/parallel_scan.h"

#include "kahypar/datastructure/fast_reset_flag_array.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/stl/scalable_unique_ptr.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/utils/range.h"

namespace mt_kahypar {
namespace ds {

// forward declaration
class IncidentNetArray;

// Iterator over the incident nets of a vertex u
class IncidentNetIterator :
  public std::iterator<std::forward_iterator_tag,    // iterator_category
                        HyperedgeID,   // value_type
                        std::ptrdiff_t,   // difference_type
                        const HyperedgeID*,   // pointer
                        HyperedgeID> {   // reference
  public:
  IncidentNetIterator(const HypernodeID u,
                      const IncidentNetArray* incident_net_array,
                      const size_t pos,
                      const bool end);

  HyperedgeID operator* () const;

  IncidentNetIterator & operator++ ();

  IncidentNetIterator operator++ (int) {
    IncidentNetIterator copy = *this;
    operator++ ();
    return copy;
  }

  bool operator!= (const IncidentNetIterator& rhs);

  bool operator== (const IncidentNetIterator& rhs);

  private:
  void next_iterator();

  HypernodeID _u;
  HypernodeID _current_u;
  HypernodeID _current_size;
  size_t _current_pos;
  const IncidentNetArray* _incident_net_array;
  bool _end;
};

// ! Class allows in-place contraction and uncontraction of the incident net array
class IncidentNetArray {

  using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
  using ThreadLocalCounter = tbb::enumerable_thread_specific<parallel::scalable_vector<size_t>>;
  using AtomicCounter = parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>>;

  using AcquireLockFunc = std::function<void (const HypernodeID)>;
  using ReleaseLockFunc = std::function<void (const HypernodeID)>;
  using CaseOneFunc = std::function<void (const HyperedgeID)>;
  using CaseTwoFunc = std::function<void (const HyperedgeID)>;
  #define NOOP_LOCK_FUNC [] (const HypernodeID) { }

  static_assert(sizeof(char) == 1);

  // Represents one incident net of a vertex.
  // A incident net is associated with a version number. Incident nets
  // with a version number greater or equal than the version number in
  // header (see Header -> current_version) are active.
  struct Entry {
    HyperedgeID e;
    HypernodeID version;
  };

  // Header of the incident net list of a vertex. The incident net lists
  // contracted into one vertex are concatenated in a double linked list.
  struct Header {
    explicit Header(const HypernodeID u) :
      prev(u),
      next(u),
      it_prev(u),
      it_next(u),
      tail(u),
      size(0),
      degree(0),
      current_version(0),
      is_head(true) { }

    // ! Previous incident net list
    HypernodeID prev;
    // ! Next incident net list
    HypernodeID next;
    // ! Previous non-empty incident net list
    HypernodeID it_prev;
    // ! Next non-empty incident net list
    HypernodeID it_next;
    // ! If we append a vertex v to the incident net list of a vertex u, we store
    // ! the previous tail of vertex v, such that we can restore the list of v
    // ! during uncontraction
    HypernodeID tail;
    // ! All incident nets between [0,size) are active
    HypernodeID size;
    // ! Degree of the vertex
    HypernodeID degree;
    // ! Current version of the incident net list
    HypernodeID current_version;
    // ! True, if the vertex is the head of a incident net list
    bool is_head;
  };

 public:
  using const_iterator = IncidentNetIterator;

  IncidentNetArray() :
    _num_hypernodes(0),
    _size_in_bytes(0),
    _index_array(),
    _incident_net_array(nullptr) { }

  IncidentNetArray(const HypernodeID num_hypernodes,
                   const HyperedgeVector& edge_vector) :
    _num_hypernodes(num_hypernodes),
    _size_in_bytes(0),
    _index_array(),
    _incident_net_array(nullptr)  {
    construct(edge_vector);
  }

  // ! Degree of the vertex
  HypernodeID nodeDegree(const HypernodeID u) const {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return header(u)->degree;
  }

  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetIterator> incidentEdges(const HypernodeID u) const {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return IteratorRange<IncidentNetIterator>(
      IncidentNetIterator(u, this, 0UL, false),
      IncidentNetIterator(u, this, 0UL, true));
  }

  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetIterator> incidentEdges(const HypernodeID u,
                                                   const size_t pos) const {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return IteratorRange<IncidentNetIterator>(
      IncidentNetIterator(u, this, pos, false),
      IncidentNetIterator(u, this, 0UL, true));
  }

  // ! Contracts two incident list of u and v, whereby u is the representative and
  // ! v the contraction partner of the contraction. The contraction involves to remove
  // ! all incident nets shared between u and v from the incident net list of v and append
  // ! the list of v to u.
  void contract(const HypernodeID u,
                const HypernodeID v,
                const kahypar::ds::FastResetFlagArray<>& shared_hes_of_u_and_v,
                const AcquireLockFunc& acquire_lock = NOOP_LOCK_FUNC,
                const ReleaseLockFunc& release_lock = NOOP_LOCK_FUNC);

  // ! Uncontract two previously contracted vertices u and v.
  // ! Uncontraction involves to decrement the version number of all incident lists contained
  // ! in v and restore all incident nets with a version number equal to the new version.
  // ! Note, uncontraction must be done in relative contraction order
  void uncontract(const HypernodeID u,
                  const HypernodeID v,
                  const AcquireLockFunc& acquire_lock = NOOP_LOCK_FUNC,
                  const ReleaseLockFunc& release_lock = NOOP_LOCK_FUNC);

  // ! Uncontract two previously contracted vertices u and v.
  // ! Uncontraction involves to decrement the version number of all incident lists contained
  // ! in v and restore all incident nets with a version number equal to the new version.
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
  void removeIncidentNets(const HypernodeID u,
                          const kahypar::ds::FastResetFlagArray<>& hes_to_remove);

  // ! Restores all previously removed incident nets
  // ! Note, function must be called in reverse order of calls to
  // ! removeIncidentNets(...) and all uncontraction that happens
  // ! between two consecutive calls to removeIncidentNets(...) must
  // ! be processed.
  void restoreIncidentNets(const HypernodeID u);

  IncidentNetArray copy(const TaskGroupID);

  IncidentNetArray copy();

  void reset();

  size_t size_in_bytes() const {
    return _size_in_bytes + sizeof(size_t) * _index_array.size();
  }

 private:
  friend class IncidentNetIterator;

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Header* header(const HypernodeID u) const {
    ASSERT(u <= _num_hypernodes, "Hypernode" << u << "does not exist");
    return reinterpret_cast<const Header*>(_incident_net_array.get() + _index_array[u]);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Header* header(const HypernodeID u) {
    return const_cast<Header*>(static_cast<const IncidentNetArray&>(*this).header(u));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Entry* firstEntry(const HypernodeID u) const {
    ASSERT(u <= _num_hypernodes, "Hypernode" << u << "does not exist");
    return reinterpret_cast<const Entry*>(_incident_net_array.get() + _index_array[u] + sizeof(Header));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Entry* firstEntry(const HypernodeID u) {
    return const_cast<Entry*>(static_cast<const IncidentNetArray&>(*this).firstEntry(u));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Entry* lastEntry(const HypernodeID u) const {
    ASSERT(u <= _num_hypernodes, "Hypernode" << u << "does not exist");
    return reinterpret_cast<const Entry*>(_incident_net_array.get() +
      _index_array[u] + sizeof(Header) + header(u)->size * sizeof(Entry));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Entry* lastEntry(const HypernodeID u) {
    return const_cast<Entry*>(static_cast<const IncidentNetArray&>(*this).lastEntry(u));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void swap(Entry* lhs, Entry* rhs) {
    Entry tmp_lhs = *lhs;
    *lhs = *rhs;
    *rhs = tmp_lhs;
  }

  // ! Restores all previously removed incident nets
  // ! Note, function must be called in reverse order of calls to
  // ! removeIncidentNets(...) and all uncontraction that happens
  // ! between two consecutive calls to removeIncidentNets(...) must
  // ! be processed.
  void restoreIncidentNets(const HypernodeID u,
                           const CaseOneFunc& case_one_func,
                           const CaseTwoFunc& case_two_func);

  void append(const HypernodeID u, const HypernodeID v);

  void splice(const HypernodeID u, const HypernodeID v);

  void removeEmptyIncidentNetList(const HypernodeID u);

  void construct(const HyperedgeVector& edge_vector);

  bool verifyIteratorPointers(const HypernodeID u) const;

  HypernodeID _num_hypernodes;
  size_t _size_in_bytes;
  Array<size_t> _index_array;
  parallel::tbb_unique_ptr<char> _incident_net_array;
};

}  // namespace ds
}  // namespace mt_kahypar
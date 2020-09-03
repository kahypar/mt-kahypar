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

class IncidentNetArray {

  using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
  using ThreadLocalCounter = tbb::enumerable_thread_specific<parallel::scalable_vector<size_t>>;
  using AtomicCounter = parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>>;

  static_assert(sizeof(char) == 1);

  struct Entry {
    HyperedgeID e;
    HypernodeID version;
  };

  struct Header {
    explicit Header(const HypernodeID u) :
      prev(u),
      next(u),
      it_prev(u),
      it_next(u),
      tail(u),
      size(0),
      current_version(0) { }

    HypernodeID prev;
    HypernodeID next;
    HypernodeID it_prev;
    HypernodeID it_next;
    HypernodeID tail;
    HypernodeID size;
    HypernodeID current_version;
  };

  class IncidentNetIterator :
    public std::iterator<std::forward_iterator_tag,    // iterator_category
                         HyperedgeID,   // value_type
                         std::ptrdiff_t,   // difference_type
                         const HyperedgeID*,   // pointer
                         HyperedgeID> {   // reference
   public:
    IncidentNetIterator(const HypernodeID u,
                        const IncidentNetArray* incident_net_array,
                        const bool end) :
      _u(u),
      _current_u(u),
      _last_u(incident_net_array->header(u)->it_prev),
      _current_pos(0),
      _incident_net_array(incident_net_array) {
      if ( end ) {
        _current_u = _last_u;
        _current_pos = incident_net_array->header(_current_u)->size;
      }

      if ( !end && _current_pos == _incident_net_array->header(_current_u)->size ) {
        next_iterator();
      }
    }

    HyperedgeID operator* () const {
      ASSERT(_current_u != _last_u || _current_pos != _incident_net_array->header(_current_u)->size);
      return _incident_net_array->firstEntry(_current_u)[_current_pos].e;
    }

    IncidentNetIterator & operator++ () {
      ASSERT(_current_u != _last_u || _current_pos != _incident_net_array->header(_current_u)->size);
      ++_current_pos;
      if ( _current_pos == _incident_net_array->header(_current_u)->size) {
        next_iterator();
      }
      return *this;
    }

    IncidentNetIterator operator++ (int) {
      IncidentNetIterator copy = *this;
      operator++ ();
      return copy;
    }

    bool operator!= (const IncidentNetIterator& rhs) {
      return _u != rhs._u || _current_u != rhs._current_u ||
             _current_pos != rhs._current_pos;
    }

    bool operator== (const IncidentNetIterator& rhs) {
      return _u == rhs._u && _current_u == rhs._current_u &&
             _current_pos == rhs._current_pos;
    }

   private:
    void next_iterator() {
      while ( _current_pos == _incident_net_array->header(_current_u)->size ) {
        if ( _current_u == _last_u ) {
          break;
        }
        _current_u = _incident_net_array->header(_current_u)->it_next;
        _current_pos = 0;
      }
    }

    HypernodeID _u;
    HypernodeID _current_u;
    HypernodeID _last_u;
    size_t _current_pos;
    const IncidentNetArray* _incident_net_array;
  };

 public:
  IncidentNetArray() :
    _num_hypernodes(0),
    _index_array(),
    _incident_net_array(nullptr) { }

  IncidentNetArray(const HypernodeID num_hypernodes,
                   const HyperedgeVector& edge_vector) :
    _num_hypernodes(num_hypernodes),
    _index_array(),
    _incident_net_array(nullptr) {
    construct(edge_vector);
  }

  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetIterator> incidentEdges(const HypernodeID u) const {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return IteratorRange<IncidentNetIterator>(
      IncidentNetIterator(u, this, false),
      IncidentNetIterator(u, this, true));
  }

  void contract(const HypernodeID u,
                const HypernodeID v,
                const kahypar::ds::FastResetFlagArray<>& shared_hes_of_u_and_v) {
    HypernodeID current_v = v;
    do {
      Header* head = header(current_v);
      const HypernodeID new_version = ++head->current_version;
      Entry* last_entry = lastEntry(current_v);
      for ( Entry* current_entry = firstEntry(current_v); current_entry != last_entry; ++current_entry ) {
        if ( shared_hes_of_u_and_v[current_entry->e] ) {
          swap(current_entry--, --last_entry);
          ASSERT(head->size > 0);
          --head->size;
        } else {
          current_entry->version = new_version;
        }
      }

      if ( head->size == 0 && current_v != v ) {
        removeEmptyIncidentNetList(current_v);
      }
      current_v = head->next;
    } while ( current_v != v );
    append(u, v);
  }

  void uncontract(const HypernodeID v) {
    ASSERT(header(v)->prev != v);
    splice(v);
    HypernodeID current_v = v;
    HypernodeID last_non_empty_entry = kInvalidHypernode;
    HypernodeID non_empty_entry_prev_v = kInvalidHypernode;
    HypernodeID non_empty_entry_next_v = kInvalidHypernode;
    do {
      Header* head = header(current_v);
      const HypernodeID size_before = head->size;
      ASSERT(head->current_version > 0);
      const HypernodeID new_version = --head->current_version;
      const Entry* last_entry = reinterpret_cast<const Entry*>(header(current_v + 1));
      for ( Entry* current_entry = lastEntry(current_v); current_entry != last_entry; ++current_entry ) {
        if ( current_entry->version == new_version ) {
          ++head->size;
        } else {
          break;
        }
      }

      if ( head->size > 0 ) {
        const bool was_non_empty_before = size_before > 0;
        if ( was_non_empty_before ) {
          non_empty_entry_prev_v = non_empty_entry_prev_v == kInvalidHypernode ?
            head->it_prev : non_empty_entry_prev_v;
          non_empty_entry_next_v = head->it_next;
        }

        if ( last_non_empty_entry != kInvalidHypernode &&
            head->it_prev != last_non_empty_entry ) {
          header(last_non_empty_entry)->it_next = current_v;
          head->it_prev = last_non_empty_entry;
        }
        last_non_empty_entry = current_v;
      }
      current_v = head->next;
    } while ( current_v != v );

    ASSERT(header(v)->size > 0);
    ASSERT(last_non_empty_entry != kInvalidHypernode);
    header(v)->it_prev = last_non_empty_entry;
    header(last_non_empty_entry)->it_next = v;

    ASSERT((non_empty_entry_next_v == kInvalidHypernode && non_empty_entry_prev_v == kInvalidHypernode) ||
           (non_empty_entry_next_v != kInvalidHypernode && non_empty_entry_prev_v != kInvalidHypernode));
    if ( non_empty_entry_next_v != kInvalidHypernode && non_empty_entry_prev_v != kInvalidHypernode ) {
      header(non_empty_entry_prev_v)->it_next = non_empty_entry_next_v;
      header(non_empty_entry_next_v)->it_prev = non_empty_entry_prev_v;
    }
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

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void append(const HypernodeID u, const HypernodeID v) {
    const HypernodeID tail_u = header(u)->prev;
    const HypernodeID tail_v = header(v)->prev;
    header(tail_u)->next = v;
    header(u)->prev = tail_v;
    header(v)->tail = tail_v;
    header(v)->prev = tail_u;
    header(tail_v)->next = u;

    const HypernodeID it_tail_u = header(u)->it_prev;
    const HypernodeID it_tail_v = header(v)->it_prev;
    header(it_tail_u)->it_next = v;
    header(u)->it_prev = it_tail_v;
    header(v)->it_prev = it_tail_u;
    header(it_tail_v)->it_next = u;

    if ( header(v)->size == 0 ) {
      removeEmptyIncidentNetList(v);
    }
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void splice(const HypernodeID v) {
    const HypernodeID tail = header(v)->tail;
    const HypernodeID prev_v = header(v)->prev;
    const HypernodeID next_tail = header(tail)->next;
    header(v)->prev = tail;
    header(tail)->next = v;
    header(next_tail)->prev = prev_v;
    header(prev_v)->next = next_tail;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void removeEmptyIncidentNetList(const HypernodeID u) {
    ASSERT(header(u)->size == 0, V(u) << V(header(u)->size));
    Header* head = header(u);
    header(head->it_prev)->it_next = head->it_next;
    header(head->it_next)->it_prev = head->it_prev;
    head->it_next = u;
    head->it_prev = u;
  }

  void construct(const HyperedgeVector& edge_vector) {
    // Accumulate degree of each vertex thread local
    const HyperedgeID num_hyperedges = edge_vector.size();
    ThreadLocalCounter local_incident_nets_per_vertex(_num_hypernodes + 1, 0);
    AtomicCounter current_incident_net_pos;
    tbb::parallel_invoke([&] {
      tbb::parallel_for(ID(0), num_hyperedges, [&](const size_t pos) {
        parallel::scalable_vector<size_t>& num_incident_nets_per_vertex =
          local_incident_nets_per_vertex.local();
        for ( const HypernodeID& pin : edge_vector[pos] ) {
          ASSERT(pin < _num_hypernodes, V(pin) << V(_num_hypernodes));
          ++num_incident_nets_per_vertex[pin + 1];
        }
      });
    }, [&] {
      _index_array.assign(_num_hypernodes + 1, 0);
      current_incident_net_pos.assign(
        _num_hypernodes, parallel::IntegralAtomicWrapper<size_t>(0));
    });

    // We sum up the number of incident nets per vertex only thread local.
    // To obtain the global number of incident nets per vertex, we iterate
    // over each thread local counter and sum it up.
    bool first_iteration = true;
    for ( const parallel::scalable_vector<size_t>& c : local_incident_nets_per_vertex ) {
      tbb::parallel_for(ID(0), _num_hypernodes + 1, [&](const size_t pos) {
        _index_array[pos] += c[pos] * sizeof(Entry) + (first_iteration ? sizeof(Header) : 0);
      });
      first_iteration = false;
    }

    // Compute start positon of the incident nets of each vertex via a parallel prefix sum
    parallel::TBBPrefixSum<size_t> incident_net_prefix_sum(_index_array);
    tbb::parallel_scan(tbb::blocked_range<size_t>(
            0UL, UI64(_num_hypernodes + 1)), incident_net_prefix_sum);
    const size_t size_in_bytes = incident_net_prefix_sum.total_sum();
    _incident_net_array = parallel::make_unique<char>(size_in_bytes);

    // Insert incident nets into incidence array
    tbb::parallel_for(ID(0), num_hyperedges, [&](const HyperedgeID he) {
      for ( const HypernodeID& pin : edge_vector[he] ) {
        Entry* entry = firstEntry(pin) + current_incident_net_pos[pin]++;
        entry->e = he;
        entry->version = 0;
      }
    });

    // Setup Header of each vertex
    tbb::parallel_for(ID(0), _num_hypernodes, [&](const HypernodeID u) {
      Header* head = header(u);
      head->prev = u;
      head->next = u;
      head->it_prev = u;
      head->it_next = u;
      head->size = current_incident_net_pos[u].load(std::memory_order_relaxed);
      head->current_version = 0;
    });
  }

  const HypernodeID _num_hypernodes;
  parallel::scalable_vector<size_t> _index_array;
  parallel::tbb_unique_ptr<char> _incident_net_array;
};
}  // namespace ds
}  // namespace mt_kahypar
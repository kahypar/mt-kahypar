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

#include "mt-kahypar/datastructures/dynamic_adjacency_array.h"

#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/datastructures/streaming_vector.h"

namespace mt_kahypar {
namespace ds {

IncidentEdgeIterator::IncidentEdgeIterator(const HypernodeID u,
                                           const DynamicAdjacencyArray* dynamic_adjacency_array,
                                           const size_t pos,
                                           const bool end):
    _u(u),
    _current_u(u),
    _current_size(dynamic_adjacency_array->header(u)->size()),
    _current_pos(pos),
    _dynamic_adjacency_array(dynamic_adjacency_array),
    _end(end) {
  if ( end ) {
    _current_pos = _current_size;
  }

  ASSERT(pos <= dynamic_adjacency_array->nodeDegree(u));
  traverse_headers();
}

HyperedgeID IncidentEdgeIterator::operator* () const {
  return _dynamic_adjacency_array->firstActiveEdge(_current_u) + _current_pos;
}

IncidentEdgeIterator & IncidentEdgeIterator::operator++ () {
  ASSERT(!_end);
  ++_current_pos;
  traverse_headers();
  return *this;
}

bool IncidentEdgeIterator::operator!= (const IncidentEdgeIterator& rhs) {
  return !(*this == rhs);
}

bool IncidentEdgeIterator::operator== (const IncidentEdgeIterator& rhs) {
  return _u == rhs._u && _end == rhs._end;
}

void IncidentEdgeIterator::traverse_headers() {
  while ( _current_pos >= _current_size ) {
    const HypernodeID last_u = _current_u;
    _current_u = _dynamic_adjacency_array->header(last_u)->it_next;
    _current_pos -= _current_size;
    _current_size = _dynamic_adjacency_array->header(_current_u)->size();
    // It can happen that due to a contraction the current vertex
    // we iterate over becomes empty or the head of the current vertex
    // changes. Therefore, we set the end flag if we reach the current
    // head of the list or it_next is equal with the current vertex (means
    // that list becomes empty due to a contraction)
    if ( _dynamic_adjacency_array->header(_current_u)->is_head ||
         last_u == _current_u ) {
      _end = true;
      break;
    }
  }
}

EdgeIterator::EdgeIterator(const HypernodeID u,
                           const DynamicAdjacencyArray* dynamic_adjacency_array,
                           std::function<bool (const HypernodeID)> filter):
    _current_u(u),
    _current_id(dynamic_adjacency_array->firstActiveEdge(u)),
    _current_last_id(dynamic_adjacency_array->firstInactiveEdge(u)),
    _dynamic_adjacency_array(dynamic_adjacency_array),
    _filter(filter) {
  traverse_headers();
}

HyperedgeID EdgeIterator::operator* () const {
  return _current_id;
}

EdgeIterator & EdgeIterator::operator++ () {
  ++_current_id;
  traverse_headers();
  return *this;
}

bool EdgeIterator::operator!= (const EdgeIterator& rhs) {
  return !(*this == rhs);
}

bool EdgeIterator::operator== (const EdgeIterator& rhs) {
  return _current_id == rhs._current_id;
}

void EdgeIterator::traverse_headers() {
  while ((!_filter(_current_u) || _current_id == _current_last_id)
            && _current_u < _dynamic_adjacency_array->_num_nodes) {
    ++_current_u;
    _current_id = _dynamic_adjacency_array->firstActiveEdge(_current_u);
    _current_last_id = _dynamic_adjacency_array->firstInactiveEdge(_current_u);
  }
}

void DynamicAdjacencyArray::construct(const EdgeVector& edge_vector, const HyperedgeWeight* edge_weight) {
  // TODO(maas): edge weights
  // Accumulate degree of each vertex thread local
  const HyperedgeID num_hyperedges = edge_vector.size();
  ThreadLocalCounter local_incident_nets_per_vertex(_num_nodes + 1, 0);
  AtomicCounter current_incident_net_pos;
  tbb::parallel_invoke([&] {
    tbb::parallel_for(ID(0), num_hyperedges, [&](const size_t pos) {
      parallel::scalable_vector<size_t>& num_incident_nets_per_vertex =
        local_incident_nets_per_vertex.local();
        ++num_incident_nets_per_vertex[edge_vector[pos].first + 1];
        ++num_incident_nets_per_vertex[edge_vector[pos].second + 1];
    });
  }, [&] {
    _index_array.assign(_num_nodes + 1, sizeof(Header) / sizeof(Edge));
    _index_array[0] = 0;
    current_incident_net_pos.assign(
      _num_nodes, parallel::IntegralAtomicWrapper<size_t>(0));
  });

  // We sum up the number of incident nets per vertex only thread local.
  // To obtain the global number of incident nets per vertex, we iterate
  // over each thread local counter and sum it up.
  for ( const parallel::scalable_vector<size_t>& c : local_incident_nets_per_vertex ) {
    tbb::parallel_for(ID(0), _num_nodes + 1, [&](const size_t pos) {
      _index_array[pos] += c[pos];
    });
  }

  // Compute start positon of the incident nets of each vertex via a parallel prefix sum
  parallel::TBBPrefixSum<HyperedgeID, Array> incident_net_prefix_sum(_index_array);
  tbb::parallel_scan(tbb::blocked_range<size_t>(
          ID(0), ID(_num_nodes + 1)), incident_net_prefix_sum);
  _size_in_bytes = incident_net_prefix_sum.total_sum() * sizeof(Edge);
  _data = parallel::make_unique<Edge>(_size_in_bytes / sizeof(Edge));

  // Insert incident nets into incidence array
  tbb::parallel_for(ID(0), num_hyperedges, [&](const HyperedgeID he) {
    HypernodeID source = edge_vector[he].first;
    HypernodeID target = edge_vector[he].second;
    const HyperedgeWeight weight = edge_weight == nullptr ? 1 : edge_weight[he];
    Edge& e1 = edge(firstEdge(source) + current_incident_net_pos[source].fetch_add(1));
    e1.source = source;
    e1.target = target;
    e1.weight = weight;
    e1.version = 0;
    e1.original_target = e1.target;
    Edge& e2 = edge(firstEdge(target) + current_incident_net_pos[target].fetch_add(1));
    e2.source = target;
    e2.target = source;
    e2.weight = weight;
    e2.version = 0;
    e2.original_target = e2.target;
  });

  // Setup Header of each vertex
  tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID u) {
    Header* head = header(u);
    head->prev = u;
    head->next = u;
    head->it_prev = u;
    head->it_next = u;
    head->degree = current_incident_net_pos[u].load(std::memory_order_relaxed);
    head->first_active = 0;
    head->first_inactive = head->degree;
    head->current_version = 0;
    head->is_head = true;
  });
}

void DynamicAdjacencyArray::contract(const HypernodeID u,
                                     const HypernodeID v,
                                     const AcquireLockFunc& acquire_lock,
                                     const ReleaseLockFunc& release_lock) {
  // iterate over edges of v and update them
  Header* head_v = header(v);
  for (HypernodeID current_v: headers(v)) {
    Header* head = header(current_v);
    const HypernodeID new_version = ++head->current_version;
    for ( HyperedgeID curr_edge = firstActiveEdge(current_v); curr_edge < firstInactiveEdge(current_v); ) {
      Edge& e = edge(curr_edge);
      if ( e.target == v ) {
        ASSERT(e.isSinglePin());
        swap_to_back(current_v, curr_edge);
        --head_v->degree;
      } else {
        e.version = new_version;
        e.source = u;
        const HyperedgeID backwardsEdge = findBackwardsEdge(e, v);
        edge(backwardsEdge).target = u;
        ++curr_edge;
      }
    }

    if ( head->size() == 0 && current_v != v ) {
      removeEmptyIncidentEdgeList(current_v);
    }
  }

  acquire_lock(u);
  // Concatenate double-linked list of u and v
  append(u, v);
  header(u)->degree += head_v->degree;
  ASSERT(verifyIteratorPointers(u), "Iterator pointers of vertex" << u << "are corrupted");
  release_lock(u);
}

void DynamicAdjacencyArray::uncontract(const HypernodeID u,
                                       const HypernodeID v,
                                       const AcquireLockFunc& acquire_lock,
                                       const ReleaseLockFunc& release_lock) {
  uncontract(u, v, [](HyperedgeID) {}, [](HyperedgeID) {}, acquire_lock, release_lock);
}

void DynamicAdjacencyArray::uncontract(const HypernodeID u,
                                       const HypernodeID v,
                                       const CaseOneFunc& case_one_func,
                                       const CaseTwoFunc& case_two_func,
                                       const AcquireLockFunc& acquire_lock,
                                       const ReleaseLockFunc& release_lock) {
  ASSERT(header(v)->prev != v);
  Header* head_u = header(u);
  Header* head_v = header(v);
  acquire_lock(u);
  // Restores the incident list of v to the time before it was appended
  // to the double-linked list of u.
  splice(u, v);
  ASSERT(verifyIteratorPointers(u), "Iterator pointers of vertex" << u << "are corrupted");
  ASSERT(head_u->degree >= head_v->degree, V(head_u->degree) << V(head_v->degree));
  head_u->degree -= head_v->degree;
  release_lock(u);

  // iterate over edges of v, update backwards edges and restore removed edges
  HypernodeID last_non_empty_v = v;
  for (HypernodeID current_v: headers(v)) {
    Header* head = header(current_v);
    ASSERT(head->current_version > 0);
    const HypernodeID new_version = --head->current_version;
    const HyperedgeID first_inactive = firstInactiveEdge(current_v);
    for (HyperedgeID curr_edge = firstActiveEdge(current_v); curr_edge < first_inactive; ++curr_edge) {
      Edge& e = edge(curr_edge);
      e.version = new_version;
      e.source = v;
      const HyperedgeID backwardsEdge = findBackwardsEdge(e, u);
      edge(backwardsEdge).target = v;
      if (e.target == u) {
        case_one_func(curr_edge);
      } else {
        case_two_func(curr_edge);
      }
    }

    const HyperedgeID last_edge = lastEdge(current_v);
    for (HyperedgeID curr_edge = first_inactive; curr_edge < last_edge; ++curr_edge) {
      Edge& e = edge(curr_edge);
      if (e.version != new_version) {
        break;
      }
      ++head->first_inactive;
      ++head_v->degree;
      ASSERT(e.isSinglePin() && e.target == v);
    }

    if (head->size() > 0) {
      restoreItLink(v, last_non_empty_v, current_v);
      last_non_empty_v = current_v;
    }
  }
  ASSERT(verifyIteratorPointers(v), "Iterator pointers of vertex" << v << "are corrupted");
}

parallel::scalable_vector<DynamicAdjacencyArray::RemovedEdgesOrWeight> DynamicAdjacencyArray::removeParallelEdges() {
  StreamingVector<RemovedEdgesOrWeight> tmp_removed_edges;
  // TODO(maas): special case for high degree nodes?
  tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID u) {
    if (header(u)->is_head) {
      vec<ParallelEdgeInformation>& local_vec = _thread_local_vec.local();
      local_vec.clear();

      // we sort all incident edges
      for (HypernodeID current_u: headers(u)) {
        const HyperedgeID first_inactive = firstInactiveEdge(current_u);
        for (HyperedgeID e = firstActiveEdge(current_u); e < first_inactive; ++e) {
          local_vec.push_back(ParallelEdgeInformation(edge(e).target, e, current_u));
        }
      }
      std::sort(local_vec.begin(), local_vec.end(), [](const auto& e1, const auto& e2) {
        return e1.target < e2.target;
      });

      // scan and swap all duplicates to front
      HyperedgeID num_duplicates = 0;
      HyperedgeWeight current_weight = 0;
      for (size_t i = 0; i + 1 < local_vec.size(); ++i) {
        const ParallelEdgeInformation& e1 = local_vec[i];
        const ParallelEdgeInformation& e2 = local_vec[i + 1];
        if (e1.target == e2.target) {
          // local_vec[i + 1] might also be parall to e1 and e2,
          // thus we remove e1 and keep e2
          ASSERT(e1.edge_id >= firstActiveEdge(e1.header_id)
                && e1.edge_id < firstInactiveEdge(e1.header_id),
                V(firstActiveEdge(e1.header_id)) << V(e1.edge_id) << V(firstInactiveEdge(e1.header_id)));
          ASSERT(e2.edge_id >= firstActiveEdge(e2.header_id)
                && e2.edge_id < firstInactiveEdge(e2.header_id),
                V(firstActiveEdge(e2.header_id)) << V(e2.edge_id) << V(firstInactiveEdge(e2.header_id)));
          current_weight += edge(e1.edge_id).weight;
          std::swap(local_vec[num_duplicates], local_vec[i]);
          ++num_duplicates;
        } else if (current_weight > 0) {
          tmp_removed_edges.stream(RemovedEdgesOrWeight { e1.header_id, edge(e1.edge_id).weight, e1.target });
          edge(e1.edge_id).weight += current_weight;
          current_weight = 0;
        }
      }
      if (current_weight > 0) {
        const ParallelEdgeInformation& last = local_vec[local_vec.size() - 1];
        tmp_removed_edges.stream(RemovedEdgesOrWeight { last.header_id, edge(last.edge_id).weight, last.target });
        edge(last.edge_id).weight += current_weight;
      }

      if (num_duplicates > 0) {
        // sort again based on edge id, so the edges can be processed in sorted order
        local_vec.resize(num_duplicates);
        std::sort(local_vec.begin(), local_vec.end(), [](const auto& e1, const auto& e2) {
          return e1.edge_id < e2.edge_id;
        });

        // swap edges and collect output
        HypernodeID current_u = local_vec[0].header_id;
        HyperedgeID current_count = 0;
        for (const ParallelEdgeInformation& e: local_vec) {
          if (current_u != e.header_id) {
            tmp_removed_edges.stream(RemovedEdgesOrWeight { current_u, current_count });
            current_u = e.header_id;
            current_count = 0;
          }
          swap_to_front(e.header_id, e.edge_id);
          ++current_count;
        }
        tmp_removed_edges.stream(RemovedEdgesOrWeight { current_u, current_count });
        header(u)->degree -= num_duplicates;
      }
    }
  });
  auto removed_edges = tmp_removed_edges.copy_parallel();
  tmp_removed_edges.clear_parallel();
  return removed_edges;
}

void DynamicAdjacencyArray::restoreParallelEdges(
      const parallel::scalable_vector<DynamicAdjacencyArray::RemovedEdgesOrWeight>& edges_to_restore) {
  // _degree_diffs does not need to be thread safe because all headers are handled separately
  _degree_diffs.assign(_num_nodes, 0);
  tbb::parallel_for(0UL, edges_to_restore.size(), [&](const size_t i) {
    const RemovedEdgesOrWeight& removed = edges_to_restore[i];
    if (removed.is_weight) {
      // restore edge weight
      const HyperedgeID first_inactive = firstInactiveEdge(removed.header);
      const HyperedgeID first = firstEdge(removed.header);
      ASSERT(first > 0); // holds because of the structure of the adjacency array
      bool found = false;
      for (HyperedgeID e = first_inactive - 1; e >= first; --e) {
        if (edge(e).target == removed.target()) {
          edge(e).weight = removed.weight();
          found = true;
          break;
        }
      }
      ASSERT(found);
    } else {
      // restore parallel edges, which have been swapped to the front
      _degree_diffs[removed.header] = removed.num_removed();
      header(removed.header)->first_active -= removed.num_removed();
      ASSERT(header(removed.header)->first_active <= header(removed.header)->first_inactive);
    }
  });

  // update node degrees
  tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID u) {
    if (header(u)->is_head) {
      for (HypernodeID current_u: headers(u)) {
        header(u)->degree += _degree_diffs[current_u];
      }
    }
  });
}

void DynamicAdjacencyArray::sortIncidentEdges() {
  tbb::parallel_for(ID(0), ID(_index_array.size()), [&](HypernodeID u) {
    std::sort(_data.get() + firstActiveEdge(u), _data.get() + firstInactiveEdge(u),
      [&](const auto& e1, const auto& e2) {
        return e1.target < e2.target;
      }
    );
  });
}

DynamicAdjacencyArray DynamicAdjacencyArray::copy(parallel_tag_t) {
  DynamicAdjacencyArray adjacency_array;
  adjacency_array._num_nodes = _num_nodes;
  adjacency_array._size_in_bytes = _size_in_bytes;

  tbb::parallel_invoke([&] {
    adjacency_array._index_array.resize(_index_array.size());
    memcpy(adjacency_array._index_array.data(), _index_array.data(),
      sizeof(HyperedgeID) * _index_array.size());
  }, [&] {
    adjacency_array._data = parallel::make_unique<Edge>(_size_in_bytes / sizeof(Edge));
    memcpy(adjacency_array._data.get(), _data.get(), _size_in_bytes);
  });

  return adjacency_array;
}

DynamicAdjacencyArray DynamicAdjacencyArray::copy() {
  DynamicAdjacencyArray adjacency_array;
  adjacency_array._num_nodes = _num_nodes;
  adjacency_array._size_in_bytes = _size_in_bytes;

  adjacency_array._index_array.resize(_index_array.size());
  memcpy(adjacency_array._index_array.data(), _index_array.data(),
    sizeof(HyperedgeID) * _index_array.size());
  adjacency_array._data = parallel::make_unique<Edge>(_size_in_bytes / sizeof(Edge));
  memcpy(adjacency_array._data.get(), _data.get(), _size_in_bytes);
  return adjacency_array;
}

HyperedgeID DynamicAdjacencyArray::findBackwardsEdge(const Edge& forward, HypernodeID source) const {
  const HypernodeID current_u = forward.original_target;
  const HyperedgeID first_inactive = firstInactiveEdge(current_u);
  for (HyperedgeID e = firstActiveEdge(current_u); e < first_inactive; ++e) {
    if (edge(e).target == source) {
      return e;
    }
  }
  ASSERT(false, "Hypernode" << current_u << "has no outgoing edge with target" << source);
  return kInvalidHyperedge;
}

void DynamicAdjacencyArray::append(const HypernodeID u, const HypernodeID v) {
  const HypernodeID tail_u = header(u)->prev;
  const HypernodeID tail_v = header(v)->prev;
  header(tail_u)->next = v;
  header(v)->prev = tail_u;
  header(tail_v)->next = u;
  header(u)->prev = tail_v;

  const HypernodeID it_tail_u = header(u)->it_prev;
  const HypernodeID it_tail_v = header(v)->it_prev;
  header(it_tail_u)->it_next = v;
  header(v)->it_prev = it_tail_u;
  header(it_tail_v)->it_next = u;
  header(u)->it_prev = it_tail_v;

  header(v)->tail = tail_v;
  header(v)->is_head = false;

  if ( header(v)->size() == 0 ) {
    removeEmptyIncidentEdgeList(v);
  }
}

void DynamicAdjacencyArray::splice(const HypernodeID u, const HypernodeID v) {
  // Restore the iterator double-linked list of u such that it does not contain
  // any incident net list of v
  const HypernodeID tail = header(v)->tail;
  HypernodeID non_empty_entry_prev_v = v;
  HypernodeID non_empty_entry_next_tail = tail;
  while ( ( non_empty_entry_prev_v == v ||
          header(non_empty_entry_prev_v)->size() == 0 ) &&
          non_empty_entry_prev_v != u ) {
    non_empty_entry_prev_v = header(non_empty_entry_prev_v)->prev;
  }
  while ( ( non_empty_entry_next_tail == tail ||
          header(non_empty_entry_next_tail)->size() == 0 ) &&
          non_empty_entry_next_tail != u ) {
    non_empty_entry_next_tail = header(non_empty_entry_next_tail)->next;
  }
  header(non_empty_entry_prev_v)->it_next = non_empty_entry_next_tail;
  header(non_empty_entry_next_tail)->it_prev = non_empty_entry_prev_v;

  // Cut out incident list of v
  const HypernodeID prev_v = header(v)->prev;
  const HypernodeID next_tail = header(tail)->next;
  header(v)->prev = tail;
  header(tail)->next = v;
  header(next_tail)->prev = prev_v;
  header(prev_v)->next = next_tail;
  header(v)->is_head = true;
}

void DynamicAdjacencyArray::removeEmptyIncidentEdgeList(const HypernodeID u) {
  ASSERT(!header(u)->is_head);
  ASSERT(header(u)->size() == 0, V(u) << V(header(u)->size()));
  Header* head = header(u);
  header(head->it_prev)->it_next = head->it_next;
  header(head->it_next)->it_prev = head->it_prev;
  head->it_next = u;
  head->it_prev = u;
}

void DynamicAdjacencyArray::restoreIteratorPointers(const HypernodeID u) {
  HypernodeID last_non_empty_u = u;
  for (HypernodeID current_u: headers(u)) {
    if (header(current_u)->size() > 0) {
      restoreItLink(u, last_non_empty_u, current_u);
      last_non_empty_u = current_u;
    }
  }
}

void DynamicAdjacencyArray::restoreItLink(const HypernodeID u, const HypernodeID prev, const HypernodeID current) {
  header(prev)->it_next = current;
  header(current)->it_prev = prev;
  header(current)->it_next = u;
  header(u)->it_prev = current;
}

bool DynamicAdjacencyArray::verifyIteratorPointers(const HypernodeID u) const {
  HypernodeID current_u = u;
  HypernodeID last_non_empty_entry = kInvalidHypernode;
  do {
    if ( header(current_u)->size() > 0 || current_u == u ) {
      if ( last_non_empty_entry != kInvalidHypernode ) {
        if ( header(current_u)->it_prev != last_non_empty_entry ) {
          return false;
        } else if ( header(last_non_empty_entry)->it_next != current_u ) {
          return false;
        }
      }
      last_non_empty_entry = current_u;
    } else {
      if ( header(current_u)->it_next != current_u ) {
        return false;
      } else if ( header(current_u)->it_prev != current_u ) {
        return false;
      }
    }

    current_u = header(current_u)->next;
  } while(current_u != u);

  if ( header(u)->it_prev != last_non_empty_entry ) {
    return false;
  } else if ( header(last_non_empty_entry)->it_next != u ) {
    return false;
  }

  return true;
}

}  // namespace ds
}  // namespace mt_kahypar

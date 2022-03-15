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

namespace mt_kahypar {
namespace ds {

IncidentEdgeIterator::IncidentEdgeIterator(const HypernodeID u,
                                           const DynamicAdjacencyArray* dynamic_adjacency_array,
                                           const size_t pos,
                                           const bool end):
    _u(u),
    _current_u(u),
    _current_size(dynamic_adjacency_array->header(u).size()),
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
  skip_invalid();
  while ( _current_pos >= _current_size ) {
    const HypernodeID last_u = _current_u;
    _current_u = _dynamic_adjacency_array->header(last_u).it_next;
    _current_pos -= _current_size;
    _current_size = _dynamic_adjacency_array->header(_current_u).size();
    // It can happen that due to a contraction the current vertex
    // we iterate over becomes empty or the head of the current vertex
    // changes. Therefore, we set the end flag if we reach the current
    // head of the list or it_next is equal with the current vertex (means
    // that list becomes empty due to a contraction)
    if ( _dynamic_adjacency_array->header(_current_u).is_head ||
         last_u == _current_u ) {
      _end = true;
      break;
    }
    skip_invalid();
  }
}

void IncidentEdgeIterator::skip_invalid() {
  while (_current_pos < _current_size &&
         !_dynamic_adjacency_array->edge(**this).isValid()) {
    ++_current_pos;
  }
}

EdgeIterator::EdgeIterator(const HypernodeID u,
                           const DynamicAdjacencyArray* dynamic_adjacency_array):
    _current_u(u),
    _current_id(dynamic_adjacency_array->firstActiveEdge(u)),
    _current_last_id(dynamic_adjacency_array->firstInactiveEdge(u)),
    _dynamic_adjacency_array(dynamic_adjacency_array) {
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
  skip_invalid();
  while (_current_id == _current_last_id && _current_u < _dynamic_adjacency_array->_num_nodes) {
    ++_current_u;
    _current_id = _dynamic_adjacency_array->firstActiveEdge(_current_u);
    _current_last_id = _dynamic_adjacency_array->firstInactiveEdge(_current_u);
    skip_invalid();
  }
}

void EdgeIterator::skip_invalid() {
  while (_current_id < _current_last_id &&
         !_dynamic_adjacency_array->edge(**this).isValid()) {
    ++_current_id;
  }
}

void DynamicAdjacencyArray::construct(const EdgeVector& edge_vector, const HyperedgeWeight* edge_weight) {
  // Accumulate degree of each vertex thread local
  const HyperedgeID num_edges = edge_vector.size();
  ThreadLocalCounter local_incident_nets_per_vertex(_num_nodes + 1, 0);
  Array<HyperedgeID> node_degrees;
  AtomicCounter current_incident_net_pos;
  tbb::parallel_invoke([&] {
    tbb::parallel_for(ID(0), num_edges, [&](const size_t pos) {
      parallel::scalable_vector<size_t>& num_incident_nets_per_vertex =
        local_incident_nets_per_vertex.local();
        ++num_incident_nets_per_vertex[edge_vector[pos].first];
        ++num_incident_nets_per_vertex[edge_vector[pos].second];
    });
  }, [&] {
    _header_array.resize(_num_nodes + 1);
  }, [&] {
    _edges.resize(2 * num_edges);
  }, [&] {
    _edge_mapping.resize(2 * num_edges);
  }, [&] {
    _degree_diffs.resize(_num_nodes);
  }, [&] {
    node_degrees.resize(_num_nodes);
  }, [&] {
    current_incident_net_pos.assign(
      _num_nodes, parallel::IntegralAtomicWrapper<size_t>(0));
  });

  // We sum up the number of incident nets per vertex only thread local.
  // To obtain the global number of incident nets per vertex, we iterate
  // over each thread local counter and sum it up.
  for ( const parallel::scalable_vector<size_t>& c : local_incident_nets_per_vertex ) {
    tbb::parallel_for(ID(0), _num_nodes, [&](const size_t pos) {
      node_degrees[pos] += c[pos];
    });
  }

  // Compute start positon of the incident nets of each vertex via a parallel prefix sum
  parallel::TBBPrefixSum<HyperedgeID, Array> incident_net_prefix_sum(node_degrees);
  tbb::parallel_scan(tbb::blocked_range<size_t>(
          ID(0), ID(_num_nodes)), incident_net_prefix_sum);

  // Setup Header of each vertex
  tbb::parallel_for(ID(0), _num_nodes + 1, [&](const HypernodeID u) {
    Header& head = header(u);
    head.prev = u;
    head.next = u;
    head.it_prev = u;
    head.it_next = u;
    head.degree = (u == _num_nodes) ? 0 : incident_net_prefix_sum[u + 1] - incident_net_prefix_sum[u];
    head.first = incident_net_prefix_sum[u];
    head.first_active = head.first;
    head.first_inactive = head.first + head.degree;
    head.is_head = true;
  });

  // Insert incident nets into incidence array
  tbb::parallel_for(ID(0), num_edges, [&](const HyperedgeID he) {
    HypernodeID source = edge_vector[he].first;
    HypernodeID target = edge_vector[he].second;
    const HyperedgeWeight weight = edge_weight == nullptr ? 1 : edge_weight[he];
    HyperedgeID id1 = firstEdge(source) + current_incident_net_pos[source].fetch_add(1);
    HyperedgeID id2 = firstEdge(target) + current_incident_net_pos[target].fetch_add(1);
    Edge& e1 = edge(id1);
    e1.source = source;
    e1.target = target;
    e1.weight = weight;
    e1.back_edge = id2;
    e1.original_source = e1.source;
    e1.unique_id = he;
    Edge& e2 = edge(id2);
    e2.source = target;
    e2.target = source;
    e2.weight = weight;
    e2.back_edge = id1;
    e2.original_source = e2.source;
    e2.unique_id = he;
  });

  // TODO(maas): is this the appropriate way to check this?
  ASSERT([&]() {
    auto removed_edges = removeSinglePinAndParallelEdges();
    return removed_edges.empty();
  }(), "Graph contains single pin edge or parallel edges!");
}

void DynamicAdjacencyArray::contract(const HypernodeID u,
                                     const HypernodeID v,
                                     const AcquireLockFunc& acquire_lock,
                                     const ReleaseLockFunc& release_lock) {
  // iterate over edges of v and update them
  Header& head_v = header(v);
  for (HypernodeID current_v: headers(v)) {
    const HyperedgeID last = firstInactiveEdge(current_v);
    for ( HyperedgeID curr_edge = firstActiveEdge(current_v); curr_edge < last; ++curr_edge ) {
      Edge& e = edge(curr_edge);
      if (e.isValid() && e.isSinglePin()) {
        ASSERT(e.source == v);
        e.setValid(false);
        --head_v.degree;
      } else if (e.isValid()) {
        e.source = u;
        edge(e.back_edge).target = u;
      }
    }
  }

  acquire_lock(u);
  // Concatenate double-linked list of u and v
  append(u, v);
  header(u).degree += head_v.degree;
  ASSERT(verifyIteratorPointers(u), "Iterator pointers of vertex" << u << "are corrupted");
  release_lock(u);
}

void DynamicAdjacencyArray::uncontract(const HypernodeID u,
                                       const HypernodeID v,
                                       const AcquireLockFunc& acquire_lock,
                                       const ReleaseLockFunc& release_lock) {
  uncontract(u, v, [](HyperedgeID) { return false; }, [](HyperedgeID) {}, [](HyperedgeID) {},
             acquire_lock, release_lock);
}

void DynamicAdjacencyArray::uncontract(const HypernodeID u,
                                       const HypernodeID v,
                                       const MarkEdgeFunc& mark_edge,
                                       const CaseOneFunc& case_one_func,
                                       const CaseTwoFunc& case_two_func,
                                       const AcquireLockFunc& acquire_lock,
                                       const ReleaseLockFunc& release_lock) {
  ASSERT(header(v).prev != v);
  Header& head_u = header(u);
  Header& head_v = header(v);
  acquire_lock(u);
  // Restores the incident list of v to the time before it was appended
  // to the double-linked list of u.
  splice(u, v);
  ASSERT(verifyIteratorPointers(u), "Iterator pointers of vertex" << u << "are corrupted");
  ASSERT(head_u.degree >= head_v.degree, V(head_u.degree) << V(head_v.degree));
  head_u.degree -= head_v.degree;
  release_lock(u);

  // iterate over edges of v, update backwards edges and restore removed edges
  HypernodeID last_non_empty_v = v;
  for (HypernodeID current_v: headers(v)) {
    const HyperedgeID first_inactive = firstInactiveEdge(current_v);
    for (HyperedgeID curr_edge = firstActiveEdge(current_v); curr_edge < first_inactive; ++curr_edge) {
      Edge& e = edge(curr_edge);
      ASSERT(e.source == u || !e.isValid());
      if (e.source == u) {
        // TODO(maas): reading e.target might be a bit of a race condition
        bool singlePin = false;
        if (e.target == u) {
          // the edge is not truly single pin if already marked
          singlePin = !mark_edge(curr_edge);
        }
        e.source = v;
        edge(e.back_edge).target = v;
        if (singlePin) {
          case_one_func(curr_edge);
        } else {
          case_two_func(curr_edge);
        }
      } else if (e.source == v) {
        e.setValid(true);
        ++head_v.degree;
      }
    }


    if (header(current_v).size() > 0) {
      restoreItLink(v, last_non_empty_v, current_v);
      last_non_empty_v = current_v;
    }
  }
  ASSERT(verifyIteratorPointers(v), "Iterator pointers of vertex" << v << "are corrupted");
}

void DynamicAdjacencyArray::streamWeight(StreamingVector<RemovedEdgesOrWeight>& tmp_removed_edges,
                                         const ParallelEdgeInformation& e, HyperedgeWeight w) {
  tmp_removed_edges.stream(RemovedEdgesOrWeight::asWeight(e.header_id, e.unique_id, edge(e.edge_id).weight));
  edge(e.edge_id).weight += w;
}

parallel::scalable_vector<DynamicAdjacencyArray::RemovedEdgesOrWeight> DynamicAdjacencyArray::removeSinglePinAndParallelEdges() {
  StreamingVector<RemovedEdgesOrWeight> tmp_removed_edges;
  initializeEdgeMapping(_edge_mapping);
  // TODO(maas): special case for high degree nodes?
  tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID u) {
    if (header(u).is_head) {
      vec<ParallelEdgeInformation>& local_vec = _thread_local_vec.local();
      local_vec.clear();

      // we sort all incident edges
      for (HypernodeID current_u: headers(u)) {
        const HyperedgeID first_inactive = firstInactiveEdge(current_u);
        for (HyperedgeID e = firstActiveEdge(current_u); e < first_inactive; ++e) {
          local_vec.emplace_back(edge(e).target, e, edge(e).unique_id, current_u);
        }
      }
      std::sort(local_vec.begin(), local_vec.end(), [](const auto& e1, const auto& e2) {
        // we need a symmetric order on edges and backedges to ensure that the
        // kept forward and backward edge are actually the same edge
        return e1.target < e2.target || (e1.target == e2.target && e1.unique_id < e2.unique_id);
      });

      // scan and swap all duplicate/single pin/invalid edges to front
      HyperedgeID num_duplicates = 0;
      HyperedgeWeight current_weight = 0;
      for (size_t i = 0; i + 1 < local_vec.size(); ++i) {
        const ParallelEdgeInformation& e1 = local_vec[i];
        const ParallelEdgeInformation& e2 = local_vec[i + 1];
        if (e1.target != e2.target && current_weight > 0) {
          streamWeight(tmp_removed_edges, e1, current_weight);
          current_weight = 0;
        }
        if (e1.target == u || e1.target == kInvalidHypernode) {
          std::swap(local_vec[num_duplicates], local_vec[i]);
          ++num_duplicates;
        } else if (e1.target == e2.target) {
          // local_vec[i + 1] might also be parall to e1 and e2,
          // thus we remove e1 and keep e2
          current_weight += edge(e1.edge_id).weight;
          std::swap(local_vec[num_duplicates], local_vec[i]);
          ++num_duplicates;
        }
      }
      if (!local_vec.empty()) {
        const ParallelEdgeInformation& last = local_vec[local_vec.size() - 1];
        if (current_weight > 0) {
          streamWeight(tmp_removed_edges, last, current_weight);
        }
        if (last.target == u || last.target == kInvalidHypernode) {
          std::swap(local_vec[num_duplicates], local_vec[local_vec.size() - 1]);
          ++num_duplicates;
        }
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
        HyperedgeID current_degree_diff = 0;
        HyperedgeID total_degree_diff = 0;
        for (const ParallelEdgeInformation& e: local_vec) {
          Header& head = header(e.header_id);
          if (current_u != e.header_id) {
            tmp_removed_edges.stream(RemovedEdgesOrWeight::asEdges(current_u, current_count, current_degree_diff));
            current_u = e.header_id;
            current_count = 0;
            current_degree_diff = 0;
          }
          ++current_count;
          if (e.target != kInvalidHypernode) {
            ++current_degree_diff;
            ++total_degree_diff;
          }

          const HyperedgeID new_id = firstActiveEdge(e.header_id);
          HyperedgeID permutation_source = new_id;
          while (_edge_mapping[permutation_source] < new_id) {
            permutation_source = _edge_mapping[permutation_source];
          }
          ASSERT(_edge_mapping[e.edge_id] == e.edge_id && _edge_mapping[permutation_source] == new_id);
          std::swap(edge(e.edge_id), edge(new_id));
          ++head.first_active;
          _edge_mapping[e.edge_id] = new_id;
          _edge_mapping[permutation_source] = e.edge_id;

          if (head.size() == 0 && current_u != u) {
            removeEmptyIncidentEdgeList(e.header_id);
          }
        }
        tmp_removed_edges.stream(RemovedEdgesOrWeight::asEdges(current_u, current_count, current_degree_diff));
        header(u).degree -= total_degree_diff;
      }
    }
  });

  applyEdgeMapping(_edge_mapping);

  HEAVY_COARSENING_ASSERT([&]() {
    for (HyperedgeID e = 0; e < _edges.size(); ++e) {
      if (edge(edge(e).back_edge).back_edge != e) {
        return false;
      }
    }
    return true;
  }());

  auto removed_edges = tmp_removed_edges.copy_parallel();
  tmp_removed_edges.clear_parallel();
  return removed_edges;
}

void DynamicAdjacencyArray::restoreSinglePinAndParallelEdges(
      const parallel::scalable_vector<DynamicAdjacencyArray::RemovedEdgesOrWeight>& edges_to_restore) {
  // _degree_diffs does not need to be thread safe because all headers are handled separately
  _degree_diffs.assign(_num_nodes, 0);
  tbb::parallel_for(0UL, edges_to_restore.size(), [&](const size_t i) {
    const RemovedEdgesOrWeight& removed = edges_to_restore[i];
    Header& head = header(removed.header);
    if (removed.is_weight) {
      // restore edge weight
      const HyperedgeID first_inactive = firstInactiveEdge(removed.header);
      bool found = false;
      for (HyperedgeID e = firstEdge(removed.header); e < first_inactive; ++e) {
        if (edge(e).unique_id == removed.edgeID()) {
          edge(e).weight = removed.weight();
          found = true;
          break;
        }
      }
      ASSERT(found);
      unused(found);
    } else {
      // restore parallel edges, which have been swapped to the front
      _degree_diffs[removed.header] = removed.degreeDiff();
      if (removed.numRemoved() > 0 && head.size() == 0 && !head.is_head) {
        // we use a negative sign to mark previously empty headers
        _degree_diffs[removed.header] = -_degree_diffs[removed.header] - 1;
      }
      head.first_active -= removed.numRemoved();
    }
  });

  // update node degrees
  tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID u) {
    if (header(u).is_head) {
      bool restore_it = false;
      for (HypernodeID current_u: headers(u)) {
        header(u).degree += _degree_diffs[current_u] >= 0 ?
                            _degree_diffs[current_u] : -_degree_diffs[current_u] - 1;
        restore_it |= _degree_diffs[current_u] < 0;
      }
      if (restore_it) {
        restoreIteratorPointers(u);
      }
    }
  });
}

void DynamicAdjacencyArray::reset() {
  // Nothing to do here
}

void DynamicAdjacencyArray::sortIncidentEdges() {
  // this is a bit complicated because we need to update the back edges
  Array<HyperedgeID> edge_permutation;
  edge_permutation.resize(_edges.size());
  initializeEdgeMapping(edge_permutation);

  tbb::parallel_for(ID(0), ID(_header_array.size()), [&](HypernodeID u) {
    // sort mapped indizes
    const HyperedgeID start = firstActiveEdge(u);
    const HyperedgeID end = firstInactiveEdge(u);
    std::sort(edge_permutation.data() + start, edge_permutation.data() + end,
      [&](const auto& e1, const auto& e2) {
        return edge(e1).target < edge(e2).target;
      }
    );

    // apply permutation
    for (size_t i = start; i < end; ++i) {
      HyperedgeID target = edge_permutation[i];
      while (target < i) {
        target = edge_permutation[target];
      }
      std::swap(_edges[i], _edges[target]);
    }
  });

  // we need the reversed permutation for the back edges
  tbb::parallel_for(ID(0), ID(edge_permutation.size()), [&](const HyperedgeID e) {
    _edge_mapping[edge_permutation[e]] = e;
  });
  applyEdgeMapping(_edge_mapping);

  HEAVY_PREPROCESSING_ASSERT([&]() {
    for (HyperedgeID e = 0; e < _edges.size(); ++e) {
      if (edge(edge(e).back_edge).back_edge != e) {
        return false;
      }
    }
    return true;
  }());
}

DynamicAdjacencyArray DynamicAdjacencyArray::copy(parallel_tag_t) const {
  DynamicAdjacencyArray adjacency_array;
  adjacency_array._num_nodes = _num_nodes;

  tbb::parallel_invoke([&] {
    adjacency_array._header_array.resize(_header_array.size());
    memcpy(adjacency_array._header_array.data(), _header_array.data(),
      sizeof(Header) * _header_array.size());
  }, [&] {
    adjacency_array._edges.resize(_edges.size());
    memcpy(adjacency_array._edges.data(), _edges.data(), sizeof(Edge) * _edges.size());
  }, [&] {
    adjacency_array._edge_mapping.resize(_edge_mapping.size());
  }, [&] {
    adjacency_array._degree_diffs.resize(_degree_diffs.size());
  });

  return adjacency_array;
}

DynamicAdjacencyArray DynamicAdjacencyArray::copy() const {
  DynamicAdjacencyArray adjacency_array;
  adjacency_array._num_nodes = _num_nodes;

  adjacency_array._header_array.resize(_header_array.size());
  memcpy(adjacency_array._header_array.data(), _header_array.data(),
    sizeof(Header) * _header_array.size());
  adjacency_array._edges.resize(_edges.size());
  memcpy(adjacency_array._edges.data(), _edges.data(), sizeof(Edge) * _edges.size());
  adjacency_array._edge_mapping.resize(_edge_mapping.size());
  adjacency_array._degree_diffs.resize(_degree_diffs.size());
  return adjacency_array;
}

void DynamicAdjacencyArray::append(const HypernodeID u, const HypernodeID v) {
  const HypernodeID tail_u = header(u).prev;
  const HypernodeID tail_v = header(v).prev;
  header(tail_u).next = v;
  header(v).prev = tail_u;
  header(tail_v).next = u;
  header(u).prev = tail_v;

  const HypernodeID it_tail_u = header(u).it_prev;
  const HypernodeID it_tail_v = header(v).it_prev;
  header(it_tail_u).it_next = v;
  header(v).it_prev = it_tail_u;
  header(it_tail_v).it_next = u;
  header(u).it_prev = it_tail_v;

  header(v).tail = tail_v;
  header(v).is_head = false;

  if ( header(v).size() == 0 ) {
    removeEmptyIncidentEdgeList(v);
  }
}

void DynamicAdjacencyArray::splice(const HypernodeID u, const HypernodeID v) {
  // Restore the iterator double-linked list of u such that it does not contain
  // any incident net list of v
  const HypernodeID tail = header(v).tail;
  HypernodeID non_empty_entry_prev_v = v;
  HypernodeID non_empty_entry_next_tail = tail;
  while ( ( non_empty_entry_prev_v == v ||
          header(non_empty_entry_prev_v).size() == 0 ) &&
          non_empty_entry_prev_v != u ) {
    non_empty_entry_prev_v = header(non_empty_entry_prev_v).prev;
  }
  while ( ( non_empty_entry_next_tail == tail ||
          header(non_empty_entry_next_tail).size() == 0 ) &&
          non_empty_entry_next_tail != u ) {
    non_empty_entry_next_tail = header(non_empty_entry_next_tail).next;
  }
  header(non_empty_entry_prev_v).it_next = non_empty_entry_next_tail;
  header(non_empty_entry_next_tail).it_prev = non_empty_entry_prev_v;

  // Cut out incident list of v
  const HypernodeID prev_v = header(v).prev;
  const HypernodeID next_tail = header(tail).next;
  header(v).prev = tail;
  header(tail).next = v;
  header(next_tail).prev = prev_v;
  header(prev_v).next = next_tail;
  header(v).is_head = true;
}

void DynamicAdjacencyArray::removeEmptyIncidentEdgeList(const HypernodeID u) {
  ASSERT(!header(u).is_head);
  ASSERT(header(u).size() == 0, V(u) << V(header(u).size()));
  Header& head = header(u);
  header(head.it_prev).it_next = head.it_next;
  header(head.it_next).it_prev = head.it_prev;
  head.it_next = u;
  head.it_prev = u;
}

void DynamicAdjacencyArray::restoreIteratorPointers(const HypernodeID u) {
  ASSERT(header(u).is_head);
  HypernodeID last_non_empty_u = u;
  for (HypernodeID current_u: headers(u)) {
    if (header(current_u).size() > 0) {
      restoreItLink(u, last_non_empty_u, current_u);
      last_non_empty_u = current_u;
    }
  }
}

void DynamicAdjacencyArray::restoreItLink(const HypernodeID u, const HypernodeID prev, const HypernodeID current) {
  header(prev).it_next = current;
  header(current).it_prev = prev;
  header(current).it_next = u;
  header(u).it_prev = current;
}

bool DynamicAdjacencyArray::verifyIteratorPointers(const HypernodeID u) const {
  HypernodeID current_u = u;
  HypernodeID last_non_empty_entry = kInvalidHypernode;
  do {
    if ( header(current_u).size() > 0 || current_u == u ) {
      if ( last_non_empty_entry != kInvalidHypernode ) {
        if ( header(current_u).it_prev != last_non_empty_entry ) {
          return false;
        } else if ( header(last_non_empty_entry).it_next != current_u ) {
          return false;
        }
      }
      last_non_empty_entry = current_u;
    } else {
      if ( header(current_u).it_next != current_u ) {
        return false;
      } else if ( header(current_u).it_prev != current_u ) {
        return false;
      }
    }

    current_u = header(current_u).next;
  } while(current_u != u);

  if ( header(u).it_prev != last_non_empty_entry ) {
    return false;
  } else if ( header(last_non_empty_entry).it_next != u ) {
    return false;
  }

  return true;
}

}  // namespace ds
}  // namespace mt_kahypar

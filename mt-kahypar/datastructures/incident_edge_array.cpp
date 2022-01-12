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

#include "mt-kahypar/datastructures/incident_edge_array.h"

#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"

namespace mt_kahypar {
namespace ds {

void print_header(const IncidentEdgeArray::Header& header) {
  LOG << V(header.prev) << ", " << V(header.next) << ", " << V(header.it_prev) << ", " << V(header.it_next);
  LOG << V(header.tail) << ", " << V(header.first_active) << ", " << V(header.first_inactive) << V(header.degree);
  LOG << V(header.current_version) << V(header.is_head);
}

IncidentEdgeIterator::IncidentEdgeIterator(const HypernodeID u,
                                           const IncidentEdgeArray* incident_edge_array,
                                           const size_t pos,
                                           const bool end):
    _u(u),
    _current_u(u),
    _current_size(incident_edge_array->header(u)->size()),
    _current_pos(pos),
    _incident_edge_array(incident_edge_array),
    _end(end) {
  if ( end ) {
    _current_pos = _current_size;
  }

  ASSERT(pos <= incident_edge_array->nodeDegree(u));
  traverse_headers();
}

HyperedgeID IncidentEdgeIterator::operator* () const {
  return _incident_edge_array->firstActiveEdge(_current_u) + _current_pos;
}

IncidentEdgeIterator & IncidentEdgeIterator::operator++ () {
  ASSERT(!_end);
  ++_current_pos;
  traverse_headers();
  return *this;
}

bool IncidentEdgeIterator::operator!= (const IncidentEdgeIterator& rhs) {
  return _u != rhs._u || _end != rhs._end;
}

bool IncidentEdgeIterator::operator== (const IncidentEdgeIterator& rhs) {
  return _u == rhs._u && _end == rhs._end;
}

void IncidentEdgeIterator::traverse_headers() {
  while ( _current_pos >= _current_size ) {
    const HypernodeID last_u = _current_u;
    _current_u = _incident_edge_array->header(last_u)->it_next;
    _current_pos -= _current_size;
    _current_size = _incident_edge_array->header(_current_u)->size();
    // It can happen that due to a contraction the current vertex
    // we iterate over becomes empty or the head of the current vertex
    // changes. Therefore, we set the end flag if we reach the current
    // head of the list or it_next is equal with the current vertex (means
    // that list becomes empty due to a contraction)
    if ( _incident_edge_array->header(_current_u)->is_head ||
         last_u == _current_u ) {
      _end = true;
      break;
    }
  }
}

void IncidentEdgeArray::construct(const EdgeVector& edge_vector, const HyperedgeWeight* edge_weight) {
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
          0UL, UI64(_num_nodes + 1)), incident_net_prefix_sum);
  _size_in_bytes = incident_net_prefix_sum.total_sum() * sizeof(Edge);
  _incident_edge_array = parallel::make_unique<Edge>(_size_in_bytes / sizeof(Edge));

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
    Edge& e2 = edge(firstEdge(target) + current_incident_net_pos[target].fetch_add(1));
    e2.source = target;
    e2.target = source;
    e2.weight = weight;
    e2.version = 0;
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

}  // namespace ds
}  // namespace mt_kahypar

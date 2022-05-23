/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2022 Nikolai Maas <nikolai.maas@student.kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "separated_nodes.h"

#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

namespace mt_kahypar::ds {


void SeparatedNodes::addNodes(const vec<std::pair<HyperedgeID, HypernodeWeight>>& nodes,
                              const vec<Edge>& edges) {
  ASSERT(_num_nodes + 1 == _nodes.size());
  ASSERT(_num_edges == _outward_edges.size() && _num_edges == _inward_edges.size());

  auto update_nodes = [&] {
    _nodes.resize(_num_nodes + nodes.size() + 1);
    tbb::parallel_for(0UL, nodes.size(), [&](const size_t& pos) {
      auto [begin, weight] = nodes[pos];
      _nodes[_num_nodes + pos] = Node(_num_edges + begin, weight);
    });
    _nodes[_num_nodes + nodes.size()] = Node(_num_edges + edges.size(), 0);
  };

  auto update_inward_edges = [&] {
    _inward_edges.resize(_num_edges + edges.size());
    tbb::parallel_for(0UL, edges.size(), [&](const size_t& pos) {
      _inward_edges[_num_edges + pos] = edges[pos];
    });
  };

  auto update_outward_edges = [&] {
    Array<parallel::IntegralAtomicWrapper<HyperedgeID>> num_incident_edges;
    Array<parallel::IntegralAtomicWrapper<HyperedgeID>> incident_edges_pos;
    Array<Edge> new_outward_edges;
    tbb::parallel_invoke([&] {
      num_incident_edges.resize(_num_graph_nodes);
    }, [&] {
      incident_edges_pos.resize(_num_graph_nodes);
    }, [&] {
      new_outward_edges.resize(_num_edges + edges.size());
    });

    tbb::parallel_for(0UL, num_incident_edges.size(), [&](const size_t& pos) {
      num_incident_edges[pos] = outwardDegree(pos);
    });

    // calculate new degrees
    tbb::parallel_for(0UL, edges.size(), [&](const size_t& pos) {
      const Edge& e = edges[pos];
      ASSERT(e.target < _num_graph_nodes);
      num_incident_edges[e.target].fetch_add(1);
    });
    parallel::TBBPrefixSum<parallel::IntegralAtomicWrapper<HyperedgeID>, Array>
            incident_edges_prefix_sum(num_incident_edges);
    tbb::parallel_scan(tbb::blocked_range<size_t>(ID(0), num_incident_edges.size()), incident_edges_prefix_sum);

    // write new edges
    tbb::parallel_for(0UL, nodes.size(), [&](const size_t& pos) {
      const size_t end = (pos + 1 == nodes.size()) ? edges.size() : nodes[pos + 1].first;
      for (size_t i = nodes[pos].first; i < end; ++i) {
        const Edge& e = edges[i];
        const auto added = incident_edges_pos[e.target].fetch_add(1);
        const HyperedgeID index = incident_edges_prefix_sum[e.target] + outwardDegree(e.target) +
                                  added;
        new_outward_edges[index] = Edge(_num_nodes + pos, e.weight);
      }
    });

    // copy old edges and calculate incident weight
    tbb::parallel_for(0UL, static_cast<size_t>(_num_graph_nodes), [&](const size_t& pos) {
      for (HyperedgeID i = 0; i < outwardDegree(pos); ++i) {
        const HyperedgeID index = incident_edges_prefix_sum[pos] + i;
        new_outward_edges[index] = _outward_edges[_graph_nodes[pos].begin + i];
      }
      HyperedgeWeight incident_weight = _graph_nodes[pos].incident_weight;
      for (HyperedgeID i = incident_edges_prefix_sum[pos] + outwardDegree(pos);
           i < incident_edges_prefix_sum[pos + 1]; ++i) {
        incident_weight += new_outward_edges[i].weight;
      }
      _graph_nodes[pos].incident_weight = incident_weight;
    });
    tbb::parallel_for(0UL, _graph_nodes.size(), [&](const size_t& pos) {
      _graph_nodes[pos].begin = incident_edges_prefix_sum[pos];
    });

    std::swap(_outward_edges, new_outward_edges);
    parallel::parallel_free(num_incident_edges, incident_edges_pos, new_outward_edges);
  };

  tbb::parallel_invoke(update_nodes, update_inward_edges, update_outward_edges);

  _num_nodes += nodes.size();
  _num_edges += edges.size();
}

// ! Copy in parallel
SeparatedNodes SeparatedNodes::copy(parallel_tag_t) const {
  SeparatedNodes sep_nodes(_num_graph_nodes);

  sep_nodes._num_nodes = _num_nodes;
  sep_nodes._num_graph_nodes = _num_graph_nodes;
  sep_nodes._num_edges = _num_edges;
  sep_nodes._total_weight = _total_weight;

  tbb::parallel_invoke([&] {
    sep_nodes._nodes.resize(_nodes.size());
    memcpy(sep_nodes._nodes.data(), _nodes.data(),
            sizeof(Node) * _nodes.size());
  }, [&] {
    sep_nodes._graph_nodes.resize(_graph_nodes.size());
    memcpy(sep_nodes._graph_nodes.data(), _graph_nodes.data(),
            sizeof(GraphNode) * _graph_nodes.size());
  }, [&] {
    sep_nodes._inward_edges.resize(_inward_edges.size());
    memcpy(sep_nodes._inward_edges.data(), _inward_edges.data(),
            sizeof(Edge) * _inward_edges.size());
  }, [&] {
    sep_nodes._outward_edges.resize(_outward_edges.size());
    memcpy(sep_nodes._outward_edges.data(), _outward_edges.data(),
            sizeof(Edge) * _outward_edges.size());
  });
  return sep_nodes;
}

// ! Copy sequential
SeparatedNodes SeparatedNodes::copy() const {
  SeparatedNodes sep_nodes(_num_graph_nodes);

  sep_nodes._num_nodes = _num_nodes;
  sep_nodes._num_graph_nodes = _num_graph_nodes;
  sep_nodes._num_edges = _num_edges;
  sep_nodes._total_weight = _total_weight;

  sep_nodes._nodes.resize(_nodes.size());
  memcpy(sep_nodes._nodes.data(), _nodes.data(),
          sizeof(Node) * _nodes.size());

  sep_nodes._graph_nodes.resize(_graph_nodes.size());
  memcpy(sep_nodes._graph_nodes.data(), _graph_nodes.data(),
          sizeof(GraphNode) * _graph_nodes.size());

  sep_nodes._inward_edges.resize(_inward_edges.size());
  memcpy(sep_nodes._inward_edges.data(), _inward_edges.data(),
          sizeof(Edge) * _inward_edges.size());

  sep_nodes._outward_edges.resize(_outward_edges.size());
  memcpy(sep_nodes._outward_edges.data(), _outward_edges.data(),
          sizeof(Edge) * _outward_edges.size());

  return sep_nodes;
}

} // namespace

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

namespace mt_kahypar::ds {


void SeparatedNodes::addNodes(const vec<std::pair<HyperedgeID, HypernodeWeight>>& nodes,
                              const vec<Edge>& edges) {
  ASSERT(_num_nodes + 1 == _nodes.size());
  ASSERT(_num_edges == _inward_edges.size());
  ASSERT(_graph_nodes_begin.empty());

  auto update_nodes = [&] {
    _nodes.resize(_num_nodes + nodes.size() + 1);
    tbb::parallel_for(0UL, nodes.size(), [&](const size_t& pos) {
      auto [begin, weight] = nodes[pos];
      _nodes[_num_nodes + pos] = Node(_num_edges + begin, weight);
    });
    _nodes[_num_nodes + nodes.size()] = Node(_num_edges + edges.size(), 0);
  };

  auto update_inward_edges_and_incident_weight = [&] {
    _inward_edges.resize(_num_edges + edges.size());
    tbb::parallel_for(0UL, edges.size(), [&](const size_t& pos) {
      const Edge& e = edges[pos];
      _inward_edges[_num_edges + pos] = e;
      _outward_incident_weight[e.target].fetch_add(e.weight);
    });
  };

  tbb::parallel_invoke(update_nodes, update_inward_edges_and_incident_weight);

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
    sep_nodes._outward_incident_weight.resize(_outward_incident_weight.size());
    for (size_t i = 0; i < _outward_incident_weight.size(); ++i) {
      sep_nodes._outward_incident_weight[i] = _outward_incident_weight[i];
    }
  }, [&] {
    sep_nodes._graph_nodes_begin.resize(_graph_nodes_begin.size());
    memcpy(sep_nodes._graph_nodes_begin.data(), _graph_nodes_begin.data(),
            sizeof(HyperedgeID) * _graph_nodes_begin.size());
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

  sep_nodes._outward_incident_weight.resize(_outward_incident_weight.size());
  for (size_t i = 0; i < _outward_incident_weight.size(); ++i) {
    sep_nodes._outward_incident_weight[i] = _outward_incident_weight[i];
  }

  sep_nodes._graph_nodes_begin.resize(_graph_nodes_begin.size());
  memcpy(sep_nodes._graph_nodes_begin.data(), _graph_nodes_begin.data(),
          sizeof(HyperedgeID) * _graph_nodes_begin.size());

  sep_nodes._inward_edges.resize(_inward_edges.size());
  memcpy(sep_nodes._inward_edges.data(), _inward_edges.data(),
          sizeof(Edge) * _inward_edges.size());

  sep_nodes._outward_edges.resize(_outward_edges.size());
  memcpy(sep_nodes._outward_edges.data(), _outward_edges.data(),
          sizeof(Edge) * _outward_edges.size());

  return sep_nodes;
}

} // namespace

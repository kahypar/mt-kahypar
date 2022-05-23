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
  ASSERT(_num_graph_nodes == _outward_incident_weight.size());
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

void SeparatedNodes::contract(const vec<HypernodeID>& communities, const HypernodeID& num_coarsened_graph_nodes) {
  ASSERT(_num_nodes + 1 == _nodes.size());
  ASSERT(_num_edges == _inward_edges.size());
  ASSERT(_num_graph_nodes == _outward_incident_weight.size() && _num_graph_nodes == communities.size());
  ASSERT(_graph_nodes_begin.empty());

  auto update_incident_weight = [&] {
    vec<parallel::IntegralAtomicWrapper<HyperedgeWeight>> new_incident_weight;
    new_incident_weight.assign(num_coarsened_graph_nodes, parallel::IntegralAtomicWrapper<HyperedgeWeight>(0));
    tbb::parallel_for(ID(0), _num_graph_nodes, [&](const HypernodeID& pos) {
      if (communities[pos] != kInvalidHypernode) {
        new_incident_weight[communities[pos]].fetch_add(_outward_incident_weight[pos]);
      }
    });

    std::swap(_outward_incident_weight, new_incident_weight);
    _num_graph_nodes = num_coarsened_graph_nodes;
  };

  auto update_inward_edges = [&] {
    Array<HyperedgeID> node_sizes;

    tbb::parallel_invoke([&] {
      // update edge targets
      tbb::parallel_for(ID(0), _num_edges, [&](const HyperedgeID& pos) {
        Edge& e = _inward_edges[pos];
        e.target = communities[e.target];
      });
    }, [&] {
      node_sizes.resize(_num_nodes);
    });

    // TODO: deduplicate code
    // deduplicate edges
    tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID& pos) {
      const HyperedgeID edges_start = _nodes[pos].begin;
      const HyperedgeID edges_end = _nodes[pos + 1].begin;
      std::sort(_inward_edges.begin() + edges_start, _inward_edges.begin() + edges_end,
                [](const Edge& e1, const Edge& e2) { return e1.target < e2.target; });

        // Deduplicate, aggregate weights and calculate minimum unique id
        //
        // <-- deduplicated --> <-- already processed --> <-- to be processed --> <-- invalid edges -->
        //                    ^                         ^
        // valid_edge_index ---        tmp_edge_index ---
        size_t valid_edge_index = edges_start;
        size_t tmp_edge_index = edges_start + 1;
        while (tmp_edge_index < edges_end && _inward_edges[tmp_edge_index].target != kInvalidHypernode) {
          HEAVY_COARSENING_ASSERT(
            [&](){
              size_t i = edges_start;
              for (; i <= valid_edge_index; ++i) {
                if (_inward_edges[i].target == kInvalidHypernode) {
                  return false;
                } else if ((i + 1 <= valid_edge_index) &&
                  _inward_edges[i].target >= _inward_edges[i + 1].target) {
                  return false;
                }
              }
              for (; i < tmp_edge_index; ++i) {
                if (_inward_edges[i].target != kInvalidHypernode) {
                  return false;
                }
              }
              return true;
            }(),
            "Invariant violated while deduplicating incident edges!"
          );

          Edge& valid_edge = _inward_edges[valid_edge_index];
          Edge& next_edge = _inward_edges[tmp_edge_index];
          if (next_edge.target != kInvalidHypernode) {
            if (valid_edge.target == next_edge.target) {
              valid_edge.weight += next_edge.weight;
              next_edge.target = kInvalidHypernode;
            } else {
              std::swap(_inward_edges[++valid_edge_index], next_edge);
            }
            ++tmp_edge_index;
          }
        }
        const bool is_non_empty = (edges_start < edges_end) &&
                                  (_inward_edges[valid_edge_index].target != kInvalidHypernode);
        const HyperedgeID contracted_size = is_non_empty ? (valid_edge_index - edges_start + 1) : 0;
        node_sizes[pos] = contracted_size;
    });

    parallel::TBBPrefixSum<HyperedgeID, Array> degree_mapping(node_sizes);
    tbb::parallel_scan(tbb::blocked_range<size_t>(
            ID(0), static_cast<size_t>(_num_nodes)), degree_mapping);
    const HyperedgeID coarsened_num_edges = degree_mapping.total_sum();

    // copy edges
    vec<Edge> new_edges;
    new_edges.resize(coarsened_num_edges);
    tbb::parallel_for(ID(0), _num_nodes, [&](const HyperedgeID& pos) {
      const HyperedgeID old_edges_start = _nodes[pos].begin;
      const HyperedgeID new_edges_start = degree_mapping[pos];
      for (size_t index = 0; index < degree_mapping.value(pos); ++index) {
        ASSERT(old_edges_start + index < _inward_edges.size() && new_edges_start + index < new_edges.size());
        new_edges[new_edges_start + index] = _inward_edges[old_edges_start + index];
      }
    });

    // update nodes
    tbb::parallel_for(0UL, _nodes.size(), [&](const size_t& pos) {
      _nodes[pos].begin = degree_mapping[pos];
    });

    std::swap(_inward_edges, new_edges);
    _num_edges = coarsened_num_edges;
  };

  tbb::parallel_invoke(update_incident_weight, update_inward_edges);
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

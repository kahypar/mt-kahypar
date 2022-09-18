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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"

#include <tbb/parallel_reduce.h>
#include <tbb/parallel_sort.h>

namespace mt_kahypar::ds {

// helper function
template<typename T>
T parallel_max(const vec<T>& data, const T& invalid) {
    return tbb::parallel_reduce(
            tbb::blocked_range<size_t>(0UL, data.size()), std::numeric_limits<T>::min(),
            [&](tbb::blocked_range<size_t>& range, T init) -> T {
            T tmp_max = init;
            for ( size_t i = range.begin(); i < range.end(); ++i ) {
              if (data[i] != invalid) {
                tmp_max = std::max(tmp_max, data[i]);
              }
            }
            return tmp_max;
            },
            [](const T& a, const T& b) { return std::max(a, b); });
}

HypernodeID SeparatedNodes::popBatch() {
  const auto [index, weight, num_internal_edges] = _batch_indices_and_weights[_batch_indices_and_weights.size() - 2];
  ASSERT(num_internal_edges != kInvalidHyperedge, V(_batch_indices_and_weights.size()));
  _batch_indices_and_weights.pop_back();
  _num_nodes = index;
  _nodes.resize(_num_nodes + 1);
  _internal_edges.resize(num_internal_edges);
  _num_edges = _nodes.back().begin;
  _inward_edges.resize(_num_edges);
  _nodes[_num_nodes] = Node(kInvalidHypernode, _num_edges, 0);
  _total_weight = weight;

  // TODO: outward edges? incident weight?!
  ASSERT(_batch_indices_and_weights.size() >= 2);
  return index;
}

void SeparatedNodes::cleanBatchState() {
  _batch_indices_and_weights.assign(2, {_num_nodes, _total_weight, _internal_edges.size()});
}

void SeparatedNodes::revealNextBatch() {
  ++_hidden_nodes_batch_index;
  ASSERT(_hidden_nodes_batch_index < _batch_indices_and_weights.size());
}

void SeparatedNodes::setSavepoint() {
  vec<Edge> copied_edges;
  vec<HyperedgeID> edge_indices;
  vec<InternalEdge> copied_internal_edges;
  tbb::parallel_invoke([&] {
    copied_edges = _inward_edges;
  }, [&] {
    edge_indices.resize(_nodes.size());
    tbb::parallel_for(0UL, edge_indices.size(), [&](const size_t& pos) {
      edge_indices[pos] = _nodes[pos].begin;
    });
  }, [&] {
    copied_internal_edges = _internal_edges;
  });
  _savepoints.push_back({std::move(copied_edges), std::move(edge_indices),
                         std::move(copied_internal_edges), _num_graph_nodes});
}

void SeparatedNodes::restoreSavepoint() {
  ASSERT(!_savepoints.empty());
  Memento& m = _savepoints.back();
  _inward_edges = std::move(m.edges);
  _internal_edges = std::move(m.internal_edges);
  ASSERT(_nodes.size() == m.edge_indices.size());
  tbb::parallel_for(0UL, _nodes.size(), [&](const size_t& pos) {
    _nodes[pos].begin = m.edge_indices[pos];
  });
  _savepoints.pop_back();
  _num_edges = _inward_edges.size();
  _num_graph_nodes = m.num_graph_nodes;
  _outward_incident_weight.clear(); // not correct anymore
}

SeparatedNodes SeparatedNodes::createCopyFromSavepoint() {
  ASSERT(!_savepoints.empty());
  Memento& m = _savepoints.back();
  ASSERT(_nodes.size() >= m.edge_indices.size() && m.edge_indices.back() == m.edges.size());
  const HypernodeID num_nodes = m.edge_indices.size() - 1;
  _outward_incident_weight.clear(); // not correct anymore

  SeparatedNodes sep_nodes(m.num_graph_nodes);
  sep_nodes._num_nodes = num_nodes;
  sep_nodes._num_graph_nodes = m.num_graph_nodes;
  sep_nodes._num_edges = m.edges.size();
  sep_nodes._total_weight = _total_weight;
  sep_nodes._inward_edges = std::move(m.edges);
  sep_nodes._internal_edges = std::move(m.internal_edges);

  tbb::parallel_invoke([&] {
    sep_nodes._nodes.resize(num_nodes + 1);
    memcpy(sep_nodes._nodes.data(), _nodes.data(), sizeof(Node) * num_nodes);
    sep_nodes._nodes.back() = Node(kInvalidHypernode, sep_nodes._num_edges, 0);
  }, [&] {
    sep_nodes._batch_indices_and_weights.clear();
    for (size_t i = 0; i < _batch_indices_and_weights.size() && std::get<0>(_batch_indices_and_weights[i]) <= num_nodes; ++i) {
      sep_nodes._hidden_nodes_batch_index = i;
      sep_nodes._batch_indices_and_weights.push_back(_batch_indices_and_weights[i]);
    }
    ASSERT(sep_nodes._hidden_nodes_batch_index > 0);
  });
  tbb::parallel_for(0UL, sep_nodes._nodes.size(), [&](const size_t& pos) {
    sep_nodes._nodes[pos].begin = m.edge_indices[pos];
  });

  _savepoints.pop_back();
  return sep_nodes;
}

void SeparatedNodes::addNodes(const Hypergraph& hg,
                              const vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>>& nodes,
                              const vec<Edge>& edges) {
  ASSERT(hg.initialNumNodes() == _num_graph_nodes);
  ASSERT(_num_nodes + 1 == _nodes.size());
  ASSERT(_num_edges == _inward_edges.size());
  // ASSERT(_num_graph_nodes == _outward_incident_weight.size());
  ASSERT(_graph_nodes_begin.empty());
  ASSERT(nodes.empty() || std::get<1>(nodes[0]) == 0);

  auto update_nodes = [&] {
    _nodes.resize(_num_nodes + nodes.size() + 1);
    tbb::parallel_for(0UL, nodes.size(), [&](const size_t& pos) {
      auto [original_node, begin, weight] = nodes[pos];
      _nodes[_num_nodes + pos] = Node(original_node, _num_edges + begin, weight);
    });
    _nodes[_num_nodes + nodes.size()] = Node(kInvalidHypernode, _num_edges + edges.size(), 0);
  };

  auto update_inward_edges_and_incident_weight = [&] {
    _inward_edges.resize(_num_edges + edges.size());
    tbb::parallel_for(0UL, edges.size(), [&](const size_t& pos) {
      const Edge& e = edges[pos];
      _inward_edges[_num_edges + pos] = e;
      if (!_outward_incident_weight.empty()) {
        _outward_incident_weight[e.target].fetch_add(e.weight);
      }
    });
  };

  StreamingVector<InternalEdge> edges_stream;
  Array<HypernodeID> graph_node_to_sep_node;
  auto add_intern_edges = [&] {
    graph_node_to_sep_node.resize(_num_graph_nodes, kInvalidHypernode);
    // initialize node set
    tbb::parallel_for(0UL, nodes.size(), [&](const size_t& pos) {
      const HypernodeID original_node = std::get<0>(nodes[pos]);
      if (original_node != kInvalidHypernode) {
        graph_node_to_sep_node[original_node] = pos + _num_nodes;
      }
    });
    // edges between added nodes
    tbb::parallel_for(0UL, nodes.size(), [&](const size_t& pos) {
      const HypernodeID original_node = std::get<0>(nodes[pos]);
      if (original_node != kInvalidHypernode) {
        for (const HyperedgeID& e: hg.incidentEdges(original_node)) {
          if (graph_node_to_sep_node[hg.edgeTarget(e)] != kInvalidHypernode) {
            const HypernodeID target = graph_node_to_sep_node[hg.edgeTarget(e)];
            if (pos + _num_nodes < target) {
              edges_stream.stream(pos + _num_nodes, target, hg.edgeWeight(e));
            }
          }
        }
      }
    });
  };

  auto update_total_weight = [&] {
    _total_weight += tbb::parallel_reduce(
            tbb::blocked_range<HypernodeID>(0UL, nodes.size()), 0,
              [&](const tbb::blocked_range<HypernodeID>& range, HypernodeWeight init) {
                HypernodeWeight weight = init;
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    weight += std::get<2>(nodes[i]);
                }
                return weight;
              }, std::plus<>());
  };

  tbb::parallel_invoke(update_nodes, update_inward_edges_and_incident_weight, update_total_weight, add_intern_edges);

  // edges from already separated nodes
  tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID& node) {
    for (const Edge& e: inwardEdges(node)) {
      if (graph_node_to_sep_node[e.target] != kInvalidHypernode) {
        edges_stream.stream(node, graph_node_to_sep_node[e.target], e.weight);
      }
    }
  });
  edges_stream.append_parallel(_internal_edges);

  _num_nodes += nodes.size();
  _num_edges += edges.size();
  _batch_indices_and_weights.push_back({_num_nodes, _total_weight, _internal_edges.size()});
}

void SeparatedNodes::contract(const vec<HypernodeID>& communities, const HypernodeID& num_coarsened_graph_nodes) {
  ASSERT(_num_nodes + 1 == _nodes.size());
  ASSERT(_num_edges == _inward_edges.size());
  ASSERT(/*_num_graph_nodes == _outward_incident_weight.size() && */_num_graph_nodes == communities.size());
  ASSERT(_graph_nodes_begin.empty());

  auto update_incident_weight = [&] {
    vec<parallel::IntegralAtomicWrapper<HyperedgeWeight>> new_incident_weight;
    new_incident_weight.assign(num_coarsened_graph_nodes, parallel::IntegralAtomicWrapper<HyperedgeWeight>(0));
    tbb::parallel_for(ID(0), static_cast<HypernodeID>(_outward_incident_weight.size()), [&](const HypernodeID& pos) {
      if (communities[pos] != kInvalidHypernode) {
        new_incident_weight[communities[pos]].fetch_add(_outward_incident_weight[pos]);
      }
    });

    std::swap(_outward_incident_weight, new_incident_weight);
    _num_graph_nodes = num_coarsened_graph_nodes;
  };

  auto update_inward_edges = [&] {
    // update edge targets
    tbb::parallel_for(ID(0), _num_edges, [&](const HyperedgeID& pos) {
      Edge& e = _inward_edges[pos];
      e.target = communities[e.target];
    });

    deduplicateEdges();
  };

  tbb::parallel_invoke(update_incident_weight, update_inward_edges);

  ASSERT(
    [&] {
      for (const Edge& e: _inward_edges) {
        if (e.target == kInvalidHypernode) { return false; }
      }
      return true;
    }()
  );
}

SeparatedNodes SeparatedNodes::coarsen(vec<HypernodeID>& communities) const {
  ASSERT(communities.size() == numVisibleNodes());
  ASSERT(_hidden_nodes_batch_index > 0);
  SeparatedNodes other(_num_graph_nodes);

  Array<HypernodeID> mapping;

  // Compute vertex ids of coarse graph with a parallel prefix sum
  mapping.assign(numVisibleNodes(), 0);
  tbb::parallel_for(ID(0), numVisibleNodes(), [&](const HyperedgeID& node) {
    ASSERT(communities[node] < mapping.size());
    mapping[communities[node]] = 1;
  });

  // Prefix sum determines vertex ids in coarse graph
  parallel::TBBPrefixSum<HyperedgeID, Array> mapping_prefix_sum(mapping);
  tbb::parallel_scan(tbb::blocked_range<size_t>(ID(0), numVisibleNodes()), mapping_prefix_sum);
  const HypernodeID coarsened_num_nodes = mapping_prefix_sum.total_sum();

  Array<parallel::IntegralAtomicWrapper<HypernodeWeight>> node_weights;
  node_weights.assign(coarsened_num_nodes, parallel::IntegralAtomicWrapper<HypernodeWeight>(0));
  Array<parallel::IntegralAtomicWrapper<HyperedgeID>> tmp_num_incident_edges;
  // sentinel for first element
  tmp_num_incident_edges.assign(coarsened_num_nodes + 1, parallel::IntegralAtomicWrapper<HyperedgeID>(0));

  tbb::parallel_for(ID(0), numVisibleNodes(), [&](const HyperedgeID& node) {
    const HypernodeID coarse_node = mapping_prefix_sum[communities[node]];
    communities[node] = coarse_node;
    node_weights[coarse_node].fetch_add(nodeWeight(node));
    tmp_num_incident_edges[coarse_node + 1].fetch_add(inwardDegree(node));
  });

  parallel::TBBPrefixSum<parallel::IntegralAtomicWrapper<HyperedgeID>, Array>
          tmp_incident_edges_prefix_sum(tmp_num_incident_edges);
  tbb::parallel_scan(tbb::blocked_range<size_t>(
          ID(0), static_cast<size_t>(coarsened_num_nodes + 1)), tmp_incident_edges_prefix_sum);
  ASSERT(tmp_incident_edges_prefix_sum.total_sum() == _nodes[numVisibleNodes()].begin);
  const HypernodeID total_coarsened_nodes = coarsened_num_nodes + _num_nodes - numVisibleNodes();

  tbb::parallel_invoke([&] {
    other._nodes.resize(total_coarsened_nodes);
  }, [&] {
    other._inward_edges.resize(_num_edges);
  }, [&] {
    other._outward_incident_weight = _outward_incident_weight;
  }, [&] {
    other._internal_edges = _internal_edges;
  });

  tbb::parallel_for(ID(0), total_coarsened_nodes, [&](const HyperedgeID& node) {
    if (node < coarsened_num_nodes) {
      other._nodes[node] = Node(kInvalidHypernode, tmp_num_incident_edges[node], node_weights[node]);
    } else {
      other._nodes[node] = _nodes[node + numVisibleNodes() - coarsened_num_nodes];
    }
  });
  other._nodes.emplace_back(kInvalidHypernode, _num_edges, 0);

  tbb::parallel_for(ID(0), _num_nodes, [&](const HyperedgeID& node) {
    HypernodeID coarse_node;
    HyperedgeID start;
    if (node < numVisibleNodes()) {
      coarse_node = communities[node];
      // modify the underlying data of the prefix sum
      start = tmp_num_incident_edges[coarse_node].fetch_add(inwardDegree(node));
    } else {
      coarse_node = node - numVisibleNodes() + coarsened_num_nodes;
      start = _nodes[coarse_node].begin;
    }
    for (HyperedgeID i = 0; i < inwardDegree(node); ++i) {
      other._inward_edges[start + i] = _inward_edges[_nodes[node].begin + i];
    }
  });

  other._num_nodes = total_coarsened_nodes;
  other._num_edges = _num_edges;
  other._total_weight = _total_weight;

  auto internal_edge_visible = [&](const InternalEdge& e) {
    return e.pin0 < coarsened_num_nodes && e.pin1 < coarsened_num_nodes;
  };

  tbb::parallel_invoke([&] {
    other.deduplicateEdges();
  }, [&] {
    auto map_node = [&](const HypernodeID& node) {
      if (node < numVisibleNodes()) {
        return communities[node];
      } else {
        return node - numVisibleNodes() + coarsened_num_nodes;
      }
    };
    tbb::parallel_for(0UL, other._internal_edges.size(), [&](const size_t& pos) {
      InternalEdge& e = other._internal_edges[pos];
      ASSERT(e.pin0 != kInvalidHypernode && e.pin1 != kInvalidHypernode);
      e.pin0 = map_node(e.pin0);
      e.pin1 = map_node(e.pin1);
      ASSERT(e.pin0 < other._num_nodes && e.pin1 < other._num_nodes);
    });
    other.deduplicateInternalEdges();
    tbb::parallel_sort(other._internal_edges.begin(), other._internal_edges.end(),
      [&](const InternalEdge& l, const InternalEdge& r) {
        return internal_edge_visible(l) && !internal_edge_visible(r);
      });
  });
  HyperedgeID visible_internal_edges = 0;
  tbb::parallel_for(ID(0), other.numInternalEdges(), [&](const HyperedgeID& e) {
    if (internal_edge_visible(other.internalEdge(e))
        && (e + 1 == other.numInternalEdges() || !internal_edge_visible(other.internalEdge(e + 1)))) {
      visible_internal_edges = e + 1;
    }
  });

  other._batch_indices_and_weights = { {coarsened_num_nodes,
        std::get<1>(_batch_indices_and_weights.at(_hidden_nodes_batch_index)), visible_internal_edges} };
  other._batch_indices_and_weights.push_back(other._batch_indices_and_weights[0]);
  for (size_t i = _hidden_nodes_batch_index + 1; i < _batch_indices_and_weights.size(); ++i) {
    auto batch = _batch_indices_and_weights[i];
    std::get<0>(batch) -= (numVisibleNodes() - coarsened_num_nodes);
    std::get<2>(batch) = kInvalidHyperedge;
    other._batch_indices_and_weights.push_back(batch);
  }
  return other;
}

void SeparatedNodes::initializeOutwardEdges() {
  ASSERT(_graph_nodes_begin.empty() && _outward_edges.empty());

  Array<parallel::IntegralAtomicWrapper<HyperedgeID>> tmp_node_degree;

  tbb::parallel_invoke([&] {
    // we use index 0 for the edge position later
    tmp_node_degree.resize(_num_graph_nodes + 1);
  }, [&] {
    _graph_nodes_begin.resize(_num_graph_nodes + 1);
  }, [&] {
    _outward_edges.resize(_num_edges);
  });

  tbb::parallel_for(ID(0), _num_edges, [&](const HyperedgeID& pos) {
    const Edge& e = _inward_edges[pos];
    tmp_node_degree[e.target + 1].fetch_add(1);
  });

  parallel::TBBPrefixSum<parallel::IntegralAtomicWrapper<HyperedgeID>, Array> 
          degree_mapping(tmp_node_degree);
  tbb::parallel_scan(tbb::blocked_range<size_t>(
          ID(0), static_cast<size_t>(_num_graph_nodes + 1)), degree_mapping);
  ASSERT(degree_mapping.total_sum() == _num_edges);

  tbb::parallel_for(0UL, _graph_nodes_begin.size(), [&](const size_t& pos) {
    _graph_nodes_begin[pos] = tmp_node_degree[pos].load();
  });

  tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID& node) {
    const HyperedgeID edges_start = _nodes[node].begin;
    const HyperedgeID edges_end = _nodes[node + 1].begin;
    for (HyperedgeID pos = edges_start; pos < edges_end; ++pos) {
      const Edge& e = _inward_edges[pos];
      const HyperedgeID index = tmp_node_degree[e.target].fetch_add(1);
      _outward_edges[index] = Edge(node, _inward_edges[pos].weight);
    }
  });

  // ASSERT( // TODO(maas): failing assertion here
  //   [&] {
  //     for (HypernodeID u = 0; u < _num_graph_nodes; ++u) {
  //       HyperedgeWeight incident_weight = 0;
  //       for (const auto& [target, weight]: outwardEdges(u)) {
  //         incident_weight += weight;
  //       }
  //       if (incident_weight != outwardIncidentWeight(u)) {
  //         return false;
  //       }
  //     }
  //     return true;
  //   }()
  // );
}

SeparatedNodes SeparatedNodes::extract(PartitionID block, const vec<HypernodeID>& graph_node_mapping,
                                       const vec<CAtomic<PartitionID>>& part_ids) const {
  ASSERT(/*_num_graph_nodes == _outward_incident_weight.size() &&*/ _num_graph_nodes == graph_node_mapping.size());
  ASSERT(part_ids.size() >= _num_nodes);
  // ASSERT(_graph_nodes_begin.empty());

  auto get_part_id = [&](HypernodeID node) {
    return part_ids[node].load(std::memory_order_relaxed);
  };

  // TODO: parallelize with tmp_node_degree
  Array<HypernodeID> sep_nodes_active;
  sep_nodes_active.resize(_num_nodes, 0);
  tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID& node) {
    if (get_part_id(node) == block) {
      sep_nodes_active[node] = 1;
    }
  });

  parallel::TBBPrefixSum<HypernodeID, Array> sep_node_mapping(sep_nodes_active);
  tbb::parallel_scan(tbb::blocked_range<size_t>(
          ID(0), static_cast<size_t>(_num_nodes)), sep_node_mapping);

  Array<HyperedgeID> tmp_node_degree;
  tmp_node_degree.resize(_num_nodes, 0);
  tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID& node) {
    if (get_part_id(node) == block) {
      tmp_node_degree[node] = inwardDegree(node);
    }
  });

  parallel::TBBPrefixSum<HyperedgeID, Array> degree_mapping(tmp_node_degree);
  tbb::parallel_scan(tbb::blocked_range<size_t>(
          ID(0), static_cast<size_t>(_num_nodes)), degree_mapping);


  SeparatedNodes other(_num_graph_nodes);
  other._num_nodes = sep_node_mapping.total_sum();
  other._num_graph_nodes = _num_graph_nodes;
  other._num_edges = degree_mapping.total_sum();

  auto set_nodes = [&] {
    vec<Node> tmp_nodes;
    tmp_nodes.resize(sep_node_mapping.total_sum());
    tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID& node) {
      if (get_part_id(node) == block) {
        Node& new_node = tmp_nodes[sep_node_mapping[node]];
        new_node = _nodes[node];
        new_node.begin = degree_mapping[node];
        new_node.original_node = kInvalidHypernode; // doesn't make sense in extracted graph
      }
    });
    tmp_nodes.emplace_back(kInvalidHypernode, degree_mapping.total_sum(), 0);
    other._nodes = std::move(tmp_nodes);
  };

  auto set_edges = [&] {
    // copy the weight vector
    vec<parallel::IntegralAtomicWrapper<HyperedgeWeight>> tmp_outward_weight = _outward_incident_weight;
    vec<Edge> tmp_inward_edges;
    tmp_inward_edges.resize(degree_mapping.total_sum());
    tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID& node) {
      if (get_part_id(node) == block) {
        const HyperedgeID start_new = degree_mapping[node];
        const HyperedgeID start_old = _nodes[node].begin;
        for (HyperedgeID i = 0; start_new + i < degree_mapping[node + 1]; ++i) {
          tmp_inward_edges[start_new + i] = _inward_edges[start_old + i];
        }
      } else {
        for (HyperedgeID i = _nodes[node].begin; i < _nodes[node + 1].begin; ++i) {
          const Edge& e = _inward_edges[i];
          if (!tmp_outward_weight.empty()) {
            tmp_outward_weight[e.target].fetch_sub(e.weight);
          }
        }
      }
    });
    other._outward_incident_weight = std::move(tmp_outward_weight);
    other._inward_edges = std::move(tmp_inward_edges);
  };

  auto set_internal_edges = [&] {
    StreamingVector<InternalEdge> internal_edges;
    tbb::parallel_for(0UL, _internal_edges.size(), [&](const size_t& pos) {
      const InternalEdge& e = _internal_edges[pos];
      if (get_part_id(e.pin0) == block && get_part_id(e.pin1) == block) {
        internal_edges.stream(sep_node_mapping[e.pin0], sep_node_mapping[e.pin1], e.weight);
      }
    });
    other._internal_edges = internal_edges.copy_parallel();
  };

  auto set_total_weight = [&] {
    other._total_weight = tbb::parallel_reduce(
            tbb::blocked_range<HypernodeID>(0UL, _num_nodes), 0,
              [&](const tbb::blocked_range<HypernodeID>& range, HypernodeWeight init) {
                HypernodeWeight weight = init;
                for (size_t i = range.begin(); i < range.end(); ++i) {
                  if (get_part_id(i) == block) {
                    weight += node(i).weight;
                  }
                }
                return weight;
              }, std::plus<>());
  };

  tbb::parallel_invoke(set_nodes, set_edges, set_internal_edges, set_total_weight);
  const HypernodeID max_node_id = parallel_max(graph_node_mapping, kInvalidHypernode);
  other.contract(graph_node_mapping, max_node_id + 1);
  other.cleanBatchState();
  return other;
}

void SeparatedNodes::deduplicateEdges() {
  Array<HyperedgeID> node_sizes;
  node_sizes.resize(_num_nodes);

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
}

void SeparatedNodes::deduplicateInternalEdges() {
  if (!_internal_edges.empty()) {
    const HyperedgeID edges_end = _internal_edges.size();
    tbb::parallel_sort(_internal_edges.begin(), _internal_edges.end(), [](const InternalEdge& e1, const InternalEdge& e2) {
      return e1.pin0 < e2.pin0 || (e1.pin0 == e2.pin0 && e1.pin1 < e2.pin1);
    });

    // Deduplicate and aggregate weights
    //
    // <-- deduplicated --> <-- already processed --> <-- to be processed -->
    //                    ^                         ^
    // valid_edge_index ---        tmp_edge_index ---
    size_t valid_edge_index = 0;
    size_t tmp_edge_index = 1;
    while (tmp_edge_index < edges_end) {
      InternalEdge& valid_edge = _internal_edges[valid_edge_index];
      InternalEdge& next_edge = _internal_edges[tmp_edge_index];
      ASSERT(next_edge.pin0 != kInvalidHypernode && next_edge.pin1 != kInvalidHypernode);
      if (valid_edge.pin0 == next_edge.pin0 && valid_edge.pin1 == next_edge.pin1) {
        valid_edge.weight += next_edge.weight;
      } else if (next_edge.pin0 != next_edge.pin1) {
        std::swap(_internal_edges[++valid_edge_index], next_edge);
      }
      ++tmp_edge_index;
    }

    size_t new_size = valid_edge_index + 1;
    // special case: first edge is single pin
    InternalEdge& first = _internal_edges[0];
    if (first.pin0 == first.pin1) {
      std::swap(first, _internal_edges[valid_edge_index]);
      --new_size;
    }
    _internal_edges.resize(new_size);
  }
}

// ! Copy in parallel
SeparatedNodes SeparatedNodes::copy(parallel_tag_t) const {
  SeparatedNodes sep_nodes(_num_graph_nodes);

  sep_nodes._num_nodes = _num_nodes;
  sep_nodes._num_graph_nodes = _num_graph_nodes;
  sep_nodes._num_edges = _num_edges;
  sep_nodes._total_weight = _total_weight;
  sep_nodes._hidden_nodes_batch_index = _hidden_nodes_batch_index;

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
    if (!_graph_nodes_begin.empty()) {
      sep_nodes._graph_nodes_begin.resize(_graph_nodes_begin.size());
      memcpy(sep_nodes._graph_nodes_begin.data(), _graph_nodes_begin.data(),
              sizeof(HyperedgeID) * _graph_nodes_begin.size());
    }
  }, [&] {
    sep_nodes._inward_edges.resize(_inward_edges.size());
    memcpy(sep_nodes._inward_edges.data(), _inward_edges.data(),
            sizeof(Edge) * _inward_edges.size());
  }, [&] {
    if (!_outward_edges.empty()) {
      sep_nodes._outward_edges.resize(_outward_edges.size());
      memcpy(sep_nodes._outward_edges.data(), _outward_edges.data(),
              sizeof(Edge) * _outward_edges.size());
    }
  }, [&] {
    sep_nodes._internal_edges.resize(_internal_edges.size());
    memcpy(sep_nodes._internal_edges.data(), _internal_edges.data(),
            sizeof(InternalEdge) * _internal_edges.size());
  }, [&] {
    sep_nodes._batch_indices_and_weights.resize(_batch_indices_and_weights.size());
    for (size_t i = 0; i < _batch_indices_and_weights.size(); ++i) {
      sep_nodes._batch_indices_and_weights[i] = _batch_indices_and_weights[i];
    }
  });
  return sep_nodes;
}

// ! Copy in parallel
SeparatedNodes SeparatedNodes::copy_first_batch(parallel_tag_t) const {
  SeparatedNodes sep_nodes(_num_graph_nodes);

  const HypernodeID num_nodes = std::get<0>(_batch_indices_and_weights[0]);
  const HyperedgeID num_edges = _nodes[num_nodes].begin;
  const HyperedgeID num_internal_edges = std::get<2>(_batch_indices_and_weights[0]);
  sep_nodes._num_nodes = num_nodes;
  sep_nodes._num_graph_nodes = _num_graph_nodes;
  sep_nodes._num_edges = num_edges;
  sep_nodes._total_weight = std::get<1>(_batch_indices_and_weights[0]);
  sep_nodes._batch_indices_and_weights = {_batch_indices_and_weights[0], _batch_indices_and_weights[1]};

  tbb::parallel_invoke([&] {
    sep_nodes._nodes.resize(num_nodes + 1);
    memcpy(sep_nodes._nodes.data(), _nodes.data(), sizeof(Node) * num_nodes);
    sep_nodes._nodes.push_back(Node(kInvalidHypernode, num_edges, 0));
  }, [&] {
    sep_nodes._outward_incident_weight.resize(_outward_incident_weight.size());
    for (size_t i = 0; i < _outward_incident_weight.size(); ++i) {
      sep_nodes._outward_incident_weight[i] = _outward_incident_weight[i];
    }
  }, [&] {
    sep_nodes._inward_edges.resize(num_edges);
    memcpy(sep_nodes._inward_edges.data(), _inward_edges.data(),
            sizeof(Edge) * num_edges);
  }, [&] {
    sep_nodes._internal_edges.resize(num_internal_edges);
    memcpy(sep_nodes._internal_edges.data(), _internal_edges.data(),
            sizeof(InternalEdge) * num_internal_edges);
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
  sep_nodes._hidden_nodes_batch_index = _hidden_nodes_batch_index;

  sep_nodes._nodes.resize(_nodes.size());
  memcpy(sep_nodes._nodes.data(), _nodes.data(),
          sizeof(Node) * _nodes.size());

  sep_nodes._outward_incident_weight.resize(_outward_incident_weight.size());
  for (size_t i = 0; i < _outward_incident_weight.size(); ++i) {
    sep_nodes._outward_incident_weight[i] = _outward_incident_weight[i];
  }

  if (!_graph_nodes_begin.empty()) {
    sep_nodes._graph_nodes_begin.resize(_graph_nodes_begin.size());
    memcpy(sep_nodes._graph_nodes_begin.data(), _graph_nodes_begin.data(),
            sizeof(HyperedgeID) * _graph_nodes_begin.size());
  }

  sep_nodes._inward_edges.resize(_inward_edges.size());
  memcpy(sep_nodes._inward_edges.data(), _inward_edges.data(),
          sizeof(Edge) * _inward_edges.size());

  if (!_outward_edges.empty()) {
    sep_nodes._outward_edges.resize(_outward_edges.size());
    memcpy(sep_nodes._outward_edges.data(), _outward_edges.data(),
            sizeof(Edge) * _outward_edges.size());
  }

  sep_nodes._internal_edges.resize(_internal_edges.size());
  memcpy(sep_nodes._internal_edges.data(), _internal_edges.data(),
          sizeof(InternalEdge) * _internal_edges.size());

  sep_nodes._batch_indices_and_weights.resize(_batch_indices_and_weights.size());
  for (size_t i = 0; i < _batch_indices_and_weights.size(); ++i) {
    sep_nodes._batch_indices_and_weights[i] = _batch_indices_and_weights[i];
  }

  return sep_nodes;
}


void SepNodesStack::contractToNLevels(size_t n) {
  ASSERT(n <= numLevels() && n > 1);
  vec<HypernodeID>& mapping = _mappings[n - 2];
  tbb::parallel_for(0UL, mapping.size(), [&](const size_t& pos) {
    HypernodeID mapped_node = mapping[pos];
    for (size_t i = n - 1; i < _mappings.size(); ++i) {
      mapped_node = _mappings[i][mapped_node];
    }
    mapping[pos] = mapped_node;
  });

  std::unique_ptr<SeparatedNodes> last_nodes = std::move(_data.back());
  _data.pop_back();
  while (_data.size() >= n) {
    _data.pop_back();
    _mappings.pop_back();
  }
  _data.push_back(std::move(last_nodes));
}

} // namespace

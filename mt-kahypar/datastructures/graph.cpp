/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "graph.h"


#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_invoke.h>
#include <tbb/enumerable_thread_specific.h>

#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/parallel/parallel_counting_sort.h"

namespace mt_kahypar::ds {

  Graph::Graph(Hypergraph& hypergraph, const LouvainEdgeWeight edge_weight_type, bool is_graph) :
          _num_nodes(0),
          _num_arcs(0),
          _total_volume(0),
          _total_weight(0),
          _max_degree(0),
          _indices(),
          _arcs(),
          _node_volumes(),
          _node_weights(),
          _tmp_graph_buffer(nullptr) {

    switch( edge_weight_type ) {
      case LouvainEdgeWeight::uniform:
        construct(hypergraph, is_graph,
                  [&](const HyperedgeWeight edge_weight,
                      const HypernodeID,
                      const HyperedgeID) {
                    return static_cast<ArcWeight>(edge_weight);
                  });
        break;
      case LouvainEdgeWeight::non_uniform:
        construct(hypergraph, is_graph,
                  [&](const HyperedgeWeight edge_weight,
                      const HypernodeID edge_size,
                      const HyperedgeID) {
                    return static_cast<ArcWeight>(edge_weight) /
                           static_cast<ArcWeight>(edge_size);
                  });
        break;
      case LouvainEdgeWeight::degree:
        construct(hypergraph, is_graph,
                  [&](const HyperedgeWeight edge_weight,
                      const HypernodeID edge_size,
                      const HyperedgeID node_degree) {
                    return static_cast<ArcWeight>(edge_weight) *
                           (static_cast<ArcWeight>(node_degree) /
                            static_cast<ArcWeight>(edge_size));
                  });
        break;
      case LouvainEdgeWeight::hybrid:
      case LouvainEdgeWeight::UNDEFINED:
      ERROR("No valid louvain edge weight");
    }
  }

  Graph::Graph(Graph&& other) :
    _num_nodes(other._num_nodes),
    _num_arcs(other._num_arcs),
    _total_volume(other._total_volume),
    _total_weight(other._total_weight),
    _max_degree(other._max_degree),
    _indices(std::move(other._indices)),
    _arcs(std::move(other._arcs)),
    _node_volumes(std::move(other._node_volumes)),
    _node_weights(std::move(other._node_weights)),
    _isolated_nodes_bitset(std::move(other._isolated_nodes_bitset)),
    _tmp_graph_buffer(other._tmp_graph_buffer) {
    other._num_nodes = 0;
    other._num_arcs = 0;
    other._total_volume = 0;
    other._total_weight = 0;
    other._max_degree = 0;
    other._tmp_graph_buffer = nullptr;
  }

  Graph& Graph::operator= (Graph&& other) {
    _num_nodes = other._num_nodes;
    _num_arcs = other._num_arcs;
    _total_volume = other._total_volume;
    _total_weight = other._total_weight;
    _max_degree = other._max_degree;
    _indices = std::move(other._indices);
    _arcs = std::move(other._arcs);
    _node_volumes = std::move(other._node_volumes);
    _node_weights = std::move(other._node_weights);
    _isolated_nodes_bitset = std::move(other._isolated_nodes_bitset);
    _tmp_graph_buffer = std::move(other._tmp_graph_buffer);
    other._num_nodes = 0;
    other._num_arcs = 0;
    other._total_volume = 0;
    other._total_weight = 0;
    other._max_degree = 0;
    other._tmp_graph_buffer = nullptr;
    return *this;
  }

  Graph::~Graph() {
    if ( _tmp_graph_buffer ) {
      delete(_tmp_graph_buffer);
    }
  }

  Graph Graph::contract_low_memory(Clustering& communities) {
    // map cluster IDs to consecutive range
    vec<NodeID> mapping(numNodes(), 0);   // TODO use memory pool?
    tbb::parallel_for(0UL, numNodes(), [&](NodeID u) { mapping[communities[u]] = 1; });
    parallel_prefix_sum(mapping.begin(), mapping.begin() + numNodes(), mapping.begin(), std::plus<>(), 0);
    NodeID num_coarse_nodes = mapping[numNodes() - 1];
    // apply mapping to cluster IDs. subtract one because prefix sum is inclusive
    tbb::parallel_for(0UL, numNodes(), [&](NodeID u) { communities[u] = mapping[communities[u]] - 1; });

    // sort nodes by cluster
    auto get_cluster = [&](NodeID u) { assert(u < communities.size()); return communities[u]; };
    vec<NodeID> nodes_sorted_by_cluster(std::move(mapping));    // reuse memory from mapping since it's no longer needed
    auto cluster_bounds = parallel::counting_sort(nodes(), nodes_sorted_by_cluster, num_coarse_nodes,
                                                  get_cluster, TBBInitializer::instance().total_number_of_threads());

    Graph coarse_graph;
    coarse_graph._num_nodes = num_coarse_nodes;
    coarse_graph._indices.resize(num_coarse_nodes + 1);
    coarse_graph._node_volumes.resize(num_coarse_nodes);
    coarse_graph._total_volume = totalVolume();
    coarse_graph._total_weight = totalWeight();

    struct ClearList {
      vec<NodeID> used;
      vec<ArcWeight> values;

      ClearList(size_t n) : values(n, 0.0) { }
    };
    tbb::enumerable_thread_specific<ClearList> clear_lists(num_coarse_nodes);
    tbb::enumerable_thread_specific<size_t> local_max_degree(0);

    // first pass generating unique coarse arcs to determine coarse node degrees
    tbb::parallel_for(0U, num_coarse_nodes, [&](NodeID cu) {
      auto& clear_list = clear_lists.local();
      ArcWeight volume_cu = 0.0;
      for (auto i = cluster_bounds[cu]; i < cluster_bounds[cu + 1]; ++i) {
        NodeID fu = nodes_sorted_by_cluster[i];
        volume_cu += nodeVolume(fu);
        for (const Arc& arc : arcsOf(fu)) {
          NodeID cv = get_cluster(arc.head);
          if (cv != cu && clear_list.values[cv] == 0.0) {
            clear_list.used.push_back(cv);
            clear_list.values[cv] = 1.0;
          }
        }
      }
      coarse_graph._indices[cu + 1] = clear_list.used.size();
      local_max_degree.local() = std::max(local_max_degree.local(), clear_list.used.size());
      for (const NodeID cv : clear_list.used) {
        clear_list.values[cv] = 0.0;
      }
      clear_list.used.clear();
      coarse_graph._node_volumes[cu] = volume_cu;
    });

    // prefix sum coarse node degrees for offsets to write the coarse arcs in second pass
    parallel_prefix_sum(coarse_graph._indices.begin(), coarse_graph._indices.end(), coarse_graph._indices.begin(), std::plus<>(), 0UL);
    size_t num_coarse_arcs = coarse_graph._indices.back();
    coarse_graph._arcs.resize(num_coarse_arcs);
    coarse_graph._num_arcs = num_coarse_arcs;
    coarse_graph._max_degree = local_max_degree.combine([](size_t lhs, size_t rhs) { return std::max(lhs, rhs); });

    // second pass generating unique coarse arcs
    tbb::parallel_for(0U, num_coarse_nodes, [&](NodeID cu) {
      auto& clear_list = clear_lists.local();
      for (auto i = cluster_bounds[cu]; i < cluster_bounds[cu+1]; ++i) {
        for (const Arc& arc : arcsOf(nodes_sorted_by_cluster[i])) {
          NodeID cv = get_cluster(arc.head);
          if (cv != cu) {
            if (clear_list.values[cv] == 0.0) {
              clear_list.used.push_back(cv);
            }
            clear_list.values[cv] += arc.weight;
          }
        }
      }
      size_t pos = coarse_graph._indices[cu];
      for (const NodeID cv : clear_list.used) {
        coarse_graph._arcs[pos++] = Arc(cv, clear_list.values[cv]);
        clear_list.values[cv] = 0.0;
      }
      clear_list.used.clear();
    });

    return coarse_graph;
  }


  /*!
 * Contracts the graph based on the community structure passed as argument.
 * In the first step the community ids are compactified (via parallel prefix sum)
 * which also determines the node ids in the coarse graph. Afterwards, we create
 * a temporary graph which contains all arcs that will not form a selfloop in the
 * coarse graph. Finally, the weights of each multiedge in that temporary graph
 * are aggregated and the result is written to the final contracted graph.
 */
  Graph Graph::contract(Clustering& communities, bool low_memory) {
    if (low_memory) {
      LOG << "LOW MEMORY!!!";
      return contract_low_memory(communities);
    }
    ASSERT(canBeUsed());
    ASSERT(_num_nodes == communities.size());
    if ( !_tmp_graph_buffer ) {
      allocateContractionBuffers();
    }
    Graph coarse_graph;
    coarse_graph._total_volume = _total_volume;
    coarse_graph._total_weight = _total_weight;

    // #################### STAGE 1 ####################
    // Compute node ids of coarse graph with a parallel prefix sum
    parallel::scalable_vector<size_t> mapping(_num_nodes, 0UL);
    ds::Array<parallel::IntegralAtomicWrapper<size_t>>& tmp_pos = _tmp_graph_buffer->tmp_pos;
    ds::Array<parallel::IntegralAtomicWrapper<size_t>>& tmp_indices = _tmp_graph_buffer->tmp_indices;
    ds::Array<parallel::AtomicWrapper<ArcWeight>>& coarse_node_volumes = _tmp_graph_buffer->tmp_node_volumes;
    ds::Array<parallel::AtomicWrapper<HypernodeWeight>>& coarse_node_weights = _tmp_graph_buffer->tmp_node_weights;
    ds::Array<char>& tmp_isolated_nodes = _tmp_graph_buffer->tmp_isolated;
    tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
      ASSERT(static_cast<size_t>(communities[u]) < _num_nodes);
      mapping[communities[u]] = 1UL;
      tmp_pos[u] = 0;
      tmp_indices[u] = 0;
      tmp_isolated_nodes[u] = 0;
      coarse_node_volumes[u].store(0.0);
      coarse_node_weights[u].store(0);
    });

    // Prefix sum determines node ids in coarse graph
    parallel::TBBPrefixSum<size_t> mapping_prefix_sum(mapping);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, _num_nodes), mapping_prefix_sum);

    // Remap community ids
    coarse_graph._num_nodes = mapping_prefix_sum.total_sum();
    tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
      communities[u] = mapping_prefix_sum[communities[u]];
    });

    // #################### STAGE 2 ####################
    // Write all arcs, that will not form a selfloop in the coarse graph, into a tmp
    // adjacence array. For that, we compute a prefix sum over the sum of all arcs
    // in each community (which are no selfloop) and write them in parallel to
    // the tmp adjacence array.
    // Compute number of arcs in tmp adjacence array with parallel prefix sum
    ASSERT(coarse_graph._num_nodes <= coarse_node_volumes.size());
    ASSERT(coarse_graph._num_nodes <= coarse_node_weights.size());
    tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
      const NodeID coarse_u = communities[u];
      ASSERT(static_cast<size_t>(coarse_u) < coarse_graph._num_nodes);
      coarse_node_volumes[coarse_u] += nodeVolume(u);     // not deterministic!
      coarse_node_weights[coarse_u] += nodeWeight(u);
      if (_isolated_nodes_bitset[u] != 0) {
        tmp_isolated_nodes[coarse_u] = 1;
      }
      for ( const Arc& arc : arcsOf(u) ) {
        const NodeID coarse_v = communities[arc.head];
        if ( coarse_u != coarse_v ) {
          ++tmp_indices[coarse_u];
        }
      }
    });

    parallel::TBBPrefixSum<parallel::IntegralAtomicWrapper<size_t>, ds::Array> tmp_indices_prefix_sum(tmp_indices);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, _num_nodes), tmp_indices_prefix_sum);

    // Write all arcs into corresponding tmp adjacence array blocks
    ds::Array<Arc>& tmp_arcs = _tmp_graph_buffer->tmp_arcs;
    ds::Array<size_t>& valid_arcs = _tmp_graph_buffer->valid_arcs;
    tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
      const NodeID coarse_u = communities[u];
      ASSERT(static_cast<size_t>(coarse_u) < coarse_graph._num_nodes);
      for ( const Arc& arc : arcsOf(u) ) {
        const NodeID coarse_v = communities[arc.head];
        if ( coarse_u != coarse_v ) {
          const size_t tmp_arcs_pos = tmp_indices_prefix_sum[coarse_u] + tmp_pos[coarse_u]++;
          ASSERT(tmp_arcs_pos < tmp_indices_prefix_sum[coarse_u + 1]);
          tmp_arcs[tmp_arcs_pos] = Arc { coarse_v, arc.weight };
          valid_arcs[tmp_arcs_pos] = 1UL;
        }
      }
    });

    // #################### STAGE 3 ####################
    // Aggregate weights of arcs that are equal in each community.
    // Therefore, we sort the arcs according to their endpoints
    // and aggregate weight of arcs with equal endpoints.
    tbb::enumerable_thread_specific<size_t> local_max_degree(0);
    tbb::parallel_for(0U, static_cast<NodeID>(coarse_graph._num_nodes), [&](const NodeID u) {
      const size_t tmp_arc_start = tmp_indices_prefix_sum[u];
      const size_t tmp_arc_end = tmp_indices_prefix_sum[u + 1];
      // commented out comparison is needed for deterministic arc weights
      // auto comp = [](const Arc& lhs, const Arc& rhs) { return std::tie(lhs.head, lhs.weight) < std::tie(rhs.head, rhs.weight); };
      auto comp = [](const Arc& lhs, const Arc& rhs) { return lhs.head < rhs.head; };
      std::sort(tmp_arcs.begin() + tmp_arc_start, tmp_arcs.begin() + tmp_arc_end, comp);

      size_t arc_rep = tmp_arc_start;
      size_t degree = tmp_arc_start < tmp_arc_end ? 1 : 0;
      for ( size_t pos = tmp_arc_start + 1; pos < tmp_arc_end; ++pos ) {
        if ( tmp_arcs[arc_rep].head == tmp_arcs[pos].head ) {
          tmp_arcs[arc_rep].weight += tmp_arcs[pos].weight;
          valid_arcs[pos] = 0UL;
        } else {
          arc_rep = pos;
          ++degree;
        }
      }
      local_max_degree.local() = std::max(local_max_degree.local(), degree);
    });
    coarse_graph._max_degree = local_max_degree.combine(
            [&](const size_t& lhs, const size_t& rhs) {
              return std::max(lhs, rhs);
            });

    // Write all arcs to coarse graph
    parallel::TBBPrefixSum<size_t, ds::Array> valid_arcs_prefix_sum(valid_arcs);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL,
                                                  tmp_indices_prefix_sum.total_sum()), valid_arcs_prefix_sum);
    coarse_graph._num_arcs = valid_arcs_prefix_sum.total_sum();

    // Move memory down to coarse graph
    coarse_graph._indices = std::move(_indices);
    coarse_graph._arcs = std::move(_arcs);
    coarse_graph._node_volumes = std::move(_node_volumes);
    coarse_graph._node_weights = std::move(_node_weights);
    coarse_graph._isolated_nodes_bitset = std::move(_isolated_nodes_bitset);

    tbb::parallel_invoke([&] {
      const size_t tmp_num_arcs = tmp_indices_prefix_sum.total_sum();
      tbb::parallel_for(0UL, tmp_num_arcs, [&](const size_t i) {
        if ( valid_arcs_prefix_sum.value(i) ) {
          const size_t pos = valid_arcs_prefix_sum[i];
          ASSERT(pos < coarse_graph._num_arcs);
          coarse_graph._arcs[pos] = tmp_arcs[i];
        }
      });
    }, [&] {
      tbb::parallel_for(0U, static_cast<NodeID>(coarse_graph._num_nodes), [&](const NodeID u) {
        const size_t start_index_pos = valid_arcs_prefix_sum[tmp_indices_prefix_sum[u]];
        ASSERT(start_index_pos <= coarse_graph._num_arcs);
        coarse_graph._indices[u] = start_index_pos;
        coarse_graph._node_volumes[u] = coarse_node_volumes[u];
        coarse_graph._node_weights[u] = coarse_node_weights[u];
        coarse_graph._isolated_nodes_bitset[u] = tmp_isolated_nodes[u];
      });
      coarse_graph._indices[coarse_graph._num_nodes] = coarse_graph._num_arcs;
    });
    coarse_graph._tmp_graph_buffer = _tmp_graph_buffer;
    _tmp_graph_buffer = nullptr;

    return coarse_graph;
  }


  Graph::Graph() :
          _num_nodes(0),
          _num_arcs(0),
          _total_volume(0),
          _total_weight(0),
          _max_degree(0),
          _indices(),
          _arcs(),
          _node_volumes(),
          _node_weights(),
          _isolated_nodes_bitset(),
          _tmp_graph_buffer(nullptr) {

  }




  /*!
   * Constructs a graph from a given hypergraph.
   */
  template<typename F>
  void Graph::construct(const Hypergraph& hypergraph,
                        const bool is_graph,
                        const F& edge_weight_func) {
    if ( is_graph ) {
      ASSERT(hypergraph.maxEdgeSize() == 2);
      _num_nodes = hypergraph.initialNumNodes();
      _num_arcs = hypergraph.initialNumPins();
      constructGraph(hypergraph, edge_weight_func);
    } else {
      _num_nodes = hypergraph.initialNumNodes() + hypergraph.initialNumEdges();
      _num_arcs = 2 * hypergraph.initialNumPins();
      constructBipartiteGraph(hypergraph, edge_weight_func);
    }

    // deterministic reduce of node volumes since double addition is not commutative or associative
    // node volumes are computed in for loop because deterministic reduce does not have dynamic load balancing
    // whereas for loop does. this important since each node incurs O(degree) time
    tbb::parallel_for(0U, NodeID(numNodes()), [&](NodeID u) { computeNodeVolume(u); });

    auto aggregate_volume = [&](const tbb::blocked_range<NodeID>& r, ArcWeight partial_volume) -> ArcWeight {
      for (NodeID u = r.begin(); u < r.end(); ++u) {
        partial_volume += nodeVolume(u);
      }
      return partial_volume;
    };
    auto r = tbb::blocked_range<NodeID>(0U, numNodes(), 1000);
    _total_volume = tbb::parallel_deterministic_reduce(r, 0.0, aggregate_volume, std::plus<>());
    _total_weight = hypergraph.totalWeight();
  }

  template<typename F>
  void Graph::constructBipartiteGraph(const Hypergraph& hypergraph,
                               F& edge_weight_func) {
    _indices.resize("Preprocessing", "indices", _num_nodes + 1);
    _arcs.resize("Preprocessing", "arcs", _num_arcs);
    _node_volumes.resize("Preprocessing", "node_volumes", _num_nodes);
    _isolated_nodes_bitset.resize("Preprocessing", "isolated_nodes", _num_nodes);

    // Initialize data structure
    const HypernodeID num_hypernodes = hypergraph.initialNumNodes();
    const HypernodeID num_hyperedges = hypergraph.initialNumEdges();
    tbb::parallel_invoke([&] {
      tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID u) {
        ASSERT(u + 1 < _indices.size());
        _indices[u + 1] = hypergraph.nodeDegree(u);
      });
    }, [&] {
      tbb::parallel_for(num_hypernodes, num_hypernodes + num_hyperedges, [&](const HyperedgeID u) {
        ASSERT(u + 1 < _indices.size());
        const HyperedgeID he = u - num_hypernodes;
        _indices[u + 1] = hypergraph.edgeSize(he);
      });
    });

    parallel::TBBPrefixSum<size_t, ds::Array> indices_prefix_sum(_indices);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, _indices.size()), indices_prefix_sum);

    tbb::enumerable_thread_specific<size_t> local_max_degree(0);
    tbb::parallel_invoke([&] {
      tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID u) {
        ASSERT(u + 1 < _indices.size());
        size_t pos = _indices[u];
        const HypernodeID hn = u;
        const HyperedgeID node_degree = hypergraph.nodeDegree(hn);
        local_max_degree.local() = std::max(
                local_max_degree.local(), static_cast<size_t>(node_degree));
        for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
          const NodeID v = he + num_hypernodes;
          const HyperedgeWeight edge_weight = hypergraph.edgeWeight(he);
          const HypernodeID edge_size = hypergraph.edgeSize(he);
          ASSERT(pos < _indices[u + 1]);
          _arcs[pos++] = Arc(v, edge_weight_func(edge_weight, edge_size, node_degree));
        }
      });
    }, [&] {
      tbb::parallel_for(num_hypernodes, num_hypernodes + num_hyperedges, [&](const HyperedgeID u) {
        ASSERT(u + 1 < _indices.size());
        size_t pos = _indices[u];
        const HyperedgeID he = u - num_hypernodes;
        const HyperedgeWeight edge_weight = hypergraph.edgeWeight(he);
        const HypernodeID edge_size = hypergraph.edgeSize(he);
        local_max_degree.local() = std::max(
                local_max_degree.local(), static_cast<size_t>(edge_size));
        for ( const HypernodeID& pin : hypergraph.pins(he) ) {
          const NodeID v = pin;
          const HyperedgeID node_degree = hypergraph.nodeDegree(pin);
          ASSERT(pos < _indices[u + 1]);
          _arcs[pos++] = Arc(v, edge_weight_func(edge_weight, edge_size, node_degree));
        }
      });
    });
    _max_degree = local_max_degree.combine([&](const size_t& lhs, const size_t& rhs) {
      return std::max(lhs, rhs);
    });
  }

  template<typename F>
  void Graph::constructGraph(const Hypergraph& hypergraph, const F& edge_weight_func) {
    _indices.resize("Preprocessing", "indices", _num_nodes + 1);
    _arcs.resize("Preprocessing", "arcs", _num_arcs);
    _node_volumes.resize("Preprocessing", "node_volumes", _num_nodes);
    _node_weights.resize("Preprocessing", "node_weights", _num_nodes);
    _isolated_nodes_bitset.resize("Preprocessing", "isolated_nodes", _num_nodes);

    // Initialize data structure
    const HypernodeID num_hypernodes = hypergraph.initialNumNodes();
    tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID u) {
      ASSERT(u + 1 < _indices.size());
      _indices[u + 1] = hypergraph.nodeDegree(u);
      _node_weights[u] = hypergraph.nodeWeight(u);
    });

    parallel::TBBPrefixSum<size_t, ds::Array> indices_prefix_sum(_indices);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, num_hypernodes + 1), indices_prefix_sum);

    tbb::enumerable_thread_specific<size_t> local_max_degree(0);
    tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID u) {
      ASSERT(u + 1 < _indices.size());
      size_t pos = _indices[u];
      const HypernodeID hn = u;
      const HyperedgeID node_degree = hypergraph.nodeDegree(hn);
      local_max_degree.local() = std::max(
              local_max_degree.local(), static_cast<size_t>(node_degree));
      for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
        const HyperedgeWeight edge_weight = hypergraph.edgeWeight(he);
        NodeID v = std::numeric_limits<NodeID>::max();
        for ( const HypernodeID& pin : hypergraph.pins(he) ) {
          if ( pin != hn ) {
            v = pin;
            break;
          }
        }
        ASSERT(v != std::numeric_limits<NodeID>::max());
        ASSERT(pos < _indices[u + 1]);
        _arcs[pos++] = Arc(v, edge_weight_func(edge_weight, ID(2), node_degree));
      }
    });
    _max_degree = local_max_degree.combine([&](const size_t& lhs, const size_t& rhs) {
      return std::max(lhs, rhs);
    });
  }

  bool Graph::canBeUsed(const bool verbose) const {
    const bool result = _indices.size() >= numNodes() + 1 && _arcs.size() >= numArcs() && _node_volumes.size() >= numNodes();
    if (verbose && !result) {
      LOG << "Some of the graph's members were stolen. For example the contract function does this. "
             "Make sure you're calling functions with a fresh graph or catch this condition and reinitialize."
             "If you do reinitialize, feel free to silence this warning by passing false to the canBeUsed function";
    }
    return result;
  }

} // namespace mt_kahypar::ds
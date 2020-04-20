/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include <functional>
#include <cmath>

#include <boost/range/irange.hpp>

#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_invoke.h>
#include <tbb/enumerable_thread_specific.h>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/clustering.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
namespace ds {

/*!
 * Scalable CSR Graph Data Structure
 * In our graph data structure the nodes are partitioned into equal-sized blocks.
 * For each block, we construct an adjacence array separately. Idea behind this is that,
 * expensive allocations are scattered among several blocks during construction and
 * contraction and can be therefore implemented in a scalable way. Furthermore, it still
 * has the advantages of a traditional adjacence array graph representation in terms
 * of cache-efficiency and memory-consumption.
 */
template<typename HyperGraph>
class GraphT {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  using TmpGraphBuffer = typename HyperGraph::TmpGraphBuffer;

 public:
  using AdjacenceIterator = typename parallel::scalable_vector<Arc>::const_iterator;

 public:
  explicit GraphT(HyperGraph& hypergraph,
                  const LouvainEdgeWeight edge_weight_type) :
    _num_nodes(0),
    _num_arcs(0),
    _total_volume(0),
    _max_degree(0),
    _indices(hypergraph.tmpGraphBuffer()->indices),
    _arcs(hypergraph.tmpGraphBuffer()->arcs),
    _node_volumes(hypergraph.tmpGraphBuffer()->node_volumes),
    _tmp_graph_buffer(hypergraph.tmpGraphBuffer()) {
    ASSERT(_tmp_graph_buffer && _tmp_graph_buffer->is_initialized);

    switch( edge_weight_type ) {
      case LouvainEdgeWeight::uniform:
        construct(hypergraph,
          [&](const HyperedgeWeight edge_weight,
              const HypernodeID,
              const HyperedgeID) {
              return static_cast<ArcWeight>(edge_weight);
            });
        break;
      case LouvainEdgeWeight::non_uniform:
        construct(hypergraph,
          [&](const HyperedgeWeight edge_weight,
              const HypernodeID edge_size,
              const HyperedgeID) {
              return static_cast<ArcWeight>(edge_weight) /
                static_cast<ArcWeight>(edge_size);
            });
        break;
      case LouvainEdgeWeight::degree:
        construct(hypergraph,
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

  // ! Number of nodes in the graph
  size_t numNodes() const {
    return _num_nodes;
  }

  // ! Number of arcs in the graph
  size_t numArcs() const {
    return _num_arcs;
  }

  // ! Iterator over all nodes of the graph
  auto nodes() const {
    return boost::irange<NodeID>(0, static_cast<NodeID>(numNodes()));
  }

    // ! Iterator over all adjacent vertices of u
  IteratorRange<AdjacenceIterator> arcsOf(const NodeID u) const {
    ASSERT(u < _num_nodes);
    return IteratorRange<AdjacenceIterator>(
      _arcs.cbegin() + _indices[u], _arcs.cbegin() + _indices[u + 1]);
  }

  // ! Degree of vertex u
  size_t degree(const NodeID u) const {
    ASSERT(u < _num_nodes);
    return _indices[u + 1] - _indices[u];
  }

  // ! Maximum degree of a vertex
  size_t max_degree() const {
    return _max_degree;
  }

  // ! Total Volume of the graph
  ArcWeight totalVolume() const {
    return _total_volume;
  }

  // ! Node volume of vertex u
  ArcWeight nodeVolume(const NodeID u) const {
    ASSERT(u < _num_nodes);
    return _node_volumes[u];
  }

  /*!
   * Contracts the graph based on the community structure passed as argument.
   * In the first step the community ids are compactified (via parallel prefix sum)
   * which also determines the node ids in the coarse graph. Afterwards, we create
   * a temporary graph which contains all arcs that will not form a selfloop in the
   * coarse graph. Finally, the weights of each multiedge in that temporary graph
   * are aggregated and the result is written to the final contracted graph.
   */
  GraphT contract(Clustering& communities) {
    ASSERT(_num_nodes == communities.size());
    ASSERT(_tmp_graph_buffer && _tmp_graph_buffer->is_initialized);
    GraphT coarse_graph(_tmp_graph_buffer);
    coarse_graph._total_volume = _total_volume;

    // #################### STAGE 1 ####################
    // Compute node ids of coarse graph with a parallel prefix sum
    utils::Timer::instance().start_timer("compute_cluster_mapping", "Compute Cluster Mapping");
    parallel::scalable_vector<size_t> mapping(_num_nodes, 0UL);
    parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>>& tmp_pos = _tmp_graph_buffer->tmp_pos;
    parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>>& tmp_indices = _tmp_graph_buffer->tmp_indices;
    parallel::scalable_vector<parallel::AtomicWrapper<ArcWeight>>& coarse_node_volumes = _tmp_graph_buffer->tmp_node_volumes;
    tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
      ASSERT(static_cast<size_t>(communities[u]) < _num_nodes);
      mapping[communities[u]] = 1UL;
      tmp_pos[u] = 0;
      tmp_indices[u] = 0;
      coarse_node_volumes[u].store(0.0);
    });

    // Prefix sum determines node ids in coarse graph
    parallel::TBBPrefixSum<size_t> mapping_prefix_sum(mapping);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, _num_nodes), mapping_prefix_sum);

    // Remap community ids
    coarse_graph._num_nodes = mapping_prefix_sum.total_sum();
    tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
      communities[u] = mapping_prefix_sum[communities[u]];
    });
    utils::Timer::instance().stop_timer("compute_cluster_mapping");

    // #################### STAGE 2 ####################
    // Write all arcs, that will not form a selfloop in the coarse graph, into a tmp
    // adjacence array. For that, we compute a prefix sum over the sum of all arcs
    // in each community (which are no selfloop) and write them in parallel to
    // the tmp adjacence array.
    utils::Timer::instance().start_timer("construct_tmp_adjacent_array", "Construct Tmp Adjacent Array");
    // Compute number of arcs in tmp adjacence array with parallel prefix sum
    ASSERT(coarse_graph._num_nodes <= coarse_node_volumes.size());
    tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
      const NodeID coarse_u = communities[u];
      ASSERT(static_cast<size_t>(coarse_u) < coarse_graph._num_nodes);
      coarse_node_volumes[coarse_u] += nodeVolume(u);
      for ( const Arc& arc : arcsOf(u) ) {
        const NodeID coarse_v = communities[arc.head];
        if ( coarse_u != coarse_v ) {
          ++tmp_indices[coarse_u];
        }
      }
    });

    parallel::TBBPrefixSum<parallel::IntegralAtomicWrapper<size_t>> tmp_indices_prefix_sum(tmp_indices);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, _num_nodes), tmp_indices_prefix_sum);

    // Write all arcs into corresponding tmp adjacence array blocks
    parallel::scalable_vector<Arc>& tmp_arcs = _tmp_graph_buffer->tmp_arcs;
    parallel::scalable_vector<size_t>& valid_arcs = _tmp_graph_buffer->valid_arcs;
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
    utils::Timer::instance().stop_timer("construct_tmp_adjacent_array");

    // #################### STAGE 3 ####################
    // Aggregate weights of arcs that are equal in each community.
    // Therefore, we sort the arcs according to their endpoints
    // and aggregate weight of arcs with equal endpoints.
    utils::Timer::instance().start_timer("contract_arcs", "Contract Arcs");
    tbb::enumerable_thread_specific<size_t> local_max_degree(0);
    tbb::parallel_for(0U, static_cast<NodeID>(coarse_graph._num_nodes), [&](const NodeID u) {
      const size_t tmp_arc_start = tmp_indices_prefix_sum[u];
      const size_t tmp_arc_end = tmp_indices_prefix_sum[u + 1];
      std::sort(tmp_arcs.begin() + tmp_arc_start, tmp_arcs.begin() + tmp_arc_end,
        [&](const Arc& lhs, const Arc& rhs) {
          return lhs.head < rhs.head;
        });

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
    parallel::TBBPrefixSum<size_t> valid_arcs_prefix_sum(valid_arcs);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL,
      tmp_indices_prefix_sum.total_sum()), valid_arcs_prefix_sum);
    coarse_graph._num_arcs = valid_arcs_prefix_sum.total_sum();

    tbb::parallel_invoke([&] {
      const size_t tmp_num_arcs = tmp_indices_prefix_sum.total_sum();
      tbb::parallel_for(0UL, tmp_num_arcs, [&](const size_t i) {
        if ( valid_arcs_prefix_sum.value(i) ) {
          const size_t pos = valid_arcs_prefix_sum[i];
          ASSERT(pos < coarse_graph._num_arcs);
          coarse_graph._arcs[pos] = std::move(tmp_arcs[i]);
        }
      });
    }, [&] {
      tbb::parallel_for(0U, static_cast<NodeID>(coarse_graph._num_nodes), [&](const NodeID u) {
        const size_t start_index_pos = valid_arcs_prefix_sum[tmp_indices_prefix_sum[u]];
        ASSERT(start_index_pos <= coarse_graph._num_arcs);
        coarse_graph._indices[u] = start_index_pos;
        coarse_graph._node_volumes[u] = coarse_node_volumes[u];
      });
      coarse_graph._indices[coarse_graph._num_nodes] = coarse_graph._num_arcs;
    });
    coarse_graph._tmp_graph_buffer = _tmp_graph_buffer;
    utils::Timer::instance().stop_timer("contract_arcs");

    return coarse_graph;
  }

 private:
  GraphT(TmpGraphBuffer* tmp_graph_buffer) :
    _num_nodes(0),
    _num_arcs(0),
    _total_volume(0),
    _max_degree(0),
    _indices(tmp_graph_buffer->indices),
    _arcs(tmp_graph_buffer->arcs),
    _node_volumes(tmp_graph_buffer->node_volumes),
    _tmp_graph_buffer(tmp_graph_buffer) { }

  /*!
   * Constructs a graph from a given hypergraph.
   */
  template<typename F>
  void construct(const HyperGraph& hypergraph,
                 const F& edge_weight_func) {
    // Test, if hypergraph is actually a graph
    const bool is_graph = tbb::parallel_reduce(tbb::blocked_range<HyperedgeID>(
      ID(0), hypergraph.initialNumEdges()), true, [&](const tbb::blocked_range<HyperedgeID>& range, bool isGraph) {
        if ( isGraph ) {
          bool tmp_is_graph = isGraph;
          for (HyperedgeID id = range.begin(); id < range.end(); ++id) {
            const HyperedgeID he = hypergraph.globalEdgeID(id);
            if ( hypergraph.edgeIsEnabled(he) ) {
              tmp_is_graph &= (hypergraph.edgeSize(he) == 2);
            }
          }
          return tmp_is_graph;
        }
        return false;
      }, [&](const bool lhs, const bool rhs) {
        return lhs && rhs;
      });

    if ( is_graph ) {
      _num_nodes = hypergraph.initialNumNodes();
      _num_arcs = 2 * hypergraph.initialNumEdges();
      constructGraph(hypergraph, edge_weight_func);
    } else {
      _num_nodes = hypergraph.initialNumNodes() + hypergraph.initialNumEdges();
      _num_arcs = 2 * hypergraph.initialNumPins();
      constructBipartiteGraph(hypergraph, edge_weight_func);
    }

    // Compute node volumes and total volume
    utils::Timer::instance().start_timer("compute_node_volumes", "Compute Node Volumes");
    _total_volume = 0.0;
    tbb::enumerable_thread_specific<ArcWeight> local_total_volume(0.0);
    tbb::parallel_for(0U, static_cast<NodeID>(numNodes()), [&](const NodeID u) {
      local_total_volume.local() += computeNodeVolume(u);
    });
    _total_volume = local_total_volume.combine(std::plus<ArcWeight>());
    utils::Timer::instance().stop_timer("compute_node_volumes");
  }

  template<typename F>
  void constructBipartiteGraph(const HyperGraph& hypergraph,
                               F& edge_weight_func) {
    // Initialize data structure
    utils::Timer::instance().start_timer("compute_node_degrees", "Compute Node Degrees");
    const HypernodeID num_hypernodes = hypergraph.initialNumNodes();
    const HypernodeID num_hyperedges = hypergraph.initialNumEdges();
    tbb::parallel_invoke([&] {
      tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID u) {
        ASSERT(u + 1 < _indices.size());
        const HypernodeID hn = hypergraph.globalNodeID(u);
        _indices[u + 1] = hypergraph.nodeDegree(hn);
      });
    }, [&] {
      tbb::parallel_for(num_hypernodes, num_hypernodes + num_hyperedges, [&](const HyperedgeID u) {
        ASSERT(u + 1 < _indices.size());
        const HyperedgeID he = hypergraph.globalEdgeID(u - num_hypernodes);
        _indices[u + 1] = hypergraph.edgeSize(he);
      });
    });

    parallel::TBBPrefixSum<size_t> indices_prefix_sum(_indices);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, _indices.size()), indices_prefix_sum);
    utils::Timer::instance().stop_timer("compute_node_degrees");

    utils::Timer::instance().start_timer("construct_arcs", "Construct Arcs");
    tbb::enumerable_thread_specific<size_t> local_max_degree(0);
    tbb::parallel_invoke([&] {
      tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID u) {
        ASSERT(u + 1 < _indices.size());
        size_t pos = _indices[u];
        const HypernodeID hn = hypergraph.globalNodeID(u);
        const HyperedgeID node_degree = hypergraph.nodeDegree(hn);
        local_max_degree.local() = std::max(
          local_max_degree.local(), static_cast<size_t>(node_degree));
        for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
          const NodeID v = hypergraph.originalEdgeID(he) + num_hypernodes;
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
        const HyperedgeID he = hypergraph.globalEdgeID(u - num_hypernodes);
        const HyperedgeWeight edge_weight = hypergraph.edgeWeight(he);
        const HypernodeID edge_size = hypergraph.edgeSize(he);
        local_max_degree.local() = std::max(
          local_max_degree.local(), static_cast<size_t>(edge_size));
        for ( const HypernodeID& pin : hypergraph.pins(he) ) {
          const NodeID v = hypergraph.originalNodeID(pin);
          const HyperedgeID node_degree = hypergraph.nodeDegree(pin);
          ASSERT(pos < _indices[u + 1]);
          _arcs[pos++] = Arc(v, edge_weight_func(edge_weight, edge_size, node_degree));
        }
      });
    });
    _max_degree = local_max_degree.combine([&](const size_t& lhs, const size_t& rhs) {
      return std::max(lhs, rhs);
    });
    utils::Timer::instance().stop_timer("construct_arcs");
  }

  template<typename F>
  void constructGraph(const HyperGraph& hypergraph,
                      const F& edge_weight_func) {
    // Initialize data structure
    utils::Timer::instance().start_timer("compute_node_degrees", "Compute Node Degrees");
    const HypernodeID num_hypernodes = hypergraph.initialNumNodes();
    tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID u) {
      ASSERT(u + 1 < _indices.size());
      const HypernodeID hn = hypergraph.globalNodeID(u);
      _indices[u + 1] = hypergraph.nodeDegree(hn);
    });

    parallel::TBBPrefixSum<size_t> indices_prefix_sum(_indices);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, num_hypernodes + 1), indices_prefix_sum);
    utils::Timer::instance().stop_timer("compute_node_degrees");

    utils::Timer::instance().start_timer("construct_arcs", "Construct Arcs");
    tbb::enumerable_thread_specific<size_t> local_max_degree(0);
    tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID u) {
      ASSERT(u + 1 < _indices.size());
      size_t pos = _indices[u];
      const HypernodeID hn = hypergraph.globalNodeID(u);
      const HyperedgeID node_degree = hypergraph.nodeDegree(hn);
      local_max_degree.local() = std::max(
        local_max_degree.local(), static_cast<size_t>(node_degree));
      for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
        const HyperedgeWeight edge_weight = hypergraph.edgeWeight(he);
        NodeID v = std::numeric_limits<NodeID>::max();
        for ( const HypernodeID& pin : hypergraph.pins(he) ) {
          if ( pin != hn ) {
            v = hypergraph.originalNodeID(pin);
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
    utils::Timer::instance().stop_timer("construct_arcs");
  }

  ArcWeight computeNodeVolume(const NodeID u) {
    ASSERT(u < _num_nodes);
    for (const Arc& arc : arcsOf(u)) {
      _node_volumes[u] += arc.weight;
    }
    return _node_volumes[u];
  }

  // ! Number of nodes
  size_t _num_nodes;
  // ! Number of arcs
  size_t _num_arcs;
  // ! Total volume of the graph (= sum of arc weights)
  ArcWeight _total_volume;
  // ! Maximum degree of a node
  size_t _max_degree;

  // ! Index Vector
  parallel::scalable_vector<size_t>& _indices;
  // ! Arcs
  parallel::scalable_vector<Arc>& _arcs;
  // ! Node Volumes (= sum of arc weights for each node)
  parallel::scalable_vector<ArcWeight>& _node_volumes;
  // ! Data that is reused throughout the louvain method
  // ! to construct and contract a graph and to prevent expensive allocations
  TmpGraphBuffer* _tmp_graph_buffer;
};

}  // namespace ds
}  // namespace mt_kahypar

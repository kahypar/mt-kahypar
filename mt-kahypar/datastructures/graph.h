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

#include <boost/range/irange.hpp>

#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>
#include <tbb/enumerable_thread_specific.h>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/clustering.h"
#include "mt-kahypar/datastructures/streaming_map.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
namespace ds {

template<typename HyperGraph>
class GraphT {

 public:
  using ArcWeight = double;

  struct Arc {
    NodeID head;
    ArcWeight weight;

    Arc() :
      head(0),
      weight(0) { }

    Arc(NodeID head, ArcWeight weight) :
      head(head),
      weight(weight) { }
  };

  using AdjacenceIterator = typename parallel::scalable_vector<Arc>::const_iterator;

  explicit GraphT(const HyperGraph& hypergraph,
                  const LouvainEdgeWeight edge_weight_type,
                  const TaskGroupID task_group_id) :
    _num_nodes(hypergraph.initialNumNodes() + hypergraph.initialNumEdges()),
    _num_arcs(2 * hypergraph.initialNumPins()),
    _total_volume(0),
    _indices(),
    _arcs(),
    _node_volumes() {
    switch( edge_weight_type ) {
      case LouvainEdgeWeight::uniform:
        construct(hypergraph, task_group_id,
          [&](const HyperedgeWeight edge_weight,
              const HypernodeID,
              const HyperedgeID) {
              return static_cast<ArcWeight>(edge_weight);
            });
        break;
      case LouvainEdgeWeight::non_uniform:
        construct(hypergraph, task_group_id,
          [&](const HyperedgeWeight edge_weight,
              const HypernodeID edge_size,
              const HyperedgeID) {
              return static_cast<ArcWeight>(edge_weight) /
                static_cast<ArcWeight>(edge_size);
            });
        break;
      case LouvainEdgeWeight::degree:
        construct(hypergraph, task_group_id,
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

  ~GraphT() {
    parallel::parallel_free(_indices, _arcs);
  }

  size_t numNodes() const {
    return _num_nodes;
  }

  size_t numArcs() const {
    return _num_arcs;
  }

  auto nodes() const {
    return boost::irange<NodeID>(0, static_cast<NodeID>(numNodes()));
  }

  IteratorRange<AdjacenceIterator> arcsOf(const NodeID u) const {
    ASSERT(u < _num_nodes);
    return IteratorRange<AdjacenceIterator>(
      _arcs.cbegin() + _indices[u],
      _arcs.cbegin() + _indices[u + 1]);
  }

  size_t degree(const NodeID u) const {
    ASSERT(u < _num_nodes);
    return _indices[u + 1] - _indices[u];
  }

  ArcWeight totalVolume() const {
    return _total_volume;
  }

  ArcWeight nodeVolume(const NodeID u) const {
    ASSERT(u < _num_nodes);
    return _node_volumes[u];
  }

  GraphT contract(Clustering& communities) {
    ASSERT(_num_nodes == communities.size());
    GraphT coarse_graph;
    coarse_graph._total_volume = _total_volume;

    // #################### STAGE 1 ####################
    // Compute node ids of coarse graph with a parallel prefix sum
    utils::Timer::instance().start_timer("compute_cluster_mapping", "Compute Cluster Mapping");
    parallel::scalable_vector<size_t> mapping(_num_nodes, 0);
    tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
      ASSERT(static_cast<size_t>(communities[u]) < mapping.size());
      mapping[communities[u]] = 1UL;
    });

    // Prefix sum determines node ids in coarse graph
    parallel::TBBPrefixSum<size_t> mapping_prefix_sum(mapping);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, mapping.size()), mapping_prefix_sum);

    // Remap community ids
    parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>> tmp_adjacent_nodes;
    parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>> tmp_pos;
    coarse_graph._num_nodes = mapping_prefix_sum.total_sum();
    tbb::parallel_invoke([&] {
      tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
        communities[u] = mapping_prefix_sum[communities[u]];
      });
    }, [&] {
      tmp_adjacent_nodes.assign(coarse_graph._num_nodes,
        parallel::IntegralAtomicWrapper<size_t>(0));
      tmp_pos.assign(coarse_graph._num_nodes,
        parallel::IntegralAtomicWrapper<size_t>(0));
      coarse_graph._indices.resize(coarse_graph._num_nodes + 1);
      coarse_graph._node_volumes.resize(coarse_graph._num_nodes);
    });
    utils::Timer::instance().stop_timer("compute_cluster_mapping");

    // #################### STAGE 2 ####################
    // Write all arcs, that are not a selfloop in the coarse graph, into a tmp
    // adjacence array. For that, we compute a prefix sum over the sum of all arcs
    // in each community (which are no selfloop) and write them in parallel to
    // the tmp adjacence array.
    utils::Timer::instance().start_timer("construct_tmp_adjacent_array", "Construct Tmp Adjacent Array");
    parallel::scalable_vector<parallel::AtomicWrapper<ArcWeight>>
      coarse_node_volumes(coarse_graph._num_nodes);
    tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
      const NodeID coarse_u = communities[u];
      coarse_node_volumes[coarse_u] += _node_volumes[u];
      for ( const Arc& arc : arcsOf(u) ) {
        const NodeID coarse_v = communities[arc.head];
        if ( coarse_u != coarse_v ) {
          ASSERT(static_cast<size_t>(coarse_u) < tmp_adjacent_nodes.size());
          ++tmp_adjacent_nodes[coarse_u];
        }
      }
    });

    parallel::TBBPrefixSum<parallel::IntegralAtomicWrapper<size_t>>
      tmp_adjacent_nodes_prefix_sum(tmp_adjacent_nodes);
    tbb::parallel_scan(tbb::blocked_range<size_t>(
      0UL, tmp_adjacent_nodes.size()), tmp_adjacent_nodes_prefix_sum);

    parallel::scalable_vector<Arc> tmp_arcs(tmp_adjacent_nodes_prefix_sum.total_sum());
    parallel::scalable_vector<size_t> valid_arcs;
    tbb::parallel_invoke([&] {
      tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
        const NodeID coarse_u = communities[u];
        ASSERT(static_cast<size_t>(coarse_u) < coarse_graph._num_nodes);
        for ( const Arc& arc : arcsOf(u) ) {
          const NodeID coarse_v = communities[arc.head];
          if ( coarse_u != coarse_v ) {
            const size_t tmp_arcs_pos = tmp_adjacent_nodes_prefix_sum[coarse_u]
              + tmp_pos[coarse_u]++;
            ASSERT(tmp_arcs_pos < tmp_adjacent_nodes_prefix_sum[coarse_u + 1]);
            tmp_arcs[tmp_arcs_pos] = Arc { coarse_v, arc.weight };
          }
        }
      });
    }, [&] {
      valid_arcs.assign(tmp_adjacent_nodes_prefix_sum.total_sum(), 1UL);
    });
    utils::Timer::instance().stop_timer("construct_tmp_adjacent_array");

    // #################### STAGE 3 ####################
    // Aggregate weights of arcs that are equal in each community.
    // Therefore, we sort the arcs according to their endpoints
    // and aggregate weight of arcs with equal endpoints.
    utils::Timer::instance().start_timer("contract_arcs", "Contract Arcs");
    tbb::parallel_for(0U, static_cast<NodeID>(coarse_graph._num_nodes), [&](const NodeID u) {
      const size_t tmp_arc_start = tmp_adjacent_nodes_prefix_sum[u];
      const size_t tmp_arc_end = tmp_adjacent_nodes_prefix_sum[u + 1];
      std::sort(tmp_arcs.begin() + tmp_arc_start, tmp_arcs.begin() + tmp_arc_end,
        [&](const Arc& lhs, const Arc& rhs) {
          return lhs.head < rhs.head;
        });

      size_t arc_rep = tmp_arc_start;
      for ( size_t pos = tmp_arc_start + 1; pos < tmp_arc_end; ++pos ) {
        if ( tmp_arcs[arc_rep].head == tmp_arcs[pos].head ) {
          tmp_arcs[arc_rep].weight += tmp_arcs[pos].weight;
          valid_arcs[pos] = 0UL;
        } else {
          arc_rep = pos;
        }
      }
    });
    utils::Timer::instance().stop_timer("contract_arcs");


    // #################### STAGE 4 ####################
    // Write aggregated arcs to coarse graph. Therefore, we compute
    // a prefix sum over all valid arcs in the tmp adjacence array which
    // determines their position in the coarse graph adjacence array.
    utils::Timer::instance().start_timer("construct_coarse_graph", "Construct Coarse Graph");
    parallel::TBBPrefixSum<size_t> indices_prefix_sum(valid_arcs);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, valid_arcs.size()), indices_prefix_sum);

    coarse_graph._num_arcs = indices_prefix_sum.total_sum();
    coarse_graph._arcs.resize(coarse_graph._num_arcs);
    tbb::parallel_invoke([&] {
      tbb::parallel_for(0UL, tmp_arcs.size(), [&](const size_t i) {
        if ( indices_prefix_sum.value(i) ) {
          const size_t pos = indices_prefix_sum[i];
          ASSERT(pos < coarse_graph._arcs.size());
          coarse_graph._arcs[pos] = std::move(tmp_arcs[i]);
        }
      });
    }, [&] {
      tbb::parallel_for(0U, static_cast<NodeID>(coarse_graph._num_nodes + 1), [&](const NodeID u) {
        const size_t start = indices_prefix_sum[tmp_adjacent_nodes_prefix_sum[u]];
        ASSERT(start <= coarse_graph._arcs.size());
        coarse_graph._indices[u] = start;
        if ( u < coarse_graph._num_nodes ) {
          coarse_graph._node_volumes[u] = coarse_node_volumes[u];
        }
      });
    });

    utils::Timer::instance().stop_timer("construct_coarse_graph");

    // Free local data parallel
    parallel::parallel_free(tmp_adjacent_nodes, tmp_pos,
      mapping, tmp_arcs, valid_arcs, coarse_node_volumes);

    return coarse_graph;
  }

 private:
  GraphT() :
    _num_nodes(0),
    _num_arcs(0),
    _total_volume(0),
    _indices(),
    _arcs(),
    _node_volumes() { }

  template<typename F>
  void construct(const HyperGraph& hypergraph,
                 const TaskGroupID task_group_id,
                 const F& edge_weight_func) {
    // Initialize data structure
    tbb::parallel_invoke([&] {
      _indices.resize(_num_nodes + 1);
    }, [&] {
      _arcs.resize(_num_arcs);
    }, [&] {
      _node_volumes.resize(_num_nodes);
    });

    // Compute index vector with a prefix sum
    _indices[0] = 0; // Sentinel
    tbb::parallel_invoke([&] {
      hypergraph.doParallelForAllNodes(task_group_id, [&](const HypernodeID hn) {
        const HypernodeID original_id = hypergraph.originalNodeID(hn);
        _indices[original_id + 1] = hypergraph.nodeDegree(hn);
      });
    }, [&] {
      hypergraph.doParallelForAllEdges(task_group_id, [&](const HyperedgeID he) {
        const HyperedgeID original_id = hypergraph.originalEdgeID(he);
        _indices[original_id + hypergraph.initialNumNodes() + 1] =
          hypergraph.edgeSize(he);
      });
    });

    parallel::TBBPrefixSum<size_t> _indices_prefix_sum(_indices);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, _indices.size()), _indices_prefix_sum);
    ASSERT(_indices.back() == _num_arcs);

    // Insert Arcs
    tbb::parallel_invoke([&] {
      hypergraph.doParallelForAllNodes(task_group_id, [&](const HypernodeID hn) {
        const NodeID u = hypergraph.originalNodeID(hn);
        const HyperedgeID node_degree = hypergraph.nodeDegree(hn);
        ASSERT(u < _num_nodes);
        size_t pos = _indices[u];
        for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
          const NodeID v = hypergraph.originalEdgeID(he) + hypergraph.initialNumNodes();
          const HyperedgeWeight edge_weight = hypergraph.edgeWeight(he);
          const HypernodeID edge_size = hypergraph.edgeSize(he);
          ASSERT(pos < _indices[u + 1]);
          _arcs[pos++] = Arc(v, edge_weight_func(edge_weight, edge_size, node_degree));
        }
      });
    }, [&] {
      hypergraph.doParallelForAllEdges(task_group_id, [&](const HyperedgeID he) {
        const NodeID u = hypergraph.originalEdgeID(he) + hypergraph.initialNumNodes();
        const HyperedgeWeight edge_weight = hypergraph.edgeWeight(he);
        const HypernodeID edge_size = hypergraph.edgeSize(he);
        ASSERT(u < _num_nodes);
        size_t pos = _indices[u];
        for ( const HypernodeID& hn : hypergraph.pins(he) ) {
          const NodeID v = hypergraph.originalNodeID(hn);
          const HyperedgeID node_degree = hypergraph.nodeDegree(hn);
          ASSERT(pos < _indices[u + 1]);
          _arcs[pos++] = Arc(v, edge_weight_func(edge_weight, edge_size, node_degree));
        }
      });
    });


    // Compute node volumes and total volume
    _total_volume = 0.0;
    tbb::enumerable_thread_specific<ArcWeight> local_total_volume(0.0);
    tbb::parallel_for(0U, static_cast<NodeID>(numNodes()), [&](const NodeID u) {
      local_total_volume.local() += computeNodeVolume(u);
    });
    _total_volume = local_total_volume.combine(std::plus<ArcWeight>());
  }

  ArcWeight computeNodeVolume(const NodeID u) {
    for (const Arc& arc : arcsOf(u)) {
      _node_volumes[u] += arc.weight;
    }
    return _node_volumes[u];
  }

  size_t hash(const NodeID u, const NodeID v) const {
    return static_cast<size_t>(u) ^ ( static_cast<size_t>(v) << 8 );
  }

  size_t _num_nodes;
  size_t _num_arcs;
  ArcWeight _total_volume;

  parallel::scalable_vector<size_t> _indices;
  parallel::scalable_vector<Arc> _arcs;
  parallel::scalable_vector<ArcWeight> _node_volumes;
};

}  // namespace ds
}  // namespace mt_kahypar

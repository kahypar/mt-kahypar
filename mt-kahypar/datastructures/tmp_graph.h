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
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/utils/range.h"

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
    _indices(),
    _arcs() {
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

 private:
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

  const size_t _num_nodes;
  const size_t _num_arcs;
  ArcWeight _total_volume;

  parallel::scalable_vector<size_t> _indices;
  parallel::scalable_vector<Arc> _arcs;
  parallel::scalable_vector<ArcWeight> _node_volumes;
};

}  // namespace ds
}  // namespace mt_kahypar

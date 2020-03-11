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

 private:
  struct ContractedArc {
    NodeID u;
    NodeID v;
    ArcWeight weight;
    bool valid;
  };

 public:
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
    parallel::TBBPrefixSum<size_t> mapping_prefix_sum(mapping);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, mapping.size()), mapping_prefix_sum);
    // Remap community ids
    tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
      communities[u] = mapping_prefix_sum[communities[u]];
    });
    utils::Timer::instance().stop_timer("compute_cluster_mapping");

    // #################### STAGE 2 ####################
    // For each node u, we contract its adjacent arcs and store them
    // in a bucket of a streaming map. The bucket is determined by a
    // hash over the two endpoints of a arc. Therefore, arcs with equal
    // endpoints are later in the same bucket in our streaming map
    // and can be processed together sequentially.
    utils::Timer::instance().start_timer("local_arc_contraction", "Local Arc Contraction");
    coarse_graph._num_nodes = mapping_prefix_sum.total_sum();
    StreamingMap<size_t, ContractedArc> arc_hash_to_contracted_arc;
    tbb::enumerable_thread_specific<parallel::scalable_vector<ArcWeight>> tmp_node_volumes(coarse_graph._num_nodes);
    tbb::parallel_invoke([&] {
      tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
        if ( degree(u) > 0 ) {
          const NodeID coarse_u = communities[u];
          tmp_node_volumes.local()[coarse_u] += _node_volumes[u];
          std::sort(_arcs.begin() + _indices[u], _arcs.begin() + _indices[u + 1],
            [&](const Arc& lhs, const Arc& rhs) {
              return communities[lhs.head] < communities[rhs.head];
            });

          size_t pos = _indices[u];
          ContractedArc c_arc = ContractedArc { coarse_u,
                                                static_cast<NodeID>(communities[_arcs[pos].head]),
                                                _arcs[pos++].weight, true };
          for ( ; pos < _indices[u + 1]; ++pos ) {
            const NodeID coarse_v = communities[_arcs[pos].head];
            if ( c_arc.v == coarse_v ) {
              c_arc.weight += _arcs[pos].weight;
            } else {
              arc_hash_to_contracted_arc.stream(hash(c_arc.u, c_arc.v), std::move(c_arc));
              c_arc = ContractedArc { coarse_u, coarse_v, _arcs[pos].weight, true };
            }
          }
          arc_hash_to_contracted_arc.stream(hash(c_arc.u, c_arc.v), std::move(c_arc));
        }
      });
    }, [&] {
      coarse_graph._indices.resize(coarse_graph._num_nodes + 1);
      coarse_graph._node_volumes.resize(coarse_graph._num_nodes);
    });

    using ArcMap = parallel::scalable_vector<parallel::scalable_vector<ContractedArc>>;
    ArcMap arc_buckets(arc_hash_to_contracted_arc.size());
    arc_hash_to_contracted_arc.copy(arc_buckets, [&](const size_t key) {
      return key % arc_buckets.size();
    });
    utils::Timer::instance().stop_timer("local_arc_contraction");

    // #################### STAGE 3 ####################
    // Arcs with same endpoints are in the same bucket now and weights can be
    // aggregated by simply sorting a bucket after endpoints.
    utils::Timer::instance().start_timer("global_arc_contraction", "Global Arc Contraction");
    parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>> tmp_indices(
      coarse_graph._num_nodes, parallel::IntegralAtomicWrapper<size_t>(0));
    parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>> arc_pos;
    tbb::parallel_invoke([&] {
      tbb::parallel_for(0UL, arc_buckets.size(), [&](const size_t bucket) {
        parallel::scalable_vector<ContractedArc>& arc_bucket = arc_buckets[bucket];
        if ( arc_bucket.size() > 0 ) {
          std::sort(arc_bucket.begin(), arc_bucket.end(),
            [&](const ContractedArc& lhs, const ContractedArc& rhs) {
              return lhs.u < rhs.u || (lhs.u == rhs.u && lhs.v < rhs.v);
            });

          for ( size_t i = 0; i < arc_bucket.size(); ++i ) {
            ContractedArc& lhs = arc_bucket[i];
            if ( lhs.valid ) {
              for ( size_t j = i + 1; j < arc_bucket.size(); ++j ) {
                ContractedArc& rhs = arc_bucket[j];
                ASSERT(rhs.valid);
                if ( lhs.u == rhs.u && lhs.v == rhs.v ) {
                  lhs.weight += rhs.weight;
                  rhs.valid = false;
                } else {
                  break;
                }
              }

              if ( lhs.u == lhs.v ) {
                // Selfloops are not part of the coarse graph
                lhs.valid = false;
              } else {
                ASSERT(static_cast<size_t>(lhs.u) < tmp_indices.size());
                ++tmp_indices[lhs.u];
              }
            }
          }
        }
      });
    }, [&] {
      arc_pos.assign(coarse_graph._num_nodes, parallel::IntegralAtomicWrapper<size_t>(0));
    });

    parallel::TBBPrefixSum<parallel::IntegralAtomicWrapper<size_t>> index_prefix_sum(tmp_indices);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, tmp_indices.size()), index_prefix_sum);
    coarse_graph._num_arcs = index_prefix_sum.total_sum();
    utils::Timer::instance().stop_timer("global_arc_contraction");

    // #################### STAGE 4 ####################
    // Construct coarse graph
    utils::Timer::instance().start_timer("construct_coarse_graph", "Contract Coarse Graph");
    coarse_graph._arcs.resize(coarse_graph._num_arcs);
    tbb::parallel_invoke([&] {
      tbb::parallel_for(0UL, arc_buckets.size(), [&](const size_t bucket) {
        parallel::scalable_vector<ContractedArc>& arc_bucket = arc_buckets[bucket];
        for ( size_t i = 0; i < arc_bucket.size(); ++i ) {
          ContractedArc& lhs = arc_bucket[i];
          if ( lhs.valid ) {
            const size_t pos = index_prefix_sum[lhs.u] + arc_pos[lhs.u]++;
            ASSERT(pos < index_prefix_sum[lhs.u + 1]);
            coarse_graph._arcs[pos] = Arc { lhs.v, lhs.weight };
          }
        }
        parallel::free(arc_bucket);
      });
    }, [&] {
      tbb::parallel_for(0U, static_cast<NodeID>(coarse_graph._num_nodes + 1), [&](const NodeID u) {
        coarse_graph._indices[u] = index_prefix_sum[u];
        for ( const parallel::scalable_vector<ArcWeight>& local_node_volumes : tmp_node_volumes ) {
          coarse_graph._node_volumes[u] += local_node_volumes[u];
        }
      });
    });
    utils::Timer::instance().stop_timer("construct_coarse_graph");

    parallel::parallel_free(mapping, tmp_indices, arc_pos);

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

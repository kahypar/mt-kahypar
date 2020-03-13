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

  static size_t NUM_BLOCKS;
  static size_t MIN_BLOCK_SIZE;

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
  template<typename T>
  using BlockScalableVector = parallel::scalable_vector<parallel::scalable_vector<T>>;

  struct AdjacenceArrayBlock {
    explicit AdjacenceArrayBlock() :
      _num_nodes(0),
      _indices(),
      _arcs(),
      _node_volumes() { }

    explicit AdjacenceArrayBlock(const size_t num_nodes) :
      _num_nodes(num_nodes),
      _indices(num_nodes + 1),
      _arcs(),
      _node_volumes(num_nodes) { }

    inline IteratorRange<AdjacenceIterator> arcsOf(const NodeID u) const {
      const size_t local_u = u % _num_nodes;
      return IteratorRange<AdjacenceIterator>(
        _arcs.cbegin() + _indices[local_u],
        _arcs.cbegin() + _indices[local_u + 1]);
    }

    inline size_t degree(const NodeID u) const {
      const size_t local_u = u % _num_nodes;
      return _indices[local_u + 1] - _indices[local_u];
    }

    inline ArcWeight nodeVolume(const NodeID u) const {
      const size_t local_u = u % _num_nodes;
      return _node_volumes[local_u];
    }

    inline ArcWeight computeNodeVolume(const NodeID u) {
      const size_t local_u = u % _num_nodes;
      for (const Arc& arc : arcsOf(u)) {
        _node_volumes[local_u] += arc.weight;
      }
      return _node_volumes[local_u];
    }

    size_t _num_nodes;
    parallel::scalable_vector<size_t> _indices;
    parallel::scalable_vector<Arc> _arcs;
    parallel::scalable_vector<ArcWeight> _node_volumes;
  };

 public:
  explicit GraphT(const HyperGraph& hypergraph,
                  const LouvainEdgeWeight edge_weight_type) :
    _num_nodes(hypergraph.initialNumNodes() + hypergraph.initialNumEdges()),
    _num_arcs(2 * hypergraph.initialNumPins()),
    _block_size(0),
    _total_volume(0),
    _blocks() {
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
    const size_t block_u = u / _block_size;
    ASSERT(block_u < _blocks.size());
    return _blocks[block_u].arcsOf(u);
  }

  size_t degree(const NodeID u) const {
    ASSERT(u < _num_nodes);
    const size_t block_u = u / _block_size;
    ASSERT(block_u < _blocks.size());
    return _blocks[block_u].degree(u);
  }

  ArcWeight totalVolume() const {
    return _total_volume;
  }

  ArcWeight nodeVolume(const NodeID u) const {
    ASSERT(u < _num_nodes);
    const size_t block_u = u / _block_size;
    ASSERT(block_u < _blocks.size());
    return _blocks[block_u].nodeVolume(u);
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
    coarse_graph._num_nodes = mapping_prefix_sum.total_sum();
    tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
      communities[u] = mapping_prefix_sum[communities[u]];
    });
    utils::Timer::instance().stop_timer("compute_cluster_mapping");

    // #################### STAGE 2 ####################
    // Write all arcs, that are not a selfloop in the coarse graph, into a tmp
    // adjacence array. For that, we compute a prefix sum over the sum of all arcs
    // in each community (which are no selfloop) and write them in parallel to
    // the tmp adjacence array. Note, our graph data structure is splitted into
    // adjacence array blocks, thus we create for all blocks tmp blocks and
    // process them in parallel.
    utils::Timer::instance().start_timer("construct_tmp_adjacent_array", "Construct Tmp Adjacent Array");
    // Compute number of adjacence array blocks in coarse graph
    size_t num_blocks = 0;
    if ( coarse_graph._num_nodes < NUM_BLOCKS ) {
      coarse_graph._block_size = coarse_graph._num_nodes;
      num_blocks = 1;
    } else {
      coarse_graph._block_size = std::max( coarse_graph._num_nodes / NUM_BLOCKS +
        (coarse_graph._num_nodes % NUM_BLOCKS != 0), MIN_BLOCK_SIZE );
      num_blocks = coarse_graph._num_nodes / coarse_graph._block_size +
        (coarse_graph._num_nodes % coarse_graph._block_size != 0);
    }

    // Initialize tmp adjacence array blocks in parallel
    coarse_graph._blocks.resize(num_blocks);
    parallel::scalable_vector<AdjacenceArrayBlock> tmp_blocks(num_blocks);
    BlockScalableVector<parallel::IntegralAtomicWrapper<size_t>> tmp_pos(num_blocks);
    BlockScalableVector<parallel::IntegralAtomicWrapper<size_t>> tmp_adjacent_nodes(num_blocks);
    tbb::parallel_for(0UL, num_blocks, [&](const size_t i) {
      const NodeID start = i * coarse_graph._block_size;
      const NodeID end = std::min((i + 1) * coarse_graph._block_size, coarse_graph._num_nodes);
      ASSERT(start < end);
      const size_t num_nodes = end - start;
      coarse_graph._blocks[i]._num_nodes = num_nodes;
      tmp_blocks[i]._num_nodes = num_nodes;
      tmp_pos[i].assign(num_nodes, parallel::IntegralAtomicWrapper<size_t>(0));
      tmp_adjacent_nodes[i].assign(num_nodes, parallel::IntegralAtomicWrapper<size_t>(0));
    });

    // Compute number of arcs in tmp adjacence array block with parallel prefix sum
    parallel::scalable_vector<parallel::AtomicWrapper<ArcWeight>> coarse_node_volumes(coarse_graph._num_nodes);
    tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
      const NodeID coarse_u = communities[u];
      const NodeID coarse_block_u = coarse_u / coarse_graph._block_size;
      const NodeID coarse_local_u = coarse_u % coarse_graph._block_size;
      ASSERT(static_cast<size_t>(coarse_u) < coarse_graph._num_nodes);
      ASSERT(static_cast<size_t>(coarse_block_u) < tmp_adjacent_nodes.size());
      ASSERT(static_cast<size_t>(coarse_local_u) < tmp_adjacent_nodes[coarse_block_u].size());
      coarse_node_volumes[coarse_u] += nodeVolume(u);
      for ( const Arc& arc : arcsOf(u) ) {
        const NodeID coarse_v = communities[arc.head];
        if ( coarse_u != coarse_v ) {
          ++tmp_adjacent_nodes[coarse_block_u][coarse_local_u];
        }
      }
    });
    parallel::scalable_vector<parallel::TBBPrefixSum<
      parallel::IntegralAtomicWrapper<size_t>>> tmp_adjacent_nodes_prefix_sum;
    for ( size_t i = 0; i < num_blocks; ++i ) {
      tmp_adjacent_nodes_prefix_sum.emplace_back(tmp_adjacent_nodes[i]);
    }
    tbb::parallel_for(0UL, num_blocks, [&](const size_t i) {
      tbb::parallel_scan(tbb::blocked_range<size_t>(
        0UL, tmp_adjacent_nodes[i].size()), tmp_adjacent_nodes_prefix_sum[i]);
      tmp_blocks[i]._arcs.resize(tmp_adjacent_nodes_prefix_sum[i].total_sum());
    });

    // Write all arcs into corresponding tmp adjacence array blocks
    tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
      const NodeID coarse_u = communities[u];
      const NodeID coarse_block_u = coarse_u / coarse_graph._block_size;
      const NodeID coarse_local_u = coarse_u % coarse_graph._block_size;
      ASSERT(static_cast<size_t>(coarse_u) < coarse_graph._num_nodes);
      ASSERT(static_cast<size_t>(coarse_block_u) < tmp_adjacent_nodes.size());
      ASSERT(static_cast<size_t>(coarse_local_u) < tmp_adjacent_nodes[coarse_block_u].size());
      for ( const Arc& arc : arcsOf(u) ) {
        const NodeID coarse_v = communities[arc.head];
        if ( coarse_u != coarse_v ) {
          const size_t tmp_arcs_pos =
            tmp_adjacent_nodes_prefix_sum[coarse_block_u][coarse_local_u] +
            tmp_pos[coarse_block_u][coarse_local_u]++;
          ASSERT(tmp_arcs_pos < tmp_adjacent_nodes_prefix_sum[coarse_block_u][coarse_local_u + 1]);
          tmp_blocks[coarse_block_u]._arcs[tmp_arcs_pos] = Arc { coarse_v, arc.weight };
        }
      }
    });
    utils::Timer::instance().stop_timer("construct_tmp_adjacent_array");


    // #################### STAGE 3 ####################
    // Aggregate weights of arcs that are equal in each community.
    // Therefore, we sort the arcs according to their endpoints
    // and aggregate weight of arcs with equal endpoints.
    tbb::enumerable_thread_specific<size_t> local_num_arcs(0);
    utils::Timer::instance().start_timer("contract_arcs", "Contract Arcs");
    tbb::parallel_for(0UL, num_blocks, [&](const size_t i) {
      parallel::scalable_vector<size_t> valid_arcs(tmp_adjacent_nodes_prefix_sum[i].total_sum(), 1UL);
      tbb::parallel_for(0U, static_cast<NodeID>(tmp_blocks[i]._num_nodes), [&, i](const NodeID u) {
        AdjacenceArrayBlock& block = tmp_blocks[i];
        const size_t tmp_arc_start = tmp_adjacent_nodes_prefix_sum[i][u];
        const size_t tmp_arc_end = tmp_adjacent_nodes_prefix_sum[i][u + 1];
        std::sort(block._arcs.begin() + tmp_arc_start, block._arcs.begin() + tmp_arc_end,
          [&](const Arc& lhs, const Arc& rhs) {
            return lhs.head < rhs.head;
          });

        size_t arc_rep = tmp_arc_start;
        for ( size_t pos = tmp_arc_start + 1; pos < tmp_arc_end; ++pos ) {
          if ( block._arcs[arc_rep].head == block._arcs[pos].head ) {
            block._arcs[arc_rep].weight += block._arcs[pos].weight;
            valid_arcs[pos] = 0UL;
          } else {
            arc_rep = pos;
          }
        }
      });

      // Write all arcs to coarse graph
      parallel::TBBPrefixSum<size_t> indices_prefix_sum(valid_arcs);
      tbb::parallel_scan(tbb::blocked_range<size_t>(0UL,
        valid_arcs.size()), indices_prefix_sum);
      local_num_arcs.local() += indices_prefix_sum.total_sum();
      AdjacenceArrayBlock& block = coarse_graph._blocks[i];
      AdjacenceArrayBlock& tmp_block = tmp_blocks[i];
      block._indices.resize(tmp_block._num_nodes + 1);
      block._arcs.resize(indices_prefix_sum.total_sum());
      block._node_volumes.resize(tmp_block._num_nodes);
      tbb::parallel_invoke([&] {
        tbb::parallel_for(0UL, tmp_block._arcs.size(), [&](const size_t j) {
          if ( indices_prefix_sum.value(j) ) {
            const size_t pos = indices_prefix_sum[j];
            ASSERT(pos < block._arcs.size());
            block._arcs[pos] = std::move(tmp_block._arcs[j]);
          }
        });
      }, [&] {
        const NodeID start_u = i * coarse_graph._block_size;
        tbb::parallel_for(0U, static_cast<NodeID>(block._num_nodes + 1), [&, i](const NodeID u) {
          const size_t start = indices_prefix_sum[tmp_adjacent_nodes_prefix_sum[i][u]];
          ASSERT(start <= block._arcs.size());
          block._indices[u] = start;
          if ( u < block._num_nodes ) {
            block._node_volumes[u] = coarse_node_volumes[start_u + u];
          }
        });
      });

      // Free local data parallel
      parallel::parallel_free(tmp_block._arcs, tmp_adjacent_nodes[i], tmp_pos[i]);
    });
    coarse_graph._num_arcs += local_num_arcs.combine(std::plus<size_t>());
    utils::Timer::instance().stop_timer("contract_arcs");

    // Free local data parallel
    parallel::parallel_free(mapping, coarse_node_volumes);

    return coarse_graph;
  }

 private:
  GraphT() :
    _num_nodes(0),
    _num_arcs(0),
    _block_size(0),
    _total_volume(0),
    _blocks() { }

  template<typename F>
  void construct(const HyperGraph& hypergraph,
                 const F& edge_weight_func) {
    // Initialize data structure
    size_t num_blocks = 0;
    if ( _num_nodes < NUM_BLOCKS ) {
      _block_size = _num_nodes;
      num_blocks = 1;
    } else {
      _block_size = std::max( _num_nodes / NUM_BLOCKS +
        (_num_nodes % NUM_BLOCKS != 0), MIN_BLOCK_SIZE );
      num_blocks = _num_nodes / _block_size + (_num_nodes % _block_size != 0);
    }

    utils::Timer::instance().start_timer("construct_arcs", "Construct Arcs");
    _blocks.resize(num_blocks);
    tbb::parallel_for(0UL, num_blocks, [&](const size_t i) {
      const NodeID start = i * _block_size;
      const NodeID end = std::min((i + 1) * _block_size, _num_nodes);
      ASSERT(start < end);
      const size_t num_nodes = end - start;

      AdjacenceArrayBlock& block = _blocks[i];
      block._num_nodes = num_nodes;
      block._indices.resize(num_nodes + 1);
      block._node_volumes.resize(num_nodes);
      const HypernodeID num_hypernodes = hypergraph.initialNumNodes();
      block._indices[0] = 0;
      tbb::parallel_for(start, end, [&](const NodeID u) {
        const NodeID local_u = u % num_nodes;
        if ( u < num_hypernodes ) {
          const HypernodeID hn = hypergraph.globalNodeID(u);
          block._indices[local_u + 1] = hypergraph.nodeDegree(hn);
        } else {
          const HyperedgeID he = hypergraph.globalEdgeID(u - num_hypernodes);
          block._indices[local_u + 1] = hypergraph.edgeSize(he);
        }
      });

      parallel::TBBPrefixSum<size_t> indices_prefix_sum(block._indices);
      tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, block._indices.size()), indices_prefix_sum);

      // Insert Arcs
      block._arcs.resize(indices_prefix_sum.total_sum());
      tbb::parallel_for(start, end, [&](const NodeID u) {
        const NodeID local_u = u % num_nodes;
        size_t pos = block._indices[local_u];
        if ( u < num_hypernodes ) {
          const HypernodeID hn = hypergraph.globalNodeID(u);
          const HyperedgeID node_degree = hypergraph.nodeDegree(hn);
          for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
            const NodeID v = hypergraph.originalEdgeID(he) + num_hypernodes;
            const HyperedgeWeight edge_weight = hypergraph.edgeWeight(he);
            const HypernodeID edge_size = hypergraph.edgeSize(he);
            ASSERT(pos < block._indices[local_u + 1]);
            block._arcs[pos++] = Arc(v, edge_weight_func(edge_weight, edge_size, node_degree));
          }
        } else {
          const HyperedgeID he = hypergraph.globalEdgeID(u - num_hypernodes);
          const HyperedgeWeight edge_weight = hypergraph.edgeWeight(he);
          const HypernodeID edge_size = hypergraph.edgeSize(he);
          for ( const HypernodeID& pin : hypergraph.pins(he) ) {
            const NodeID v = hypergraph.originalNodeID(pin);
            const HyperedgeID node_degree = hypergraph.nodeDegree(pin);
            ASSERT(pos < block._indices[local_u + 1]);
            block._arcs[pos++] = Arc(v, edge_weight_func(edge_weight, edge_size, node_degree));
          }
        }
      });
    });
    utils::Timer::instance().stop_timer("construct_arcs");

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

  ArcWeight computeNodeVolume(const NodeID u) {
    const size_t block_u = u / _block_size;
    ASSERT(block_u < _blocks.size());
    return _blocks[block_u].computeNodeVolume(u);
  }

  size_t _num_nodes;
  size_t _num_arcs;
  size_t _block_size;
  ArcWeight _total_volume;

  parallel::scalable_vector<AdjacenceArrayBlock> _blocks;
};

template <typename HyperGraph>
size_t GraphT<HyperGraph>::NUM_BLOCKS = 512;
template <typename HyperGraph>
size_t GraphT<HyperGraph>::MIN_BLOCK_SIZE = 32;

}  // namespace ds
}  // namespace mt_kahypar

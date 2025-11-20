/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#pragma once

#include <atomic>

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/weight/hypernode_weight_common.h"

namespace mt_kahypar {
namespace ds {

template<class Hypergraph>
class FixedVertexSupport {

  static constexpr bool debug = false;

  struct FixedVertexData {
    // ! Fixed vertex block ID
    PartitionID block;
    // ! Number of fixed vertices contracted onto this node
    HypernodeID fixed_vertex_contraction_cnt;
    // ! Weight at the time it becomes fixed
    // HypernodeWeight fixed_vertex_weight;
    // ! Spin lock to syncronize contractions
    SpinLock sync;
  };

 public:
  FixedVertexSupport() :
    _num_nodes(0),
    _k(kInvalidPartition),
    _hg(nullptr),
    _total_fixed_vertex_weight(),
    _fixed_vertex_block_weights(),
    _max_block_weights(),
    _fixed_vertex_data() { }

  FixedVertexSupport(const HypernodeID num_nodes,
                     const Dimension dimension,
                     const PartitionID k) :
    _num_nodes(num_nodes),
    _k(k),
    _hg(nullptr),
    _total_fixed_vertex_weight(dimension, 0),
    _fixed_vertex_block_weights(k, dimension, 0, false),
    _max_block_weights(k, dimension, std::numeric_limits<HNWeightScalar>::max(), false),
    _fixed_vertex_data(num_nodes, FixedVertexData { kInvalidPartition, 0, SpinLock() }),
    _fixed_vertex_hn_weights(num_nodes, dimension, 0, false) { }

  FixedVertexSupport(const FixedVertexSupport&) = delete;
  FixedVertexSupport & operator= (const FixedVertexSupport &) = delete;

  FixedVertexSupport(FixedVertexSupport&&) = default;
  FixedVertexSupport & operator= (FixedVertexSupport &&) = default;

  void setHypergraph(const Hypergraph* hg) {
    _hg = hg;
  }

  void setMaxBlockWeight(const HypernodeWeightArray& max_block_weights) {
    if ( hasFixedVertices() ) {
      ASSERT(max_block_weights.size() >= static_cast<size_t>(_k));
      _max_block_weights = max_block_weights.copy();
    }
  }

  PartitionID numBlocks() const {
    return _k;
  }

  // ####################### Fixed Vertex Block Weights #######################

  bool hasFixedVertices() const {
    return !weight::isZero(_total_fixed_vertex_weight.load(std::memory_order_relaxed));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightAtomicCRef totalFixedVertexWeight() const {
    return _total_fixed_vertex_weight.get().load(std::memory_order_relaxed);
  }

  // ! Returns the weight of all fixed vertices assigned to the corresponding block
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightAtomicCRef fixedVertexBlockWeight(const PartitionID block) const {
    ASSERT(block != kInvalidPartition && block < _k);
    return _fixed_vertex_block_weights[block].load(std::memory_order_relaxed);
  }

  // ####################### Fixed Vertex Information #######################

  // ! Fixes a node to a block
  void fixToBlock(const HypernodeID hn, const PartitionID block) {
    ASSERT(_hg);
    ASSERT(hn < _num_nodes);
    ASSERT(block != kInvalidPartition && block < _k);
    PartitionID expected = kInvalidPartition;
    PartitionID desired = block;
    if ( __atomic_compare_exchange_n(&_fixed_vertex_data[hn].block,
           &expected, desired, false, __ATOMIC_ACQ_REL, __ATOMIC_RELAXED) ) {
      const HNWeightConstRef weight_of_hn = _hg->nodeWeight(hn);
      _fixed_vertex_data[hn].fixed_vertex_contraction_cnt = 1;
      _fixed_vertex_hn_weights[hn] = weight_of_hn;
      weight::eval(_fixed_vertex_block_weights[block].fetch_add(
        weight_of_hn, std::memory_order_relaxed));
      weight::eval(_total_fixed_vertex_weight.get().fetch_add(
        weight_of_hn, std::memory_order_relaxed));
    } else {
      ASSERT(_fixed_vertex_data[hn].block == block,
        "Try to fix hypernode" << hn << "to block" << block
        << ", but it is already fixed to block" << _fixed_vertex_data[hn].block);
    }
  }

  // ! Returns whether or not the node is fixed to a block
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool isFixed(const HypernodeID hn) const {
    return hn < _num_nodes && fixedVertexBlock(hn) != kInvalidPartition;
  }

  // ! Returns the fixed vertex block of the node
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE PartitionID fixedVertexBlock(const HypernodeID hn) const {
    ASSERT(hn < _num_nodes);
    return __atomic_load_n(&_fixed_vertex_data[hn].block, __ATOMIC_RELAXED);
  }

  // ####################### (Un)contractions #######################

  // ! Contracts v onto u. If v is a fixed vertex than u becomes also an fixed vertex.
  // ! If u and v are fixed vertices, then both must be assigned to same block
  // ! The function returns false, if u and v are fixed and are assigned to different blocks
  bool contract(const HypernodeID u, const HypernodeID v);

  // ! Contracts v onto u. Unlike `contract`, this assumes than any contractions of other nodes
  // ! onto v form a new cluster without v instead of joining the same cluster (see deterministic coarsening).
  bool contractWithoutChains(const HypernodeID u, const HypernodeID v);

  // ! Uncontract v from u. This reverts the corresponding contraction operation of v onto u.
  void uncontract(const HypernodeID u, const HypernodeID v);

  // ####################### Miscellaneous #######################

  // ! Only for testing
  bool verifyClustering(const vec<HypernodeID>& cluster_ids) const;

  FixedVertexSupport<Hypergraph> copy() const {
    FixedVertexSupport<Hypergraph> cpy;
    cpy._num_nodes = _num_nodes;
    cpy._k = _k;
    cpy._hg = _hg;
    cpy._total_fixed_vertex_weight = _total_fixed_vertex_weight.copy();
    cpy._fixed_vertex_block_weights = _fixed_vertex_block_weights.copy();
    cpy._max_block_weights = _max_block_weights.copy();
    cpy._fixed_vertex_data = _fixed_vertex_data;
    return cpy;
  }

  size_t size_in_bytes() const {
    const auto dimension = _total_fixed_vertex_weight.dimension();
    return ( 2 * sizeof(HNWeightScalar) * dimension) * _k +
      (sizeof(FixedVertexData) + sizeof(HNWeightScalar) * dimension) * _num_nodes +
      sizeof(HNWeightScalar) * dimension;
  }

  Dimension dimension() const {
    return _total_fixed_vertex_weight.dimension();
  }

 private:
  bool contractImpl(const HypernodeID u, const HypernodeID v, bool ignore_v);

  // ! Number of nodes
  HypernodeID _num_nodes;

  // ! Number of blocks
  PartitionID _k;

  // ! Underlying hypergraph
  const Hypergraph* _hg;

  // ! Total weight of all fixed vertices
  AllocatedHNWeight _total_fixed_vertex_weight;

  // ! Weight of all vertices fixed to a block
  HypernodeWeightArray _fixed_vertex_block_weights;

  // ! Maximum allowed fixed vertex block weight
  HypernodeWeightArray _max_block_weights;

  // ! Fixed vertex block IDs of each node
  vec<FixedVertexData> _fixed_vertex_data;

  // ! Fixed vertex weights of each node
  HypernodeWeightArray _fixed_vertex_hn_weights;
};

}  // namespace ds
}  // namespace mt_kahypar

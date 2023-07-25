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
    // ! Spin lock to syncronize contractions
    SpinLock sync;
  };

 public:
  FixedVertexSupport() :
    _num_nodes(0),
    _k(kInvalidPartition),
    _hg(nullptr),
    _total_fixed_vertex_weight(0),
    _fixed_vertex_block_weights(),
    _fixed_vertex_data() { }

  FixedVertexSupport(const HypernodeID num_nodes,
                     const PartitionID k) :
    _num_nodes(num_nodes),
    _k(k),
    _hg(nullptr),
    _total_fixed_vertex_weight(0),
    _fixed_vertex_block_weights(k, CAtomic<HypernodeWeight>(0) ),
    _fixed_vertex_data(num_nodes, FixedVertexData { kInvalidPartition, 0, SpinLock() }) { }

  FixedVertexSupport(const FixedVertexSupport&) = delete;
  FixedVertexSupport & operator= (const FixedVertexSupport &) = delete;

  FixedVertexSupport(FixedVertexSupport&&) = default;
  FixedVertexSupport & operator= (FixedVertexSupport &&) = default;

  void setHypergraph(const Hypergraph* hg) {
    _hg = hg;
  }

  // ####################### Fixed Vertex Block Weights #######################

  bool hasFixedVertices() const {
    return _total_fixed_vertex_weight.load(std::memory_order_relaxed) > 0;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeWeight totalFixedVertexWeight() const {
    return _total_fixed_vertex_weight.load(std::memory_order_relaxed);
  }

  // ! Returns the weight of all fixed vertices assigned to the corresponding block
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeWeight fixedVertexBlockWeight(const PartitionID block) const {
    ASSERT(block != kInvalidPartition && block < _k);
    return _fixed_vertex_block_weights[block].load(std::memory_order_relaxed);
  }

  // ####################### Fixed Vertex Information #######################

  // ! Fixes a node to a block
  void fixToBlock(const HypernodeID hn, const PartitionID block) {
    ASSERT(_hg);
    ASSERT(hn < _num_nodes);
    ASSERT(block != kInvalidPartition && block < _k);
    ASSERT(_fixed_vertex_data[hn].block == kInvalidPartition,
      "Hypernode" << hn << "already fixed to a block");
    _fixed_vertex_data[hn].block = block;
    _fixed_vertex_data[hn].fixed_vertex_contraction_cnt = 1;
    _fixed_vertex_block_weights[block].fetch_add(
      _hg->nodeWeight(hn), std::memory_order_relaxed);
    _total_fixed_vertex_weight.fetch_add(
      _hg->nodeWeight(hn), std::memory_order_relaxed);
  }

  // ! Returns whether or not the node is fixed to a block
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool isFixed(const HypernodeID hn) const {
    ASSERT(hn < _num_nodes);
    return fixedVertexBlock(hn) != kInvalidPartition;
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
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool contract(const HypernodeID u, const HypernodeID v) {
    ASSERT(_hg);
    ASSERT(u < _num_nodes && v < _num_nodes);
    bool success = true;
    bool u_becomes_fixed = false;
    bool v_becomes_fixed = false;
    _fixed_vertex_data[v].sync.lock();
    _fixed_vertex_data[u].sync.lock();
    // Note that this does not produce a deadlock since our clustering and
    // contraction algorithm ensures that there are no cyclic dependencies
    const bool is_fixed_u = isFixed(u);
    const bool is_fixed_v = isFixed(v);
    const bool both_fixed = is_fixed_u && is_fixed_v;
    if ( !is_fixed_u && is_fixed_v ) {
      // u becomes a fixed vertex
      _fixed_vertex_data[u].block = fixedVertexBlock(v);
      _fixed_vertex_data[u].fixed_vertex_contraction_cnt = 1;
      u_becomes_fixed = true;
    } else if ( is_fixed_u && ! is_fixed_v ) {
      // v is not fixed, but it is contracted onto a fixed vertex
      // => we have to add the node weight of v to the total and
      // fixed vertex block weight of u
      _fixed_vertex_data[v].block = fixedVertexBlock(u);
      v_becomes_fixed = true;
    } else if ( both_fixed ) {
      if ( fixedVertexBlock(u) == fixedVertexBlock(v) ) {
        ++_fixed_vertex_data[u].fixed_vertex_contraction_cnt;
      } else {
        // Both nodes are fixed vertices, but are assigned to different blocks
        // => contraction is not allowed
        success = false;
      }
    }
    _fixed_vertex_data[u].sync.unlock();
    _fixed_vertex_data[v].sync.unlock();

    if ( u_becomes_fixed || v_becomes_fixed ) {
      // u becomes a fixed vertex or v is contracted onto a fixed vertex
      // => add node weight to total and block weight
      const PartitionID fixed_vertex_block_of_u = fixedVertexBlock(u);
      const HypernodeWeight weight_of_u =
        u_becomes_fixed * _hg->nodeWeight(u) + v_becomes_fixed * _hg->nodeWeight(v);
      _fixed_vertex_block_weights[fixed_vertex_block_of_u].fetch_add(
        weight_of_u, std::memory_order_relaxed);
      _total_fixed_vertex_weight.fetch_add(
        weight_of_u, std::memory_order_relaxed);
    }
    return success;
  }

  // ! Uncontract v from u. This reverts the corresponding contraction operation of v onto u.
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void uncontract(const HypernodeID u, const HypernodeID v) {
    ASSERT(_hg);
    ASSERT(u < _num_nodes && v < _num_nodes);
    if ( isFixed(u) && isFixed(v) ) {
      if ( _fixed_vertex_data[v].fixed_vertex_contraction_cnt > 0 ) {
        // v was fixed before the contraction
        _fixed_vertex_data[u].sync.lock();
        HypernodeID& current_contraction_cnt_of_u = _fixed_vertex_data[u].fixed_vertex_contraction_cnt;
        const bool was_u_fixed_before = current_contraction_cnt_of_u > 0;
        const HypernodeID contraction_cnt_of_u_after = was_u_fixed_before ?
          --current_contraction_cnt_of_u : current_contraction_cnt_of_u;
        _fixed_vertex_data[u].sync.unlock();
        if ( contraction_cnt_of_u_after == 0 ) {
          // u was not fixed before the contraction
          const PartitionID fixed_vertex_block_of_u = _fixed_vertex_data[u].block;
          const HypernodeWeight weight_of_u = _hg->nodeWeight(u);
          _fixed_vertex_block_weights[fixed_vertex_block_of_u].fetch_sub(
            weight_of_u, std::memory_order_relaxed);
          _total_fixed_vertex_weight.fetch_sub(
            weight_of_u, std::memory_order_relaxed);
          // Make u a not fixed vertex again
          _fixed_vertex_data[u].block = kInvalidPartition;
        }
      } else {
        // v was not fixed before the contraction
        const PartitionID fixed_vertex_block_of_v = _fixed_vertex_data[v].block;
        const HypernodeWeight weight_of_v = _hg->nodeWeight(v);
        _fixed_vertex_block_weights[fixed_vertex_block_of_v].fetch_sub(
          weight_of_v, std::memory_order_relaxed);
        _total_fixed_vertex_weight.fetch_sub(
          weight_of_v, std::memory_order_relaxed);
        // Make v a not fixed vertex again
        _fixed_vertex_data[v].block = kInvalidPartition;
      }
    }

  }

 private:
  // ! Number of nodes
  HypernodeID _num_nodes;

  // ! Number of blocks
  PartitionID _k;

  // ! Underlying hypergraph
  const Hypergraph* _hg;

  // ! Total weight of all fixed vertices
  CAtomic<HypernodeWeight> _total_fixed_vertex_weight;

  // ! Weight of all vertices fixed to a block
  vec< CAtomic<HypernodeWeight> > _fixed_vertex_block_weights;

  // ! Fixed vertex block IDs of each node
  vec<FixedVertexData> _fixed_vertex_data;
};

}  // namespace ds
}  // namespace mt_kahypar

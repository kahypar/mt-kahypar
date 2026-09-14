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

#include "mt-kahypar/datastructures/fixed_vertex_support.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/macros.h"

namespace mt_kahypar {
namespace ds {

FixedVertexSupport::FixedVertexSupport() :
  _num_nodes(0),
  _k(kInvalidPartition),
  _total_fixed_vertex_weight(0),
  _fixed_vertex_block_weights(),
  _max_block_weights(),
  _fixed_vertex_data() { }

FixedVertexSupport::FixedVertexSupport(const HypernodeID num_nodes, const PartitionID k) :
  _num_nodes(num_nodes),
  _k(k),
  _total_fixed_vertex_weight(0),
  _fixed_vertex_block_weights(k, CAtomic<HypernodeWeight>(0) ),
  _max_block_weights(k, std::numeric_limits<HypernodeWeight>::max()),
  _fixed_vertex_data(num_nodes, FixedVertexData { kInvalidPartition, 0, 0, SpinLock() }) { }

void FixedVertexSupport::setMaxBlockWeight(const std::vector<HypernodeWeight>& max_block_weights) {
  if ( hasFixedVertices() ) {
    ASSERT(max_block_weights.size() >= static_cast<size_t>(_k));
    _max_block_weights = max_block_weights;
  }
}

template<class Hypergraph>
void FixedVertexSupport::fixToBlock(const Hypergraph& hg, const HypernodeID hn, const PartitionID block) {
  ASSERT(hn < _num_nodes);
  ASSERT(block != kInvalidPartition && block < _k);
  PartitionID expected = kInvalidPartition;
  PartitionID desired = block;
  if (std::atomic_ref(_fixed_vertex_data[hn].block)
      .compare_exchange_strong(expected, desired, std::memory_order::acq_rel, std::memory_order::relaxed)) {
    const HypernodeWeight weight_of_hn = hg.nodeWeight(hn);
    _fixed_vertex_data[hn].fixed_vertex_contraction_cnt = 1;
    _fixed_vertex_data[hn].fixed_vertex_weight = weight_of_hn;
    _fixed_vertex_block_weights[block].fetch_add(
      weight_of_hn, std::memory_order_relaxed);
    _total_fixed_vertex_weight.fetch_add(
      weight_of_hn, std::memory_order_relaxed);
  } else {
    ASSERT(_fixed_vertex_data[hn].block == block,
      "Try to fix hypernode" << hn << "to block" << block
      << ", but it is already fixed to block" << _fixed_vertex_data[hn].block);
  }
}

template<typename Hypergraph>
bool FixedVertexSupport::contract(const Hypergraph& hg, const HypernodeID u, const HypernodeID v) {
  return contractImpl(hg, u, v, false);
}

template<typename Hypergraph>
bool FixedVertexSupport::contractWithoutChains(const Hypergraph& hg, const HypernodeID u, const HypernodeID v) {
  return contractImpl(hg, u, v, true);
}

template<typename Hypergraph>
bool FixedVertexSupport::contractImpl(const Hypergraph& hg, const HypernodeID u, const HypernodeID v, bool ignore_v) {
  ASSERT(u < _num_nodes && v < _num_nodes);
  bool success = true;
  bool u_becomes_fixed = false;
  bool v_becomes_fixed = false;
  const bool is_fixed_v = isFixed(v);
  const HypernodeWeight weight_of_u = hg.nodeWeight(u);
  const HypernodeWeight weight_of_v = hg.nodeWeight(v);
  PartitionID fixed_vertex_block = kInvalidPartition;
  _fixed_vertex_data[u].sync.lock();
  // If we contract a node v onto another node u, all contractions onto v are completed
  // => we therefore do not have to lock v
  const bool is_fixed_u = isFixed(u);
  const bool both_fixed = is_fixed_u && is_fixed_v;
  if ( !is_fixed_u && is_fixed_v ) {
    // u becomes a fixed vertex since v is a fixed vertex
    fixed_vertex_block = fixedVertexBlock(v);
    u_becomes_fixed = true;
  } else if ( is_fixed_u && !is_fixed_v ) {
    // v becomes a fixed vertex since it is contracted onto a fixed vertex
    fixed_vertex_block = fixedVertexBlock(u);
    v_becomes_fixed = true;
  } else if ( both_fixed ) {
    if ( fixedVertexBlock(u) == fixedVertexBlock(v) ) {
      ASSERT(_fixed_vertex_data[u].fixed_vertex_contraction_cnt > 0);
      ASSERT(_fixed_vertex_data[v].fixed_vertex_contraction_cnt > 0);
      ++_fixed_vertex_data[u].fixed_vertex_contraction_cnt;
    } else {
      // Both nodes are fixed vertices, but are assigned to different blocks
      // => contraction is not allowed
      success = false;
    }
  }

  if ( success && ( u_becomes_fixed || v_becomes_fixed ) ) {
    ASSERT(fixed_vertex_block != kInvalidPartition);
    ASSERT(!(u_becomes_fixed && v_becomes_fixed));
    // Either u or v becomes a fixed vertex. Therefore, the fixed vertex block weight changes.
    // To guarantee that we find a feasible initial partition, we ensure that the new block weight
    // is smaller than the maximum allowed block weight.
    const HypernodeWeight delta_weight =
      u_becomes_fixed * weight_of_u + v_becomes_fixed * weight_of_v;
    const HypernodeWeight block_weight_after =
      _fixed_vertex_block_weights[fixed_vertex_block].add_fetch(
        delta_weight, std::memory_order_relaxed);
    if ( likely( block_weight_after <= _max_block_weights[fixed_vertex_block] ) ) {
      _total_fixed_vertex_weight.fetch_add(delta_weight, std::memory_order_relaxed);
      if ( u_becomes_fixed ) {
        ASSERT(isFixed(v));
        ASSERT(_fixed_vertex_data[u].fixed_vertex_contraction_cnt == 0);
        // Block weight update was successful => set fixed vertex block of u
        _fixed_vertex_data[u].block = fixedVertexBlock(v);
        _fixed_vertex_data[u].fixed_vertex_contraction_cnt = 1;
        _fixed_vertex_data[u].fixed_vertex_weight = weight_of_u;
      }
    } else {
      // The new fixed vertex block weight is larger than the maximum allowed bock weight
      // => revert block weight update and forbid contraction
      _fixed_vertex_block_weights[fixed_vertex_block].sub_fetch(
        delta_weight, std::memory_order_relaxed);
      v_becomes_fixed = false;
      success = false;
    }
  }
  _fixed_vertex_data[u].sync.unlock();

  if ( !ignore_v && v_becomes_fixed ) {
    // Our contraction algorithm ensures that there are no concurrent contractions onto v
    // if v is contracted onto another node. We therefore can set the fixed vertex block of
    // v outside the lock
    _fixed_vertex_data[v].block = fixed_vertex_block;
    _fixed_vertex_data[v].fixed_vertex_weight = weight_of_v;
  }
  return success;
}

void FixedVertexSupport::uncontract(const HypernodeID u, const HypernodeID v) {
  ASSERT(u < _num_nodes && v < _num_nodes);
  if ( isFixed(v) ) {
    if ( _fixed_vertex_data[v].fixed_vertex_contraction_cnt > 0 ) {
      // v was fixed before the contraction
      _fixed_vertex_data[u].sync.lock();
      ASSERT(_fixed_vertex_data[u].fixed_vertex_contraction_cnt > 0);
      const HypernodeID contraction_cnt_of_u_after =
        --_fixed_vertex_data[u].fixed_vertex_contraction_cnt;
      _fixed_vertex_data[u].sync.unlock();
      if ( contraction_cnt_of_u_after == 0 ) {
        // u was not fixed before the contraction
        const PartitionID fixed_vertex_block_of_u = _fixed_vertex_data[u].block;
        const HypernodeWeight weight_of_u = _fixed_vertex_data[u].fixed_vertex_weight;
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
      const HypernodeWeight weight_of_v = _fixed_vertex_data[v].fixed_vertex_weight;
      _fixed_vertex_block_weights[fixed_vertex_block_of_v].fetch_sub(
        weight_of_v, std::memory_order_relaxed);
      _total_fixed_vertex_weight.fetch_sub(
        weight_of_v, std::memory_order_relaxed);
      // Make v a not fixed vertex again
      _fixed_vertex_data[v].block = kInvalidPartition;
    }
  }
}

FixedVertexSupport FixedVertexSupport::copy() const {
  FixedVertexSupport cpy;
  cpy._num_nodes = _num_nodes;
  cpy._k = _k;
  cpy._total_fixed_vertex_weight = _total_fixed_vertex_weight;
  cpy._fixed_vertex_block_weights = _fixed_vertex_block_weights;
  cpy._max_block_weights = _max_block_weights;
  cpy._fixed_vertex_data = _fixed_vertex_data;
  return cpy;
}

template<typename Hypergraph>
bool FixedVertexSupport::verifyClustering(const Hypergraph& hg, const vec<HypernodeID>& cluster_ids) const {
  vec<PartitionID> fixed_vertex_blocks(hg.initialNumNodes(), kInvalidPartition);
  for ( const HypernodeID& hn : hg.nodes() ) {
    if ( hg.isFixed(hn) ) {
      if ( fixed_vertex_blocks[cluster_ids[hn]] != kInvalidPartition &&
            fixed_vertex_blocks[cluster_ids[hn]] != hg.fixedVertexBlock(hn)) {
        LOG << "There are two nodes assigned to same cluster that belong to different fixed vertex blocks";
        return false;
      }
      fixed_vertex_blocks[cluster_ids[hn]] = hg.fixedVertexBlock(hn);
    }
  }

  vec<HypernodeWeight> expected_block_weights(_k, 0);
  for ( const HypernodeID& hn : hg.nodes() ) {
    if ( fixed_vertex_blocks[cluster_ids[hn]] != kInvalidPartition ) {
      if ( !isFixed(cluster_ids[hn]) ) {
        LOG << "Cluster" << cluster_ids[hn] << "should be fixed to block"
            << fixed_vertex_blocks[cluster_ids[hn]];
        return false;
      }
      expected_block_weights[fixed_vertex_blocks[cluster_ids[hn]]] += hg.nodeWeight(hn);
    }
  }

  for ( PartitionID block = 0; block < _k; ++block ) {
    if ( fixedVertexBlockWeight(block) != expected_block_weights[block] ) {
      LOG << "Fixed vertex block" << block << "should have weight" << expected_block_weights[block]
          << ", but it is" << fixedVertexBlockWeight(block);
      return false;
    }
  }
  return true;
}


namespace {
  #define FIX_TO_BLOCK(X) void FixedVertexSupport::fixToBlock(const X& hg, const HypernodeID hn, const PartitionID block);
  #define CONTRACT(X) bool FixedVertexSupport::contract(const X& hg, const HypernodeID u, const HypernodeID v)
  #define CONTRACT_WITHOUT_CHAINS(X) bool FixedVertexSupport::contractWithoutChains(const X& hg, const HypernodeID u, const HypernodeID v)
  #define CONTRACT_IMPL(X) bool FixedVertexSupport::contractImpl(const X& hg, const HypernodeID u, const HypernodeID v, bool ignore_v)
  #define VERIFY_CLUSTERING(X) bool FixedVertexSupport::verifyClustering(const X& hg, const vec<HypernodeID>& cluster_ids) const
}

INSTANTIATE_FUNC_WITH_HYPERGRAPHS(FIX_TO_BLOCK)
INSTANTIATE_FUNC_WITH_HYPERGRAPHS(CONTRACT)
INSTANTIATE_FUNC_WITH_HYPERGRAPHS(CONTRACT_WITHOUT_CHAINS)
INSTANTIATE_FUNC_WITH_HYPERGRAPHS(CONTRACT_IMPL)
INSTANTIATE_FUNC_WITH_HYPERGRAPHS(VERIFY_CLUSTERING)

} // namespace ds
} // namespace mt_kahypar

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

#include "mt-kahypar/datastructures/dynamic_graph.h"
#include "mt-kahypar/datastructures/dynamic_graph_factory.h"
#include "mt-kahypar/macros.h"

namespace mt_kahypar {
namespace ds {

template<typename Hypergraph>
FixedVertexSupport<Hypergraph>::FixedVertexSupport() :
  _num_nodes(0),
  _k(kInvalidPartition),
  _hg(nullptr),
  _total_fixed_vertex_weight(0),
  _fixed_vertex_block_weights(),
  _max_block_weights(),
  _fixed_vertex_data(),
  _constraint_graph(nullptr),
  _hypergraph_id_to_graph_id(
        std::make_unique<ds::FixedSizeSparseMap<HypernodeID, HypernodeID>>(0)
    ) { }

template<typename Hypergraph>
FixedVertexSupport<Hypergraph>::FixedVertexSupport(const HypernodeID num_nodes,
                                                   const PartitionID k) :
  _num_nodes(num_nodes),
  _k(k),
  _hg(nullptr),
  _total_fixed_vertex_weight(0),
  _fixed_vertex_block_weights(k, CAtomic<HypernodeWeight>(0) ),
  _max_block_weights(k, std::numeric_limits<HypernodeWeight>::max()),
  _fixed_vertex_data(num_nodes, FixedVertexData { kInvalidPartition, 0, 0, SpinLock() }),
  _constraint_graph(nullptr),
  _hypergraph_id_to_graph_id(
        std::make_unique<ds::FixedSizeSparseMap<HypernodeID, HypernodeID>>(0)
    ) { }

// the following definitions are necessary to avoid issues with the unique_ptr deleter
template<typename Hypergraph>
FixedVertexSupport<Hypergraph>::FixedVertexSupport(FixedVertexSupport&&) = default;
template<typename Hypergraph>
FixedVertexSupport<Hypergraph>& FixedVertexSupport<Hypergraph>::operator=(FixedVertexSupport<Hypergraph> &&) = default;
template<typename Hypergraph>
FixedVertexSupport<Hypergraph>::~FixedVertexSupport() = default;


template<typename Hypergraph>
bool FixedVertexSupport<Hypergraph>::contract(const HypernodeID u, const HypernodeID v) {
  return contractImpl(u, v, false);
}

template<typename Hypergraph>
bool FixedVertexSupport<Hypergraph>::contractWithoutChains(const HypernodeID u, const HypernodeID v) {
  return contractImpl(u, v, true);
}

template<typename Hypergraph>
bool FixedVertexSupport<Hypergraph>::contractImpl(const HypernodeID u, const HypernodeID v, bool ignore_v) {
  ASSERT(_hg);
  ASSERT(u < _num_nodes && v < _num_nodes);
  bool success = true;
  bool u_becomes_fixed = false;
  bool v_becomes_fixed = false;
  const bool is_fixed_v = isFixed(v);
  const HypernodeWeight weight_of_u = _hg->nodeWeight(u);
  const HypernodeWeight weight_of_v = _hg->nodeWeight(v);
  PartitionID fixed_vertex_block = kInvalidPartition;

  // If we contract a node v onto another node u, all contractions onto v are completed
  // => we therefore do not have to lock v
  const bool is_fixed_u = isFixed(u);
  const bool both_fixed = is_fixed_u && is_fixed_v;
  if ( !is_fixed_u && is_fixed_v ) {
    _fixed_vertex_data[u].sync.lock();
    // u becomes a fixed vertex since v is a fixed vertex
    fixed_vertex_block = fixedVertexBlock(v);
    u_becomes_fixed = true;
  } else if ( is_fixed_u && !is_fixed_v ) {
    _fixed_vertex_data[u].sync.lock();
    // v becomes a fixed vertex since it is contracted onto a fixed vertex
    fixed_vertex_block = fixedVertexBlock(u);
    v_becomes_fixed = true;
  } else if ( both_fixed ) {
    _fixed_vertex_data[u].sync.lock();
    if ( fixedVertexBlock(u) == fixedVertexBlock(v) ) {
      ASSERT(_fixed_vertex_data[u].fixed_vertex_contraction_cnt > 0);
      ASSERT(_fixed_vertex_data[v].fixed_vertex_contraction_cnt > 0);
      ++_fixed_vertex_data[u].fixed_vertex_contraction_cnt;
    } else {
      // Both nodes are fixed vertices, but are assigned to different blocks
      // => contraction is not allowed
      success = false;
    }
  } else if (hasNegativeConstraints()) {
    HypernodeID u1 = u;
    HypernodeID v1 = v;
    // if v is in constraints and u not, dont contract it in u for now.
    // TODO: contract them but change all occurencies of v to u in constraints
    if (getConstraintIdFromHypergraphId(v, v1) && !getConstraintIdFromHypergraphId(u, u1)) {
      success = false;
    } else if (getConstraintIdFromHypergraphId(u, u1) && getConstraintIdFromHypergraphId(v, v1)) {
      // if anti constraints should be checked this was wrong locked
      bool locked = false;
      while (!locked) {
        if (u1 == v1) {
          success = false;
          break;
        } else if(u1 < v1) {
          // these are the real locks of nodes in fixed vertex support
          _fixed_vertex_data[u1].sync.lock();
          _fixed_vertex_data[v1].sync.lock();
        } else {
          _fixed_vertex_data[v1].sync.lock();
          _fixed_vertex_data[u1].sync.lock();
        }
        HypernodeID u2 = invalidNode;
        HypernodeID v2 = invalidNode;
        getConstraintIdFromHypergraphId(u, u2);
        getConstraintIdFromHypergraphId(v, v2);
        if (u1 == u2 && v1 == v2) {
          // no contraction happened between getting nodes and locking them
          locked = true;
        } else {
          // contraction happend, get new nodes and lock them
          _fixed_vertex_data[u1].sync.unlock();
          _fixed_vertex_data[v1].sync.unlock();
          u1 = u2;
          v1 = v2;
        }
      }
      if(success) {
        // if both nodes are in the constraint graph and no neighbors
        if (constraintExistsForPair(u,v)) {
          success = false;
        } else {
          success = _constraint_graph->registerContraction(u1, v1);
          if (success) {
            // node weight to node id in hypergraph mapping gets broken here
            // because new weight of u1 is u1 + v1
            size_t num_contractions = _constraint_graph->contract(v1, std::numeric_limits<HypernodeWeight>::max());
            if (num_contractions == 0) {
              success = false;
              ASSERT(false, "could not contract nodes" << u1 << v1);
            } else {
              // v1 gets contracted in u1
              // both hg nodes point to the same new contracted node u1
              _hypergraph_id_to_graph_id->operator[](v) = u1;
            }
          } else {
            ASSERT(false, "could not register contraction beteween nodes" << u1 << v1);
          }
        }
        _fixed_vertex_data[u1].sync.unlock();
        _fixed_vertex_data[v1].sync.unlock();
      }
    }
    return success;
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

template<typename Hypergraph>
void FixedVertexSupport<Hypergraph>::uncontract(const HypernodeID u, const HypernodeID v) {
  ASSERT(_hg);
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

vec<std::pair<HypernodeID, HypernodeID>> trasform_node_vector(const vec<std::pair<HypernodeID,HypernodeID>>& node_vector,
                                                          std::unique_ptr<ds::FixedSizeSparseMap<HypernodeID, HypernodeID>>& hypergraph_id_to_graph_id,
                                                          vec<HypernodeWeight>& node_weight, 
                                                          HypernodeID& num_nodes) {
  /**
   * Transforms random nodeIDs to ongoing ids for and graph.
   * Get Hypernode count and set the nodeWeights accordingly to the read nodeIDs
   */
  HypernodeID node_count = 0;
  vec<std::pair<HypernodeID, HypernodeID>> new_node_vector;
  new_node_vector.reserve(node_vector.size());
  for (const auto& node_pair : node_vector) {
    if (!hypergraph_id_to_graph_id->contains(node_pair.first)) {
      hypergraph_id_to_graph_id->operator[](node_pair.first) = node_count;
      node_weight.push_back(node_pair.first);
      node_count++;
    }
    if (!hypergraph_id_to_graph_id->contains(node_pair.second)) {
      hypergraph_id_to_graph_id->operator[](node_pair.second) = node_count;
      node_weight.push_back(node_pair.second);
      node_count++;
    }
    new_node_vector.push_back(std::make_pair( hypergraph_id_to_graph_id->get(node_pair.first), 
                                              hypergraph_id_to_graph_id->get(node_pair.second)));
  }
  num_nodes = node_count;
  return new_node_vector;
}

template<typename Hypergraph>
void FixedVertexSupport<Hypergraph>::setNegativeConstraints(const vec<std::pair<HypernodeID, HypernodeID>>& constraints) {
  _constraints = constraints;
  _hypergraph_id_to_graph_id->setMaxSize(constraints.size() * 3);
  vec<HypernodeWeight> node_weight;
  HypernodeID num_nodes;
  vec<std::pair<HypernodeID, HypernodeID>> transformed_constraints = trasform_node_vector(constraints, 
                                                                                          _hypergraph_id_to_graph_id,
                                                                                          node_weight, 
                                                                                          num_nodes);
  vec<HyperedgeWeight> edge_weight(transformed_constraints.size(), HyperedgeWeight(1));
  _constraint_graph = std::make_unique<DynamicGraph>(
    DynamicGraphFactory::construct_from_graph_edges(num_nodes,
                                                    transformed_constraints.size(), 
                                                    transformed_constraints,
                                                    edge_weight.data(),
                                                    node_weight.data(),
                                                    true));
}

template<typename Hypergraph>
void FixedVertexSupport<Hypergraph>::setNegativeConstraints(const FixedVertexSupport<Hypergraph>& fixed_vertices) {
    _constraint_graph = std::make_unique<DynamicGraph>(fixed_vertices._constraint_graph->copy());
    _constraints = fixed_vertices._constraints;
    _hypergraph_id_to_graph_id = std::make_unique<ds::FixedSizeSparseMap<HypernodeID, HypernodeID>>(fixed_vertices._hypergraph_id_to_graph_id->copy());
  }

template<typename Hypergraph>
bool FixedVertexSupport<Hypergraph>::constraintExistsForPair(const HypernodeID u, const HypernodeID v) const {
    HypernodeID u1;
    HypernodeID v1;
    return (getConstraintIdFromHypergraphId(u, u1) && getConstraintIdFromHypergraphId(v, v1) && _constraint_graph->isIncidentTo(u1, v1));
  }

template<typename Hypergraph>
FixedVertexSupport<Hypergraph> FixedVertexSupport<Hypergraph>::copy() const {
  FixedVertexSupport<Hypergraph> cpy;
  cpy._num_nodes = _num_nodes;
  cpy._k = _k;
  cpy._hg = _hg;
  cpy._total_fixed_vertex_weight = _total_fixed_vertex_weight;
  cpy._fixed_vertex_block_weights = _fixed_vertex_block_weights;
  cpy._max_block_weights = _max_block_weights;
  cpy._fixed_vertex_data = _fixed_vertex_data;
  if (_constraint_graph != nullptr) {
    cpy._constraint_graph = std::make_unique<DynamicGraph>(_constraint_graph->copy());
    cpy._constraints = _constraints;
    cpy._hypergraph_id_to_graph_id = std::make_unique<ds::FixedSizeSparseMap<HypernodeID, HypernodeID>>(_hypergraph_id_to_graph_id->copy());
  }
  return cpy;
}

} // namespace ds
} // namespace mt_kahypar

#include "mt-kahypar/datastructures/static_graph.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/dynamic_graph.h"
#include "mt-kahypar/datastructures/dynamic_hypergraph.h"

template class mt_kahypar::ds::FixedVertexSupport<mt_kahypar::ds::StaticHypergraph>;
template class mt_kahypar::ds::FixedVertexSupport<mt_kahypar::ds::StaticGraph>;
template class mt_kahypar::ds::FixedVertexSupport<mt_kahypar::ds::DynamicHypergraph>;
template class mt_kahypar::ds::FixedVertexSupport<mt_kahypar::ds::DynamicGraph>;
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

#pragma once

#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace ds {

/*!
 * Special hypergraph data structure that allows to concurrently modifiy hypergraph.
 * Following operations are supported:
 *  1.) Contract Vertices
 *  2.) Remove Hyperedges
 *  3.) Replace Hyperedges
 * Note, that several vertices and hyperedges can be modified concurrently, but not the
 * same vertex or hyperedge.
 */
template<typename Hypergraph,
         typename HypergraphFactory>
class SparsifierHypergraph {

  using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
  using AtomicNodeDegreeVector = parallel::scalable_vector<parallel::IntegralAtomicWrapper<HyperedgeID>>;
  using PinIterator = parallel::scalable_vector<HypernodeID>::const_iterator;

 public:
  SparsifierHypergraph(const Hypergraph& hypergraph) :
    _num_nodes(hypergraph.initialNumNodes()),
    _num_removed_nodes(0),
    _num_edges(hypergraph.initialNumEdges()),
    _edge_vector(),
    _hyperedge_weight(),
    _he_included_in_sparsified_hg(),
    _mapping(),
    _hypernode_weight(),
    _node_degrees(),
    _hn_included_in_sparsified_hg() {
    initialize(hypergraph);
  }

  SparsifierHypergraph(const HypernodeID num_hypernodes,
                       const HyperedgeID num_hyperedges) :
    _num_nodes(num_hypernodes),
    _num_removed_nodes(0),
    _num_edges(num_hyperedges),
    _edge_vector(),
    _hyperedge_weight(),
    _he_included_in_sparsified_hg(),
    _mapping(),
    _hypernode_weight(),
    _node_degrees(),
    _hn_included_in_sparsified_hg() {
    tbb::parallel_invoke([&] {
      _edge_vector.resize(num_hyperedges);
    }, [&] {
      _hyperedge_weight.resize(num_hyperedges);
    }, [&] {
      _hypernode_weight.resize(num_hypernodes);
    });
  }

  ~SparsifierHypergraph() {
    tbb::parallel_invoke([&] {
      parallel::parallel_free(_edge_vector);
    }, [&] {
      parallel::parallel_free(_hyperedge_weight,
        _he_included_in_sparsified_hg,
        _hypernode_weight, _node_degrees,
        _hn_included_in_sparsified_hg);
    });
  }

  // ####################### General Hypergraph Stats #######################

  HypernodeID numNodes() const {
    return _num_nodes;
  }

  HypernodeID numRemovedNodes() {
    return _num_removed_nodes.combine(std::plus<HypernodeID>());
  }

  HyperedgeID numEdges() const {
    return _num_edges;
  }

  // ####################### Iterators #######################

  const parallel::scalable_vector<HypernodeID>& pins(const HyperedgeID e) const {
    ASSERT(e < _num_edges);
    return _edge_vector[e];
  }

  // ####################### Hypernode Information #######################

  HypernodeWeight nodeWeight(const HypernodeID u) const {
    ASSERT(u < _num_nodes, V(u) << V(_num_nodes));
    return _hypernode_weight[u];
  }

  HyperedgeID nodeDegree(const HypernodeID u) const {
    ASSERT(u < _num_nodes);
    return _node_degrees[u];
  }

  bool nodeIsEnabled(const HypernodeID u) const {
    ASSERT(u < _num_nodes);
    return _hn_included_in_sparsified_hg[u];
  }

  // ####################### Hyperedge Information #######################

  HyperedgeWeight edgeWeight(const HyperedgeID e) const {
    ASSERT(e < _num_edges);
    return _hyperedge_weight[e];
  }

  bool edgeIsEnabled(const HyperedgeID e) const {
    ASSERT(e < _num_edges);
    return _he_included_in_sparsified_hg[e];
  }

  // ####################### Modification Functions #######################

  // ! Contract two vertices
  // ! Note, that contraction is not applied to hyperedges. The contraction
  // ! partner v is still present in its incident hyperedges. Contraction
  // ! is applied once function sparsify() is called. However, keep in mind
  // ! that parallel hyperedges and single-pin hyperedges are not removed.
  // ! Furthermore, contracting two already contracted vertices results in
  // ! undefined behaviour.
  void contract(const HypernodeID u, const HypernodeID v) {
    ASSERT(u != v);
    ASSERT(u < _num_nodes);
    ASSERT(v < _num_nodes);
    _mapping[v] = u;
    _hypernode_weight[u] += _hypernode_weight[v];
    _hypernode_weight[v] = 0;
    _hn_included_in_sparsified_hg[v] = 0;
    ++_num_removed_nodes.local();
  }

  void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
    ASSERT(e < _num_edges);
    _hyperedge_weight[e] = weight;
  }

  void remove(const HyperedgeID e) {
    ASSERT(e < _num_edges);
    for ( const HypernodeID& pin : _edge_vector[e] ) {
      --_node_degrees[pin];
    }
    _edge_vector[e].clear();
    _hyperedge_weight[e] = 0;
    _he_included_in_sparsified_hg[e] = 0;
  }

  void replace(const HyperedgeID e, parallel::scalable_vector<HyperedgeID>&& hyperedge) {
    ASSERT(e < _num_edges);
    for ( const HypernodeID& pin : _edge_vector[e]) {
      --_node_degrees[pin];
    }
    _edge_vector[e] = std::move(hyperedge);
    for ( const HypernodeID& pin : _edge_vector[e]) {
      ++_node_degrees[pin];
    }
  }

  // ####################### Sparsification Functions #######################

  // ! Sparsifies the hypergraph by applying all modifications.
  // ! Note, class is not useable any more afterwards.
  Hypergraph sparsify() {
    // Compute number of nodes and edges in sparsified hypergraph
    parallel::TBBPrefixSum<HyperedgeID> he_prefix_sum(_he_included_in_sparsified_hg);
    parallel::TBBPrefixSum<HyperedgeID> hn_prefix_sum(_hn_included_in_sparsified_hg);
    tbb::parallel_invoke([&] {
      tbb::parallel_scan(tbb::blocked_range<size_t>(ID(0), _num_edges), he_prefix_sum);
    }, [&] {
      tbb::parallel_scan(tbb::blocked_range<size_t>(ID(0), _num_nodes), hn_prefix_sum);
    });

    // Apply vertex id of sparsified hypergraph to mapping
    tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID u) {
      _mapping[u] = hn_prefix_sum[_mapping[u]];
    });

    // Apply mapping to all enabled hyperedges
    tbb::parallel_for(ID(0), _num_edges, [&](const HyperedgeID& e) {
      if ( edgeIsEnabled(e) ) {
        parallel::scalable_vector<HypernodeID>& hyperedge = _edge_vector[e];
        for ( HypernodeID& pin : hyperedge ) {
          ASSERT(pin < _mapping.size());
          pin = _mapping[pin];
        }
        std::sort(hyperedge.begin(), hyperedge.end());
        hyperedge.erase(std::unique(hyperedge.begin(), hyperedge.end()), hyperedge.end());
      }
    });

    const HypernodeID num_hypernodes = hn_prefix_sum.total_sum();
    const HyperedgeID num_hyperedges = he_prefix_sum.total_sum();
    SparsifierHypergraph sparsified_hypergraph(num_hypernodes, num_hyperedges);
    tbb::parallel_invoke([&] {
      tbb::parallel_for(ID(0), _num_edges, [&](const HyperedgeID e) {
        if ( he_prefix_sum.value(e) ) {
          const HyperedgeID sparsified_id = he_prefix_sum[e];
          sparsified_hypergraph._edge_vector[sparsified_id] = std::move(_edge_vector[e]);
          sparsified_hypergraph._hyperedge_weight[sparsified_id] = _hyperedge_weight[e];
        }
      });
    }, [&] {
      tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID u) {
        if ( hn_prefix_sum.value(u) ) {
          const HyperedgeID sparsified_id = hn_prefix_sum[u];
          sparsified_hypergraph._hypernode_weight[sparsified_id] = _hypernode_weight[u];
        }
      });
    });

    return HypergraphFactory::construct(
      num_hypernodes, num_hyperedges,
      sparsified_hypergraph._edge_vector,
      sparsified_hypergraph._hyperedge_weight.data(),
      sparsified_hypergraph._hypernode_weight.data());
  }

  // ! Returns a mapping from the sparsified hypergraph
  parallel::scalable_vector<HypernodeID>&& getMapping() {
    return std::move(_mapping);
  }

 private:
  void initialize(const Hypergraph& hypergraph) {
    tbb::parallel_invoke([&] {
      tbb::parallel_invoke([&] {
        _edge_vector.resize(_num_edges);
      }, [&] {
        _hyperedge_weight.resize(_num_edges);
      }, [&] {
        _he_included_in_sparsified_hg.assign(_num_edges, 1);
      }, [&] {
        _node_degrees.assign(_num_nodes, parallel::IntegralAtomicWrapper<HyperedgeID>(0));
      });

      tbb::parallel_for(ID(0), hypergraph.initialNumEdges(), [&](const HyperedgeID he) {
        if ( hypergraph.edgeIsEnabled(he) ) {
          _hyperedge_weight[he] = hypergraph.edgeWeight(he);
          _edge_vector[he].reserve(hypergraph.edgeSize(he));
          for ( const HypernodeID& pin : hypergraph.pins(he) ) {
            _edge_vector[he].push_back(pin);
            ++_node_degrees[pin];
          }
          std::sort(_edge_vector[he].begin(), _edge_vector[he].end());
        } else {
          _he_included_in_sparsified_hg[he] = 0;
        }
      });
    }, [&] {
      tbb::parallel_invoke([&] {
        _mapping.resize(_num_nodes);
      }, [&] {
        _hypernode_weight.resize(_num_nodes);
      }, [&] {
        _hn_included_in_sparsified_hg.assign(_num_nodes, 1);
      });

      tbb::parallel_for(ID(0), hypergraph.initialNumNodes(), [&](const HypernodeID hn) {
        _mapping[hn] = hn;
        if ( hypergraph.nodeIsEnabled(hn) ) {
          _hypernode_weight[hn] = hypergraph.nodeWeight(hn);
        } else {
          _hn_included_in_sparsified_hg[hn] = 0;
          ++_num_removed_nodes.local();
        }
      });
    });
  }

  // Stats
  const HypernodeID _num_nodes;
  tbb::enumerable_thread_specific<HypernodeID> _num_removed_nodes;
  const HyperedgeID _num_edges;

  // Hyperedges
  HyperedgeVector _edge_vector;
  parallel::scalable_vector<HyperedgeWeight> _hyperedge_weight;
  parallel::scalable_vector<HyperedgeID> _he_included_in_sparsified_hg;

  // Hypernodes
  parallel::scalable_vector<HypernodeID> _mapping;
  parallel::scalable_vector<HypernodeWeight> _hypernode_weight;
  AtomicNodeDegreeVector _node_degrees;
  parallel::scalable_vector<HypernodeID> _hn_included_in_sparsified_hg;
};

} // namespace ds
} // namespace mt_kahypar
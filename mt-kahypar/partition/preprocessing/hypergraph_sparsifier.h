/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <vector>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/clustering.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
template<typename TypeTraits>
class HypergraphSparsifierT {
 public:

  using HyperGraph = typename TypeTraits::HyperGraph;
  using PartitionedHyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
  using HyperGraphFactory = typename TypeTraits::HyperGraphFactory;
  using TBB = typename TypeTraits::TBB;
  using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;

  HypergraphSparsifierT(const Context& context,
                        const TaskGroupID task_group_id) :
    _context(context),
    _task_group_id(task_group_id),
    _is_init(false),
    _sparsified_hg(),
    _sparsified_partitioned_hg(),
    _removed_hes(),
    _removed_hns(),
    _mapping() { }

  HypergraphSparsifierT(const HypergraphSparsifierT&) = delete;
  HypergraphSparsifierT & operator= (const HypergraphSparsifierT &) = delete;

  HypergraphSparsifierT(HypergraphSparsifierT&&) = delete;
  HypergraphSparsifierT & operator= (HypergraphSparsifierT &&) = delete;


  // ####################### Sparsification Functions #######################

  bool isInitialized() const {
    return _is_init;
  }

  HyperGraph& sparsifiedHypergraph() {
    ASSERT(_is_init);
    return _sparsified_hg;
  }

  PartitionedHyperGraph& sparsifiedPartitionedHypergraph() {
    ASSERT(_is_init);
    return _sparsified_partitioned_hg;
  }

  void sparsify(const HyperGraph& hypergraph) {
    ASSERT(!_is_init);
    ASSERT(_mapping.size() == 0, "contractDegreeZeroHypernodes was triggered before");
    HyperedgeVector edge_vector;
    parallel::scalable_vector<int> vertices_to_numa_node;
    parallel::scalable_vector<HyperedgeWeight> hyperedge_weight;
    parallel::scalable_vector<HypernodeWeight> hypernode_weight;

    tbb::parallel_invoke([&] {
      _mapping.resize(hypergraph.initialNumNodes());
    }, [&]{
      hypernode_weight.resize(hypergraph.initialNumNodes());
    }, [&] {
      vertices_to_numa_node.resize(hypergraph.initialNumNodes());
    });

    // #################### STAGE 1 ####################
    // Perform Degree-Zero Contractions
    // Degree-Zero hypernodes are contracted to supervertices such that
    // each supervertex has a weight smaller than the maximum allowed
    // node weight.
    utils::Timer::instance().start_timer("degree_zero_contraction", "Degree-Zero Contractions");
    const HypernodeID num_hypernodes = degreeZeroSparsification(
      hypergraph, vertices_to_numa_node, hypernode_weight);
    utils::Timer::instance().stop_timer("degree_zero_contraction");

    // #################### STAGE 2 ####################
    // Heavy Hyperedge Removal
    // If the weight of all pins of a hyperedge is greater than a
    // certain threshold, we remove them from the hypergraph
    utils::Timer::instance().start_timer("heavy_hyperedge_removal", "Heavy HE Removal");
    const HyperedgeID num_hyperedges = heavyHyperedgeRemovalSparsification(
      hypergraph, edge_vector, hyperedge_weight);
    utils::Timer::instance().stop_timer("heavy_hyperedge_removal");

    // #################### STAGE 3 ####################
    // High Degree Vertex Removal
    // Removes vertices that have a degree greater than a certain threshold
    // from all hyperedges. High-degree vertices are removed from hyperedges
    // as long as their are more than two pins left.
    utils::Timer::instance().start_timer("high_degree_vertex_removal", "High-Degree HN Removal");
    highDegreeVertexSparsification(num_hypernodes, edge_vector, hyperedge_weight);
    utils::Timer::instance().stop_timer("high_degree_vertex_removal");

    // #################### STAGE 4 ####################
    // Construct sparsified hypergraph
    utils::Timer::instance().start_timer("construct_sparsified_hypergraph", "Construct Sparsified HG");
    if ( TBB::instance().num_used_numa_nodes() == 1 ) {
      _sparsified_hg = HyperGraphFactory::construct(
        _task_group_id, num_hypernodes, num_hyperedges,
        edge_vector, hyperedge_weight.data(), hypernode_weight.data());
    } else {
      _sparsified_hg = HyperGraphFactory::construct(
        _task_group_id, num_hypernodes, num_hyperedges,
        edge_vector, std::move(vertices_to_numa_node),
        hyperedge_weight.data(), hypernode_weight.data());
    }
    _sparsified_partitioned_hg = PartitionedHyperGraph(
      _context.partition.k, _task_group_id, _sparsified_hg);
    utils::Timer::instance().stop_timer("construct_sparsified_hypergraph");

    tbb::parallel_invoke([&] {
      parallel::parallel_free(edge_vector);
    }, [&] {
      parallel::parallel_free(hyperedge_weight, hypernode_weight);
    });

    _is_init = true;
  }


  // ! Removes single-pin hyperedges since they do not contribute
  // ! to cut or km1 metric.
  HyperedgeID removeSingleNodeHyperedges(HyperGraph& hypergraph) {
    ASSERT(!_is_init);
    HyperedgeID num_removed_single_node_hes = 0;
    for (const HyperedgeID& he : hypergraph.edges()) {
      if (hypergraph.edgeSize(he) == 1) {
        ++num_removed_single_node_hes;
        hypergraph.removeEdge(he);
        _removed_hes.push_back(he);
      }
    }
    return num_removed_single_node_hes;
  }

  // ! Contracts degree-zero vertices to degree-zero supervertices
  // ! We contract sets of degree-zero vertices such that the weight of
  // ! each supervertex is less than or equal than the maximum allowed
  // ! node weight for a vertex during coarsening.
  HypernodeID contractDegreeZeroHypernodes(HyperGraph& hypergraph) {
    ASSERT(!_is_init);
    _mapping.assign(hypergraph.initialNumNodes(), kInvalidHypernode);
    HypernodeID num_removed_degree_zero_hypernodes = 0;
    HypernodeID last_degree_zero_representative = kInvalidHypernode;
    HypernodeWeight last_degree_zero_weight = 0;
    for (const HypernodeID& hn : hypergraph.nodes()) {
      if ( hypergraph.nodeDegree(hn) == 0 ) {
        bool was_removed = false;
        if ( last_degree_zero_representative != kInvalidHypernode ) {
          const HypernodeWeight weight = hypergraph.nodeWeight(hn);
          if ( last_degree_zero_weight + weight <= _context.coarsening.max_allowed_node_weight ) {
            // Remove vertex and aggregate its weight in its represenative supervertex
            ++num_removed_degree_zero_hypernodes;
            hypergraph.removeHypernode(hn);
            _removed_hns.push_back(hn);
            was_removed = true;
            _mapping[hypergraph.originalNodeID(hn)] = last_degree_zero_representative;
            last_degree_zero_weight += weight;
            hypergraph.setNodeWeight(last_degree_zero_representative, last_degree_zero_weight);
          }
        }

        if ( !was_removed ) {
          last_degree_zero_representative = hn;
          last_degree_zero_weight = hypergraph.nodeWeight(hn);
        }
      }
    }
    return num_removed_degree_zero_hypernodes;
  }

  // ####################### Restore Operations #######################

  // ! Restore single-pin hyperedges
  void restoreSingleNodeHyperedges(PartitionedHyperGraph& hypergraph) {
    for (const HyperedgeID& he : _removed_hes) {
      hypergraph.restoreSinglePinHyperedge(he);
    }
  }

  // ! Restore degree-zero vertices
  // ! Each removed degree-zero vertex is assigned to the block of its supervertex.
  void restoreDegreeZeroHypernodes(PartitionedHyperGraph& hypergraph) {
    for ( const HypernodeID& hn : _removed_hns ) {
      const HypernodeID original_id = hypergraph.originalNodeID(hn);
      ASSERT(original_id < _mapping.size());
      const HypernodeID representative = _mapping[original_id];
      ASSERT(representative != kInvalidHypernode);
      // Restore degree-zero vertex and assign it to the block
      // of its supervertex
      hypergraph.enableHypernode(hn);
      hypergraph.setNodeWeight(representative,
        hypergraph.nodeWeight(representative) - hypergraph.nodeWeight(hn));
      hypergraph.setNodePart(hn, hypergraph.partID(representative));
    }
  }

  void undoSparsification(PartitionedHyperGraph& hypergraph) {
    ASSERT(_is_init);
    hypergraph.doParallelForAllNodes(_task_group_id, [&](const HypernodeID hn) {
      const HypernodeID original_id = hypergraph.originalNodeID(hn);
      ASSERT(original_id < _mapping.size());
      const HypernodeID original_sparsified_id = _mapping[original_id];
      const HypernodeID sparsified_hn = _sparsified_partitioned_hg.globalNodeID(original_sparsified_id);
      ASSERT(_sparsified_partitioned_hg.nodeIsEnabled(sparsified_hn));
      hypergraph.setNodePart(hn, _sparsified_partitioned_hg.partID(sparsified_hn));
    });
    hypergraph.initializeNumCutHyperedges(_task_group_id);
  }

  // ####################### Other #######################

  void assignAllDegreeZeroHypernodesToSameCommunity(HyperGraph& hypergraph, ds::Clustering& clustering) {
    ASSERT(hypergraph.initialNumNodes() == clustering.size());
    PartitionID community_id = -1;
    for ( const HypernodeID& hn : hypergraph.nodes() ) {
      if ( hypergraph.nodeDegree(hn) == 0 ) {
        if ( community_id >= 0 ) {
          clustering[hypergraph.originalNodeID(hn)] = community_id;
        } else {
          community_id = clustering[hypergraph.originalNodeID(hn)];
        }
      }
    }
  }

 private:
  // ! Similiar to contractDegreeZeroHypernodes, but contraction is not applied
  // ! directly to hypergraph. Instead a mapping is computed that maps each vertex
  // ! of the original hypergraph to its supervertex and the weight of each
  // ! supervertex is aggregated in the hypernode weight vector.
  HypernodeID degreeZeroSparsification(const HyperGraph& hypergraph,
                                       parallel::scalable_vector<int>& vertices_to_numa_node,
                                       parallel::scalable_vector<HypernodeWeight>& hypernode_weight) {
    ASSERT(hypergraph.initialNumNodes() == vertices_to_numa_node.size());
    ASSERT(hypergraph.initialNumNodes() == hypernode_weight.size());
    ASSERT(hypergraph.initialNumNodes() == _mapping.size());

    HypernodeID num_hypernodes = 0;
    HypernodeID last_degree_zero_representative = kInvalidHypernode;
    for (const HypernodeID& hn : hypergraph.nodes()) {
      const HypernodeID original_id = hypergraph.originalNodeID(hn);
      if ( hypergraph.nodeDegree(hn) == 0 ) {
        bool is_removed = false;
        if ( last_degree_zero_representative != kInvalidHypernode ) {
          const HypernodeWeight weight = hypergraph.nodeWeight(hn);
          if ( hypernode_weight[last_degree_zero_representative] + weight <=
               _context.coarsening.max_allowed_node_weight ) {
            // Aggregate weight of each degree-zero vertex in
            // its represenative supervertex
            _mapping[original_id] = last_degree_zero_representative;
            hypernode_weight[last_degree_zero_representative] += weight;
            is_removed = true;
          }
        }

        if ( !is_removed ) {
          hypernode_weight[num_hypernodes] = hypergraph.nodeWeight(hn);
          vertices_to_numa_node[num_hypernodes] = common::get_numa_node_of_vertex(hn);
          last_degree_zero_representative = num_hypernodes;
          _mapping[original_id] = num_hypernodes++;
        }
      } else {
        hypernode_weight[num_hypernodes] = hypergraph.nodeWeight(hn);
        vertices_to_numa_node[num_hypernodes] = common::get_numa_node_of_vertex(hn);
        _mapping[original_id] = num_hypernodes++;
      }
    }
    hypernode_weight.resize(num_hypernodes);
    vertices_to_numa_node.resize(num_hypernodes);

    return num_hypernodes;
  }

  // ! Removes hyperedges where the weight of all pins is greater
  // ! than a certain threshold. The threshold is specified in
  // ! '_context.initial_partitioning.max_hyperedge_pin_weight'
  HyperedgeID heavyHyperedgeRemovalSparsification(const HyperGraph& hypergraph,
                                                  HyperedgeVector& edge_vector,
                                                  parallel::scalable_vector<HyperedgeWeight>& hyperedge_weight) {
    ASSERT(hypergraph.initialNumNodes() == _mapping.size());

    utils::Timer::instance().start_timer("compute_included_hes", "Calc Included HEs");
    parallel::scalable_vector<HyperedgeID> include_in_sparsified_hg(hypergraph.initialNumEdges(), 0UL);
    hypergraph.doParallelForAllEdges(_task_group_id, [&](const HyperedgeID he) {
      HypernodeWeight pin_weight = 0;
      for ( const HypernodeID& pin : hypergraph.pins(he) ) {
        pin_weight += hypergraph.nodeWeight(pin);
      }
      // Hyperedge will be include in sparsified hypergraph if its weight of
      // all pins is less than a predefined upper bound
      if ( pin_weight <= _context.initial_partitioning.sparsification.max_hyperedge_pin_weight ) {
        include_in_sparsified_hg[hypergraph.originalEdgeID(he)] = 1UL;
      }
    });
    utils::Timer::instance().stop_timer("compute_included_hes");

    utils::Timer::instance().start_timer("sparsified_prefix_sum", "Sparsified HE Prefix Sum");
    parallel::TBBPrefixSum<HyperedgeID> sparsified_hg_prefix_sum(include_in_sparsified_hg);
    tbb::parallel_scan(tbb::blocked_range<size_t>(
      0UL, hypergraph.initialNumEdges()), sparsified_hg_prefix_sum);
    utils::Timer::instance().stop_timer("sparsified_prefix_sum");

    utils::Timer::instance().start_timer("add_sparsified_hes", "Add Sparsified HEs");
    tbb::parallel_invoke([&] {
      edge_vector.resize(sparsified_hg_prefix_sum.total_sum());
    }, [&]{
      hyperedge_weight.resize(sparsified_hg_prefix_sum.total_sum());
    });
    hypergraph.doParallelForAllEdges(_task_group_id, [&](const HyperedgeID he) {
      const HyperedgeID original_id = hypergraph.originalEdgeID(he);
      const bool included_in_sparsified_hg = sparsified_hg_prefix_sum.value(original_id);
      if ( included_in_sparsified_hg ) {
        const HyperedgeID pos = sparsified_hg_prefix_sum[original_id];
        hyperedge_weight[pos] = hypergraph.edgeWeight(he);
        for ( const HypernodeID& pin : hypergraph.pins(he) ) {
          edge_vector[pos].push_back(_mapping[hypergraph.originalNodeID(pin)]);
        }
      }
    });
    utils::Timer::instance().stop_timer("add_sparsified_hes");

    return edge_vector.size();
  }

  // ! Removes vertices that have a degree greater than a certain threshold
  // ! from all hyperedges. High-degree vertices are removed from hyperedges
  // ! as long as their are more than two pins left.
  void highDegreeVertexSparsification(const HypernodeID num_hypernodes,
                                      HyperedgeVector& edge_vector,
                                      const parallel::scalable_vector<HyperedgeWeight>& hyperedge_weight) {
    ASSERT(edge_vector.size() == hyperedge_weight.size());
    const HyperedgeID num_hyperedges = edge_vector.size();
    const HyperedgeID high_degree_threshold = static_cast<double>(num_hyperedges) *
      _context.initial_partitioning.sparsification.high_degree_threshold_factor;

    utils::Timer::instance().start_timer("calculate_vertex_degrees", "Calc. Vertex Degrees");
    parallel::scalable_vector<parallel::IntegralAtomicWrapper<HyperedgeID>> vertex_degree(
      num_hypernodes, parallel::IntegralAtomicWrapper<HyperedgeID>(0));
    tbb::parallel_for(0UL, num_hyperedges, [&](const size_t i) {
      for ( const HypernodeID& pin : edge_vector[i] ) {
        ++vertex_degree[pin];
      }
    });
    utils::Timer::instance().stop_timer("calculate_vertex_degrees");

    utils::Timer::instance().start_timer("remove_high_degree_vertices_from_hes", "Remove High-Degree HNs from HEs");
    tbb::parallel_for(0UL, num_hyperedges, [&](const size_t i) {
      if ( hyperedge_weight[i] <= _context.initial_partitioning.sparsification.high_degree_hyperedge_weight_threshold ) {
        parallel::scalable_vector<HypernodeID>& hyperedge = edge_vector[i];
        std::sort(hyperedge.begin(), hyperedge.end(),
          [&](const HypernodeID lhs, const HypernodeID rhs) {
            return vertex_degree[lhs] < vertex_degree[rhs];
          });
        while ( hyperedge.size() > 2 ) {
          const HypernodeID pin = hyperedge.back();
          if ( vertex_degree[pin] > high_degree_threshold ) {
            hyperedge.pop_back();
          } else {
            break;
          }
        }
      }
    });
    utils::Timer::instance().stop_timer("remove_high_degree_vertices_from_hes");
  }

  const Context& _context;
  const TaskGroupID _task_group_id;

  bool _is_init;
  HyperGraph _sparsified_hg;
  PartitionedHyperGraph _sparsified_partitioned_hg;
  parallel::scalable_vector<HyperedgeID> _removed_hes;
  parallel::scalable_vector<HypernodeID> _removed_hns;
  parallel::scalable_vector<HypernodeID> _mapping;
};
}  // namespace mt_kahypar

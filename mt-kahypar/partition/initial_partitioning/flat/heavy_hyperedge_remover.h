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

#include "tbb/parallel_invoke.h"
#include "tbb/parallel_scan.h"
#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"

namespace mt_kahypar {
template<typename TypeTraits>
class HeavyHyperedgeRemoverT {
  using HyperGraph = typename TypeTraits::HyperGraph;
  using PartitionedHyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
  using HyperGraphFactory = typename TypeTraits::HyperGraphFactory;
  using TBB = typename TypeTraits::TBB;
  using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;

 public:
  HeavyHyperedgeRemoverT(const Context& context,
                         const TaskGroupID task_group_id) :
    _is_init(false),
    _sparsified_hg(),
    _sparsified_partitioned_hg(),
    _context(context),
    _task_group_id(task_group_id) { }

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
    HyperedgeVector edge_vector;
    parallel::scalable_vector<int> vertices_to_numa_node;
    parallel::scalable_vector<HyperedgeWeight> hyperedge_weight;
    parallel::scalable_vector<HypernodeWeight> hypernode_weight;

    tbb::parallel_invoke([&] {
      utils::Timer::instance().start_timer("setup_hyperedges", "Setup Hyperedges", true);
      utils::Timer::instance().start_timer("compute_included_hes", "Calc Included HEs", true);
      parallel::scalable_vector<HyperedgeID> include_in_sparsified_hg(hypergraph.initialNumEdges(), 0UL);
      hypergraph.doParallelForAllEdges(_task_group_id, [&](const HyperedgeID he) {
        HypernodeWeight pin_weight = 0;
        for ( const HypernodeID& pin : hypergraph.pins(he) ) {
          pin_weight += hypergraph.nodeWeight(pin);
        }
        // Hyperedge will be include in sparsified hypergraph if its weight of
        // all pins is less than a predefined upper bound
        if ( pin_weight <= _context.initial_partitioning.max_hyperedge_pin_weight ) {
          include_in_sparsified_hg[hypergraph.originalEdgeID(he)] = 1UL;
        }
      });
      utils::Timer::instance().stop_timer("compute_included_hes");

      utils::Timer::instance().start_timer("sparsified_prefix_sum", "Sparsified HE Prefix Sum", true);
      parallel::TBBPrefixSum<HyperedgeID> sparsified_hg_prefix_sum(include_in_sparsified_hg);
      tbb::parallel_scan(tbb::blocked_range<size_t>(
        0UL, hypergraph.initialNumEdges()), sparsified_hg_prefix_sum);
      utils::Timer::instance().stop_timer("sparsified_prefix_sum");

      utils::Timer::instance().start_timer("add_sparsified_hes", "Add Sparsified HEs", true);
      edge_vector.resize(sparsified_hg_prefix_sum.total_sum());
      hyperedge_weight.resize(sparsified_hg_prefix_sum.total_sum());
      hypergraph.doParallelForAllEdges(_task_group_id, [&](const HyperedgeID he) {
        const HyperedgeID original_id = hypergraph.originalEdgeID(he);
        const bool included_in_sparsified_hg = sparsified_hg_prefix_sum.value(original_id);
        if ( included_in_sparsified_hg ) {
          const HyperedgeID pos = sparsified_hg_prefix_sum[original_id];
          hyperedge_weight[pos] = hypergraph.edgeWeight(he);
          for ( const HypernodeID& pin : hypergraph.pins(he) ) {
            edge_vector[pos].push_back(hypergraph.originalNodeID(pin));
          }
        }
      });
      utils::Timer::instance().stop_timer("add_sparsified_hes");
      utils::Timer::instance().stop_timer("setup_hyperedges");
    }, [&] {
      utils::Timer::instance().start_timer("setup_hypernodes", "Setup Hypernodes", true);
      hypernode_weight.resize(hypergraph.initialNumNodes());
      vertices_to_numa_node.resize(hypergraph.initialNumNodes());
      hypergraph.doParallelForAllNodes(_task_group_id, [&](const HypernodeID hn) {
        const HypernodeID original_id = hypergraph.originalNodeID(hn);
        hypernode_weight[original_id] = hypergraph.nodeWeight(hn);
        vertices_to_numa_node[original_id] = common::get_numa_node_of_vertex(hn);
      });
      utils::Timer::instance().stop_timer("setup_hypernodes");
    });

    utils::Timer::instance().start_timer("construct_sparsified_hypergraph", "Construct Sparsified HG");
    const HypernodeID num_hypernodes = hypergraph.initialNumNodes();
    const HyperedgeID num_hyperedges = edge_vector.size();
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

  void applyPartition(PartitionedHyperGraph& partitioned_hypergraph) {
    ASSERT(_is_init);
    ASSERT(partitioned_hypergraph.initNumNodes() == _sparsified_partitioned_hg.initNumNodes());
    partitioned_hypergraph.doParallelForAllNodes(_task_group_id, [&](const HypernodeID hn) {
      ASSERT(_sparsified_partitioned_hg.nodeIsEnabled(hn));
      partitioned_hypergraph.setNodePart(hn, _sparsified_partitioned_hg.partID(hn));
    });
    partitioned_hypergraph.initializeNumCutHyperedges(_task_group_id);
  }

 private:
  bool _is_init;
  HyperGraph _sparsified_hg;
  PartitionedHyperGraph _sparsified_partitioned_hg;
  const Context& _context;
  const TaskGroupID _task_group_id;

};
}  // namespace mt_kahypar

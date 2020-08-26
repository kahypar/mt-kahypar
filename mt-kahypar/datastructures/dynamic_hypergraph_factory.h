/*******************************************************************************
 * This file is part of KaHyPar.
 *
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

#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_invoke.h"
#include "tbb/parallel_scan.h"

#include "kahypar/utils/math.h"

#include "mt-kahypar/datastructures/dynamic_hypergraph.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
namespace ds {

class DynamicHypergraphFactory {

  using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
  using Counter = parallel::scalable_vector<size_t>;
  using AtomicCounter = parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>>;
  using ThreadLocalCounter = tbb::enumerable_thread_specific<Counter>;
  using ThreadLocalBitset = tbb::enumerable_thread_specific<parallel::scalable_vector<bool>>;

 public:
  static DynamicHypergraph construct(const TaskGroupID task_group_id,
                                    const HypernodeID num_hypernodes,
                                    const HyperedgeID num_hyperedges,
                                    const HyperedgeVector& edge_vector,
                                    const HyperedgeWeight* hyperedge_weight = nullptr,
                                    const HypernodeWeight* hypernode_weight = nullptr,
                                    const bool stable_construction_of_incident_edges = false) {
    DynamicHypergraph hypergraph;
    hypergraph._num_hypernodes = num_hypernodes;
    hypergraph._num_hyperedges = num_hyperedges;
    tbb::parallel_invoke([&] {
      hypergraph._hypernodes.resize(num_hypernodes);
    }, [&] {
      hypergraph._incident_nets.resize(num_hypernodes);
    }, [&] {
      hypergraph._hyperedges.resize(num_hyperedges + 1);
    }, [&] {
      hypergraph._removable_single_pin_and_parallel_nets =
        kahypar::ds::FastResetFlagArray<>(num_hyperedges);
    });
    hypergraph._removable_incident_nets = ThreadLocalBitset(num_hyperedges, false);

    ASSERT(edge_vector.size() == num_hyperedges);

    // Compute number of pins per hyperedge and number
    // of incident nets per vertex
    utils::Timer::instance().start_timer("compute_ds_sizes", "Precompute DS Size", true);
    Counter num_pins_per_hyperedge(num_hyperedges, 0);
    ThreadLocalCounter local_incident_nets_per_vertex(num_hypernodes, 0);
    tbb::enumerable_thread_specific<size_t> local_max_edge_size(0UL);
    tbb::parallel_for(ID(0), num_hyperedges, [&](const size_t pos) {
      Counter& num_incident_nets_per_vertex = local_incident_nets_per_vertex.local();
      num_pins_per_hyperedge[pos] = edge_vector[pos].size();
      local_max_edge_size.local() = std::max(
        local_max_edge_size.local(), edge_vector[pos].size());
      for ( const HypernodeID& pin : edge_vector[pos] ) {
        ++num_incident_nets_per_vertex[pin];
      }
    });
    hypergraph._max_edge_size = local_max_edge_size.combine(
      [&](const size_t lhs, const size_t rhs) {
        return std::max(lhs, rhs);
      });

    // We sum up the number of incident nets per vertex only thread local.
    // To obtain the global number of incident nets per vertex, we iterate
    // over each thread local counter and sum it up.
    size_t cnt = 0;
    AtomicCounter num_incident_nets_per_vertex(
      num_hypernodes, parallel::IntegralAtomicWrapper<size_t>(0));
    for ( Counter& c : local_incident_nets_per_vertex ) {
      ++cnt;
      const bool allocate_incident_nets = cnt == local_incident_nets_per_vertex.size();
      tbb::parallel_for(ID(0), num_hypernodes, [&](const size_t pos) {
        num_incident_nets_per_vertex[pos] += c[pos];
        if ( allocate_incident_nets ) {
          hypergraph._incident_nets[pos].resize(num_incident_nets_per_vertex[pos]);
        }
      });
    }
    utils::Timer::instance().stop_timer("compute_ds_sizes");

    // Compute prefix sum over the number of pins per hyperedge and the.
    // The prefix sum is used than as
    // start position for each hyperedge in the incidence array.
    utils::Timer::instance().start_timer("compute_incidence_array_prefix_sum", "Compute Incidence Array PS", true);
    parallel::TBBPrefixSum<size_t> pin_prefix_sum(num_pins_per_hyperedge);
    tbb::parallel_scan(tbb::blocked_range<size_t>(
      0UL, UI64(num_hyperedges)), pin_prefix_sum);
    utils::Timer::instance().stop_timer("compute_incidence_array_prefix_sum");

    utils::Timer::instance().start_timer("setup_hypergraph", "Setup hypergraph", true);
    hypergraph._num_pins = pin_prefix_sum.total_sum();
    hypergraph._total_degree = pin_prefix_sum.total_sum();
    hypergraph._incidence_array.resize(hypergraph._num_pins);

    tbb::parallel_invoke([&] {
      hypergraph._acquired_hes.assign(
        num_hyperedges, parallel::IntegralAtomicWrapper<bool>(false));
      tbb::parallel_for(ID(0), num_hyperedges, [&](const size_t pos) {
        // Setup hyperedges
        DynamicHypergraph::Hyperedge& hyperedge = hypergraph._hyperedges[pos];
        hyperedge.enable();
        hyperedge.setFirstEntry(pin_prefix_sum[pos]);
        hyperedge.setSize(pin_prefix_sum.value(pos));
        if ( hyperedge_weight ) {
          hyperedge.setWeight(hyperedge_weight[pos]);
        }

        const HyperedgeID he = pos;
        size_t incidence_array_pos = hyperedge.firstEntry();
        size_t hash = kEdgeHashSeed;
        for ( const HypernodeID& pin : edge_vector[pos] ) {
          ASSERT(incidence_array_pos < hyperedge.firstInvalidEntry());
          ASSERT(pin < num_hypernodes);
          // Compute hash of hyperedge
          hash += kahypar::math::hash(pin);
          // Add pin to incidence array
          hypergraph._incidence_array[incidence_array_pos++] = pin;
          // Add hyperedge he as a incident net to pin
          const size_t incident_nets_pos = --num_incident_nets_per_vertex[pin];
          ASSERT(incident_nets_pos < hypergraph._incident_nets[pin].size());
          hypergraph._incident_nets[pin][incident_nets_pos] = he;
        }
        hyperedge.hash() = hash;
      });
      // Sentinel
      hypergraph._hyperedges[num_hyperedges].enable();
      hypergraph._hyperedges[num_hyperedges].setFirstEntry(hypergraph._num_pins);
    }, [&] {
      tbb::parallel_invoke([&] {
        hypergraph._acquired_hns.assign(
          num_hypernodes, parallel::IntegralAtomicWrapper<bool>(false));
      }, [&] {
        hypergraph._contraction_tree.initialize(num_hypernodes);
      });
      tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID hn) {
        // Setup hypernodes
        DynamicHypergraph::Hypernode& hypernode = hypergraph._hypernodes[hn];
        hypernode.enable();
        if ( hypernode_weight ) {
          hypernode.setWeight(hypernode_weight[hn]);
        }
      });
    }, [&] {
      // graph edge ID mapping
      hypergraph._num_graph_edges_up_to.resize(num_hyperedges + 1);
      tbb::parallel_for(0U, num_hyperedges, [&](const HyperedgeID e) {
        const size_t edge_size = edge_vector[e].size();   // hypergraph.edgeSize(e) is not yet constructed
        hypergraph._num_graph_edges_up_to[e+1] = static_cast<HyperedgeID>(edge_size == 2);
      }, tbb::static_partitioner());
      hypergraph._num_graph_edges_up_to[0] = 0;

      parallel::TBBPrefixSum<HyperedgeID, Array> scan_graph_edges(hypergraph._num_graph_edges_up_to);
      tbb::parallel_scan(tbb::blocked_range<size_t>(0, num_hyperedges + 1), scan_graph_edges);
      hypergraph._num_graph_edges = scan_graph_edges.total_sum();
    });

    if (stable_construction_of_incident_edges) {
      // sort incident hyperedges of each node, so their ordering is independent of scheduling (and the same as a typical sequential implementation)
      tbb::parallel_for(ID(0), num_hypernodes, [&](HypernodeID u) {
        std::sort(hypergraph._incident_nets[u].begin(), hypergraph._incident_nets[u].end());
      });
    }

    // Compute total weight of hypergraph
    hypergraph.updateTotalWeight(task_group_id);
    utils::Timer::instance().stop_timer("setup_hypergraph");
    return hypergraph;
  }

  /**
   * Compactifies a given hypergraph such that it only contains enabled vertices and hyperedges within
   * a consecutive range of IDs.
   */
  static std::pair<DynamicHypergraph,
                   parallel::scalable_vector<HypernodeID> > compactify(const TaskGroupID task_group_id,
                                                                       const DynamicHypergraph& hypergraph) {
    HypernodeID num_hypernodes = 0;
    HyperedgeID num_hyperedges = 0;
    parallel::scalable_vector<HypernodeID> hn_mapping;
    parallel::scalable_vector<HyperedgeID> he_mapping;
    // Computes a mapping for vertices and hyperedges to a consecutive range of IDs
    // in the compactified hypergraph via a parallel prefix sum
    utils::Timer::instance().start_timer("compactify_hn_and_he_ids", "Compactify HN and HE IDs");
    tbb::parallel_invoke([&] {
      hn_mapping.assign(hypergraph._num_hypernodes + 1, 0);
      hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
        hn_mapping[hn + 1] = ID(1);
      });

      parallel::TBBPrefixSum<HypernodeID, parallel::scalable_vector> hn_mapping_prefix_sum(hn_mapping);
      tbb::parallel_scan(tbb::blocked_range<size_t>(
        0UL, hypergraph._num_hypernodes + 1), hn_mapping_prefix_sum);
      num_hypernodes = hn_mapping_prefix_sum.total_sum();
      hn_mapping.resize(hypergraph._num_hypernodes);
    }, [&] {
      he_mapping.assign(hypergraph._num_hyperedges + 1, 0);
      hypergraph.doParallelForAllEdges([&](const HyperedgeID& he) {
        he_mapping[he + 1] = ID(1);
      });

      parallel::TBBPrefixSum<HyperedgeID, parallel::scalable_vector> he_mapping_prefix_sum(he_mapping);
      tbb::parallel_scan(tbb::blocked_range<size_t>(
        0UL, hypergraph._num_hyperedges + 1), he_mapping_prefix_sum);
      num_hyperedges = he_mapping_prefix_sum.total_sum();
      he_mapping.resize(hypergraph._num_hyperedges);
    });
    utils::Timer::instance().stop_timer("compactify_hn_and_he_ids");

    // Remap pins of each hyperedge
    utils::Timer::instance().start_timer("remap_pins", "Remap Pins");
    using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
    HyperedgeVector edge_vector;
    parallel::scalable_vector<HyperedgeWeight> hyperedge_weights;
    parallel::scalable_vector<HypernodeWeight> hypernode_weights;
    tbb::parallel_invoke([&] {
      hypernode_weights.resize(num_hypernodes);
      hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
        const HypernodeID mapped_hn = hn_mapping[hn];
        ASSERT(mapped_hn < num_hypernodes);
        hypernode_weights[mapped_hn] = hypergraph.nodeWeight(hn);
      });
    }, [&] {
      edge_vector.resize(num_hyperedges);
      hyperedge_weights.resize(num_hyperedges);
      hypergraph.doParallelForAllEdges([&](const HyperedgeID he) {
        const HyperedgeID mapped_he = he_mapping[he];
        ASSERT(mapped_he < num_hyperedges);
        hyperedge_weights[mapped_he] = hypergraph.edgeWeight(he);
        for ( const HypernodeID pin : hypergraph.pins(he) ) {
          edge_vector[mapped_he].push_back(hn_mapping[pin]);
        }
      });
    });
    utils::Timer::instance().stop_timer("remap_pins");

    // Construct compactified hypergraph
    utils::Timer::instance().start_timer("construct_compactified_hypergraph", "Construct Compactified Hypergraph");
    DynamicHypergraph compactified_hypergraph = DynamicHypergraphFactory::construct(
      task_group_id, num_hypernodes, num_hyperedges,
      edge_vector, hyperedge_weights.data(), hypernode_weights.data());
    compactified_hypergraph._removed_degree_zero_hn_weight = hypergraph._removed_degree_zero_hn_weight;
    compactified_hypergraph._total_weight += hypergraph._removed_degree_zero_hn_weight;
    utils::Timer::instance().stop_timer("construct_compactified_hypergraph");

    // Set community ids
    utils::Timer::instance().start_timer("initialize_communities", "Initialized Communities");
    hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
      const HypernodeID mapped_hn = hn_mapping[hn];
      compactified_hypergraph.setCommunityID(mapped_hn, hypergraph.communityID(hn));
    });
    compactified_hypergraph.initializeCommunities();
    utils::Timer::instance().stop_timer("initialize_communities");

    tbb::parallel_invoke([&] {
      parallel::parallel_free(he_mapping,
        hyperedge_weights, hypernode_weights);
    }, [&] {
      parallel::parallel_free(edge_vector);
    });

    return std::make_pair(std::move(compactified_hypergraph), std::move(hn_mapping));
  }

 private:
  DynamicHypergraphFactory() { }
};

} // namespace ds
} // namespace mt_kahypar
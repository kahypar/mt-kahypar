/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "static_hypergraph_factory.h"

#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>

#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar::ds {

  StaticHypergraph StaticHypergraphFactory::construct(
          const TaskGroupID task_group_id,
          const HypernodeID num_hypernodes,
          const HyperedgeID num_hyperedges,
          const HyperedgeVector& edge_vector,
          const HyperedgeWeight* hyperedge_weight,
          const HypernodeWeight* hypernode_weight,
          const bool stable_construction_of_incident_edges) {
    StaticHypergraph hypergraph;
    hypergraph._num_hypernodes = num_hypernodes;
    hypergraph._num_hyperedges = num_hyperedges;
    hypergraph._hypernodes.resize(num_hypernodes + 1);
    hypergraph._hyperedges.resize(num_hyperedges + 1);

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
        ASSERT(pin < num_hypernodes, V(pin) << V(num_hypernodes));
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
    Counter num_incident_nets_per_vertex(num_hypernodes, 0);
    for ( Counter& c : local_incident_nets_per_vertex ) {
      tbb::parallel_for(ID(0), num_hypernodes, [&](const size_t pos) {
        num_incident_nets_per_vertex[pos] += c[pos];
      });
    }
    utils::Timer::instance().stop_timer("compute_ds_sizes");

    // Compute prefix sum over the number of pins per hyperedge and the
    // number of incident nets per vertex. The prefix sum is used than as
    // start position for each hyperedge resp. hypernode in the incidence
    // resp. incident nets array.
    utils::Timer::instance().start_timer("compute_prefix_sums", "Compute Prefix Sums", true);
    parallel::TBBPrefixSum<size_t> pin_prefix_sum(num_pins_per_hyperedge);
    parallel::TBBPrefixSum<size_t> incident_net_prefix_sum(num_incident_nets_per_vertex);
    tbb::parallel_invoke([&] {
      tbb::parallel_scan(tbb::blocked_range<size_t>(
              0UL, UI64(num_hyperedges)), pin_prefix_sum);
    }, [&] {
      tbb::parallel_scan(tbb::blocked_range<size_t>(
              0UL, UI64(num_hypernodes)), incident_net_prefix_sum);
    });
    utils::Timer::instance().stop_timer("compute_prefix_sums");

    utils::Timer::instance().start_timer("setup_hypergraph", "Setup hypergraph", true);
    ASSERT(pin_prefix_sum.total_sum() == incident_net_prefix_sum.total_sum());
    hypergraph._num_pins = pin_prefix_sum.total_sum();
    hypergraph._total_degree = incident_net_prefix_sum.total_sum();
    hypergraph._incident_nets.resize(hypergraph._num_pins);
    hypergraph._incidence_array.resize(hypergraph._num_pins);

    AtomicCounter incident_nets_position(num_hypernodes,
                                         parallel::IntegralAtomicWrapper<size_t>(0));

    auto setup_hyperedges = [&] {
      tbb::parallel_for(ID(0), num_hyperedges, [&](const size_t pos) {
        StaticHypergraph::Hyperedge& hyperedge = hypergraph._hyperedges[pos];
        hyperedge.enable();
        hyperedge.setFirstEntry(pin_prefix_sum[pos]);
        hyperedge.setSize(pin_prefix_sum.value(pos));
        if ( hyperedge_weight ) {
          hyperedge.setWeight(hyperedge_weight[pos]);
        }

        const HyperedgeID he = pos;
        size_t incidence_array_pos = hyperedge.firstEntry();
        for ( const HypernodeID& pin : edge_vector[pos] ) {
          ASSERT(incidence_array_pos < hyperedge.firstInvalidEntry());
          ASSERT(pin < num_hypernodes);
          // Add pin to incidence array
          hypergraph._incidence_array[incidence_array_pos++] = pin;
          // Add hyperedge he as a incident net to pin
          const size_t incident_nets_pos = incident_net_prefix_sum[pin] + incident_nets_position[pin]++;
          ASSERT(incident_nets_pos < incident_net_prefix_sum[pin + 1]);
          hypergraph._incident_nets[incident_nets_pos] = he;
        }
      });
    };

    auto setup_hypernodes = [&] {
      tbb::parallel_for(ID(0), num_hypernodes, [&](const size_t pos) {
        StaticHypergraph::Hypernode& hypernode = hypergraph._hypernodes[pos];
        hypernode.enable();
        hypernode.setFirstEntry(incident_net_prefix_sum[pos]);
        hypernode.setSize(incident_net_prefix_sum.value(pos));
        if ( hypernode_weight ) {
          hypernode.setWeight(hypernode_weight[pos]);
        }
      });
    };

    auto init_communities = [&] {
      hypergraph._community_ids.resize(num_hypernodes, 0);
    };

    tbb::parallel_invoke(setup_hyperedges, setup_hypernodes, init_communities);

    if (stable_construction_of_incident_edges) {
      // sort incident hyperedges of each node, so their ordering is independent of scheduling (and the same as a typical sequential implementation)
      tbb::parallel_for(ID(0), num_hypernodes, [&](HypernodeID u) {
        auto b = hypergraph._incident_nets.begin() + hypergraph.hypernode(u).firstEntry();
        auto e = hypergraph._incident_nets.begin() + hypergraph.hypernode(u).firstInvalidEntry();
        std::sort(b, e);
      });
    }

    // Add Sentinels
    hypergraph._hypernodes.back() = StaticHypergraph::Hypernode(hypergraph._incident_nets.size());
    hypergraph._hyperedges.back() = StaticHypergraph::Hyperedge(hypergraph._incidence_array.size());

    hypergraph.computeAndSetTotalNodeWeight(task_group_id);

    utils::Timer::instance().stop_timer("setup_hypergraph");
    return hypergraph;
  }

}
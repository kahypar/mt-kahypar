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

#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

namespace mt_kahypar {
namespace ds {

class StaticHypergraphFactory {

  using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
  using Counter = parallel::scalable_vector<size_t>;
  using AtomicCounter = parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>>;
  using ThreadLocalCounter = tbb::enumerable_thread_specific<Counter>;

 public:
  static StaticHypergraph construct(const HypernodeID num_hypernodes,
                                    const HyperedgeID num_hyperedges,
                                    const HyperedgeVector& edge_vector,
                                    const HyperedgeWeight* hyperedge_weight = nullptr,
                                    const HypernodeWeight* hypernode_weight = nullptr) {
    StaticHypergraph hypergraph;
    hypergraph._num_hypernodes = num_hypernodes;
    hypergraph._num_hyperedges = num_hyperedges;
    hypergraph._hypernodes.resize(num_hypernodes);
    hypergraph._hyperedges.resize(num_hyperedges);

    ASSERT(edge_vector.size() == num_hyperedges);

    // Compute number of pins per hyperedge and number
    // of incident nets per vertex
    Counter num_pins_per_hyperedge(num_hyperedges, 0);
    ThreadLocalCounter local_incident_nets_per_vertex(num_hypernodes, 0);
    tbb::parallel_for(0UL, num_hyperedges, [&](const size_t pos) {
      Counter& num_incident_nets_per_vertex = local_incident_nets_per_vertex.local();
      num_pins_per_hyperedge[pos] = edge_vector[pos].size();
      for ( const HypernodeID& pin : edge_vector[pos] ) {
        ASSERT(pin < num_hypernodes);
        ++num_incident_nets_per_vertex[pin];
      }
    });

    // We sum up the number of incident nets per vertex only thread local.
    // To obtain the global number of incident nets per vertex, we iterate
    // over each thread local counter and sum it up.
    Counter num_incident_nets_per_vertex(num_hypernodes, 0);
    for ( Counter& c : local_incident_nets_per_vertex ) {
      tbb::parallel_for(0UL, num_hypernodes, [&](const size_t pos) {
        num_incident_nets_per_vertex[pos] += c[pos];
      });
    }

    // Compute prefix sum over the number of pins per hyperedge and the
    // number of incident nets per vertex. The prefix sum is used than as
    // start position for each hyperedge resp. hypernode in the incidence
    // resp. incident nets array.
    parallel::TBBPrefixSum<size_t> pin_prefix_sum(num_pins_per_hyperedge);
    parallel::TBBPrefixSum<size_t> incident_net_prefix_sum(num_incident_nets_per_vertex);
    tbb::parallel_invoke([&] {
      tbb::parallel_scan(tbb::blocked_range<size_t>(
        0UL, num_hyperedges), pin_prefix_sum);
    }, [&] {
      tbb::parallel_scan(tbb::blocked_range<size_t>(
        0UL, num_hypernodes), incident_net_prefix_sum);
    });

    ASSERT(pin_prefix_sum.total_sum() == incident_net_prefix_sum.total_sum());
    hypergraph._num_pins = pin_prefix_sum.total_sum();
    hypergraph._incident_nets.resize(hypergraph._num_pins);
    hypergraph._incidence_array.resize(hypergraph._num_pins);

    AtomicCounter incident_nets_position(num_hypernodes,
      parallel::IntegralAtomicWrapper<size_t>(0));
    tbb::parallel_invoke([&] {
      tbb::parallel_for(0UL, num_hyperedges, [&](const size_t pos) {
        // Setup hyperedges
        StaticHypergraph::Hyperedge& hyperedge = hypergraph._hyperedges[pos];
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
          const size_t incident_nets_pos = incident_net_prefix_sum[pin] +
            incident_nets_position[pin]++;
          ASSERT(incident_nets_pos < incident_net_prefix_sum[pin + 1]);
          hypergraph._incident_nets[incident_nets_pos] = he;
        }
        hyperedge.hash() = hash;
      });
    }, [&] {
      tbb::parallel_for(0UL, num_hypernodes, [&](const size_t pos) {
        // Setup hypernodes
        StaticHypergraph::Hypernode& hypernode = hypergraph._hypernodes[pos];
        hypernode.enable();
        hypernode.setFirstEntry(incident_net_prefix_sum[pos]);
        hypernode.setSize(incident_net_prefix_sum.value(pos));
        if ( hypernode_weight ) {
          hypernode.setWeight(hypernode_weight[pos]);
        }
      });
    });

    // Compute total weight of hypergraph
    hypergraph._total_weight = tbb::parallel_reduce(
      tbb::blocked_range<HypernodeID>(0UL, num_hypernodes), 0,
      [&](const tbb::blocked_range<HypernodeID>& range, HypernodeWeight init) {
        HypernodeWeight weight = init;
        for (HypernodeID hn = range.begin(); hn < range.end(); ++hn) {
          weight += hypergraph._hypernodes[hn].weight();
        }
        return weight;
      }, std::plus<HypernodeWeight>());

    return hypergraph;
  }


 private:
  StaticHypergraphFactory() { }

};

} // namespace ds
} // namespace mt_kahypar
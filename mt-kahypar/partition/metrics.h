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

#include <cmath>

#include <algorithm>
#include <vector>

#include "mt-kahypar/definitions.h"

namespace mt_kahypar {
namespace metrics {

/**
 * Counts the number of pins which refers to an other numa node than the numa
 * node which its corresponding hyperedge belongs to
 */
static inline HyperedgeWeight remotePinCount(const Hypergraph& hypergraph) {
  int used_numa_nodes = TBBNumaArena::instance().num_used_numa_nodes();
  HyperedgeWeight remote_pin_count = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    int he_node = StreamingHypergraph::get_numa_node_of_hyperedge(he);
    std::vector<size_t> pin_count_on_node(used_numa_nodes, 0);
    for ( const HypernodeID& pin : hypergraph.pins(he) ) {
      int hn_node = StreamingHypergraph::get_numa_node_of_vertex(pin);
      ASSERT(hn_node < used_numa_nodes);
      ++pin_count_on_node[hn_node];
    }

    HyperedgeWeight he_remote_pin_count = 0;
    for ( int node = 0; node < used_numa_nodes; ++node ) {
      if ( he_node != node ) {
        he_remote_pin_count += pin_count_on_node[node];
      }
    }

    remote_pin_count += he_remote_pin_count;
  }
  return remote_pin_count;
}

static inline double avgHyperedgeDegree(const Hypergraph& hypergraph) {
  return static_cast<double>(hypergraph.currentNumPins()) / hypergraph.currentNumEdges();
}

static inline double avgHypernodeDegree(const Hypergraph& hypergraph) {
  return static_cast<double>(hypergraph.currentNumPins()) / hypergraph.currentNumNodes();
}

static inline HyperedgeID hypernodeDegreeRank(const Hypergraph& hypergraph, const size_t rank) {
  std::vector<HyperedgeID> hn_degrees;
  hn_degrees.reserve(hypergraph.currentNumNodes());
  for (const auto& hn : hypergraph.nodes()) {
    hn_degrees.push_back(hypergraph.nodeDegree(hn));
  }
  ASSERT(!hn_degrees.empty(), "Hypergraph does not contain any hypernodes");
  std::sort(hn_degrees.begin(), hn_degrees.end());

  ASSERT(rank < hn_degrees.size());
  return hn_degrees[rank];
}

} // namespace metrics
} // namespace mt_kahypar
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

static inline HyperedgeWeight communicationVolume(const Hypergraph& hypergraph) {
  int used_numa_nodes = TBBNumaArena::instance().num_used_numa_nodes();
  HyperedgeWeight communication_volume = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    std::vector<size_t> pin_count_on_node(used_numa_nodes, 0);
    for ( const HypernodeID& pin : hypergraph.pins(he) ) {
      int node = StreamingHypergraph::get_numa_node_of_vertex(pin);
      ASSERT(node < used_numa_nodes);
      ++pin_count_on_node[node];
    }

    HyperedgeWeight connectivity = 0;
    for ( int i = 0; i < used_numa_nodes; ++i ) {
      if ( pin_count_on_node[i] > 0 ) {
        ++connectivity;
      }
    }

    ASSERT(connectivity > 0);
    communication_volume += (connectivity - 1) * hypergraph.edgeWeight(he);
  }
  return communication_volume;
}

static inline double avgHyperedgeDegree(const Hypergraph& hypergraph) {
  return static_cast<double>(hypergraph.currentNumPins()) / hypergraph.currentNumEdges();
}

static inline double avgHypernodeDegree(const Hypergraph& hypergraph) {
  return static_cast<double>(hypergraph.currentNumPins()) / hypergraph.currentNumNodes();
}

} // namespace metrics
} // namespace mt_kahypar
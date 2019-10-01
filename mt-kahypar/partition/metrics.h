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

// TODO(heuer): Replace this once connectivity is available in hypergraph

static inline PartitionID connectivity(const Hypergraph& hypergraph, const HyperedgeID he) {
  std::vector<bool> parts(hypergraph.k(), false);
  PartitionID connectivity = 0;
  for ( const HypernodeID& pin : hypergraph.pins(he) ) {
    PartitionID part_id = hypergraph.partID(pin);
    if ( !parts[part_id] ) {
      ++connectivity;
      parts[part_id] = true;
    }
  }
  return connectivity;
}

static inline HyperedgeWeight hyperedgeCut(const Hypergraph& hypergraph) {
  HyperedgeWeight cut = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    if ( connectivity(hypergraph, he) > 1 ) {
      cut += hypergraph.edgeWeight(he);
    }
  }
  return cut;
}

static inline HyperedgeWeight km1(const Hypergraph& hypergraph) {
  HyperedgeWeight km1 = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    km1 += std::max(connectivity(hypergraph, he) - 1, 0) * hypergraph.edgeWeight(he);
  }
  return km1;
}

static inline HyperedgeWeight soed(const Hypergraph& hypergraph) {
  HyperedgeWeight soed = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    PartitionID conn = connectivity(hypergraph, he);
    if ( conn > 1 ) {
      soed += conn * hypergraph.edgeWeight(he);
    }
  }
  return soed;
}

static inline double absorption(const Hypergraph& hypergraph) {
  double absorption_val = 0.0;
  for (PartitionID part = 0; part < hypergraph.k(); ++part) {
    for (const HyperedgeID& he : hypergraph.edges()) {
      HypernodeID pin_count_in_part = 0;
      for ( const HypernodeID& pin : hypergraph.pins(he) ) {
        if ( hypergraph.partID(pin) == part ) {
          ++pin_count_in_part;
        }
      }
      if (pin_count_in_part > 0 && hypergraph.edgeSize(he) > 1) {
        absorption_val += static_cast<double>((pin_count_in_part - 1)) / (hypergraph.edgeSize(he) - 1)
                          * hypergraph.edgeWeight(he);
      }
    }
  }
  return absorption_val;
}

static inline HyperedgeWeight objective(const Hypergraph& hg, const Objective& objective) {
  switch (objective) {
    case Objective::cut: return hyperedgeCut(hg);
    case Objective::km1: return km1(hg);
    default:
      LOG << "Unknown Objective";
      exit(-1);
  }
}


static inline double imbalance(const Hypergraph& hypergraph, const Context& context) {
  ASSERT(context.partition.perfect_balance_part_weights.size() == (size_t) context.partition.k);

  double max_balance = (hypergraph.partWeight(0) /
                        static_cast<double>(context.partition.perfect_balance_part_weights[0]));

  for (PartitionID i = 1; i < context.partition.k; ++i) {
    const double balance_i =
      (hypergraph.partWeight(i) /
       static_cast<double>(context.partition.perfect_balance_part_weights[i]));
    max_balance = std::max(max_balance, balance_i);
  }

  return max_balance - 1.0;
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
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

#include "mt-kahypar/partition/metrics.h"


#include <cmath>
#include <algorithm>

namespace mt_kahypar::metrics {
  HyperedgeWeight hyperedgeCut(const PartitionedHypergraph& hypergraph, const bool parallel) {
    if ( parallel ) {
      tbb::enumerable_thread_specific<HyperedgeWeight> cut(0);
      tbb::parallel_for(ID(0), hypergraph.initialNumEdges(), [&](const HyperedgeID he) {
        if (hypergraph.edgeIsEnabled(he) && hypergraph.connectivity(he) > 1) {
          cut.local() += hypergraph.edgeWeight(he);
        }
      });
      return cut.combine(std::plus<>()) / (Hypergraph::is_graph ? 2 : 1);
    } else {
      HyperedgeWeight cut = 0;
      for (const HyperedgeID& he : hypergraph.edges()) {
        if (hypergraph.connectivity(he) > 1) {
          cut += hypergraph.edgeWeight(he);
        }
      }
      return cut / (Hypergraph::is_graph ? 2 : 1);
    }
  }

  HyperedgeWeight km1(const PartitionedHypergraph& hypergraph, const bool parallel) {
    if ( parallel ) {
      tbb::enumerable_thread_specific<HyperedgeWeight> km1(0);
      tbb::parallel_for(ID(0), hypergraph.initialNumEdges(), [&](const HyperedgeID he) {
        if (hypergraph.edgeIsEnabled(he)) {
          km1.local() += std::max(hypergraph.connectivity(he) - 1, 0) * hypergraph.edgeWeight(he);
        }
      });
      return km1.combine(std::plus<>()) / (Hypergraph::is_graph ? 2 : 1);
    } else {
      HyperedgeWeight km1 = 0;
      for (const HyperedgeID& he : hypergraph.edges()) {
        km1 += std::max(hypergraph.connectivity(he) - 1, 0) * hypergraph.edgeWeight(he);
      }
      return km1 / (Hypergraph::is_graph ? 2 : 1);
    }
  }

  HyperedgeWeight soed(const PartitionedHypergraph& hypergraph, const bool parallel) {
    if ( parallel ) {
      tbb::enumerable_thread_specific<HyperedgeWeight> soed(0);
      tbb::parallel_for(ID(0), hypergraph.initialNumEdges(), [&](const HyperedgeID he) {
        if ( hypergraph.edgeIsEnabled(he) ) {
          PartitionID connectivity = hypergraph.connectivity(he);
          if (connectivity > 1) {
            soed.local() += connectivity * hypergraph.edgeWeight(he);
          }
        }
      });
      return soed.combine(std::plus<>()) / (Hypergraph::is_graph ? 2 : 1);
    } else {
      HyperedgeWeight soed = 0;
      for (const HyperedgeID& he : hypergraph.edges()) {
        PartitionID connectivity = hypergraph.connectivity(he);
        if (connectivity > 1) {
          soed += connectivity * hypergraph.edgeWeight(he);
        }
      }
      return soed / (Hypergraph::is_graph ? 2 : 1);
    }
  }

  bool isBalanced(const PartitionedHypergraph& phg, const Context& context) {
    size_t num_empty_parts = 0;
    for (PartitionID i = 0; i < context.partition.k; ++i) {
      if (phg.partWeight(i) > context.partition.max_part_weights[i]) {
        return false;
      }
      if (phg.partWeight(i) == 0) {
        num_empty_parts++;
      }
    }
    return num_empty_parts <= phg.numRemovedHypernodes();
  }

  HyperedgeWeight objective(const PartitionedHypergraph& hg, const kahypar::Objective& objective,
                                          const bool parallel) {
    switch (objective) {
      case kahypar::Objective::cut: return hyperedgeCut(hg, parallel);
      case kahypar::Objective::km1: return km1(hg, parallel);
      default:
      ERROR("Unknown Objective");
    }
  }

  double imbalance(const PartitionedHypergraph& hypergraph, const Context& context) {
    ASSERT(context.partition.perfect_balance_part_weights.size() == (size_t)context.partition.k);

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

  std::pair<HyperedgeWeight, HyperedgeWeight> minMaxLoad(const PartitionedHypergraph& hypergraph, const bool parallel) {
    vec<HyperedgeWeight> vol(hypergraph.k(), 0);
    if ( parallel ) {
      tbb::enumerable_thread_specific<vec<HyperedgeWeight>> ets_vol(hypergraph.k(), 0);
      tbb::parallel_for(ID(0), hypergraph.initialNumEdges(), [&](const HyperedgeID he) {
        if (hypergraph.edgeIsEnabled(he)) {
          for (const auto p : hypergraph.connectivitySet(he)) {
            ets_vol.local()[p] += hypergraph.edgeWeight(he);
          }
        }
      });
      tbb::parallel_for(ID(0), hypergraph.initialNumNodes(), [&](const HypernodeID hn) {
            ets_vol.local()[hypergraph.partID(hn)] += hypergraph.weightOfDisabledEdges(hn);
      });
      for (const auto &v : ets_vol) {
        for (PartitionID p = 0; p < hypergraph.k(); ++p) {
          vol[p] += v[p];
        }
      }
    } else {
      for (const HyperedgeID& he : hypergraph.edges()) {
        if (hypergraph.edgeIsEnabled(he)) {
          for (const auto p : hypergraph.connectivitySet(he)) {
            vol[p] += hypergraph.edgeWeight(he);
          }
        }
      }
      for (const HypernodeID hn : hypergraph.nodes()) {
            vol[hypergraph.partID(hn)] += hypergraph.weightOfDisabledEdges(hn);
      }
    }
    return std::make_pair(*std::min_element(vol.begin(), vol.end()), *std::max_element(vol.begin(), vol.end()));
  }

  HyperedgeWeight judiciousLoad(const PartitionedHypergraph& hypergraph, const bool parallel) {
    HyperedgeWeight max_load = minMaxLoad(hypergraph, parallel).second;
    ASSERT(std::all_of(hypergraph.nodes().begin(), hypergraph.nodes().end(), [&](const HypernodeID &hn) {
      return static_cast<HyperedgeWeight>(hypergraph.nodeDegree(hn)) <= max_load;
    }));
    return max_load;
  }

  HyperedgeWeight minLoad(const PartitionedHypergraph& hypergraph, const bool parallel) {
    return minMaxLoad(hypergraph, parallel).first;
  }

  HyperedgeID maxHnDegree(const PartitionedHypergraph& hypergraph) {
    return hypergraph.nodeDegree(*std::max_element(hypergraph.nodes().begin(), hypergraph.nodes().end(),
                                                   [&hypergraph](const HypernodeID a, const HypernodeID b) {
                                                     return hypergraph.nodeDegree(a) < hypergraph.nodeDegree(b);
                                                   }));
  }

} // namespace mt_kahypar::metrics

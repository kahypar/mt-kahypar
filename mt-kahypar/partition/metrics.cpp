/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include "mt-kahypar/partition/metrics.h"

#include <cmath>
#include <algorithm>

#include "mt-kahypar/definitions.h"

namespace mt_kahypar::metrics {

template<typename PartitionedHypergraph>
HyperedgeWeight hyperedgeCut(const PartitionedHypergraph& hypergraph, const bool parallel) {
  if ( parallel ) {
    tbb::enumerable_thread_specific<HyperedgeWeight> cut(0);
    hypergraph.doParallelForAllEdges([&](const HyperedgeID he) {
      if (hypergraph.connectivity(he) > 1) {
        cut.local() += hypergraph.edgeWeight(he);
      }
    });
    return cut.combine(std::plus<>()) / (PartitionedHypergraph::is_graph ? 2 : 1);
  } else {
    HyperedgeWeight cut = 0;
    for (const HyperedgeID& he : hypergraph.edges()) {
      if (hypergraph.connectivity(he) > 1) {
        cut += hypergraph.edgeWeight(he);
      }
    }
    return cut / (PartitionedHypergraph::is_graph ? 2 : 1);
  }
}

template<typename PartitionedHypergraph>
HyperedgeWeight km1(const PartitionedHypergraph& hypergraph, const bool parallel) {
  if ( parallel ) {
    tbb::enumerable_thread_specific<HyperedgeWeight> km1(0);
    hypergraph.doParallelForAllEdges([&](const HyperedgeID he) {
      km1.local() += std::max(hypergraph.connectivity(he) - 1, 0) * hypergraph.edgeWeight(he);
    });
    return km1.combine(std::plus<>()) / (PartitionedHypergraph::is_graph ? 2 : 1);
  } else {
    HyperedgeWeight km1 = 0;
    for (const HyperedgeID& he : hypergraph.edges()) {
      km1 += std::max(hypergraph.connectivity(he) - 1, 0) * hypergraph.edgeWeight(he);
    }
    return km1 / (PartitionedHypergraph::is_graph ? 2 : 1);
  }
}

template<typename PartitionedHypergraph>
HyperedgeWeight soed(const PartitionedHypergraph& hypergraph, const bool parallel) {
  if ( parallel ) {
    tbb::enumerable_thread_specific<HyperedgeWeight> soed(0);
    hypergraph.doParallelForAllEdges([&](const HyperedgeID he) {
      PartitionID connectivity = hypergraph.connectivity(he);
      if (connectivity > 1) {
        soed.local() += connectivity * hypergraph.edgeWeight(he);
      }
    });
    return soed.combine(std::plus<>()) / (PartitionedHypergraph::is_graph ? 2 : 1);
  } else {
    HyperedgeWeight soed = 0;
    for (const HyperedgeID& he : hypergraph.edges()) {
      PartitionID connectivity = hypergraph.connectivity(he);
      if (connectivity > 1) {
        soed += connectivity * hypergraph.edgeWeight(he);
      }
    }
    return soed / (PartitionedHypergraph::is_graph ? 2 : 1);
  }
}

template<typename PartitionedHypergraph>
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
  return context.partition.preset_type == PresetType::large_k ||
    num_empty_parts <= phg.numRemovedHypernodes();
}

template<typename PartitionedHypergraph>
HyperedgeWeight objective(const PartitionedHypergraph& hg,
                          const Objective& objective,
                          const bool parallel) {
  switch (objective) {
    case Objective::cut: return hyperedgeCut(hg, parallel);
    case Objective::km1: return km1(hg, parallel);
    default:
    ERR("Unknown Objective");
  }
}

template<typename PartitionedHypergraph>
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

namespace {
#define HYPEREDGE_CUT(X) HyperedgeWeight hyperedgeCut(const X& hypergraph, bool parallel)
#define KM1(X) HyperedgeWeight km1(const X& hypergraph, bool parallel)
#define SOED(X) HyperedgeWeight soed(const X& hypergraph, bool parallel)
#define OBJECTIVE(X) HyperedgeWeight objective(const X& hg, const Objective& objective, bool parallel)
#define IS_BALANCED(X) bool isBalanced(const X& phg, const Context& context)
#define IMBALANCE(X) double imbalance(const X& hypergraph, const Context& context)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(HYPEREDGE_CUT)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(KM1)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(SOED)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(OBJECTIVE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(IS_BALANCED)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(IMBALANCE)

} // namespace mt_kahypar::metrics
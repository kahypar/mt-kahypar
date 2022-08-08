/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {

struct Metrics {
  HyperedgeWeight km1;
  HyperedgeWeight cut;
  double imbalance;

  void updateMetric(const HyperedgeWeight value, const Mode mode, const kahypar::Objective objective) {
    if (mode == Mode::recursive_bipartitioning || objective == kahypar::Objective::cut) {
      // in recursive bisection, km1 is also optimized via the cut net metric
      cut = value;
    } else {
      ASSERT(objective == kahypar::Objective::km1);
      km1 = value;
    }
  }

  HyperedgeWeight getMetric(const Mode mode, const kahypar::Objective objective) {
    if (mode == Mode::recursive_bipartitioning || objective == kahypar::Objective::cut) {
      // in recursive bisection, km1 is also optimized via the cut net metric
      return cut;
    } else {
      ASSERT(objective == kahypar::Objective::km1);
      return km1;
    }
  }
};

namespace metrics {

HyperedgeWeight hyperedgeCut(const PartitionedHypergraph& hypergraph, bool parallel = true);

HyperedgeWeight km1(const PartitionedHypergraph& hypergraph, bool parallel = true);

HyperedgeWeight soed(const PartitionedHypergraph& hypergraph, bool parallel = true);

bool isBalanced(const PartitionedHypergraph& phg, const Context& context);

HyperedgeWeight objective(
        const PartitionedHypergraph& hg,
        const kahypar::Objective& objective,
        bool parallel = true);

double imbalance(const PartitionedHypergraph& hypergraph, const Context& context);

}  // namespace metrics
}  // namespace mt_kahypar

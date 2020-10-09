/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "kahypar/partition/metrics.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
class IRefiner {

 public:
  IRefiner(const IRefiner&) = delete;
  IRefiner(IRefiner&&) = delete;
  IRefiner & operator= (const IRefiner &) = delete;
  IRefiner & operator= (IRefiner &&) = delete;

  virtual ~IRefiner() = default;

  void initialize(PartitionedHypergraph& hypergraph) {
    initializeImpl(hypergraph);
  }

  bool refine(PartitionedHypergraph& hypergraph,
              const parallel::scalable_vector<HypernodeID>& refinement_nodes,
              kahypar::Metrics& best_metrics,
              const double time_limit) {
    return refineImpl(hypergraph, refinement_nodes, best_metrics, time_limit);
  }

 protected:
  IRefiner() = default;

 private:
  virtual void initializeImpl(PartitionedHypergraph& hypergraph) = 0;

  virtual bool refineImpl(PartitionedHypergraph& hypergraph,
                          const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                          kahypar::Metrics& best_metrics,
                          const double time_limit) = 0;
};

}  // namespace mt_kahypar

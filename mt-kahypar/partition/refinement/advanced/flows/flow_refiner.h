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

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/advanced/i_advanced_refiner.h"

namespace mt_kahypar {

class FlowRefiner final : public IAdvancedRefiner {

  static constexpr bool debug = false;

 public:
  explicit FlowRefiner(const Hypergraph&,
                       const Context& context) :
    _context(context),
    _num_threads(0) { }

  FlowRefiner(const FlowRefiner&) = delete;
  FlowRefiner(FlowRefiner&&) = delete;
  FlowRefiner & operator= (const FlowRefiner &) = delete;
  FlowRefiner & operator= (FlowRefiner &&) = delete;

  virtual ~FlowRefiner() = default;

 protected:

 private:
  void initializeImpl(const PartitionedHypergraph&) {

  }

  MoveSequence refineImpl(const PartitionedHypergraph&,
                          const vec<HypernodeID>& refinement_nodes) {
     // LOG << V(refinement_nodes.size());
     return MoveSequence { { }, 0 };
  }

  PartitionID maxNumberOfBlocksPerSearchImpl() const {
    return 2;
  }

  void setNumThreadsForSearchImpl(const size_t num_threads) {
    _num_threads = num_threads;
  }

  bool isMaximumProblemSizeReachedImpl(ProblemStats& stats) const {
    return stats.numNodes() >= 200;
  }

  const Context& _context;
  size_t _num_threads;
};
}  // namespace mt_kahypar
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


#include "kahypar/meta/registrar.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/advanced/i_advanced_refiner.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

class ILPRefiner final : public IAdvancedRefiner {

 public:
  explicit ILPRefiner(const Hypergraph&,
                      const Context& context,
                      const TaskGroupID) :
    _context(context),
    _num_threads(0),
    _max_num_blocks(0) { }

  ILPRefiner(const ILPRefiner&) = delete;
  ILPRefiner(ILPRefiner&&) = delete;
  ILPRefiner & operator= (const ILPRefiner &) = delete;
  ILPRefiner & operator= (ILPRefiner &&) = delete;

  virtual ~ILPRefiner() = default;

 protected:

 private:
  void initializeImpl(const PartitionedHypergraph&) {
    _max_num_blocks = utils::Randomize::instance().getRandomInt(
      2, _context.refinement.advanced.ilp.max_allowed_blocks, sched_getcpu());
  }

  MoveSequence refineImpl(const PartitionedHypergraph&,
                          const vec<HypernodeID>&) {
    return MoveSequence { {}, 0 };
  }

  PartitionID maxNumberOfBlocksPerSearchImpl() const {
    ASSERT(_max_num_blocks >= 2);
    ASSERT(_max_num_blocks <= _context.refinement.advanced.ilp.max_allowed_blocks);
    return _max_num_blocks;
  }

  void setNumThreadsForSearchImpl(const size_t num_threads) {
    _num_threads = num_threads;
  }

  bool isMaximumProblemSizeReachedImpl(ProblemStats& stats) const {
    return estimatedNumberOfNonZeros(stats) >=
      _context.refinement.advanced.ilp.max_non_zeros;
  }

  size_t estimatedNumberOfNonZeros(const ProblemStats& stats) const {
    const PartitionID k = stats.numContainedBlocks();
    const HypernodeID num_nodes = stats.numNodes();
    const HypernodeID num_pins = stats.numPins();
    if ( _context.partition.objective == kahypar::Objective::km1 ) {
      return k * ( 2 * num_nodes + 2 * num_pins );
    } else if ( _context.partition.objective == kahypar::Objective::cut ) {
      const HyperedgeID num_edges = stats.numEdges();
      return k * ( 2 * num_nodes + 3 * ( num_pins - num_edges ) );
    } else {
      return std::numeric_limits<size_t>::max();
    }
  }

  const Context& _context;
  size_t _num_threads;
  PartitionID _max_num_blocks;
};
}  // namespace mt_kahypar

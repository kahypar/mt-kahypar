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

#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/partition/refinement/advanced/i_advanced_refiner.h"

namespace mt_kahypar {

class AdvancedRefinerMock final : public IAdvancedRefiner {

 public:
  explicit AdvancedRefinerMock(const Hypergraph&,
                               const Context& context,
                               const TaskGroupID) :
    _context(context) { }

  AdvancedRefinerMock(const AdvancedRefinerMock&) = delete;
  AdvancedRefinerMock(AdvancedRefinerMock&&) = delete;
  AdvancedRefinerMock & operator= (const AdvancedRefinerMock &) = delete;
  AdvancedRefinerMock & operator= (AdvancedRefinerMock &&) = delete;

  virtual ~AdvancedRefinerMock() = default;

 protected:

 private:
  void initializeImpl(const PartitionedHypergraph&) {

  }

  MoveSequence refineImpl(const PartitionedHypergraph&,
                          const vec<HypernodeID>&) {
    return MoveSequence { {}, 0 };
  }

  PartitionID maxNumberOfBlocksPerSearchImpl() const {
    return 4;
  }

  void setNumThreadsForSearchImpl(const int) {

  }

  bool isMaximumProblemSizeReachedImpl(const ProblemStats&) const {
    return true;
  }

  const Context& _context;
};

#define REGISTER_ADVANCED_REFINER(id, refiner)                                                                            \
  static kahypar::meta::Registrar<AdvancedRefinementFactory> register_ ## refiner(                                            \
    id,                                                                                                           \
    [](const Hypergraph& hypergraph, const Context& context, const TaskGroupID task_group_id) -> IAdvancedRefiner* {    \
    return new refiner(hypergraph, context, task_group_id);                                                                      \
  })

REGISTER_ADVANCED_REFINER(AdvancedRefinementAlgorithm::mock, AdvancedRefinerMock);

}  // namespace mt_kahypar

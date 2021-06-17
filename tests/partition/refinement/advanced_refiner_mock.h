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

using RefineFunc = std::function<MoveSequence(const PartitionedHypergraph&, const vec<HypernodeID>&, const size_t)>;
using MaxProblemSizeFunc = std::function<bool(ProblemStats&)>;

class AdvancedRefinerMockControl {

  #define NOOP_REFINE_FUNC [] (const PartitionedHypergraph&, const vec<HypernodeID>&, const size_t) { return MoveSequence { {}, 0, 0 }; }
  #define NOOP_MAX_PROB_SIZE_FUNC [&] (ProblemStats&) { return false; }

 public:
  AdvancedRefinerMockControl(const AdvancedRefinerMockControl&) = delete;
  AdvancedRefinerMockControl & operator= (const AdvancedRefinerMockControl &) = delete;

  AdvancedRefinerMockControl(AdvancedRefinerMockControl&&) = delete;
  AdvancedRefinerMockControl & operator= (AdvancedRefinerMockControl &&) = delete;

  static AdvancedRefinerMockControl& instance() {
    static AdvancedRefinerMockControl instance;
    return instance;
  }

 private:
  explicit AdvancedRefinerMockControl() :
    max_num_blocks(2),
    refine_func(NOOP_REFINE_FUNC),
    max_prob_size_func(NOOP_MAX_PROB_SIZE_FUNC) { }

 public:

  void reset() {
    max_num_blocks = 2;
    refine_func = NOOP_REFINE_FUNC;
    max_prob_size_func = NOOP_MAX_PROB_SIZE_FUNC;
  }

  PartitionID max_num_blocks;
  RefineFunc refine_func;
  MaxProblemSizeFunc max_prob_size_func;
};

class AdvancedRefinerMock final : public IAdvancedRefiner {

 public:
  explicit AdvancedRefinerMock(const Hypergraph&,
                               const Context& context) :
    _context(context),
    _max_num_blocks(AdvancedRefinerMockControl::instance().max_num_blocks),
    _num_threads(0),
    _refine_func(AdvancedRefinerMockControl::instance().refine_func),
    _max_prob_size_func(AdvancedRefinerMockControl::instance().max_prob_size_func)  { }

  AdvancedRefinerMock(const AdvancedRefinerMock&) = delete;
  AdvancedRefinerMock(AdvancedRefinerMock&&) = delete;
  AdvancedRefinerMock & operator= (const AdvancedRefinerMock &) = delete;
  AdvancedRefinerMock & operator= (AdvancedRefinerMock &&) = delete;

  virtual ~AdvancedRefinerMock() = default;

 protected:

 private:
  void initializeImpl(const PartitionedHypergraph&) { }

  MoveSequence refineImpl(const PartitionedHypergraph& phg,
                          const vec<HypernodeID>& refinement_nodes) {
    return _refine_func(phg, refinement_nodes, _num_threads);
  }

  void setBlockPairsImpl(const vec<BlockPair>&) { }

  PartitionID maxNumberOfBlocksPerSearchImpl() const {
    return _max_num_blocks;
  }

  void setNumThreadsForSearchImpl(const size_t num_threads) {
    _num_threads = num_threads;
  }

  bool isMaximumProblemSizeReachedImpl(ProblemStats& stats) const {
    return _max_prob_size_func(stats);
  }

  const Context& _context;
  const PartitionID _max_num_blocks;
  size_t _num_threads;
  RefineFunc _refine_func;
  MaxProblemSizeFunc _max_prob_size_func;
};

#define REGISTER_ADVANCED_REFINER(id, refiner)                                          \
  static kahypar::meta::Registrar<AdvancedRefinementFactory> register_ ## refiner(      \
    id,                                                                                 \
    [](const Hypergraph& hypergraph, const Context& context) -> IAdvancedRefiner* {     \
    return new refiner(hypergraph, context);                                            \
  })

REGISTER_ADVANCED_REFINER(AdvancedRefinementAlgorithm::mock, AdvancedRefinerMock);

}  // namespace mt_kahypar

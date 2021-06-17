/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/partition/refinement/advanced/i_advanced_refiner.h"

namespace mt_kahypar {
class DoNothingAdvancedRefiner final : public IAdvancedRefiner {
 public:
  template <typename ... Args>
  explicit DoNothingAdvancedRefiner(Args&& ...) noexcept { }
  DoNothingAdvancedRefiner(const DoNothingAdvancedRefiner&) = delete;
  DoNothingAdvancedRefiner(DoNothingAdvancedRefiner&&) = delete;
  DoNothingAdvancedRefiner & operator= (const DoNothingAdvancedRefiner &) = delete;
  DoNothingAdvancedRefiner & operator= (DoNothingAdvancedRefiner &&) = delete;

 private:
  void initializeImpl(const PartitionedHypergraph&) override final { }

  MoveSequence refineImpl(const PartitionedHypergraph&,
                  const vec<HypernodeID>&) override final {
    return MoveSequence { {}, 0, 0 };
  }

  void setBlockPairsImpl(const vec<BlockPair>&) {}

  PartitionID maxNumberOfBlocksPerSearchImpl() const {
    return 2;
  }

  void setNumThreadsForSearchImpl(const size_t) {}

  bool isMaximumProblemSizeReachedImpl(ProblemStats& stats) const {
    return stats.numNodes() >= std::numeric_limits<HypernodeID>::max();
  }
};
}  // namespace kahypar

/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "mt-kahypar/partition/refinement/flows/i_flow_refiner.h"

namespace mt_kahypar {
class DoNothingFlowRefiner final : public IFlowRefiner {
 public:
  template <typename ... Args>
  explicit DoNothingFlowRefiner(Args&& ...) noexcept { }
  DoNothingFlowRefiner(const DoNothingFlowRefiner&) = delete;
  DoNothingFlowRefiner(DoNothingFlowRefiner&&) = delete;
  DoNothingFlowRefiner & operator= (const DoNothingFlowRefiner &) = delete;
  DoNothingFlowRefiner & operator= (DoNothingFlowRefiner &&) = delete;

 private:
  void initializeImpl(const PartitionedHypergraph&) override final { }

  MoveSequence refineImpl(const PartitionedHypergraph&,
                          const Subhypergraph&,
                          const HighResClockTimepoint&) override final {
    return MoveSequence { {}, 0 };
  }

  PartitionID maxNumberOfBlocksPerSearchImpl() const {
    return 2;
  }

  void setNumThreadsForSearchImpl(const size_t) {}
};
}  // namespace kahypar

/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2015 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include "mt-kahypar/partition/refinement/i_refiner.h"

namespace mt_kahypar {
class DoNothingRefiner final : public IRefiner {
 public:
  template <typename ... Args>
  explicit DoNothingRefiner(Args&& ...) noexcept { }
  DoNothingRefiner(const DoNothingRefiner&) = delete;
  DoNothingRefiner(DoNothingRefiner&&) = delete;
  DoNothingRefiner & operator= (const DoNothingRefiner &) = delete;
  DoNothingRefiner & operator= (DoNothingRefiner &&) = delete;

 private:
  void initializeImpl(PartitionedHypergraph&) final { }

  bool refineImpl(PartitionedHypergraph&,
                  const parallel::scalable_vector<HypernodeID>&,
                  kahypar::Metrics &,
                  const double) final {
    return false;
  }
};

    class DoNothingAsyncRefiner final : public IAsyncRefiner {
    public:
        template <typename ... Args>
        explicit DoNothingAsyncRefiner(Args&& ...) noexcept { }
        DoNothingAsyncRefiner(const DoNothingAsyncRefiner&) = delete;
        DoNothingAsyncRefiner(DoNothingAsyncRefiner&&) = delete;
        DoNothingAsyncRefiner & operator= (const DoNothingAsyncRefiner &) = delete;
        DoNothingAsyncRefiner & operator= (DoNothingAsyncRefiner &&) = delete;

        int64_t getNumTotalAttemptedMoves() const override {
          return 0;
        }

        int64_t getNumTotalMovedNodes() const override {
          return 0;
        }

    private:

        bool refineImpl(PartitionedHypergraph&,
                        const parallel::scalable_vector<HypernodeID>&,
                        metrics::ThreadSafeMetrics &,
                        const double) final {
            return false;
        }

        void resetForGroup(ds::ContractionGroupID) final {
        }
    };

}  // namespace kahypar

/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Nikolai Maas <nikolai.maas@kit.edu>
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

#pragma once

#include "mt-kahypar/partition/refinement/fm/fm_commons.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/unconstrained_strategy.h"

namespace mt_kahypar {

class CombinedStrategy {
public:
  static constexpr bool uses_gain_cache = true;
  static constexpr bool maintain_gain_cache_between_rounds = true;
  static constexpr bool is_unconstrained = true;

  CombinedStrategy(const Context& context,
                   FMSharedData& sharedData,
                   FMStats& runStats) :
      default_strategy(context, sharedData, runStats),
      unconstrained_strategy(context, sharedData, runStats) { }

  template<typename DispatchedStrategyApplicatorFn>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void applyWithDispatchedStrategy(size_t /*taskID*/, size_t round, DispatchedStrategyApplicatorFn applicator_fn) {
    if ((round % 2) == 1) {
      applicator_fn(static_cast<GainCacheStrategy&>(default_strategy));
    } else {
      applicator_fn(static_cast<UnconstrainedStrategy&>(unconstrained_strategy));
    }
  }

  void changeNumberOfBlocks(const PartitionID new_k) {
    default_strategy.changeNumberOfBlocks(new_k);
    unconstrained_strategy.changeNumberOfBlocks(new_k);
  }

  void memoryConsumption(utils::MemoryTreeNode *parent) const {
    default_strategy.memoryConsumption(parent);
    // TODO
  }

  static bool isUnconstrainedRound(size_t round) {
    return (round % 2) == 0;
  }

 private:
  GainCacheStrategy default_strategy;
  UnconstrainedStrategy unconstrained_strategy;
};

}
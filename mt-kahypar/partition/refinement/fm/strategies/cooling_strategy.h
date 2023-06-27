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

class CoolingStrategy {
public:
  static constexpr bool uses_gain_cache = true;
  static constexpr bool maintain_gain_cache_between_rounds = true;
  static constexpr bool is_unconstrained = true;

  CoolingStrategy(const Context& context,
                   FMSharedData& sharedData,
                   FMStats& runStats) :
      context(context),
      shared_data(sharedData),
      default_strategy(context, sharedData, runStats),
      unconstrained_strategy(context, sharedData, runStats, context.refinement.fm.imbalance_penalty_min),
      improvement_first_round(0),
      disable_unconstrained(false) { }

  template<typename DispatchedStrategyApplicatorFn>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void applyWithDispatchedStrategy(size_t /*taskID*/, const size_t round, DispatchedStrategyApplicatorFn applicator_fn) {
    size_t n_rounds = context.refinement.fm.unconstrained_rounds;
    size_t u_round = round;
    if (context.refinement.fm.activate_unconstrained_dynamically) {
      if (round == 1) {
        // first round: measure improvement of constrained FM
        applicator_fn(static_cast<GainCacheStrategy&>(default_strategy));
        return;
      } else if (round == 2) {
        // second round: measure improvement of unconstrained FM
        ASSERT(!shared_data.unconstrained.disabled);
        improvement_first_round = shared_data.previous_improvement_absolute;
        unconstrained_strategy.setPenaltyFactor(context.refinement.fm.penalty_for_activation_test);
        unconstrained_strategy.setUpperBound(context.refinement.fm.unconstrained_upper_bound);
        applicator_fn(static_cast<UnconstrainedStrategy&>(unconstrained_strategy));
        return;
      } else if (round == 3) {
        // third round: decide whether to use unconstrained FM
        if (shared_data.previous_improvement_absolute < improvement_first_round) {
          shared_data.unconstrained.disabled.store(true, std::memory_order_relaxed);
        }
      }
      u_round -= 2;
      n_rounds = std::min(n_rounds, context.refinement.fm.multitry_rounds - 2);
    }

    auto interpolate = [&](double start, double end) {
      if (u_round == 0) {
        return start;
      }
      double summed = (n_rounds - u_round - 1) * start + u_round * end;
      return summed / static_cast<double>(n_rounds - 1);
    };

    if (u_round < n_rounds && !shared_data.unconstrained.disabled) {
      ASSERT(isUnconstrainedRound(round, context));
      double penalty = interpolate(context.refinement.fm.imbalance_penalty_min,
                                  context.refinement.fm.imbalance_penalty_max);
      unconstrained_strategy.setPenaltyFactor(penalty);
      if (context.refinement.fm.unconstrained_upper_bound >= 1 && context.refinement.fm.unconstrained_upper_bound_min >= 1) {
        double upper_bound = interpolate(context.refinement.fm.unconstrained_upper_bound,
                                        context.refinement.fm.unconstrained_upper_bound_min);
        unconstrained_strategy.setUpperBound(upper_bound);
      }
      applicator_fn(static_cast<UnconstrainedStrategy&>(unconstrained_strategy));
    } else {
      applicator_fn(static_cast<GainCacheStrategy&>(default_strategy));
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

  static bool isUnconstrainedRound(size_t round, const Context& context) {
    if (context.refinement.fm.activate_unconstrained_dynamically) {
      return round == 2 || (round > 2 && round - 2 < context.refinement.fm.unconstrained_rounds);
    } else {
      return round < context.refinement.fm.unconstrained_rounds;
    }
  }

 private:
  const Context& context;
  FMSharedData& shared_data;
  GainCacheStrategy default_strategy;
  UnconstrainedStrategy unconstrained_strategy;
  Gain improvement_first_round;
  bool disable_unconstrained;
};

}
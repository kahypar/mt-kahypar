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

template<typename TypeTraits, typename GainTypes>
class CoolingStrategy: public IFMStrategy {
  using LocalFM = LocalizedKWayFM<TypeTraits, GainTypes>;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

 public:
  CoolingStrategy(const Context& context, const FMSharedData&):
      context(context),
      current_round(0),
      current_penalty(context.refinement.fm.imbalance_penalty_min),
      current_upper_bound(context.refinement.fm.unconstrained_upper_bound),
      absolute_improvement_first_round(kInvalidGain),
      unconstrained_is_enabled(true) { }

 private:
  virtual bool dispatchedFindMovesImpl(localized_k_way_fm_t local_fm, mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                       size_t task_id, size_t num_seeds, size_t round) final {
    LocalFM& my_fm = utils::cast<LocalFM>(local_fm);
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);

    if (round == 0) {
      unconstrained_is_enabled = true;
    }
    if (current_round != round && isUnconstrainedRound(round)
        && current_round.compare_exchange_strong(current_round, round)) {
      updateValues(round);
    }

    if (isUnconstrainedRound(round)) {
      LocalUnconstrainedStrategy local_strategy = my_fm.template initializeDispatchedStrategy<LocalUnconstrainedStrategy>();
      return my_fm.findMoves(local_strategy, phg, task_id, num_seeds);
    } else {
      LocalGainCacheStrategy local_strategy = my_fm.template initializeDispatchedStrategy<LocalGainCacheStrategy>();
      return my_fm.findMoves(local_strategy, phg, task_id, num_seeds);
    }
  }

  virtual bool isUnconstrainedRoundImpl(size_t round) const final {
    if (!unconstrained_is_enabled) {
      return false;
    }
    if (context.refinement.fm.activate_unconstrained_dynamically) {
      return round == 1 || (round > 1 && round - 2 < context.refinement.fm.unconstrained_rounds);
    } else {
      return round < context.refinement.fm.unconstrained_rounds;
    }
  }

  virtual bool includesUnconstrainedImpl() const final {
    return true;
  }

  virtual void reportImprovementImpl(size_t round, Gain absolute_improvement, double relative_improvement) final {
    if (round == 0) {
      absolute_improvement_first_round = absolute_improvement;
    } else if (round == 1
               && context.refinement.fm.activate_unconstrained_dynamically
               && absolute_improvement < absolute_improvement_first_round) {
        // this is the decision point whether unconstrained or constrained FM is used
        unconstrained_is_enabled = false;
    } else if (relative_improvement < context.refinement.fm.unconstrained_min_improvement) {
      unconstrained_is_enabled = false;
    }
  }

  void updateValues(size_t round) {
    ASSERT(unconstrained_is_enabled && isUnconstrainedRound(round));
    if (context.refinement.fm.activate_unconstrained_dynamically) {
      if (round == 1) {
        current_penalty = context.refinement.fm.penalty_for_activation_test;
        current_upper_bound = context.refinement.fm.unconstrained_upper_bound;
      } else if (round > 1 && isUnconstrainedRound(round)) {
        size_t n_rounds = std::min(context.refinement.fm.unconstrained_rounds, context.refinement.fm.multitry_rounds - 2);
        calculateInterpolation(round - 2, n_rounds);
      }
    } else if (isUnconstrainedRound(round)) {
      calculateInterpolation(round, context.refinement.fm.unconstrained_rounds);
    }
  }

  void calculateInterpolation(size_t round, size_t n_rounds) {
    ASSERT(unconstrained_is_enabled && round < context.refinement.fm.multitry_rounds);
    auto interpolate = [&](double start, double end) {
      if (round == 0) {
        return start;
      }
      double summed = (n_rounds - round - 1) * start + round * end;
      return summed / static_cast<double>(n_rounds - 1);
    };

    if (round < n_rounds) {
      // interpolate values for current penalty and upper bound
      current_penalty = interpolate(context.refinement.fm.imbalance_penalty_min,
                                    context.refinement.fm.imbalance_penalty_max);
      if (context.refinement.fm.unconstrained_upper_bound >= 1) {
        if (context.refinement.fm.unconstrained_upper_bound_min >= 1) {
          current_upper_bound = interpolate(context.refinement.fm.unconstrained_upper_bound,
                                            context.refinement.fm.unconstrained_upper_bound_min);
        } else {
          current_upper_bound = context.refinement.fm.unconstrained_upper_bound;
        }
      }
    }
  }

  const Context& context;
  parallel::IntegralAtomicWrapper<size_t> current_round;
  double current_penalty;
  double current_upper_bound;
  Gain absolute_improvement_first_round;
  bool unconstrained_is_enabled;
};

}
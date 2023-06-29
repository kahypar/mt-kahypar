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
#include "mt-kahypar/partition/refinement/fm/localized_kway_fm_core.h"
#include "mt-kahypar/partition/refinement/fm/strategies/i_fm_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/unconstrained_strategy.h"

namespace mt_kahypar {

  template<typename TypeTraits, typename GainTypes>
  class CombinedStrategy: public IFMStrategy {
    using Base = IFMStrategy;

  public:
    using LocalFM = LocalizedKWayFM<TypeTraits, GainTypes>;
    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

    CombinedStrategy(const Context& context, FMSharedData& sharedData):
            Base(context, sharedData) { }

    bool dispatchedFindMoves(LocalFM& local_fm, PartitionedHypergraph& phg, size_t task_id, size_t num_seeds, size_t round) {
      if (isUnconstrainedRound(round)) {
        LocalUnconstrainedStrategy local_strategy = local_fm.template initializeDispatchedStrategy<LocalUnconstrainedStrategy>();
        return local_fm.findMoves(local_strategy, phg, task_id, num_seeds);
      } else {
        LocalGainCacheStrategy local_strategy = local_fm.template initializeDispatchedStrategy<LocalGainCacheStrategy>();
        return local_fm.findMoves(local_strategy, phg, task_id, num_seeds);
      }
    }

  private:
    virtual void findMovesImpl(localized_k_way_fm_t local_fm, mt_kahypar_partitioned_hypergraph_t& phg,
                               size_t num_tasks, size_t num_seeds, size_t round,
                               ds::StreamingVector<HypernodeID>& locally_locked_vertices) final {
      Base::findMovesWithConcreteStrategy<CombinedStrategy>(
              local_fm, phg, num_tasks, num_seeds, round, locally_locked_vertices);
    }

    virtual bool isUnconstrainedRoundImpl(size_t round) const final {
      return (round % 2) == 0;
    }

    virtual bool includesUnconstrainedImpl() const final {
      return true;
    }
  };

}

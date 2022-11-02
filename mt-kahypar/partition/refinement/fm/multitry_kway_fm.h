/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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


#include "mt-kahypar/partition/context.h"

#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/fm/localized_kway_fm_core.h"
#include "mt-kahypar/partition/refinement/fm/global_rollback.h"



namespace mt_kahypar {

template<typename FMStrategy>
class MultiTryKWayFM final : public IRefiner {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;


public:

  MultiTryKWayFM(const Hypergraph& hypergraph,
                 const Context& c) :
    initial_num_nodes(hypergraph.initialNumNodes()),
    context(c),
    sharedData(hypergraph.initialNumNodes(), context),
    globalRollback(hypergraph, context),
    ets_fm([&] { return constructLocalizedKWayFMSearch(); })
  {
    if (context.refinement.fm.obey_minimal_parallelism) {
      sharedData.finishedTasksLimit = std::min(8UL, context.shared_memory.num_threads);
    }
  }

  bool refineImpl(PartitionedHypergraph& phg,
                  const vec<HypernodeID>& refinement_nodes,
                  Metrics& metrics,
                  double time_limit) final ;

  void initializeImpl(PartitionedHypergraph& phg) final ;

  void roundInitialization(PartitionedHypergraph& phg,
                           const vec<HypernodeID>& refinement_nodes);


  LocalizedKWayFM<FMStrategy> constructLocalizedKWayFMSearch() {
    return LocalizedKWayFM<FMStrategy>(context, initial_num_nodes, sharedData);
  }

  static double improvementFraction(Gain gain, HyperedgeWeight old_km1) {
    if (old_km1 == 0)
      return 0;
    else
      return static_cast<double>(gain) / static_cast<double>(old_km1);
  }

  void printMemoryConsumption();

  bool is_initialized = false;
  bool enable_light_fm = false;
  const HypernodeID initial_num_nodes;
  const Context& context;
  FMSharedData sharedData;
  GlobalRollback globalRollback;
  tbb::enumerable_thread_specific<LocalizedKWayFM<FMStrategy>> ets_fm;
};

} // namespace mt_kahypar

/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_queue.h>

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/flows/i_flow_refiner.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

namespace mt_kahypar {

template<typename TypeTraits>
class FlowRefinerAdapter {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

  struct ActiveSearch {
    HighResClockTimepoint start;
    double running_time;
    bool reaches_time_limit;
  };

public:
  explicit FlowRefinerAdapter(const HyperedgeID num_hyperedges,
                              const Context& context) :
    _num_hyperedges(num_hyperedges),
    _context(context),
    _refiner(),
    _search_lock(),
    _num_refinements(0),
    _average_running_time(0.0) {
    for ( size_t i = 0; i < _context.shared_memory.num_threads; ++i ) {
      _refiner.emplace_back(nullptr);
    }
  }

  FlowRefinerAdapter(const FlowRefinerAdapter&) = delete;
  FlowRefinerAdapter(FlowRefinerAdapter&&) = delete;

  FlowRefinerAdapter & operator= (const FlowRefinerAdapter &) = delete;
  FlowRefinerAdapter & operator= (FlowRefinerAdapter &&) = delete;

  void initialize();

  // ! Associates a refiner with a search id.
  // ! Returns true, if there is an idle refiner left.
  void registerNewSearch(const PartitionedHypergraph& phg,
                         const size_t refiner_idx);

  MoveSequence refine(const PartitionedHypergraph& phg,
                      const Subhypergraph& sub_hg,
                      HighResClockTimepoint start,
                      size_t refiner_idx);

  // ! Makes the refiner associated with the corresponding search id
  // ! available again
  void finalizeSearch(double running_time, bool reaches_time_limit);

  double timeLimit() const {
    return shouldSetTimeLimit() ?
      std::max(_context.refinement.flows.time_limit_factor *
        _average_running_time, 0.1) : std::numeric_limits<double>::max();
  }

private:
  std::unique_ptr<IFlowRefiner> initializeRefiner();

  bool shouldSetTimeLimit() const {
    return _num_refinements > static_cast<size_t>(_context.partition.k) &&
      _context.refinement.flows.time_limit_factor > 1.0;
  }

  const HyperedgeID _num_hyperedges;
  const Context& _context;

  // ! Available refiners
  vec<std::unique_ptr<IFlowRefiner>> _refiner;
  SpinLock _search_lock;

  size_t _num_refinements;
  double _average_running_time;

};

}  // namespace kahypar

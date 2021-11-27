/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_queue.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/flows/i_flow_refiner.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

namespace mt_kahypar {

class FlowRefinerAdapter {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;
  static constexpr size_t INVALID_REFINER_IDX = std::numeric_limits<size_t>::max();

  struct ActiveSearch {
    size_t refiner_idx;
    HighResClockTimepoint start;
    double running_time;
    bool reaches_time_limit;
  };

public:
  explicit FlowRefinerAdapter(const Hypergraph& hg,
                                  const Context& context) :
    _hg(hg),
    _context(context),
    _unused_refiners(),
    _refiner(),
    _search_lock(),
    _active_searches(),
    _num_used_threads_lock(),
    _num_used_threads(0),
    _num_refinements(0),
    _average_running_time(0.0) {
    for ( size_t i = 0; i < numAvailableRefiner(); ++i ) {
      _refiner.emplace_back(nullptr);
      _unused_refiners.push(i);
    }
  }

  FlowRefinerAdapter(const FlowRefinerAdapter&) = delete;
  FlowRefinerAdapter(FlowRefinerAdapter&&) = delete;

  FlowRefinerAdapter & operator= (const FlowRefinerAdapter &) = delete;
  FlowRefinerAdapter & operator= (FlowRefinerAdapter &&) = delete;

  // ! Associates a refiner with a search id.
  // ! Returns true, if there is an idle refiner left.
  bool registerNewSearch(const SearchID search_id,
                         const PartitionedHypergraph& phg);

  MoveSequence refine(const SearchID search_id,
                      const PartitionedHypergraph& phg,
                      const Subhypergraph& sub_hg);

  // ! Returns the maximum number of blocks which is allowed to be
  // ! contained in the problem of the refiner associated with
  // ! corresponding search id
  PartitionID maxNumberOfBlocks(const SearchID search_id);

  // ! Makes the refiner associated with the corresponding search id
  // ! available again
  void finalizeSearch(const SearchID search_id);

  void reset();

  size_t numAvailableRefiner() const {
    ASSERT(_context.refinement.flows.num_threads_per_search > 0);
    return _context.shared_memory.num_threads / _context.refinement.flows.num_threads_per_search
      + (_context.shared_memory.num_threads % _context.refinement.flows.num_threads_per_search != 0);
  }

  double runningTime(const SearchID search_id) const {
    ASSERT(static_cast<size_t>(search_id) < _active_searches.size());
    return _active_searches[search_id].running_time;
  }

  double timeLimit() const {
    return shouldSetTimeLimit() ?
      std::max(_context.refinement.flows.time_limit_factor *
        _average_running_time, 0.1) : std::numeric_limits<double>::max();
  }

  // ! Only for testing
  size_t numUsedThreads() const {
    return _num_used_threads.load(std::memory_order_relaxed);
  }

private:
  std::unique_ptr<IFlowRefiner> initializeRefiner();

  bool shouldSetTimeLimit() const {
    return _num_refinements > static_cast<size_t>(_context.partition.k) &&
      _context.refinement.flows.time_limit_factor > 1.0;
  }

  const Hypergraph& _hg;
  const Context& _context;

  // ! Indices of unused refiners
  tbb::concurrent_queue<size_t> _unused_refiners;
  // ! Available refiners
  vec<std::unique_ptr<IFlowRefiner>> _refiner;
  // ! Mapping from search id to refiner
  SpinLock _search_lock;
  tbb::concurrent_vector<ActiveSearch> _active_searches;

  SpinLock _num_used_threads_lock;
  // ! Number of used threads
  CAtomic<size_t> _num_used_threads;

  size_t _num_refinements;
  double _average_running_time;

};

}  // namespace kahypar

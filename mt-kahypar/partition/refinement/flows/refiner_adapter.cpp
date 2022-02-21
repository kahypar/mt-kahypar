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

#include "mt-kahypar/partition/refinement/flows/refiner_adapter.h"

#include "mt-kahypar/partition/factories.h"

namespace mt_kahypar {

namespace {
  #define NOW std::chrono::high_resolution_clock::now()
  #define RUNNING_TIME(X) std::chrono::duration<double>(NOW - X).count();
}

bool FlowRefinerAdapter::registerNewSearch(const SearchID search_id,
                                               const PartitionedHypergraph& phg) {
  bool success = true;
  size_t refiner_idx = INVALID_REFINER_IDX;
  if ( _unused_refiners.try_pop(refiner_idx) ) {
    // Note, search id are usually consecutive starting from 0.
    // However, this function is not called in increasing search id order.
    _search_lock.lock();
    while ( static_cast<size_t>(search_id) >= _active_searches.size() ) {
      _active_searches.push_back(ActiveSearch { INVALID_REFINER_IDX, NOW, 0.0, false });
    }
    _search_lock.unlock();

    if ( !_refiner[refiner_idx] ) {
      // Lazy initialization of refiner
      _refiner[refiner_idx] = initializeRefiner();
    }

    _active_searches[search_id].refiner_idx = refiner_idx;
    _active_searches[search_id].start = NOW;
    _refiner[refiner_idx]->initialize(phg);
    _refiner[refiner_idx]->updateTimeLimit(timeLimit());
  } else {
    success = false;
  }
  return success;
}

MoveSequence FlowRefinerAdapter::refine(const SearchID search_id,
                                            const PartitionedHypergraph& phg,
                                            const Subhypergraph& sub_hg) {
  ASSERT(static_cast<size_t>(search_id) < _active_searches.size());
  ASSERT(_active_searches[search_id].refiner_idx != INVALID_REFINER_IDX);

  // Perform refinement
  const size_t refiner_idx = _active_searches[search_id].refiner_idx;
  const size_t num_free_threads = _threads.acquireFreeThreads();
  _refiner[refiner_idx]->setNumThreadsForSearch(num_free_threads);
  MoveSequence moves = _refiner[refiner_idx]->refine(phg, sub_hg, _active_searches[search_id].start);
  _threads.releaseThreads(num_free_threads);
  _active_searches[search_id].reaches_time_limit = moves.state == MoveSequenceState::TIME_LIMIT;
  return moves;
}

PartitionID FlowRefinerAdapter::maxNumberOfBlocks(const SearchID search_id) {
  ASSERT(static_cast<size_t>(search_id) < _active_searches.size());
  ASSERT(_active_searches[search_id].refiner_idx != INVALID_REFINER_IDX);
  const size_t refiner_idx = _active_searches[search_id].refiner_idx;
  return _refiner[refiner_idx]->maxNumberOfBlocksPerSearch();
}

void FlowRefinerAdapter::finalizeSearch(const SearchID search_id) {
  ASSERT(static_cast<size_t>(search_id) < _active_searches.size());
  const double running_time = RUNNING_TIME(_active_searches[search_id].start);
  _active_searches[search_id].running_time = running_time;

  //Update average running time
  _search_lock.lock();
  if ( !_active_searches[search_id].reaches_time_limit ) {
    _average_running_time = (running_time + _num_refinements *
      _average_running_time) / static_cast<double>(_num_refinements + 1);
    ++_num_refinements;
  }
  _search_lock.unlock();

  // Search position of refiner associated with the search id
  if ( shouldSetTimeLimit() ) {
    for ( size_t idx = 0; idx < _refiner.size(); ++idx ) {
      if ( _refiner[idx] ) {
        _refiner[idx]->updateTimeLimit(timeLimit());
      }
    }
  }

  ASSERT(_active_searches[search_id].refiner_idx != INVALID_REFINER_IDX);
  _unused_refiners.push(_active_searches[search_id].refiner_idx);
  _active_searches[search_id].refiner_idx = INVALID_REFINER_IDX;
}

void FlowRefinerAdapter::initialize(const size_t max_parallelism) {
  _num_parallel_refiners = max_parallelism;
  _threads.num_threads = _context.shared_memory.num_threads;
  _threads.num_parallel_refiners = max_parallelism;
  _threads.num_active_refiners = 0;
  _threads.num_used_threads = 0;

  _unused_refiners.clear();
  for ( size_t i = 0; i < numAvailableRefiner(); ++i ) {
    _unused_refiners.push(i);
  }
  _active_searches.clear();
  _num_refinements = 0;
  _average_running_time = 0.0;
}

std::unique_ptr<IFlowRefiner> FlowRefinerAdapter::initializeRefiner() {
  return FlowRefinementFactory::getInstance().createObject(
    _context.refinement.flows.algorithm, _hg, _context);
}

}

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

#include "mt-kahypar/partition/refinement/advanced/refiner_adapter.h"

#include "mt-kahypar/partition/factories.h"

namespace mt_kahypar {

bool AdvancedRefinerAdapter::registerNewSearch(const SearchID search_id,
                                               const PartitionedHypergraph& phg) {
  bool success = true;
  _refiner_lock.lock();
  if ( _num_unused_refiners > 0 ) {
    // Note, search id are usually consecutive starting from 0.
    // However, this function is not called in increasing search id order.
    while ( static_cast<size_t>(search_id) >= _search_to_refiner.size() ) {
      _search_to_refiner.emplace_back(nullptr);
    }

    const size_t refiner_idx = --_num_unused_refiners;
    if ( !_refiner[refiner_idx] ) {
      // Lazy initialization of refiner
      initializeRefiner(_refiner[refiner_idx]);
    }

    // Associate refiner with search id
    ASSERT(!_search_to_refiner[search_id]);
    _search_to_refiner[search_id] = _refiner[refiner_idx].get();
    _search_to_refiner[search_id]->initialize(phg);
  } else {
    success = false;
  }
  _refiner_lock.unlock();
  return success;
}

MoveSequence AdvancedRefinerAdapter::refine(const SearchID search_id,
                                            const PartitionedHypergraph& phg,
                                            const vec<HypernodeID>& refinement_nodes) {
  ASSERT(static_cast<size_t>(search_id) < _search_to_refiner.size());
  ASSERT(_search_to_refiner[search_id]);

  // Determine number of free threads for current search
  _num_used_threads_lock.lock();
  const size_t num_free_threads = std::min(
    _context.refinement.advanced.num_threads_per_search,
    _context.shared_memory.num_threads - _num_used_threads.load(std::memory_order_relaxed));
  _num_used_threads += num_free_threads;
  _num_used_threads_lock.unlock();

  // Perform refinement
  ASSERT(num_free_threads > 0);
  _search_to_refiner[search_id]->setNumThreadsForSearch(num_free_threads);
  MoveSequence moves = _search_to_refiner[search_id]->refine(phg, refinement_nodes);
  _num_used_threads -= num_free_threads;
  return moves;
}

bool AdvancedRefinerAdapter::isMaximumProblemSizeReached(const SearchID search_id,
                                                         ProblemStats& stats) {
  ASSERT(static_cast<size_t>(search_id) < _search_to_refiner.size());
  ASSERT(_search_to_refiner[search_id]);
  return _search_to_refiner[search_id]->isMaximumProblemSizeReached(stats);
}

PartitionID AdvancedRefinerAdapter::maxNumberOfBlocks(const SearchID search_id) {
  ASSERT(static_cast<size_t>(search_id) < _search_to_refiner.size());
  ASSERT(_search_to_refiner[search_id]);
  return _search_to_refiner[search_id]->maxNumberOfBlocksPerSearch();
}

void AdvancedRefinerAdapter::finalizeSearch(const SearchID search_id) {
  ASSERT(static_cast<size_t>(search_id) < _search_to_refiner.size());
  ASSERT(_search_to_refiner[search_id]);
  _refiner_lock.lock();
  IAdvancedRefiner* refiner = _search_to_refiner[search_id];
  // Search position of refiner associated with the search id
  size_t idx = _num_unused_refiners;
  for ( ; idx < _refiner.size(); ++idx ) {
    if ( _refiner[idx].get() == refiner ) {
      break;
    }
  }
  ASSERT(idx < _refiner.size());
  // Swap refiner associated with the search id to free part
  std::swap(_refiner[_num_unused_refiners++], _refiner[idx]);
  _search_to_refiner[search_id] = nullptr;
  _refiner_lock.unlock();
}

void AdvancedRefinerAdapter::reset() {
  _num_unused_refiners = _refiner.size();
  _search_to_refiner.clear();
  _num_used_threads.store(0, std::memory_order_relaxed);
}

void AdvancedRefinerAdapter::initializeRefiner(std::unique_ptr<IAdvancedRefiner>& refiner) {
  refiner = AdvancedRefinementFactory::getInstance().createObject(
    _context.refinement.advanced.algorithm, _hg, _context);
}

}
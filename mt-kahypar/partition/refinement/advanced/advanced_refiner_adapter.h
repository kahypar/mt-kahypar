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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/advanced/i_advanced_refiner.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

namespace mt_kahypar {

class AdvancedRefinerAdapter {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

public:
  explicit AdvancedRefinerAdapter(const Hypergraph& hg,
                                  const Context& context,
                                  const TaskGroupID task_group_id) :
    _hg(hg),
    _context(context),
    _task_group_id(task_group_id),
    _refiner_lock(),
    _num_unused_refiners(0),
    _refiner(),
    _search_to_refiner(),
    _num_used_threads_lock(),
    _num_used_threads(0) {
    for ( size_t i = 0; i < numAvailableRefiner(); ++i ) {
      _refiner.emplace_back(nullptr);
    }
    _num_unused_refiners = _refiner.size();
  }

  AdvancedRefinerAdapter(const AdvancedRefinerAdapter&) = delete;
  AdvancedRefinerAdapter(AdvancedRefinerAdapter&&) = delete;

  AdvancedRefinerAdapter & operator= (const AdvancedRefinerAdapter &) = delete;
  AdvancedRefinerAdapter & operator= (AdvancedRefinerAdapter &&) = delete;

  // ! Associates a refiner with a search id.
  // ! Returns true, if there is an idle refiner left.
  bool registerNewSearch(const SearchID search_id,
                         const PartitionedHypergraph& phg);

  MoveSequence refine(const SearchID search_id,
                      const PartitionedHypergraph& phg,
                      const vec<HypernodeID>& refinement_nodes);

  // ! Returns wheather or not more nodes can be added to problem
  // ! of the refiner associated with the corresponding search id
  bool isMaximumProblemSizeReached(const SearchID search_id,
                                   const ProblemStats& stats);

  // ! Returns the maximum number of blocks which is allowed to be
  // ! contained in the problem of the refiner associated with
  // ! corresponding search id
  PartitionID maxNumberOfBlocks(const SearchID search_id);

  // ! Makes the refiner associated with the corresponding search id
  // ! available again
  void finalizeSearch(const SearchID search_id);

  void reset();

  size_t numAvailableRefiner() const {
    ASSERT(_context.refinement.advanced.num_threads_per_search > 0);
    return _context.shared_memory.num_threads / _context.refinement.advanced.num_threads_per_search
      + (_context.shared_memory.num_threads % _context.refinement.advanced.num_threads_per_search != 0);
  }

  // ! Only for testing
  size_t numUsedThreads() const {
    return _num_used_threads.load(std::memory_order_relaxed);
  }

private:
  void initializeRefiner(std::unique_ptr<IAdvancedRefiner>& refiner);

  const Hypergraph& _hg;
  const Context& _context;
  const TaskGroupID _task_group_id;

  SpinLock _refiner_lock;
  // ! Number of unused refiners
  size_t _num_unused_refiners;
  // ! Available refiners
  vec<std::unique_ptr<IAdvancedRefiner>> _refiner;
  // ! Mapping from search id to refiner
  vec<IAdvancedRefiner*> _search_to_refiner;

  SpinLock _num_used_threads_lock;
  // ! Number of used threads
  CAtomic<size_t> _num_used_threads;

};

}  // namespace kahypar

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

#include "mt-kahypar/partition/refinement/flows/refiner_adapter.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar {

template<typename TypeTraits>
void FlowRefinerAdapter<TypeTraits>::registerNewSearch(const PartitionedHypergraph& phg,
                                                       const size_t refiner_idx) {
  ALWAYS_ASSERT(refiner_idx < _refiner.size());
  if ( !_refiner[refiner_idx] ) {
    // Lazy initialization of refiner
    _refiner[refiner_idx] = initializeRefiner();
  }

  mt_kahypar_partitioned_hypergraph_const_t partitioned_hg = utils::partitioned_hg_const_cast(phg);
  _refiner[refiner_idx]->initialize(partitioned_hg);
  _refiner[refiner_idx]->updateTimeLimit(timeLimit());
}

template<typename TypeTraits>
MoveSequence FlowRefinerAdapter<TypeTraits>::refine(const PartitionedHypergraph& phg,
                                                    const Subhypergraph& sub_hg,
                                                    HighResClockTimepoint start,
                                                    size_t refiner_idx) {
  mt_kahypar_partitioned_hypergraph_const_t partitioned_hg = utils::partitioned_hg_const_cast(phg);
  MoveSequence moves = _refiner[refiner_idx]->refine(partitioned_hg, sub_hg, start);
  return moves;
}

template<typename TypeTraits>
void FlowRefinerAdapter<TypeTraits>::finalizeSearch(double running_time, bool reaches_time_limit) {
  // Update average running time
  _search_lock.lock();
  if ( !reaches_time_limit ) {
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
}

template<typename TypeTraits>
void FlowRefinerAdapter<TypeTraits>::initialize() {
  _num_refinements = 0;
  _average_running_time = 0.0;
}

template<typename TypeTraits>
std::unique_ptr<IFlowRefiner> FlowRefinerAdapter<TypeTraits>::initializeRefiner() {
  return FlowRefinementFactory::getInstance().createObject(
    _context.refinement.flows.algorithm, _num_hyperedges, _context);
}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(FlowRefinerAdapter)

}

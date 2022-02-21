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

#pragma once

#include "kahypar/partition/metrics.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/refinement/flows/flow_common.h"

namespace mt_kahypar {

class IFlowRefiner {

 public:
  IFlowRefiner(const IFlowRefiner&) = delete;
  IFlowRefiner(IFlowRefiner&&) = delete;
  IFlowRefiner & operator= (const IFlowRefiner &) = delete;
  IFlowRefiner & operator= (IFlowRefiner &&) = delete;

  virtual ~IFlowRefiner() = default;

  void initialize(const PartitionedHypergraph& hypergraph) {
    initializeImpl(hypergraph);
  }

  MoveSequence refine(const PartitionedHypergraph& hypergraph,
                      const Subhypergraph& sub_hg,
                      const HighResClockTimepoint& start) {
    return refineImpl(hypergraph, sub_hg, start);
  }

  // ! Returns the maximum number of blocks that can be refined
  // ! per search with this refinement algorithm
  PartitionID maxNumberOfBlocksPerSearch() const {
    return maxNumberOfBlocksPerSearchImpl();
  }

  // ! Set the number of threads that is used for the next search
  void setNumThreadsForSearch(const size_t num_threads) {
    setNumThreadsForSearchImpl(num_threads);
  }

  // ! Updates the time limit (in seconds)
  void updateTimeLimit(const double time_limit) {
    _time_limit = time_limit;
  }


 protected:
  IFlowRefiner() = default;

  double _time_limit = 0;

 private:
  virtual void initializeImpl(const PartitionedHypergraph& hypergraph) = 0;

  virtual MoveSequence refineImpl(const PartitionedHypergraph& hypergraph,
                                  const Subhypergraph& sub_hg,
                                  const HighResClockTimepoint& start) = 0;

  virtual PartitionID maxNumberOfBlocksPerSearchImpl() const = 0;

  virtual void setNumThreadsForSearchImpl(const size_t num_threads) = 0;
};

}  // namespace mt_kahypar

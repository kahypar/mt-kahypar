/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include <tbb/concurrent_vector.h>

#include "algorithm/hyperflowcutter.h"
#include "algorithm/sequential_push_relabel.h"
#include "algorithm/parallel_push_relabel.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/flows/i_flow_refiner.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/datastructures/thread_safe_fast_reset_flag_array.h"
#include "mt-kahypar/parallel/stl/scalable_queue.h"
#include "mt-kahypar/partition/refinement/flows/sequential_construction.h"
#include "mt-kahypar/partition/refinement/flows/parallel_construction.h"
#include "mt-kahypar/partition/refinement/flows/flow_hypergraph_builder.h"

namespace mt_kahypar {

class FlowRefiner final : public IFlowRefiner {

  static constexpr bool debug = false;

 public:
  explicit FlowRefiner(const Hypergraph& hg,
                       const Context& context) :
    _phg(nullptr),
    _context(context),
    _num_available_threads(0),
    _block_0(kInvalidPartition),
    _block_1(kInvalidPartition),
    _flow_hg(),
    _sequential_hfc(_flow_hg, context.partition.seed),
    _parallel_hfc(_flow_hg, context.partition.seed),
    _whfc_to_node(),
    _sequential_construction(hg, _flow_hg, _sequential_hfc, context),
    _parallel_construction(hg, _flow_hg, _parallel_hfc, context)
    {
      _sequential_hfc.find_most_balanced = _context.refinement.flows.find_most_balanced_cut;
      _sequential_hfc.timer.active = false;
      _sequential_hfc.forceSequential(true);
      _sequential_hfc.setBulkPiercing(context.refinement.flows.pierce_in_bulk);

      _parallel_hfc.find_most_balanced = _context.refinement.flows.find_most_balanced_cut;
      _parallel_hfc.timer.active = false;
      _parallel_hfc.forceSequential(false);
      _sequential_hfc.setBulkPiercing(context.refinement.flows.pierce_in_bulk);
  }

  FlowRefiner(const FlowRefiner&) = delete;
  FlowRefiner(FlowRefiner&&) = delete;
  FlowRefiner & operator= (const FlowRefiner &) = delete;
  FlowRefiner & operator= (FlowRefiner &&) = delete;

  virtual ~FlowRefiner() = default;

 protected:

 private:
  void initializeImpl(const PartitionedHypergraph& phg) {
    _phg = &phg;
    _time_limit = std::numeric_limits<double>::max();
    _block_0 = kInvalidPartition;
    _block_1 = kInvalidPartition;
    _flow_hg.clear();
    _whfc_to_node.clear();
  }

  MoveSequence refineImpl(const PartitionedHypergraph& phg,
                          const Subhypergraph& sub_hg,
                          const HighResClockTimepoint& start);

  bool runFlowCutter(const FlowProblem& flow_problem,
                     const HighResClockTimepoint& start,
                     bool& time_limit_reached);

  FlowProblem constructFlowHypergraph(const PartitionedHypergraph& phg,
                                      const Subhypergraph& sub_hg);

  PartitionID maxNumberOfBlocksPerSearchImpl() const {
    return 2;
  }

  void setNumThreadsForSearchImpl(const size_t num_threads) {
    _num_available_threads = num_threads;
  }

  bool canHyperedgeBeDropped(const PartitionedHypergraph& phg,
                             const HyperedgeID he) {
    return _context.partition.objective == kahypar::Objective::cut &&
      phg.pinCountInPart(he, _block_0) + phg.pinCountInPart(he, _block_1) < phg.edgeSize(he);
  }

  const PartitionedHypergraph* _phg;
  const Context& _context;
  using IFlowRefiner::_time_limit;
  size_t _num_available_threads;

  mutable PartitionID _block_0;
  mutable PartitionID _block_1;
  FlowHypergraphBuilder _flow_hg;
  whfc::HyperFlowCutter<whfc::SequentialPushRelabel> _sequential_hfc;
  whfc::HyperFlowCutter<whfc::ParallelPushRelabel> _parallel_hfc;

  vec<HypernodeID> _whfc_to_node;
  SequentialConstruction _sequential_construction;
  ParallelConstruction _parallel_construction;
};
}  // namespace mt_kahypar

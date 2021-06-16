/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include "datastructure/flow_hypergraph_builder.h"
#include "algorithm/hyperflowcutter.h"
#include "algorithm/dinic.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/advanced/i_advanced_refiner.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/parallel/stl/scalable_queue.h"

namespace mt_kahypar {

class FlowRefiner final : public IAdvancedRefiner {

  static constexpr bool debug = false;

  struct FlowProblem {
    whfc::Node source;
    whfc::Node sink;
    HyperedgeWeight total_cut;
    HyperedgeWeight non_removable_cut;
    HypernodeWeight weight_of_block_0;
    HypernodeWeight weight_of_block_1;
  };

 public:
  explicit FlowRefiner(const Hypergraph&,
                       const Context& context) :
    _phg(nullptr),
    _context(context),
    _num_threads(0),
    _scaling(1.0 + _context.refinement.advanced.flows.alpha *
      std::min(0.05, _context.partition.epsilon)),
    _block_0(kInvalidPartition),
    _block_1(kInvalidPartition),
    _flow_hg(),
    _hfc(_flow_hg, context.partition.seed),
    _node_to_whfc(),
    _visited_hns(),
    _visited_hes(),
    _tmp_pins(),
    _cut_hes() {
    _hfc.find_most_balanced =
      _context.refinement.advanced.flows.find_most_balanced_cut;
    _hfc.timer.active = false;
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
    _block_0 = kInvalidPartition;
    _block_1 = kInvalidPartition;
    _flow_hg.clear();
    _node_to_whfc.clear();
    _visited_hes.clear();
  }

  MoveSequence refineImpl(const PartitionedHypergraph& phg,
                          const vec<HypernodeID>& refinement_nodes);

  FlowProblem constructFlowHypergraph(const PartitionedHypergraph& phg,
                                      const vec<HypernodeID>& refinement_nodes);

  void determineDistanceFromCut(const PartitionedHypergraph& phg,
                                const whfc::Node source,
                                const whfc::Node sink);

  PartitionID maxNumberOfBlocksPerSearchImpl() const {
    return 2;
  }

  void setNumThreadsForSearchImpl(const size_t num_threads) {
    _num_threads = num_threads;
  }

  bool isMaximumProblemSizeReachedImpl(ProblemStats& stats) const;

  bool canHyperedgeBeDropped(const PartitionedHypergraph& phg,
                             const HyperedgeID he) {
    return _context.partition.objective == kahypar::Objective::cut &&
      phg.pinCountInPart(he, _block_0) + phg.pinCountInPart(he, _block_1) < phg.edgeSize(he);
  }

  const PartitionedHypergraph* _phg;
  const Context& _context;
  size_t _num_threads;
  double _scaling;

  mutable PartitionID _block_0;
  mutable PartitionID _block_1;
  whfc::FlowHypergraphBuilder _flow_hg;
  whfc::HyperFlowCutter<whfc::Dinic> _hfc;

  ds::DynamicSparseMap<HypernodeID, whfc::Node> _node_to_whfc;
  ds::DynamicSparseSet<HypernodeID> _visited_hns;
  ds::DynamicSparseSet<HyperedgeID> _visited_hes;
  vec<HypernodeID> _tmp_pins;
  vec<HyperedgeID> _cut_hes;
};
}  // namespace mt_kahypar
/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@student.kit.edu>
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

#include "mt-kahypar/partition/refinement/i_refiner.h"
#include <mt-kahypar/partition/refinement/judicious/pq_strategy.h>
namespace mt_kahypar {
class JudiciousRefiner final : public IRefiner {

  static constexpr bool debug = false;

public:
  explicit JudiciousRefiner(Hypergraph& hypergraph,
                            const Context& context) :
    _hypergraph(hypergraph),
    _context(context),
    _refinement_nodes(context.partition.k),
    _gain_cache(context, hypergraph.initialNumNodes()),
    _part_loads(static_cast<size_t>(context.partition.k)),
    _neighbor_deduplicator(hypergraph.initialNumNodes(), 0) {}

  JudiciousRefiner(const JudiciousRefiner&) = delete;
  JudiciousRefiner(JudiciousRefiner&&) = delete;

  JudiciousRefiner & operator= (const JudiciousRefiner &) = delete;
  JudiciousRefiner & operator= (JudiciousRefiner &&) = delete;

 private:

  bool refineImpl(PartitionedHypergraph& hypergraph,
                  const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                  kahypar::Metrics& best_metrics,
                  double) final ;

  void initializeImpl(PartitionedHypergraph& phg) final;
  void calculateRefinementNodes(const PartitionedHypergraph& phg, const PartitionID heaviest_part);
  void doRefinement(PartitionedHypergraph& phg, const PartitionID part_id);
  void revertToBestLocalPrefix(PartitionedHypergraph& phg, const size_t bestGainIndex);
  void updateNeighbors(PartitionedHypergraph& phg, const Move& move);
  void reset();
private:
  bool _is_initialized = false;
  Hypergraph& _hypergraph;
  const Context& _context;
  vec<HypernodeID> _refinement_nodes;
  JudiciousGainCache _gain_cache;
  vec<HyperedgeID> _edgesWithGainChanges;
  vec<Move> _moves;
  ds::ExclusiveHandleHeap<ds::MaxHeap<HyperedgeWeight, PartitionID>> _part_loads;
  // ! Used after a move. Stores whether a neighbor of the just moved vertex has already been updated.
  vec<HypernodeID> _neighbor_deduplicator;
  HypernodeID _deduplication_time = 1;
  Gain _total_improvement = 0;
  Gain _best_improvement = 0;
  HyperedgeWeight _last_load = 0;
  bool _reached_lower_bound = false;
};
}

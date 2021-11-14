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

#include "mt-kahypar/partition/refinement/i_refiner.h"
#include <mt-kahypar/partition/refinement/judicious_gain_cache.h>
namespace mt_kahypar {
class JudiciousRefiner final : public IRefiner {
public:
  explicit JudiciousRefiner(Hypergraph& hypergraph,
                            const Context& context) :
    _hypergraph(hypergraph),
    _context(context),
    _refinement_nodes(context.partition.k),
    _gain_cache(context, hypergraph.initialNumNodes()),
    _part_loads(static_cast<size_t>(context.partition.k)),
    _part_load_margin(context.refinement.judicious.part_load_margin),
    _min_load_ratio(context.refinement.judicious.min_load_ratio),
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
  void calculateRefinementNodes(PartitionedHypergraph& phg);
  Gain doRefinement(PartitionedHypergraph& phg, PartitionID part_id);
  void revertToBestLocalPrefix(PartitionedHypergraph& phg, size_t bestGainIndex, bool update_gain_cache = false);
  void updateNeighbors(PartitionedHypergraph& phg, const Move& move);
private:
  bool _is_initialized = false;
  Hypergraph& _hypergraph;
  const Context& _context;
  vec<vec<HypernodeID>> _refinement_nodes;
  JudiciousGainCache _gain_cache;
  vec<HyperedgeID> _edgesWithGainChanges;
  vec<Move> _moves;
  ds::ExclusiveHandleHeap<ds::MaxHeap<HyperedgeWeight, PartitionID>> _part_loads;
  const double _part_load_margin;
  const double _min_load_ratio;
  // ! Used after a move. Stores whether a neighbor of the just moved vertex has already been updated.
  vec<HypernodeID> _neighbor_deduplicator;
  HypernodeID _deduplication_time = 1;
};
}

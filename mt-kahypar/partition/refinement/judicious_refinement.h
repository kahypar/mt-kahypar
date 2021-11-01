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
    _border_nodes(context.partition.k),
    _gain_cache(context, hypergraph.initialNumNodes()),
    _part_weights(static_cast<size_t>(context.partition.k)),
    _move_status(hypergraph.initialNumNodes(), false) {}

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
  void calculateBorderNodes(PartitionedHypergraph& phg);
  Gain doRefinement(PartitionedHypergraph& phg, PartitionID part_id);
  void revertToBestLocalPrefix(PartitionedHypergraph& phg, size_t bestGainIndex);
  void updateNeighbors(PartitionedHypergraph& phg, const Move& move);
private:
  bool _is_initialized = false;
  Hypergraph& _hypergraph;
  const Context& _context;
  vec<vec<HypernodeID>> _border_nodes;
  JudiciousGainCache _gain_cache;
  vec<HyperedgeID> _edgesWithGainChanges;
  vec<Move> _moves;
  ds::ExclusiveHandleHeap<ds::MaxHeap<HypernodeWeight, PartitionID>> _part_weights;
  const double _part_weight_margin = 1.03;
  const double _min_improvement = 1.05;
  vec<bool> _move_status;
};
}

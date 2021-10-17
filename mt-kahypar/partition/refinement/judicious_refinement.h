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
namespace mt_kahypar {
class JudiciousRefiner final : public IRefiner {
public:
  explicit JudiciousRefiner(Hypergraph& hypergraph,
                            const Context& context) :
    _hypergraph(hypergraph),
    _context(context) {}

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
private:
  Hypergraph& _hypergraph;
  const Context& _context;
};
}

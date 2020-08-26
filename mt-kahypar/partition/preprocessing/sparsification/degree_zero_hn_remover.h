/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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
#include "mt-kahypar/datastructures/clustering.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {

class DegreeZeroHypernodeRemover {

 public:
  DegreeZeroHypernodeRemover(const Context& context) :
    _context(context),
    _removed_hns() { }

  DegreeZeroHypernodeRemover(const DegreeZeroHypernodeRemover&) = delete;
  DegreeZeroHypernodeRemover & operator= (const DegreeZeroHypernodeRemover &) = delete;

  DegreeZeroHypernodeRemover(DegreeZeroHypernodeRemover&&) = delete;
  DegreeZeroHypernodeRemover & operator= (DegreeZeroHypernodeRemover &&) = delete;

  // ! Contracts degree-zero vertices to degree-zero supervertices
  // ! We contract sets of degree-zero vertices such that the weight of
  // ! each supervertex is less than or equal than the maximum allowed
  // ! node weight for a vertex during coarsening.
  HypernodeID removeDegreeZeroHypernodes(Hypergraph& hypergraph) {
    HypernodeID num_removed_degree_zero_hypernodes = 0;
    for ( const HypernodeID& hn : hypergraph.nodes()  ) {
      if ( hypergraph.nodeDegree(hn) == 0 && hypergraph.nodeWeight(hn) == 1 ) {
        hypergraph.removeDegreeZeroHypernode(hn);
        _removed_hns.push_back(hn);
        ++num_removed_degree_zero_hypernodes;
      }
    }
    return num_removed_degree_zero_hypernodes;
  }

  // ! Restore degree-zero vertices
  // ! Each removed degree-zero vertex is assigned to the block of its supervertex.
  void restoreDegreeZeroHypernodes(PartitionedHypergraph& hypergraph) {
    parallel::scalable_vector<PartitionID> non_overloaded_blocks(_context.partition.k, 0);
    std::iota(non_overloaded_blocks.begin(), non_overloaded_blocks.end(), 0);
    for ( const HypernodeID& hn : _removed_hns ) {
      ASSERT(!non_overloaded_blocks.empty());
      PartitionID to = non_overloaded_blocks.back();
      while ( hypergraph.partWeight(to) >= _context.partition.max_part_weights[to] ) {
        ASSERT(non_overloaded_blocks.size() > 1);
        non_overloaded_blocks.pop_back();
        to = non_overloaded_blocks.back();
      }
      hypergraph.restoreDegreeZeroHypernode(hn);
      hypergraph.setNodePart(hn, to);
    }
  }

  void reset() {
    _removed_hns.clear();
  }

 private:
  const Context& _context;
  parallel::scalable_vector<HypernodeID> _removed_hns;
};

}  // namespace mt_kahypar

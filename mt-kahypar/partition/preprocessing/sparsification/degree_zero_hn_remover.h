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

#include "tbb/parallel_sort.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/datastructures/streaming_vector.h"

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

  // ! Remove all degree zero vertices
  HypernodeID removeDegreeZeroHypernodes(Hypergraph& hypergraph) {
    const HypernodeID current_num_nodes =
      hypergraph.initialNumNodes() - hypergraph.numRemovedHypernodes();
    HypernodeID num_removed_degree_zero_hypernodes = 0;
    for ( const HypernodeID& hn : hypergraph.nodes()  ) {
      if ( current_num_nodes - num_removed_degree_zero_hypernodes <= _context.coarsening.contraction_limit) {
        break;
      }
      if ( hypergraph.nodeDegree(hn) == 0 ) {
        hypergraph.removeDegreeZeroHypernode(hn);
        _removed_hns.push_back(hn);
        ++num_removed_degree_zero_hypernodes;
      }
    }
    return num_removed_degree_zero_hypernodes;
  }

  // ! Restore degree-zero vertices
  void restoreDegreeZeroHypernodes(PartitionedHypergraph& hypergraph) {
    // Sort degree-zero vertices in decreasing order of their weight
    tbb::parallel_sort(_removed_hns.begin(), _removed_hns.end(),
      [&](const HypernodeID& lhs, const HypernodeID& rhs) {
        return hypergraph.nodeWeight(lhs) > hypergraph.nodeWeight(rhs);
      });
    // Sort blocks of partition in increasing order of their weight
    parallel::scalable_vector<PartitionID> blocks(_context.partition.k, 0);
    std::iota(blocks.begin(), blocks.end(), 0);
    std::sort(blocks.begin(), blocks.end(),
      [&](const PartitionID& lhs, const PartitionID& rhs) {
        return hypergraph.partWeight(lhs) < hypergraph.partWeight(rhs);
      });

    // Perform Bin-Packing
    for ( const HypernodeID& hn : _removed_hns ) {
      PartitionID to = blocks.front();
      hypergraph.restoreDegreeZeroHypernode(hn, to);
      PartitionID i = 0;
      while ( i + 1 < _context.partition.k &&
              hypergraph.partWeight(blocks[i]) > hypergraph.partWeight(blocks[i + 1]) ) {
        std::swap(blocks[i], blocks[i + 1]);
        ++i;
      }
    }
    _removed_hns.clear();
  }

 private:
  const Context& _context;
  parallel::scalable_vector<HypernodeID> _removed_hns;
};

}  // namespace mt_kahypar

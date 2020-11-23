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
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {

class LargeHyperedgeRemover {

 #define LARGE_HE_THRESHOLD ID(50000)

 public:
  LargeHyperedgeRemover(const Context& context) :
    _context(context),
    _removed_hes() { }

  LargeHyperedgeRemover(const LargeHyperedgeRemover&) = delete;
  LargeHyperedgeRemover & operator= (const LargeHyperedgeRemover &) = delete;

  LargeHyperedgeRemover(LargeHyperedgeRemover&&) = delete;
  LargeHyperedgeRemover & operator= (LargeHyperedgeRemover &&) = delete;

  // ! Removes large hyperedges from the hypergraph
  // ! Returns the number of removed large hyperedges.
  HypernodeID removeLargeHyperedges(Hypergraph& hypergraph) {
    HypernodeID num_removed_large_hyperedges = 0;
    for ( const HyperedgeID& he : hypergraph.edges() ) {
      if ( hypergraph.edgeSize(he) > largeHyperedgeThreshold() ) {
        hypergraph.removeLargeEdge(he);
        _removed_hes.push_back(he);
        ++num_removed_large_hyperedges;
      }
    }
    std::reverse(_removed_hes.begin(), _removed_hes.end());
    return num_removed_large_hyperedges;
  }

  // ! Restores all previously removed large hyperedges
  void restoreLargeHyperedges(PartitionedHypergraph& hypergraph) {
    HyperedgeWeight delta = 0;
    for ( const HyperedgeID& he : _removed_hes ) {
      hypergraph.restoreLargeEdge(he);
      if ( _context.partition.objective == kahypar::Objective::cut ) {
        delta += (hypergraph.connectivity(he) > 1 ? hypergraph.edgeWeight(he) : 0);
       } else {
         delta += (hypergraph.connectivity(he) - 1) * hypergraph.edgeWeight(he);
       }
    }

    if ( _context.partition.verbose_output && delta > 0 ) {
      LOG << RED << "Restoring of" << _removed_hes.size() << "large hyperedges (|e| >"
          << largeHyperedgeThreshold() << ") increased" << _context.partition.objective
          << "by" << delta << END;
    }
  }

  HypernodeID largeHyperedgeThreshold() const {
    return std::max(_context.partition.large_hyperedge_size_threshold, LARGE_HE_THRESHOLD);
  }

  void reset() {
    _removed_hes.clear();
  }

 private:
  const Context& _context;
  parallel::scalable_vector<HypernodeID> _removed_hes;
};

}  // namespace mt_kahypar

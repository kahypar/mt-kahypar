/*******************************************************************************
 * This file is part of KaHyPar.
 *
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


#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"


namespace mt_kahypar {

/*!
 * RECURSIVE BISECTION INITIAL PARTITIONER
 * The recursive bisection initial partitioner starts by performing a parallel multilevel bisection.
 * Once the hypergraph is bisected both blocks are partitioned recursively in parallel until the
 * desired number of blocks are reached.
*/

class RecursiveBisectionInitialPartitioner : public IInitialPartitioner {
 private:
  static constexpr bool enable_heavy_assert = false;

 public:
  RecursiveBisectionInitialPartitioner(PartitionedHypergraph& hypergraph,
                                        const Context& context);

  RecursiveBisectionInitialPartitioner(const RecursiveBisectionInitialPartitioner&) = delete;
  RecursiveBisectionInitialPartitioner(RecursiveBisectionInitialPartitioner&&) = delete;
  RecursiveBisectionInitialPartitioner & operator= (const RecursiveBisectionInitialPartitioner &) = delete;
  RecursiveBisectionInitialPartitioner & operator= (RecursiveBisectionInitialPartitioner &&) = delete;

 private:
  void initialPartitionImpl() final ;

  PartitionedHypergraph& _hg;
  const Context& _context;
};


}  // namespace mt_kahypar

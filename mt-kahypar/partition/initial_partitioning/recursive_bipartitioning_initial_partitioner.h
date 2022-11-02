/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#pragma once


#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"


namespace mt_kahypar {

/*!
 * RECURSIVE BIPARTITIONING INITIAL PARTITIONER
 * The recursive bipartitioning initial partitioner starts by performing a parallel multilevel bisection.
 * Once the hypergraph is bisected both blocks are partitioned recursively in parallel until the
 * desired number of blocks are reached.
*/

class RecursiveBipartitioningInitialPartitioner : public IInitialPartitioner {
 private:
  static constexpr bool enable_heavy_assert = false;

 public:
  RecursiveBipartitioningInitialPartitioner(PartitionedHypergraph& hypergraph,
                                        const Context& context);

  RecursiveBipartitioningInitialPartitioner(const RecursiveBipartitioningInitialPartitioner&) = delete;
  RecursiveBipartitioningInitialPartitioner(RecursiveBipartitioningInitialPartitioner&&) = delete;
  RecursiveBipartitioningInitialPartitioner & operator= (const RecursiveBipartitioningInitialPartitioner &) = delete;
  RecursiveBipartitioningInitialPartitioner & operator= (RecursiveBipartitioningInitialPartitioner &&) = delete;

 private:
  void initialPartitionImpl() final ;

  PartitionedHypergraph& _hg;
  const Context& _context;
};


}  // namespace mt_kahypar

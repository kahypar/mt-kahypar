/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

namespace mt_kahypar::multilevel {

// ! Performs multilevel partitioning on the given hypergraph
// ! in TBB blocking-style.
PartitionedHypergraph partition(Hypergraph& hypergraph, const Context& context);
// ! Performs multilevel partitioning on the given hypergraph
// ! in TBB continuation-style.
// ! Note, the final partitioned hypergraph is moved into the
// ! passed partitioned hypergraph object.
void partition_async(Hypergraph& hypergraph, PartitionedHypergraph& partitioned_hypergraph,
                     const Context& context, tbb::task* parent);

// ! Performs a multilevel partitioning v-cycle on the given hypergraph
// ! in TBB blocking-style.
void partitionVCycle(Hypergraph& hypergraph, PartitionedHypergraph& partitioned_hypergraph,
                     const Context& context);

}  // namespace mt_kahypar

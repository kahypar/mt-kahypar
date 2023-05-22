/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/partitioned_graph.h"
#include "mt-kahypar/partition/process_mapping/process_graph.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {

class ProcessGraphPartitioner {

  static constexpr bool debug = false;

  using PartitionedGraph = ds::PartitionedGraph<ds::StaticGraph>;

 public:
  static void partition(ProcessGraph& process_graph,
                        const Context& context);

 private:
  ProcessGraphPartitioner() { }

  static void recursive_bisection(PartitionedGraph& graph,
                                  const Context& context,
                                  const PartitionID k0,
                                  const PartitionID k1);

  static void recursively_bisect_block(PartitionedGraph& graph,
                                       const Context& context,
                                       const PartitionID block,
                                       const PartitionID k0,
                                       const PartitionID k1);

  static void brute_force_best_bisection(PartitionedGraph& graph,
                                         const HypernodeWeight weight_block_0,
                                         const HypernodeWeight weight_block_1);
};

}  // namespace kahypar

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
#include "mt-kahypar/partition/process_mapping/target_graph.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {

template<typename CommunicationHypergraph>
class DualBipartitioning {

  static constexpr bool debug = false;

  using PartitionedGraph = ds::PartitionedGraph<ds::StaticGraph>;

 public:
  // ! This function assumes that the communication hypergraph is already
  // ! partitioned into k blocks via recursive bipartitioning. It then also
  // ! recursively bipartitions the target graph into k blocks and then
  // ! applies the partition of the target graph to the communication hypergraph.
  static void mapToTargetGraph(CommunicationHypergraph& communication_hg,
                                const TargetGraph& target_graph,
                                const Context& context);

 private:
  DualBipartitioning() { }

  // ! Recursively bisection the target graph into (k1 - k0) blocks.
  // ! The final partition corresponds to block [k0, ..., k1) of the
  // ! input target graph.
  static void recursive_bisection(PartitionedGraph& target_graph,
                                  const Context& context,
                                  const PartitionID k0,
                                  const PartitionID k1);

  // ! Recursively bisects block `block` of the target graph into (k1 - k0) blocks.
  static void recursively_bisect_block(PartitionedGraph& target_graph,
                                       const Context& context,
                                       const PartitionID block,
                                       const PartitionID k0,
                                       const PartitionID k1);

  // ! Brutes forces the optimal bisection of the target graph where
  // ! block 0 has weight equals `weight_block_0` and block 1 has weight
  // ! equals `weight_block_1`.
  static void brute_force_optimal_bisection(PartitionedGraph& target_graph,
                                            const HypernodeWeight weight_block_0,
                                            const HypernodeWeight weight_block_1);
};

}  // namespace kahypar

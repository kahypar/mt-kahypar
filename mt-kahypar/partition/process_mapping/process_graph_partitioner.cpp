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

#include "mt-kahypar/partition/process_mapping/process_graph_partitioner.h"

#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/initial_partitioning/pool_initial_partitioner.h"
#include "mt-kahypar/partition/process_mapping/set_enumerator.h"

namespace mt_kahypar {

namespace {
Context setBisectionContext(const Context& context,
                            const PartitionID k) {
  Context b_context(context);
  b_context.partition.k = 2;
  b_context.partition.epsilon = 0;
  b_context.partition.objective = Objective::cut;
  b_context.partition.gain_policy = GainPolicy::cut_for_graphs;
  b_context.partition.verbose_output = false;
  b_context.partition.instance_type = InstanceType::graph;
  b_context.partition.partition_type = MULTILEVEL_GRAPH_PARTITIONING;
  b_context.initial_partitioning.mode = Mode::direct;
  b_context.type = ContextType::initial_partitioning;
  b_context.refinement = b_context.initial_partitioning.refinement;

  const PartitionID k0 = k / 2 + k % 2;
  const PartitionID k1 = k / 2;
  b_context.partition.perfect_balance_part_weights.assign(2, 0);
  b_context.partition.max_part_weights.assign(2, 0);
  b_context.partition.perfect_balance_part_weights[0] = k0;
  b_context.partition.perfect_balance_part_weights[1] = k1;
  b_context.partition.max_part_weights[0] = k0;
  b_context.partition.max_part_weights[1] = k1;

  return b_context;
}
}

void ProcessGraphPartitioner::partition(ProcessGraph& process_graph,
                                        const Context& context) {
  // Recursively partition process graph
  const PartitionID k = process_graph.numBlocks();
  PartitionedGraph graph(k, process_graph.graph());
  recursive_bisection(graph, context, 0, k);

  // Apply partition on process graph
  for ( const HypernodeID& hn : graph.nodes() ) {
    process_graph.setPartID(hn, graph.partID(hn));
  }
}

void ProcessGraphPartitioner::recursive_bisection(PartitionedGraph& graph,
                                                  const Context& context,
                                                  const PartitionID k0,
                                                  const PartitionID k1) {
  ASSERT(k1 - k0 >= 2);
  ASSERT(k1 - k0 == graph.initialNumNodes());
  const PartitionID k = k1 - k0;

  // Compute bisection of graph
  PartitionedGraph bisection(2, graph.hypergraph());
  Context b_context = setBisectionContext(context, k);
  if ( UL(k) > context.process_mapping.bisection_brute_fore_threshold ) {
    Pool<StaticGraphTypeTraits>::bipartition(bisection, b_context, false);
  } else {
    brute_force_best_bisection(bisection,
      b_context.partition.max_part_weights[0],
      b_context.partition.max_part_weights[1]);
  }

  DBG << V(metrics::quality(bisection, b_context.partition.objective))
      << V(k0) << V(k1) << V(bisection.partWeight(0)) << V(bisection.partWeight(1));

  // Apply bisection to the input hypergraph
  const PartitionID block_0 = 0;
  const PartitionID block_1 = k / 2 + (k % 2);
  for ( const HypernodeID& hn : graph.nodes() ) {
    const PartitionID from = bisection.partID(hn);
    const PartitionID to = from == 0 ? block_0 : block_1;
    ASSERT(from != kInvalidPartition && from < graph.k());
    ASSERT(graph.partID(hn) == kInvalidPartition);
    graph.setOnlyNodePart(hn, to);
  }
  graph.initializePartition();

  ASSERT(metrics::quality(bisection, b_context) == metrics::quality(graph, b_context));
  // Recursively bisect each block
  PartitionID rb_k0 = k / 2 + k % 2;
  PartitionID rb_k1 = k / 2;
  if ( rb_k0 >= 2 ) {
    recursively_bisect_block(graph, context, block_0, 0, rb_k0);
  }
  if ( rb_k1 >= 2 ) {
    recursively_bisect_block(graph, context, block_1, rb_k0, rb_k0 + rb_k1);
  }
}

void ProcessGraphPartitioner::recursively_bisect_block(PartitionedGraph& graph,
                                                       const Context& context,
                                                       const PartitionID block,
                                                       const PartitionID k0,
                                                       const PartitionID k1) {
  // Extract block of partition
  auto extracted_block = graph.extract(block, nullptr, false, true);
  ds::StaticGraph& extracted_graph = extracted_block.hg;
  auto& hn_mapping = extracted_block.hn_mapping;

  ASSERT(extracted_graph.initialNumNodes() > 0);
  // Recursively bisect graph
  const PartitionID k = k1 - k0;
  PartitionedGraph rb_graph(k, extracted_graph);
  recursive_bisection(rb_graph, context, k0, k1);

  // Apply partition from recursion to input graph
  for ( const HypernodeID& hn : graph.nodes() ) {
    const PartitionID from = graph.partID(hn);
    if ( from == block ) {
      ASSERT(hn < hn_mapping.size());
      const PartitionID to = block + rb_graph.partID(hn_mapping[hn]);
      ASSERT(to != kInvalidPartition && to < graph.k());
      if ( from != to ) {
        graph.changeNodePart(hn, from, to);
      }
    }
  }
}

void ProcessGraphPartitioner::brute_force_best_bisection(PartitionedGraph& graph,
                                                         const HypernodeWeight weight_block_0,
                                                         const HypernodeWeight weight_block_1) {
  unused(weight_block_1);
  ASSERT(weight_block_0 + weight_block_1 == graph.initialNumNodes());
  // Enumerate all bisections
  SetEnumerator bisection_enumerator(graph.initialNumNodes(), weight_block_0);
  ds::Bitset best_bisection(graph.initialNumNodes());
  HyperedgeWeight best_objective = std::numeric_limits<HyperedgeWeight>::max();
  for ( const ds::StaticBitset& bisection : bisection_enumerator ) {
    for ( const HypernodeID& hn : graph.nodes() ) {
      graph.setOnlyNodePart(hn, bisection.isSet(hn) ? 0 : 1);
    }
    graph.initializePartition();
    const HyperedgeWeight objective = metrics::quality(graph, Objective::cut);
    if ( objective < best_objective ) {
      best_bisection = bisection.copy();
      best_objective = objective;
    }
    graph.resetPartition();
  }

  // Apply best bisection
  for ( const HypernodeID& hn : graph.nodes() ) {
    graph.setOnlyNodePart(hn, best_bisection.isSet(hn) ? 0 : 1);
  }
  graph.initializePartition();
  ASSERT(graph.partWeight(0) == weight_block_0);
  ASSERT(graph.partWeight(1) == weight_block_1);
}

}  // namespace kahypar

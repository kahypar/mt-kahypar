/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "libmtkahypargp.h"

#include "tbb/parallel_for.h"
#include "tbb/parallel_invoke.h"

#ifndef USE_GRAPH_PARTITIONER
#define USE_GRAPH_PARTITIONER
#endif

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/partition/partitioner.h"
#include "mt-kahypar/partition/registries/register_memory_pool.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/utilities.h"

namespace {
  template<typename T>
  using vec = mt_kahypar::parallel::scalable_vector<T>;
}

namespace gp {

using Graph = mt_kahypar::Hypergraph;
using PartitionedGraph = mt_kahypar::PartitionedHypergraph;

mt_kahypar_graph_t* mt_kahypar_read_graph_from_file(const char* file_name,
                                                    const mt_kahypar_context_t* context,
                                                    const mt_kahypar_file_format_type_t file_format) {
  Graph* graph = new Graph();
  const mt_kahypar::Context& c = *reinterpret_cast<const mt_kahypar::Context*>(context);
  mt_kahypar::FileFormat format = file_format == HMETIS ?
    mt_kahypar::FileFormat::hMetis : mt_kahypar::FileFormat::Metis;

  *graph = mt_kahypar::io::readInputFile(file_name, format,
    c.preprocessing.stable_construction_of_incident_edges);

  return reinterpret_cast<mt_kahypar_graph_t*>(graph);
}

mt_kahypar_graph_t* mt_kahypar_create_graph(const mt_kahypar_hypernode_id_t num_vertices,
                                            const mt_kahypar_hyperedge_id_t num_edges,
                                            const mt_kahypar_hypernode_id_t* edges,
                                            const mt_kahypar_hyperedge_weight_t* edge_weights,
                                            const mt_kahypar_hypernode_weight_t* vertex_weights) {
  // Transform adjacence array into adjacence list
  vec<std::pair<mt_kahypar::HypernodeID, mt_kahypar::HypernodeID>> edge_vector(num_edges);
  tbb::parallel_for(0UL, num_edges, [&](const mt_kahypar::HyperedgeID& he) {
    edge_vector[he] = std::make_pair(edges[2*he], edges[2*he + 1]);
  });

  Graph* graph = new Graph();
  *graph = mt_kahypar::HypergraphFactory::construct_from_graph_edges(
    num_vertices, num_edges, edge_vector, edge_weights, vertex_weights);

  return reinterpret_cast<mt_kahypar_graph_t*>(graph);
}

void mt_kahypar_free_graph(mt_kahypar_graph_t* graph) {
  if (graph == nullptr) {
    return;
  }
  delete reinterpret_cast<Graph*>(graph);
}

mt_kahypar_hypernode_id_t mt_kahypar_num_nodes(mt_kahypar_graph_t* graph) {
  return reinterpret_cast<Graph*>(graph)->initialNumNodes();
}

mt_kahypar_hyperedge_id_t mt_kahypar_num_edges(mt_kahypar_graph_t* graph) {
  return reinterpret_cast<Graph*>(graph)->initialNumEdges() / 2;
}

mt_kahypar_hypernode_id_t mt_kahypar_total_weight(mt_kahypar_graph_t* graph) {
  return reinterpret_cast<Graph*>(graph)->totalWeight();
}

void mt_kahypar_free_partitioned_graph(mt_kahypar_partitioned_graph_t* partitioned_graph) {
  if (partitioned_graph == nullptr) {
    return;
  }
  delete reinterpret_cast<PartitionedGraph*>(partitioned_graph);
}

namespace {
  void prepare_context(mt_kahypar::Context& context) {
    context.partition.mode = mt_kahypar::Mode::direct;
    context.shared_memory.num_threads = mt_kahypar::TBBInitializer::instance().total_number_of_threads();
    context.utility_id = mt_kahypar::utils::Utilities::instance().registerNewUtilityObjects();

    context.partition.perfect_balance_part_weights.clear();
    if ( !context.partition.use_individual_part_weights ) {
      context.partition.max_part_weights.clear();
    }

    if ( context.partition.objective == mt_kahypar::Objective::cut &&
         context.refinement.label_propagation.algorithm ==
         mt_kahypar::LabelPropagationAlgorithm::label_propagation_km1 ) {
      context.refinement.label_propagation.algorithm =
        mt_kahypar::LabelPropagationAlgorithm::label_propagation_cut;
    }

    if ( context.partition.objective == mt_kahypar::Objective::cut &&
         context.initial_partitioning.refinement.label_propagation.algorithm ==
         mt_kahypar::LabelPropagationAlgorithm::label_propagation_km1 ) {
      context.initial_partitioning.refinement.label_propagation.algorithm =
        mt_kahypar::LabelPropagationAlgorithm::label_propagation_cut;
    }

    if ( context.partition.objective == mt_kahypar::Objective::km1 &&
         context.refinement.label_propagation.algorithm ==
         mt_kahypar::LabelPropagationAlgorithm::label_propagation_cut ) {
      context.refinement.label_propagation.algorithm =
        mt_kahypar::LabelPropagationAlgorithm::label_propagation_km1;
    }

    if ( context.partition.objective == mt_kahypar::Objective::km1 &&
         context.initial_partitioning.refinement.label_propagation.algorithm ==
         mt_kahypar::LabelPropagationAlgorithm::label_propagation_cut ) {
      context.initial_partitioning.refinement.label_propagation.algorithm =
        mt_kahypar::LabelPropagationAlgorithm::label_propagation_km1;
    }
  }
}

mt_kahypar_partitioned_graph_t* mt_kahypar_partition(mt_kahypar_graph_t* graph,
                                                     mt_kahypar_context_t* context) {
  Graph& gr = *reinterpret_cast<Graph*>(graph);
  PartitionedGraph* p_graph = new PartitionedGraph();
  mt_kahypar::Context& c = *reinterpret_cast<mt_kahypar::Context*>(context);
  prepare_context(c);
  mt_kahypar::utils::Randomize::instance().setSeed(c.partition.seed);

  // Partition Graph
  *p_graph = mt_kahypar::partition(gr, c);

  return reinterpret_cast<mt_kahypar_partitioned_graph_t*>(p_graph);
}

void mt_kahypar_improve_partition(mt_kahypar_partitioned_graph_t* partitioned_graph,
                                  mt_kahypar_context_t* context,
                                  const size_t num_vcycles) {
  PartitionedGraph& p_graph = *reinterpret_cast<PartitionedGraph*>(partitioned_graph);
  mt_kahypar::Context& c = *reinterpret_cast<mt_kahypar::Context*>(context);
  c.partition.num_vcycles = num_vcycles;
  prepare_context(c);
  mt_kahypar::utils::Randomize::instance().setSeed(c.partition.seed);

  // Perform V-Cycle
  mt_kahypar::partitionVCycle(p_graph, c);
}

mt_kahypar_partitioned_graph_t* mt_kahypar_create_partitioned_graph(mt_kahypar_graph_t* graph,
                                                                    const mt_kahypar_partition_id_t num_blocks,
                                                                    const mt_kahypar_partition_id_t* partition) {
  Graph& gr = *reinterpret_cast<Graph*>(graph);
  PartitionedGraph* partitioned_graph = new PartitionedGraph(num_blocks, gr, mt_kahypar::parallel_tag_t { });

  const mt_kahypar::HypernodeID num_nodes = gr.initialNumNodes();
  tbb::parallel_for(ID(0), num_nodes, [&](const mt_kahypar::HypernodeID& hn) {
    partitioned_graph->setOnlyNodePart(hn, partition[hn]);
  });
  partitioned_graph->initializePartition();

  return reinterpret_cast<mt_kahypar_partitioned_graph_t*>(partitioned_graph);
}

mt_kahypar_partitioned_graph_t* mt_kahypar_read_partition_from_file(mt_kahypar_graph_t* graph,
                                                                    const mt_kahypar_partition_id_t num_blocks,
                                                                    const char* partition_file) {
  std::vector<mt_kahypar::PartitionID> partition;
  mt_kahypar::io::readPartitionFile(partition_file, partition);
  return gp::mt_kahypar_create_partitioned_graph(graph, num_blocks, partition.data());
}

void mt_kahypar_write_partition_to_file(const mt_kahypar_partitioned_graph_t* partitioned_graph,
                                        const char* partition_file) {
  const PartitionedGraph& gr = *reinterpret_cast<const PartitionedGraph*>(partitioned_graph);
  mt_kahypar::io::writePartitionFile(gr, partition_file);
}

void mt_kahypar_get_partition(const mt_kahypar_partitioned_graph_t* partitioned_graph,
                              mt_kahypar_partition_id_t* partition) {
  ASSERT(partition != nullptr);
  const PartitionedGraph* p_graph =
    reinterpret_cast<const PartitionedGraph*>(partitioned_graph);
  p_graph->doParallelForAllNodes([&](const mt_kahypar::HypernodeID& hn) {
    partition[hn] = p_graph->partID(hn);
  });
}

void mt_kahypar_get_block_weights(const mt_kahypar_partitioned_graph_t* partitioned_graph,
                                  mt_kahypar_hypernode_weight_t* block_weights) {
  ASSERT(block_weights != nullptr);
  const PartitionedGraph* p_graph =
    reinterpret_cast<const PartitionedGraph*>(partitioned_graph);
  for ( mt_kahypar::PartitionID i = 0; i < p_graph->k(); ++i ) {
    block_weights[i] = p_graph->partWeight(i);
  }
}

double mt_kahypar_imbalance(const mt_kahypar_partitioned_graph_t* partitioned_graph,
                            const mt_kahypar_context_t* context) {
  return mt_kahypar::metrics::imbalance(
    *reinterpret_cast<const PartitionedGraph*>(partitioned_graph),
    *reinterpret_cast<const mt_kahypar::Context*>(context));
}

mt_kahypar_hyperedge_weight_t mt_kahypar_cut(const mt_kahypar_partitioned_graph_t* partitioned_graph) {
  return mt_kahypar::metrics::hyperedgeCut(
    *reinterpret_cast<const PartitionedGraph*>(partitioned_graph));
}

} // namespace gp
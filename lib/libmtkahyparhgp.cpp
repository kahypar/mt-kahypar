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

#include "libmtkahyparhgp.h"

#include "tbb/parallel_for.h"
#include "tbb/parallel_invoke.h"

#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/definitions.h"
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

namespace hgp {

mt_kahypar_hypergraph_t* mt_kahypar_read_hypergraph_from_file(const char* file_name,
                                                              const mt_kahypar_context_t* context,
                                                              const mt_kahypar_file_format_type_t file_format) {
  mt_kahypar::Hypergraph* hypergraph = new mt_kahypar::Hypergraph();
  const mt_kahypar::Context& c = *reinterpret_cast<const mt_kahypar::Context*>(context);
  mt_kahypar::FileFormat format = file_format == HMETIS ?
    mt_kahypar::FileFormat::hMetis : mt_kahypar::FileFormat::Metis;

  *hypergraph = mt_kahypar::io::readInputFile(file_name, format,
    c.preprocessing.stable_construction_of_incident_edges);

  return reinterpret_cast<mt_kahypar_hypergraph_t*>(hypergraph);
}

mt_kahypar_hypergraph_t* mt_kahypar_create_hypergraph(const mt_kahypar_hypernode_id_t num_vertices,
                                                      const mt_kahypar_hyperedge_id_t num_hyperedges,
                                                      const size_t* hyperedge_indices,
                                                      const mt_kahypar_hyperedge_id_t* hyperedges,
                                                      const mt_kahypar_hyperedge_weight_t* hyperedge_weights,
                                                      const mt_kahypar_hypernode_weight_t* vertex_weights) {
  // Transform adjacence array into adjacence list
  vec<vec<mt_kahypar::HypernodeID>> edge_vector(num_hyperedges);
  tbb::parallel_for(0UL, num_hyperedges, [&](const mt_kahypar::HyperedgeID& he) {
    const size_t num_pins = hyperedge_indices[he + 1] - hyperedge_indices[he];
    edge_vector[he].resize(num_pins);
    for ( size_t i = 0; i < num_pins; ++i ) {
      edge_vector[he][i] = hyperedges[hyperedge_indices[he] + i];
    }
  });

  mt_kahypar::Hypergraph* hypergraph = new mt_kahypar::Hypergraph();
  *hypergraph = mt_kahypar::HypergraphFactory::construct(
    num_vertices, num_hyperedges, edge_vector,
    hyperedge_weights, vertex_weights);

  return reinterpret_cast<mt_kahypar_hypergraph_t*>(hypergraph);
}

void mt_kahypar_free_hypergraph(mt_kahypar_hypergraph_t* hypergraph) {
  if (hypergraph == nullptr) {
    return;
  }
  delete reinterpret_cast<mt_kahypar::Hypergraph*>(hypergraph);
}

mt_kahypar_hypernode_id_t mt_kahypar_num_nodes(mt_kahypar_hypergraph_t* hypergraph) {
  return reinterpret_cast<mt_kahypar::Hypergraph*>(hypergraph)->initialNumNodes();
}

mt_kahypar_hyperedge_id_t mt_kahypar_num_hyperedges(mt_kahypar_hypergraph_t* hypergraph) {
  return reinterpret_cast<mt_kahypar::Hypergraph*>(hypergraph)->initialNumEdges();
}

mt_kahypar_hypernode_id_t mt_kahypar_num_pins(mt_kahypar_hypergraph_t* hypergraph) {
  return reinterpret_cast<mt_kahypar::Hypergraph*>(hypergraph)->initialNumPins();
}

mt_kahypar_hypernode_id_t mt_kahypar_total_weight(mt_kahypar_hypergraph_t* hypergraph) {
  return reinterpret_cast<mt_kahypar::Hypergraph*>(hypergraph)->totalWeight();
}

void mt_kahypar_free_partitioned_hypergraph(mt_kahypar_partitioned_hypergraph_t* partitioned_hg) {
  if (partitioned_hg == nullptr) {
    return;
  }
  delete reinterpret_cast<mt_kahypar::PartitionedHypergraph*>(partitioned_hg);
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

mt_kahypar_partitioned_hypergraph_t* mt_kahypar_partition(mt_kahypar_hypergraph_t* hypergraph,
                                                          mt_kahypar_context_t* context) {
  mt_kahypar::Hypergraph& hg = *reinterpret_cast<mt_kahypar::Hypergraph*>(hypergraph);
  mt_kahypar::PartitionedHypergraph* phg = new mt_kahypar::PartitionedHypergraph();
  mt_kahypar::Context& c = *reinterpret_cast<mt_kahypar::Context*>(context);
  prepare_context(c);
  mt_kahypar::utils::Randomize::instance().setSeed(c.partition.seed);

  // Partition Hypergraph
  *phg = mt_kahypar::partition(hg, c);

  return reinterpret_cast<mt_kahypar_partitioned_hypergraph_t*>(phg);
}

void mt_kahypar_improve_partition(mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                  mt_kahypar_context_t* context,
                                  const size_t num_vcycles) {
  mt_kahypar::PartitionedHypergraph& phg = *reinterpret_cast<mt_kahypar::PartitionedHypergraph*>(partitioned_hg);
  mt_kahypar::Context& c = *reinterpret_cast<mt_kahypar::Context*>(context);
  c.partition.num_vcycles = num_vcycles;
  prepare_context(c);
  mt_kahypar::utils::Randomize::instance().setSeed(c.partition.seed);

  // Perform V-Cycle
  mt_kahypar::partitionVCycle(phg, c);
}

mt_kahypar_partitioned_hypergraph_t* mt_kahypar_create_partitioned_hypergraph(mt_kahypar_hypergraph_t* hypergraph,
                                                                              const mt_kahypar_partition_id_t num_blocks,
                                                                              const mt_kahypar_partition_id_t* partition) {
  mt_kahypar::Hypergraph& hg = *reinterpret_cast<mt_kahypar::Hypergraph*>(hypergraph);
  mt_kahypar::PartitionedHypergraph* partitioned_hg = new mt_kahypar::PartitionedHypergraph(num_blocks, hg, mt_kahypar::parallel_tag_t { });

  const mt_kahypar::HypernodeID num_nodes = hg.initialNumNodes();
  tbb::parallel_for(ID(0), num_nodes, [&](const mt_kahypar::HypernodeID& hn) {
    partitioned_hg->setOnlyNodePart(hn, partition[hn]);
  });
  partitioned_hg->initializePartition();

  return reinterpret_cast<mt_kahypar_partitioned_hypergraph_t*>(partitioned_hg);
}

mt_kahypar_partitioned_hypergraph_t* mt_kahypar_read_partition_from_file(mt_kahypar_hypergraph_t* hypergraph,
                                                                         const mt_kahypar_partition_id_t num_blocks,
                                                                         const char* partition_file) {
  std::vector<mt_kahypar::PartitionID> partition;
  mt_kahypar::io::readPartitionFile(partition_file, partition);
  return hgp::mt_kahypar_create_partitioned_hypergraph(hypergraph, num_blocks, partition.data());
}

void mt_kahypar_write_partition_to_file(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                        const char* partition_file) {
  const mt_kahypar::PartitionedHypergraph& hg = *reinterpret_cast<const mt_kahypar::PartitionedHypergraph*>(partitioned_hg);
  mt_kahypar::io::writePartitionFile(hg, partition_file);
}

void mt_kahypar_get_partition(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                              mt_kahypar_partition_id_t* partition) {
  ASSERT(partition != nullptr);
  const mt_kahypar::PartitionedHypergraph* phg =
    reinterpret_cast<const mt_kahypar::PartitionedHypergraph*>(partitioned_hg);
  phg->doParallelForAllNodes([&](const mt_kahypar::HypernodeID& hn) {
    partition[hn] = phg->partID(hn);
  });
}

void mt_kahypar_get_block_weights(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                  mt_kahypar_hypernode_weight_t* block_weights) {
  ASSERT(block_weights != nullptr);
  const mt_kahypar::PartitionedHypergraph* phg =
    reinterpret_cast<const mt_kahypar::PartitionedHypergraph*>(partitioned_hg);
  for ( mt_kahypar::PartitionID i = 0; i < phg->k(); ++i ) {
    block_weights[i] = phg->partWeight(i);
  }
}

double mt_kahypar_imbalance(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                            const mt_kahypar_context_t* context) {
  return mt_kahypar::metrics::imbalance(
    *reinterpret_cast<const mt_kahypar::PartitionedHypergraph*>(partitioned_hg),
    *reinterpret_cast<const mt_kahypar::Context*>(context));
}

mt_kahypar_hyperedge_weight_t mt_kahypar_cut(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg) {
  return mt_kahypar::metrics::hyperedgeCut(
    *reinterpret_cast<const mt_kahypar::PartitionedHypergraph*>(partitioned_hg));
}

mt_kahypar_hyperedge_weight_t mt_kahypar_km1(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg) {
  return mt_kahypar::metrics::km1(
    *reinterpret_cast<const mt_kahypar::PartitionedHypergraph*>(partitioned_hg));
}

mt_kahypar_hyperedge_weight_t mt_kahypar_soed(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg) {
  return mt_kahypar::metrics::soed(
    *reinterpret_cast<const mt_kahypar::PartitionedHypergraph*>(partitioned_hg));
}

} // namespace hgp
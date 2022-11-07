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

#include "libmtkahypar.h"

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


mt_kahypar_context_t* mt_kahypar_context_new() {
  return reinterpret_cast<mt_kahypar_context_t*>(new mt_kahypar::Context(false));
}

void mt_kahypar_free_context(mt_kahypar_context_t* context) {
  if (context == nullptr) {
    return;
  }
  delete reinterpret_cast<mt_kahypar::Context*>(context);
}

void mt_kahypar_configure_context_from_file(mt_kahypar_context_t* kahypar_context,
                                            const char* ini_file_name) {
  mt_kahypar::parseIniToContext(
    *reinterpret_cast<mt_kahypar::Context*>(kahypar_context), ini_file_name);
}

MT_KAHYPAR_API void mt_kahypar_load_preset(mt_kahypar_context_t* context,
                                           const mt_kahypar_preset_type_t preset) {
  mt_kahypar::Context& c = *reinterpret_cast<mt_kahypar::Context*>(context);
  switch(preset) {
    case DETERMINISTIC:
      c.load_deterministic_preset();
      break;
    case SPEED:
      c.load_default_preset();
      break;
    case HIGH_QUALITY:
      c.load_default_flow_preset();
      break;
  }
}

int mt_kahypar_set_context_parameter(mt_kahypar_context_t* context,
                                     const mt_kahypar_context_parameter_type_t type,
                                     const char* value) {
  mt_kahypar::Context& c = *reinterpret_cast<mt_kahypar::Context*>(context);
  switch(type) {
    case NUM_BLOCKS:
      c.partition.k = atoi(value);
      if ( c.partition.k > 0 ) return 0; /** success **/
      else return 2; /** integer conversion error **/
    case EPSILON:
      c.partition.epsilon = atof(value);
      return 0;
    case OBJECTIVE:
      {
        std::string objective(value);
        if ( objective == "km1" ) {
          c.partition.objective = mt_kahypar::Objective::km1;
          return 0;
        } else if ( objective == "cut" ) {
          c.partition.objective = mt_kahypar::Objective::cut;
          return 0;
        }
        return 3;
      }
    case SEED:
      c.partition.seed = atoi(value);
      return 0;
    case NUM_VCYCLES:
      c.partition.num_vcycles = atoi(value);
      return 0;
    case VERBOSE:
      c.partition.verbose_output = atoi(value);
      return 0;
  }
  return 1; /** no valid parameter type **/
}

void mt_kahypar_set_partitioning_parameters(mt_kahypar_context_t* context,
                                            const mt_kahypar_partition_id_t num_blocks,
                                            const double epsilon,
                                            const mt_kahypar_objective_t objective,
                                            const size_t seed) {
  mt_kahypar::Context& c = *reinterpret_cast<mt_kahypar::Context*>(context);
  c.partition.k = num_blocks;
  c.partition.epsilon = epsilon;
  c.partition.objective = objective == KM1 ?
    mt_kahypar::Objective::km1 : mt_kahypar::Objective::cut;
  c.partition.seed = seed;
}


void mt_kahypar_initialize_thread_pool(const size_t num_threads,
                                       const bool interleaved_allocations) {
  size_t P = num_threads;
  size_t num_available_cpus = mt_kahypar::HardwareTopology::instance().num_cpus();
  if ( num_available_cpus < num_threads ) {
    WARNING("There are currently only" << num_available_cpus << "cpus available."
      << "Setting number of threads from" << num_threads
      << "to" << num_available_cpus);
    P = num_available_cpus;
  }

  // Initialize TBB task arenas on numa nodes
  mt_kahypar::TBBInitializer::instance(P);

  if ( interleaved_allocations ) {
    // We set the membind policy to interleaved allocations in order to
    // distribute allocations evenly across NUMA nodes
    hwloc_cpuset_t cpuset = mt_kahypar::TBBInitializer::instance().used_cpuset();
    mt_kahypar::parallel::HardwareTopology<>::instance().activate_interleaved_membind_policy(cpuset);
    hwloc_bitmap_free(cpuset);
  }
}

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

MT_KAHYPAR_API void mt_kahypar_free_partitioned_hypergraph(mt_kahypar_partitioned_hypergraph_t* partitioned_hg) {
  if (partitioned_hg == nullptr) {
    return;
  }
  delete reinterpret_cast<mt_kahypar::PartitionedHypergraph*>(partitioned_hg);
}

mt_kahypar_partitioned_hypergraph_t* mt_kahypar_partition(mt_kahypar_hypergraph_t* hypergraph,
                                                          mt_kahypar_context_t* context) {
  mt_kahypar::Hypergraph& hg = *reinterpret_cast<mt_kahypar::Hypergraph*>(hypergraph);
  mt_kahypar::PartitionedHypergraph* phg = new mt_kahypar::PartitionedHypergraph();
  mt_kahypar::Context& c = *reinterpret_cast<mt_kahypar::Context*>(context);
  c.partition.mode = mt_kahypar::Mode::direct;
  c.shared_memory.num_threads = mt_kahypar::TBBInitializer::instance().total_number_of_threads();
  c.utility_id = mt_kahypar::utils::Utilities::instance().registerNewUtilityObjects();

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
  c.partition.mode = mt_kahypar::Mode::direct;
  c.partition.num_vcycles = num_vcycles;
  c.shared_memory.num_threads = mt_kahypar::TBBInitializer::instance().total_number_of_threads();
  c.utility_id = mt_kahypar::utils::Utilities::instance().registerNewUtilityObjects();

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
  return mt_kahypar_create_partitioned_hypergraph(hypergraph, num_blocks, partition.data());
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
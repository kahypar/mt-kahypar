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

#include <string>
#include <sstream>
#include <type_traits>

#include "mtkahypartypes.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/conversion.h"
#include "mt-kahypar/partition/partitioner_facade.h"
#include "mt-kahypar/partition/mapping/target_graph.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/registries/registry.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/exception.h"
#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/io/presets.h"


using namespace mt_kahypar;

namespace lib {
  using StaticGraph = typename StaticGraphTypeTraits::Hypergraph;
  using DynamicGraph = typename DynamicGraphTypeTraits::Hypergraph;
  using StaticHypergraph = typename StaticHypergraphTypeTraits::Hypergraph;
  using DynamicHypergraph = typename DynamicHypergraphTypeTraits::Hypergraph;

  using StaticPartitionedGraph = typename StaticGraphTypeTraits::PartitionedHypergraph;
  using DynamicPartitionedGraph = typename DynamicGraphTypeTraits::PartitionedHypergraph;
  using StaticPartitionedHypergraph = typename StaticHypergraphTypeTraits::PartitionedHypergraph;
  using DynamicPartitionedHypergraph = typename DynamicHypergraphTypeTraits::PartitionedHypergraph;
  using SparsePartitionedHypergraph = typename LargeKHypergraphTypeTraits::PartitionedHypergraph;

  using StaticHypergraphFactory = typename ds::StaticHypergraph::Factory;
  using DynamicHypergraphFactory = typename ds::DynamicHypergraph::Factory;
  using StaticGraphFactory = typename ds::StaticGraph::Factory;
  using DynamicGraphFactory = typename ds::DynamicGraph::Factory;


// ####################### General Helper Functions #######################

void initialize(const size_t num_threads, const bool interleaved_allocations, const bool print_warnings) {
  size_t P = num_threads;
  #ifndef KAHYPAR_DISABLE_HWLOC
    size_t num_available_cpus = HardwareTopology::instance().num_cpus();
    if ( num_available_cpus < num_threads ) {
      P = num_available_cpus;
      if (print_warnings) {
        WARNING("There are currently only" << num_available_cpus << "cpus available."
          << "Setting number of threads from" << num_threads << "to" << num_available_cpus);
      }
    }
  #endif

  // Initialize TBB task arenas on numa nodes
  TBBInitializer::instance(P);

  #ifndef KAHYPAR_DISABLE_HWLOC
    if ( interleaved_allocations ) {
      // We set the membind policy to interleaved allocations in order to
      // distribute allocations evenly across NUMA nodes
      hwloc_cpuset_t cpuset = TBBInitializer::instance().used_cpuset();
      parallel::HardwareTopology<>::instance().activate_interleaved_membind_policy(cpuset);
      hwloc_bitmap_free(cpuset);
    }
  #else
    unused(interleaved_allocations);
  #endif

  register_algorithms_and_policies();
}

bool is_compatible(mt_kahypar_hypergraph_t hypergraph, mt_kahypar_preset_type_t preset) {
  switch ( preset ) {
    case DEFAULT:
    case QUALITY:
    case DETERMINISTIC:
    case LARGE_K:
      return hypergraph.type == STATIC_GRAPH || hypergraph.type == STATIC_HYPERGRAPH;
    case HIGHEST_QUALITY:
      return hypergraph.type == DYNAMIC_GRAPH || hypergraph.type == DYNAMIC_HYPERGRAPH;
  }
  return false;
}

bool is_compatible(mt_kahypar_partitioned_hypergraph_t partitioned_hg, mt_kahypar_preset_type_t preset) {
  switch ( preset ) {
    case DEFAULT:
    case QUALITY:
    case DETERMINISTIC:
      return partitioned_hg.type == MULTILEVEL_GRAPH_PARTITIONING ||
             partitioned_hg.type == MULTILEVEL_HYPERGRAPH_PARTITIONING;
    case LARGE_K:
      return partitioned_hg.type == MULTILEVEL_GRAPH_PARTITIONING ||
             partitioned_hg.type == LARGE_K_PARTITIONING;
    case HIGHEST_QUALITY:
      return partitioned_hg.type == N_LEVEL_GRAPH_PARTITIONING ||
             partitioned_hg.type == N_LEVEL_HYPERGRAPH_PARTITIONING;
  }
  return false;
}

void check_if_all_relevant_parameters_are_set(Context& context) {
  bool success = true;
  auto check_parameter = [&](bool is_uninitialized, const char* warning_msg) {
    if (is_uninitialized) {
      success = false;
      if (context.partition.verbose_output) {
        WARNING(warning_msg);
      }
    }
  };

  check_parameter(context.partition.preset_type == PresetType::UNDEFINED, "Preset type not specified.");
  check_parameter(context.partition.k == std::numeric_limits<PartitionID>::max(), "Number of blocks not specified.");
  check_parameter(context.partition.epsilon == std::numeric_limits<double>::max(), "Imbalance not specified.");
  check_parameter(context.partition.objective == Objective::UNDEFINED, "Objective function not specified.");
  if (!success) {
    throw InvalidInputException("A required context parameter is not set. Required are: preset type, k, epsilon, objective");
  }
}

Context context_from_file(const char* ini_file_name) {
  Context context(false);
  parseIniToContext(context, ini_file_name, true);
  return context;
}

Context context_from_preset(PresetType preset) {
  Context context(false);
  auto preset_option_list = loadPreset(preset);
  presetToContext(context, preset_option_list, true);
  return context;
}

void prepare_context(Context& context) {
  context.shared_memory.original_num_threads = mt_kahypar::TBBInitializer::instance().total_number_of_threads();
  context.shared_memory.num_threads = mt_kahypar::TBBInitializer::instance().total_number_of_threads();
  context.utility_id = mt_kahypar::utils::Utilities::instance().registerNewUtilityObjects();

  context.partition.perfect_balance_part_weights.clear();
  if ( !context.partition.use_individual_part_weights ) {
    context.partition.max_part_weights.clear();
  }
}

InstanceType get_instance_type(mt_kahypar_hypergraph_t hypergraph) {
  switch ( hypergraph.type ) {
    case STATIC_GRAPH:
    case DYNAMIC_GRAPH:
      return InstanceType::graph;
    case STATIC_HYPERGRAPH:
    case DYNAMIC_HYPERGRAPH:
      return InstanceType::hypergraph;
    case NULLPTR_HYPERGRAPH:
      return InstanceType::UNDEFINED;
  }
  return InstanceType::UNDEFINED;
}

InstanceType get_instance_type(mt_kahypar_partitioned_hypergraph_t partitioned_hg) {
  switch ( partitioned_hg.type ) {
    case MULTILEVEL_GRAPH_PARTITIONING:
    case N_LEVEL_GRAPH_PARTITIONING:
      return InstanceType::graph;
    case MULTILEVEL_HYPERGRAPH_PARTITIONING:
    case N_LEVEL_HYPERGRAPH_PARTITIONING:
    case LARGE_K_PARTITIONING:
      return InstanceType::hypergraph;
    case NULLPTR_PARTITION:
      return InstanceType::UNDEFINED;
  }
  return InstanceType::UNDEFINED;
}

mt_kahypar_preset_type_t get_preset_c_type(const PresetType preset) {
  switch ( preset ) {
    case PresetType::default_preset: return DEFAULT;
    case PresetType::quality: return QUALITY;
    case PresetType::highest_quality: return HIGHEST_QUALITY;
    case PresetType::deterministic: return DETERMINISTIC;
    case PresetType::large_k: return LARGE_K;
    case PresetType::UNDEFINED: return DEFAULT;
  }
  return DEFAULT;
}

std::string incompatibility_description(mt_kahypar_hypergraph_t hypergraph) {
  std::stringstream ss;
  switch ( hypergraph.type ) {
    case STATIC_GRAPH:
      ss << "The hypergraph uses the static graph data structure which can be only used "
         << "in combination with the following presets: "
         << "DEFAULT, QUALITY, DETERMINISTIC and LARGE_K"; break;
    case DYNAMIC_GRAPH:
      ss << "The hypergraph uses the dynamic graph data structure which can be only used "
         << "in combination with the following preset: "
         << "HIGHEST_QUALITY"; break;
    case STATIC_HYPERGRAPH:
      ss << "The hypergraph uses the static hypergraph data structure which can be only used "
         << "in combination with the following presets: "
         << "DEFAULT, QUALITY, DETERMINISTIC and LARGE_K"; break;
    case DYNAMIC_HYPERGRAPH:
      ss << "The hypergraph uses the dynamic hypergraph data structure which can be only used "
         << "in combination with the following preset: "
         << "HIGHEST_QUALITY"; break;
    case NULLPTR_HYPERGRAPH:
      ss << "The hypergraph holds a nullptr. "
         << "Did you forgot to construct or load a hypergraph?"; break;
  }
  return ss.str();
}

void check_compatibility(mt_kahypar_hypergraph_t hypergraph,
                         mt_kahypar_preset_type_t preset) {
  if ( !is_compatible(hypergraph, preset) ) {
    throw UnsupportedOperationException(incompatibility_description(hypergraph));
  }
}

std::string incompatibility_description(mt_kahypar_partitioned_hypergraph_t partitioned_hg) {
  std::stringstream ss;
  switch ( partitioned_hg.type ) {
    case MULTILEVEL_GRAPH_PARTITIONING:
      ss << "The partitioned hypergraph uses the data structures for multilevel graph partitioning "
         << "which can be only used in combination with the following presets: "
         << "DEFAULT, QUALITY, DETERMINISTIC, and LARGE_K"; break;
    case N_LEVEL_GRAPH_PARTITIONING:
      ss << "The partitioned hypergraph uses the data structures for n-level graph partitioning "
         << "which can be only used in combination with the following preset: "
         << "HIGHEST_QUALITY"; break;
    case MULTILEVEL_HYPERGRAPH_PARTITIONING:
      ss << "The partitioned hypergraph uses the data structures for multilevel hypergraph partitioning "
         << "which can be only used in combination with the following presets: "
         << "DEFAULT, QUALITY, and DETERMINISTIC"; break;
    case N_LEVEL_HYPERGRAPH_PARTITIONING:
      ss << "The partitioned hypergraph uses the data structures for n-level hypergraph partitioning "
         << "which can be only used in combination with the following preset: "
         << "HIGHEST_QUALITY"; break;
    case LARGE_K_PARTITIONING:
      ss << "The partitioned hypergraph uses the data structures for large k hypergraph partitioning "
         << "which can be only used in combination with the following preset: "
         << "LARGE_K"; break;
    case NULLPTR_PARTITION:
      ss << "The hypergraph holds a nullptr. "
         << "Did you forgot to construct or load a hypergraph?"; break;
  }
  return ss.str();
}

void check_compatibility(mt_kahypar_partitioned_hypergraph_t partitioned_hg,
                         mt_kahypar_preset_type_t preset) {
  if ( !is_compatible(partitioned_hg, preset) ) {
    throw UnsupportedOperationException(incompatibility_description(partitioned_hg));
  }
}


mt_kahypar_hypergraph_t create_hypergraph(const Context& context,
                                          const mt_kahypar_hypernode_id_t num_vertices,
                                          const mt_kahypar_hyperedge_id_t num_hyperedges,
                                          const vec<vec<HypernodeID>>& edge_vector,
                                          const mt_kahypar_hyperedge_weight_t* hyperedge_weights,
                                          const mt_kahypar_hypernode_weight_t* vertex_weights) {
  // TODO
  // const bool stable_construction = context.partition.preset_type == PresetType::deterministic ? true : false;
  switch ( context.partition.preset_type ) {
    case PresetType::deterministic:
    case PresetType::large_k:
    case PresetType::default_preset:
    case PresetType::quality:
      return mt_kahypar_hypergraph_t {
        reinterpret_cast<mt_kahypar_hypergraph_s*>(new ds::StaticHypergraph(
          StaticHypergraphFactory::construct(num_vertices, num_hyperedges,
            edge_vector, hyperedge_weights, vertex_weights, true))), STATIC_HYPERGRAPH };
    case PresetType::highest_quality:
      return mt_kahypar_hypergraph_t {
        reinterpret_cast<mt_kahypar_hypergraph_s*>(new ds::DynamicHypergraph(
          DynamicHypergraphFactory::construct(num_vertices, num_hyperedges,
            edge_vector, hyperedge_weights, vertex_weights, true))), DYNAMIC_HYPERGRAPH };
    case PresetType::UNDEFINED:
      break;
  }
  throw InvalidParameterException("Invalid preset type.");
}

mt_kahypar_hypergraph_t create_graph(const Context& context,
                                     const mt_kahypar_hypernode_id_t num_vertices,
                                     const mt_kahypar_hyperedge_id_t num_edges,
                                     const vec<std::pair<HypernodeID, HypernodeID>>& edge_vector,
                                     const mt_kahypar_hyperedge_weight_t* edge_weights,
                                     const mt_kahypar_hypernode_weight_t* vertex_weights) {
  // TODO
  // const bool stable_construction = context.partition.preset_type == PresetType::deterministic ? true : false;
  switch ( context.partition.preset_type ) {
    case PresetType::deterministic:
    case PresetType::large_k:
    case PresetType::default_preset:
    case PresetType::quality:
      return mt_kahypar_hypergraph_t {
        reinterpret_cast<mt_kahypar_hypergraph_s*>(new ds::StaticGraph(
          StaticGraphFactory::construct_from_graph_edges(num_vertices, num_edges,
            edge_vector, edge_weights, vertex_weights, true))), STATIC_GRAPH };
    case PresetType::highest_quality:
      return mt_kahypar_hypergraph_t {
        reinterpret_cast<mt_kahypar_hypergraph_s*>(new ds::DynamicGraph(
          DynamicGraphFactory::construct_from_graph_edges(num_vertices, num_edges,
            edge_vector, edge_weights, vertex_weights, true))), DYNAMIC_GRAPH };
    case PresetType::UNDEFINED:
      break;
  }
  throw InvalidParameterException("Invalid preset type.");
}

template<typename PartitionedHypergraph, typename Hypergraph>
mt_kahypar_partitioned_hypergraph_t create_partitioned_hypergraph(Hypergraph& hg,
                                                                  const mt_kahypar_partition_id_t num_blocks,
                                                                  const mt_kahypar_partition_id_t* partition) {
  PartitionedHypergraph partitioned_hg(num_blocks, hg, parallel_tag_t { });
  const mt_kahypar::HypernodeID num_nodes = hg.initialNumNodes();
  tbb::parallel_for(ID(0), num_nodes, [&](const mt_kahypar::HypernodeID& hn) {
    partitioned_hg.setOnlyNodePart(hn, partition[hn]);
  });
  partitioned_hg.initializePartition();
  return mt_kahypar_partitioned_hypergraph_t { reinterpret_cast<mt_kahypar_partitioned_hypergraph_s*>(
    new PartitionedHypergraph(std::move(partitioned_hg))), PartitionedHypergraph::TYPE };
}

mt_kahypar_partitioned_hypergraph_t create_partitioned_hypergraph(mt_kahypar_hypergraph_t hypergraph,
                                                                  const Context& context,
                                                                  const mt_kahypar_partition_id_t num_blocks,
                                                                  const mt_kahypar_partition_id_t* partition) {
  if ( hypergraph.type == STATIC_GRAPH || hypergraph.type == DYNAMIC_GRAPH ) {
    switch ( context.partition.preset_type ) {
      case PresetType::large_k:
      case PresetType::deterministic:
      case PresetType::default_preset:
      case PresetType::quality:
        ASSERT(hypergraph.type == STATIC_GRAPH);
        return lib::create_partitioned_hypergraph<StaticPartitionedGraph>(
          utils::cast<ds::StaticGraph>(hypergraph), num_blocks, partition);
      case PresetType::highest_quality:
        ASSERT(hypergraph.type == DYNAMIC_GRAPH);
        return lib::create_partitioned_hypergraph<DynamicPartitionedGraph>(
          utils::cast<ds::DynamicGraph>(hypergraph), num_blocks, partition);
      case PresetType::UNDEFINED: break;
    }
  } else {
    switch ( context.partition.preset_type ) {
      case PresetType::large_k:
        ASSERT(hypergraph.type == STATIC_HYPERGRAPH);
        return lib::create_partitioned_hypergraph<SparsePartitionedHypergraph>(
          utils::cast<ds::StaticHypergraph>(hypergraph), num_blocks, partition);
      case PresetType::deterministic:
      case PresetType::default_preset:
      case PresetType::quality:
        ASSERT(hypergraph.type == STATIC_HYPERGRAPH);
        return lib::create_partitioned_hypergraph<StaticPartitionedHypergraph>(
          utils::cast<ds::StaticHypergraph>(hypergraph), num_blocks, partition);
      case PresetType::highest_quality:
        ASSERT(hypergraph.type == DYNAMIC_HYPERGRAPH);
        return lib::create_partitioned_hypergraph<DynamicPartitionedHypergraph>(
          utils::cast<ds::DynamicHypergraph>(hypergraph), num_blocks, partition);
      case PresetType::UNDEFINED: break;
    }
  }
  return mt_kahypar_partitioned_hypergraph_t { nullptr, NULLPTR_PARTITION };
}

template<typename PartitionedHypergraph>
void get_partition(const PartitionedHypergraph& partitioned_hg,
                   mt_kahypar_partition_id_t* partition) {
  ASSERT(partition != nullptr);
  partitioned_hg.doParallelForAllNodes([&](const HypernodeID& hn) {
    partition[hn] = partitioned_hg.partID(hn);
  });
}

template<typename PartitionedHypergraph>
void get_block_weights(const PartitionedHypergraph& partitioned_hg,
                       mt_kahypar_hypernode_weight_t* block_weights) {
  ASSERT(block_weights != nullptr);
  for ( PartitionID i = 0; i < partitioned_hg.k(); ++i ) {
    block_weights[i] = partitioned_hg.partWeight(i);
  }
}

void set_individual_block_weights(Context& context,
                                  const mt_kahypar_partition_id_t num_blocks,
                                  const mt_kahypar_hypernode_weight_t* block_weights) {
  context.partition.use_individual_part_weights = true;
  context.partition.max_part_weights.assign(num_blocks, 0);
  for ( mt_kahypar_partition_id_t i = 0; i < num_blocks; ++i ) {
    context.partition.max_part_weights[i] = block_weights[i];
  }
}

template<typename PartitionedHypergraph>
double imbalance(const PartitionedHypergraph& partitioned_graph, const Context& context) {
  Context c(context);
  c.setupPartWeights(partitioned_graph.totalWeight());
  return metrics::imbalance(partitioned_graph, c);
}

// ####################### Partitioning #######################

mt_kahypar_partitioned_hypergraph_t partitionImpl(mt_kahypar_hypergraph_t hg, Context& context, TargetGraph* target_graph) {
  check_compatibility(hg, lib::get_preset_c_type(context.partition.preset_type));
  check_if_all_relevant_parameters_are_set(context);
  context.partition.instance_type = lib::get_instance_type(hg);
  context.partition.partition_type = to_partition_c_type(context.partition.preset_type, context.partition.instance_type);
  prepare_context(context);
  context.partition.num_vcycles = 0;
  return PartitionerFacade::partition(hg, context, target_graph);
}

mt_kahypar_partitioned_hypergraph_t partition(mt_kahypar_hypergraph_t hg, const Context& context) {
  Context partition_context(context);
  return partitionImpl(hg, partition_context, nullptr);
}

mt_kahypar_partitioned_hypergraph_t map(mt_kahypar_hypergraph_t hg, const ds::StaticGraph& graph, const Context& context) {
  TargetGraph target_graph(graph.copy());
  Context partition_context(context);
  partition_context.partition.objective = Objective::steiner_tree;
  return partitionImpl(hg, partition_context, &target_graph);
}


// ####################### V-Cycles #######################

void improveImpl(mt_kahypar_partitioned_hypergraph_t phg,
                  Context& context,
                  const size_t num_vcycles,
                  TargetGraph* target_graph) {
  check_compatibility(phg, lib::get_preset_c_type(context.partition.preset_type));
  check_if_all_relevant_parameters_are_set(context);
  context.partition.instance_type = lib::get_instance_type(phg);
  context.partition.partition_type = to_partition_c_type(context.partition.preset_type, context.partition.instance_type);
  prepare_context(context);
  context.partition.num_vcycles = num_vcycles;
  PartitionerFacade::improve(phg, context, target_graph);
}

void improve(mt_kahypar_partitioned_hypergraph_t phg, const Context& context, const size_t num_vcycles) {
  Context partition_context(context);
  improveImpl(phg, partition_context, num_vcycles, nullptr);
}

void improveMapping(mt_kahypar_partitioned_hypergraph_t phg,
                    const ds::StaticGraph& graph,
                    const Context& context,
                    const size_t num_vcycles) {
  TargetGraph target_graph(graph.copy());
  Context partition_context(context);
  partition_context.partition.objective = Objective::steiner_tree;
  improveImpl(phg, partition_context, num_vcycles, &target_graph);
}


// ####################### Generic Handling of Different Graph Types #######################

struct NoReturn {};

template<typename ReturnT, bool Throwing, typename Func>
ReturnT switch_hg(mt_kahypar_hypergraph_t hg, Func f) {
  switch ( hg.type ) {
    case STATIC_GRAPH:
      return f(utils::cast<StaticGraph>(hg));
    case DYNAMIC_GRAPH:
      return f(utils::cast<DynamicGraph>(hg));
    case STATIC_HYPERGRAPH:
      return f(utils::cast<StaticHypergraph>(hg));
    case DYNAMIC_HYPERGRAPH:
      return f(utils::cast<DynamicHypergraph>(hg));
    case NULLPTR_HYPERGRAPH: break;
  }
  if constexpr (Throwing) {
    throw UnsupportedOperationException("Input is not a valid hypergraph.");
  }
  return ReturnT{};
}

template<typename ReturnT, bool Throwing, typename Func>
ReturnT switch_graph(mt_kahypar_hypergraph_t hg, Func f) {
  switch ( hg.type ) {
    case STATIC_GRAPH:
      return f(utils::cast<StaticGraph>(hg));
    case DYNAMIC_GRAPH:
      return f(utils::cast<DynamicGraph>(hg));
    case STATIC_HYPERGRAPH:
    case DYNAMIC_HYPERGRAPH:
    case NULLPTR_HYPERGRAPH:
      break;
  }
  if constexpr (Throwing) {
    throw UnsupportedOperationException("Input is not a valid hypergraph.");
  }
  return ReturnT{};
}

template<typename ReturnT, typename Func>
ReturnT switch_phg_throwing_impl(mt_kahypar_partitioned_hypergraph_t phg, Func f) {
  switch ( phg.type ) {
    case MULTILEVEL_GRAPH_PARTITIONING:
      return f(utils::cast<StaticPartitionedGraph>(phg));
    case N_LEVEL_GRAPH_PARTITIONING:
      return f(utils::cast<DynamicPartitionedGraph>(phg));
    case MULTILEVEL_HYPERGRAPH_PARTITIONING:
      return f(utils::cast<StaticPartitionedHypergraph>(phg));
    case N_LEVEL_HYPERGRAPH_PARTITIONING:
      return f(utils::cast<DynamicPartitionedHypergraph>(phg));
    case LARGE_K_PARTITIONING:
      return f(utils::cast<SparsePartitionedHypergraph>(phg));
    case NULLPTR_PARTITION:
      break;
  }
  throw UnsupportedOperationException("Input is not a valid partitioned hypergraph.");
}

template<typename ReturnT = NoReturn, typename Func>
ReturnT switch_phg_throwing(mt_kahypar_partitioned_hypergraph_t phg, Func f) {
  if constexpr ( std::is_same_v<ReturnT, NoReturn> ) {
    return switch_phg_throwing_impl<NoReturn>(phg, [=](auto& phg) {
      f(phg);
      return NoReturn{};
    });
  } else {
    return switch_phg_throwing_impl<ReturnT>(phg, f);
  }
}


// ####################### Generic Implementations #######################

template<bool Throwing>
HypernodeID num_nodes(mt_kahypar_hypergraph_t hypergraph) {
  return switch_hg<PartitionID, Throwing>(hypergraph, [](const auto& hg) {
    return hg.initialNumNodes();
  });
}

template<bool Throwing>
HyperedgeID num_edges(mt_kahypar_hypergraph_t hypergraph) {
  return switch_hg<HyperedgeID, Throwing>(hypergraph, [](const auto& hg) {
    return hg.initialNumEdges();
  });
}

template<bool Throwing>
HypernodeID num_pins(mt_kahypar_hypergraph_t hypergraph) {
  return switch_hg<HypernodeID, Throwing>(hypergraph, [](const auto& hg) {
    return hg.initialNumPins();
  });
}

template<bool Throwing>
HypernodeWeight total_weight(mt_kahypar_hypergraph_t hypergraph) {
  return switch_hg<HypernodeWeight, Throwing>(hypergraph, [](const auto& hg) {
    return hg.totalWeight();
  });
}

template<bool Throwing>
HyperedgeID node_degree(mt_kahypar_hypergraph_t hypergraph, HypernodeID hn) {
  return switch_hg<HyperedgeID, Throwing>(hypergraph, [=](const auto& hg) {
    return hg.nodeDegree(hn);
  });
}

template<bool Throwing>
HypernodeWeight node_weight(mt_kahypar_hypergraph_t hypergraph, HypernodeID hn) {
  return switch_hg<HypernodeWeight, Throwing>(hypergraph, [=](const auto& hg) {
    return hg.nodeWeight(hn);
  });
}

template<bool Throwing>
bool is_fixed(mt_kahypar_hypergraph_t hypergraph, HypernodeID hn) {
  return switch_hg<bool, Throwing>(hypergraph, [=](const auto& hg) {
    return hg.isFixed(hn);
  });
}

template<bool Throwing>
PartitionID fixed_vertex_block(mt_kahypar_hypergraph_t hypergraph, HypernodeID hn) {
  return switch_hg<PartitionID, Throwing>(hypergraph, [=](const auto& hg) {
    return hg.isFixed(hn) ? hg.fixedVertexBlock(hn) : kInvalidPartition;
  });
}

template<bool Throwing>
HypernodeID edge_size(mt_kahypar_hypergraph_t hypergraph, HyperedgeID he) {
  return switch_hg<HypernodeID, Throwing>(hypergraph, [=](const auto& hg) {
    return hg.edgeSize(he);
  });
}

template<bool Throwing>
HyperedgeWeight edge_weight(mt_kahypar_hypergraph_t hypergraph, HyperedgeID he) {
  return switch_hg<HyperedgeWeight, Throwing>(hypergraph, [=](const auto& hg) {
    return hg.edgeWeight(he);
  });
}

} // namespace lib

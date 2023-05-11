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

#include <iostream>
#include <string>

#include "include/libmtkahypartypes.h"
#include "mt-kahypar/macros.h"

namespace mt_kahypar {

enum class Type : int8_t {
  Unweighted = 0,
  EdgeWeights = 1,
  NodeWeights = 10,
  EdgeAndNodeWeights = 11,
};

enum class FileFormat : int8_t {
  hMetis = 0,
  Metis = 1,
};

enum class InstanceType : int8_t {
  ENABLE_GRAPHS(graph = 0 COMMA)
  hypergraph = 1,
  UNDEFINED = 2
};

enum class PresetType : int8_t {
  deterministic,
  ENABLE_LARGE_K(large_k COMMA)
  default_preset,
  default_flows,
  ENABLE_N_LEVEL(quality_preset COMMA)
  ENABLE_N_LEVEL(quality_flows COMMA)
  UNDEFINED
};

enum class ContextType : bool {
  main,
  initial_partitioning
};

enum class Mode : uint8_t {
  recursive_bipartitioning,
  direct,
  deep_multilevel,
  UNDEFINED
};

enum class Objective : uint8_t {
  cut,
  km1,
  soed,
  UNDEFINED
};

enum class GainPolicy : uint8_t {
  km1,
  cut,
  soed,
  ENABLE_GRAPHS(cut_for_graphs COMMA)
  none
};

enum class LouvainEdgeWeight : uint8_t {
  hybrid,
  uniform,
  non_uniform,
  degree,
  UNDEFINED
};

enum class SimiliarNetCombinerStrategy : uint8_t {
  union_nets,
  max_size,
  importance,
  UNDEFINED
};

enum class CoarseningAlgorithm : uint8_t {
  multilevel_coarsener,
  deterministic_multilevel_coarsener,
  ENABLE_N_LEVEL(nlevel_coarsener COMMA)
  UNDEFINED
};

enum class RatingFunction : uint8_t {
  heavy_edge,
  ENABLE_EXPERIMENTAL_FEATURES(sameness COMMA)
  UNDEFINED
};

enum class HeavyNodePenaltyPolicy : uint8_t {
  no_penalty,
  ENABLE_EXPERIMENTAL_FEATURES(multiplicative_penalty COMMA)
  ENABLE_EXPERIMENTAL_FEATURES(additive COMMA)
  UNDEFINED
};

enum class AcceptancePolicy : uint8_t {
  ENABLE_EXPERIMENTAL_FEATURES(best COMMA)
  best_prefer_unmatched,
  UNDEFINED
};

enum class InitialPartitioningAlgorithm : uint8_t {
  greedy_round_robin_fm = 0,
  greedy_global_fm = 1,
  greedy_sequential_fm = 2,
  random = 3,
  bfs = 4,
  label_propagation = 5,
  greedy_round_robin_max_net = 6,
  greedy_global_max_net = 7,
  greedy_sequential_max_net = 8,
  UNDEFINED = 9
};

enum class LabelPropagationAlgorithm : uint8_t {
  label_propagation,
  deterministic,
  do_nothing
};

enum class JetAlgorithm : uint8_t {
  precomputed_ordered,
  greedy_unordered,
  do_nothing
};

enum class FMAlgorithm : uint8_t {
  kway_fm,
  do_nothing
};

enum class FlowAlgorithm : uint8_t {
  flow_cutter,
  mock,
  do_nothing
};

enum class RebalancingAlgorithm : uint8_t {
  simple_rebalancer,
  do_nothing
};

std::ostream & operator<< (std::ostream& os, const Type& type);

std::ostream & operator<< (std::ostream& os, const FileFormat& type);

std::ostream & operator<< (std::ostream& os, const InstanceType& type);

std::ostream & operator<< (std::ostream& os, const PresetType& type);

std::ostream & operator<< (std::ostream& os, const mt_kahypar_partition_type_t& type);

std::ostream & operator<< (std::ostream& os, const ContextType& type);

std::ostream & operator<< (std::ostream& os, const Mode& mode);

std::ostream & operator<< (std::ostream& os, const Objective& objective);

std::ostream & operator<< (std::ostream& os, const GainPolicy& type);

std::ostream & operator<< (std::ostream& os, const LouvainEdgeWeight& type);

std::ostream & operator<< (std::ostream& os, const SimiliarNetCombinerStrategy& strategy);

std::ostream & operator<< (std::ostream& os, const CoarseningAlgorithm& algo);

std::ostream & operator<< (std::ostream& os, const HeavyNodePenaltyPolicy& heavy_hn_policy);

std::ostream & operator<< (std::ostream& os, const AcceptancePolicy& acceptance_policy);

std::ostream & operator<< (std::ostream& os, const RatingFunction& func);

std::ostream & operator<< (std::ostream& os, const InitialPartitioningAlgorithm& algo);

std::ostream & operator<< (std::ostream& os, const LabelPropagationAlgorithm& algo);

std::ostream & operator<< (std::ostream& os, const JetAlgorithm& algo);

std::ostream & operator<< (std::ostream& os, const FMAlgorithm& algo);

std::ostream & operator<< (std::ostream& os, const FlowAlgorithm& algo);

std::ostream & operator<< (std::ostream& os, const RebalancingAlgorithm& algo);

Mode modeFromString(const std::string& mode);

InstanceType instanceTypeFromString(const std::string& type);

PresetType presetTypeFromString(const std::string& type);

Objective objectiveFromString(const std::string& obj);

LouvainEdgeWeight louvainEdgeWeightFromString(const std::string& type);

SimiliarNetCombinerStrategy similiarNetCombinerStrategyFromString(const std::string& type);

CoarseningAlgorithm coarseningAlgorithmFromString(const std::string& type);

HeavyNodePenaltyPolicy heavyNodePenaltyFromString(const std::string& penalty);

AcceptancePolicy acceptanceCriterionFromString(const std::string& crit);

RatingFunction ratingFunctionFromString(const std::string& function);

InitialPartitioningAlgorithm initialPartitioningAlgorithmFromString(const std::string& algo);

LabelPropagationAlgorithm labelPropagationAlgorithmFromString(const std::string& type);

JetAlgorithm jetAlgorithmFromString(const std::string& type);

FMAlgorithm fmAlgorithmFromString(const std::string& type);

FlowAlgorithm flowAlgorithmFromString(const std::string& type);

RebalancingAlgorithm rebalancingAlgorithmFromString(const std::string& type);

}  // namesapce mt_kahypar

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
  graph = 0,
  hypergraph = 1,
  UNDEFINED = 2
};

enum class PresetType : int8_t {
  deterministic,
  default_preset,
  default_flows,
  quality_preset,
  quality_flows,
  UNDEFINED
};

enum class Paradigm : int8_t {
  multilevel,
  nlevel
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
  UNDEFINED
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
  nlevel_coarsener,
  UNDEFINED
};

enum class RatingFunction : uint8_t {
  heavy_edge,
  sameness,
  UNDEFINED
};

enum class HeavyNodePenaltyPolicy : uint8_t {
  no_penalty,
  multiplicative_penalty,
  additive,
  UNDEFINED
};

enum class AcceptancePolicy : uint8_t {
  best,
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
  label_propagation_km1,
  label_propagation_cut,
  deterministic,
  do_nothing
};

enum class FMAlgorithm : uint8_t {
  fm_gain_cache,
  fm_gain_cache_on_demand,
  fm_gain_delta,
  fm_recompute_gain,
  do_nothing
};

enum class FlowAlgorithm : uint8_t {
  flow_cutter,
  mock,
  do_nothing
};

std::ostream & operator<< (std::ostream& os, const Type& type);

std::ostream & operator<< (std::ostream& os, const FileFormat& type);

std::ostream & operator<< (std::ostream& os, const InstanceType& type);

std::ostream & operator<< (std::ostream& os, const PresetType& type);

std::ostream & operator<< (std::ostream& os, const Paradigm& paradigm);

std::ostream & operator<< (std::ostream& os, const ContextType& type);

std::ostream & operator<< (std::ostream& os, const Mode& mode);

std::ostream & operator<< (std::ostream& os, const Objective& objective);

std::ostream & operator<< (std::ostream& os, const LouvainEdgeWeight& type);

std::ostream & operator<< (std::ostream& os, const SimiliarNetCombinerStrategy& strategy);

std::ostream & operator<< (std::ostream& os, const CoarseningAlgorithm& algo);

std::ostream & operator<< (std::ostream& os, const HeavyNodePenaltyPolicy& heavy_hn_policy);

std::ostream & operator<< (std::ostream& os, const AcceptancePolicy& acceptance_policy);

std::ostream & operator<< (std::ostream& os, const RatingFunction& func);

std::ostream & operator<< (std::ostream& os, const InitialPartitioningAlgorithm& algo);

std::ostream & operator<< (std::ostream& os, const LabelPropagationAlgorithm& algo);

std::ostream & operator<< (std::ostream& os, const FMAlgorithm& algo);

std::ostream & operator<< (std::ostream& os, const FlowAlgorithm& algo);

Mode modeFromString(const std::string& mode);

InstanceType instanceTypeFromString(const std::string& type);

PresetType presetTypeFromString(const std::string& type);

LouvainEdgeWeight louvainEdgeWeightFromString(const std::string& type);

SimiliarNetCombinerStrategy similiarNetCombinerStrategyFromString(const std::string& type);

CoarseningAlgorithm coarseningAlgorithmFromString(const std::string& type);

HeavyNodePenaltyPolicy heavyNodePenaltyFromString(const std::string& penalty);

AcceptancePolicy acceptanceCriterionFromString(const std::string& crit);

RatingFunction ratingFunctionFromString(const std::string& function);

InitialPartitioningAlgorithm initialPartitioningAlgorithmFromString(const std::string& algo);

LabelPropagationAlgorithm labelPropagationAlgorithmFromString(const std::string& type);

FMAlgorithm fmAlgorithmFromString(const std::string& type);

FlowAlgorithm flowAlgorithmFromString(const std::string& type);

}  // namesapce mt_kahypar

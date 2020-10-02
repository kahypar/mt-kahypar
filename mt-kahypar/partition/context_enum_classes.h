/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
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

enum class Paradigm : int8_t {
  multilevel,
  nlevel
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

enum class InitialPartitioningMode : uint8_t {
  direct,
  recursive,
  recursive_bisection,
  UNDEFINED
};

enum class LabelPropagationAlgorithm : uint8_t {
  label_propagation_km1,
  label_propagation_cut,
  do_nothing
};

enum class FMAlgorithm : uint8_t {
  fm_gain_cache,
  fm_gain_cache_on_demand,
  fm_gain_delta,
  fm_recompute_gain,
  do_nothing
};

std::ostream & operator<< (std::ostream& os, const Type& type);

std::ostream & operator<< (std::ostream& os, const Paradigm& paradigm);

std::ostream & operator<< (std::ostream& os, const LouvainEdgeWeight& type);

std::ostream & operator<< (std::ostream& os, const SimiliarNetCombinerStrategy& strategy);

std::ostream & operator<< (std::ostream& os, const CoarseningAlgorithm& algo);

std::ostream & operator<< (std::ostream& os, const HeavyNodePenaltyPolicy& heavy_hn_policy);

std::ostream & operator<< (std::ostream& os, const AcceptancePolicy& acceptance_policy);

std::ostream & operator<< (std::ostream& os, const RatingFunction& func);

std::ostream & operator<< (std::ostream& os, const InitialPartitioningAlgorithm& algo);

std::ostream & operator<< (std::ostream& os, const InitialPartitioningMode& mode);

std::ostream & operator<< (std::ostream& os, const LabelPropagationAlgorithm& algo);

std::ostream & operator<< (std::ostream& os, const FMAlgorithm& algo);

LouvainEdgeWeight louvainEdgeWeightFromString(const std::string& type);

SimiliarNetCombinerStrategy similiarNetCombinerStrategyFromString(const std::string& type);

CoarseningAlgorithm coarseningAlgorithmFromString(const std::string& type);

HeavyNodePenaltyPolicy heavyNodePenaltyFromString(const std::string& penalty);

AcceptancePolicy acceptanceCriterionFromString(const std::string& crit);

RatingFunction ratingFunctionFromString(const std::string& function);

InitialPartitioningAlgorithm initialPartitioningAlgorithmFromString(const std::string& algo);

InitialPartitioningMode initialPartitioningModeFromString(const std::string& mode);

LabelPropagationAlgorithm labelPropagationAlgorithmFromString(const std::string& type);

FMAlgorithm fmAlgorithmFromString(const std::string& type);

}  // namesapce mt_kahypar

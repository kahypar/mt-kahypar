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

#include "mt-kahypar/macros.h"

namespace mt_kahypar {

enum class Type : int8_t {
  Unweighted = 0,
  EdgeWeights = 1,
  NodeWeights = 10,
  EdgeAndNodeWeights = 11,
};

enum class InitialHyperedgeDistribution : uint8_t {
  equally,
  random,
  all_on_one,
  UNDEFINED
};

enum class CommunityAssignmentObjective : uint8_t {
  vertex_objective,
  pin_objective,
  UNDEFINED
};

enum class CommunityAssignmentStrategy : uint8_t {
  bin_packing,
  UNDEFINED
};

enum class LouvainEdgeWeight : uint8_t {
  hybrid,
  uniform,
  non_uniform,
  degree,
  UNDEFINED
};

enum class CommunityLoadBalancingStrategy : uint8_t {
  size_constraint,
  none
};

enum class CoarseningAlgorithm : uint8_t {
  community_coarsener,
  UNDEFINED
};

enum class RatingFunction : uint8_t {
  heavy_edge,
  UNDEFINED
};

enum class HeavyNodePenaltyPolicy : uint8_t {
  no_penalty,
  multiplicative_penalty,
  edge_frequency_penalty,
  UNDEFINED
};

enum class AcceptancePolicy : uint8_t {
  best,
  best_prefer_unmatched,
  UNDEFINED
};

enum class InitialPartitioningMode : uint8_t {
  direct,
  recursive,
  UNDEFINED
};

enum class LabelPropagationAlgorithm : uint8_t {
  label_propagation_km1,
  label_propagation_cut,
  do_nothing
};

enum class ExecutionType : uint8_t {
  exponential,
  multilevel,
  constant,
  UNDEFINED
};

std::ostream& operator<< (std::ostream& os, const Type& type) {
  switch (type) {
    case Type::Unweighted: return os << "unweighted";
    case Type::EdgeWeights: return os << "edge_weights";
    case Type::NodeWeights: return os << "node_weights";
    case Type::EdgeAndNodeWeights: return os << "edge_and_node_weights";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(type);
}

std::ostream& operator<< (std::ostream& os, const InitialHyperedgeDistribution& strategy) {
  switch (strategy) {
    case InitialHyperedgeDistribution::equally: return os << "equally";
    case InitialHyperedgeDistribution::random: return os << "random";
    case InitialHyperedgeDistribution::all_on_one: return os << "all_on_one";
    case InitialHyperedgeDistribution::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(strategy);
}


std::ostream& operator<< (std::ostream& os, const CommunityAssignmentObjective& objective) {
  switch (objective) {
    case CommunityAssignmentObjective::vertex_objective: return os << "vertex_objective";
    case CommunityAssignmentObjective::pin_objective: return os << "pin_objective";
    case CommunityAssignmentObjective::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(objective);
}

std::ostream& operator<< (std::ostream& os, const CommunityAssignmentStrategy& strategy) {
  switch (strategy) {
    case CommunityAssignmentStrategy::bin_packing: return os << "bin_packing";
    case CommunityAssignmentStrategy::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(strategy);
}

std::ostream& operator<< (std::ostream& os, const LouvainEdgeWeight& type) {
  switch (type) {
    case LouvainEdgeWeight::hybrid: return os << "hybrid";
    case LouvainEdgeWeight::uniform: return os << "uniform";
    case LouvainEdgeWeight::non_uniform: return os << "non_uniform";
    case LouvainEdgeWeight::degree: return os << "degree";
    case LouvainEdgeWeight::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(type);
}

std::ostream& operator<< (std::ostream& os, const CommunityLoadBalancingStrategy& strategy) {
  switch (strategy) {
    case CommunityLoadBalancingStrategy::size_constraint: return os << "size_constraint";
    case CommunityLoadBalancingStrategy::none: return os << "none";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(strategy);
}

std::ostream& operator<< (std::ostream& os, const CoarseningAlgorithm& algo) {
  switch (algo) {
    case CoarseningAlgorithm::community_coarsener: return os << "community_coarsener";
    case CoarseningAlgorithm::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(algo);
}

std::ostream& operator<< (std::ostream& os, const HeavyNodePenaltyPolicy& heavy_hn_policy) {
  switch (heavy_hn_policy) {
    case HeavyNodePenaltyPolicy::multiplicative_penalty: return os << "multiplicative";
    case HeavyNodePenaltyPolicy::no_penalty: return os << "no_penalty";
    case HeavyNodePenaltyPolicy::edge_frequency_penalty: return os << "edge_frequency_penalty";
    case HeavyNodePenaltyPolicy::UNDEFINED: return os << "UNDEFINED";
  }
  return os << static_cast<uint8_t>(heavy_hn_policy);
}

std::ostream& operator<< (std::ostream& os, const AcceptancePolicy& acceptance_policy) {
  switch (acceptance_policy) {
    case AcceptancePolicy::best: return os << "best";
    case AcceptancePolicy::best_prefer_unmatched: return os << "best_prefer_unmatched";
    case AcceptancePolicy::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(acceptance_policy);
}

std::ostream& operator<< (std::ostream& os, const RatingFunction& func) {
  switch (func) {
    case RatingFunction::heavy_edge: return os << "heavy_edge";
    case RatingFunction::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(func);
}

std::ostream& operator<< (std::ostream& os, const InitialPartitioningMode& mode) {
  switch (mode) {
    case InitialPartitioningMode::direct: return os << "direct";
    case InitialPartitioningMode::recursive: return os << "recursive";
    case InitialPartitioningMode::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(mode);
}

std::ostream& operator<< (std::ostream& os, const LabelPropagationAlgorithm& algo) {
  switch (algo) {
    case LabelPropagationAlgorithm::label_propagation_km1: return os << "label_propagation_km1";
    case LabelPropagationAlgorithm::label_propagation_cut: return os << "label_propagation_cut";
    case LabelPropagationAlgorithm::do_nothing: return os << "do_nothing";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(algo);
}

std::ostream& operator<< (std::ostream& os, const ExecutionType& type) {
  switch (type) {
    case ExecutionType::exponential: return os << "exponential";
    case ExecutionType::multilevel: return os << "multilevel";
    case ExecutionType::constant: return os << "constant";
    case ExecutionType::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(type);
}

static InitialHyperedgeDistribution initialHyperedgeDistributionFromString(const std::string& strategy) {
  if (strategy == "equally") {
    return InitialHyperedgeDistribution::equally;
  } else if (strategy == "random") {
    return InitialHyperedgeDistribution::random;
  } else if (strategy == "all_on_one") {
    return InitialHyperedgeDistribution::all_on_one;
  }
  LOG << "No valid community assignment strategy.";
  exit(0);
  return InitialHyperedgeDistribution::UNDEFINED;
}

static CommunityAssignmentObjective communityAssignmentObjectiveFromString(const std::string& objective) {
  if (objective == "vertex_objective") {
    return CommunityAssignmentObjective::vertex_objective;
  } else if (objective == "pin_objective") {
    return CommunityAssignmentObjective::pin_objective;
  }
  LOG << "No valid community assignment objective.";
  exit(0);
  return CommunityAssignmentObjective::UNDEFINED;
}

static CommunityAssignmentStrategy communityAssignmentStrategyFromString(const std::string& strategy) {
  if (strategy == "bin_packing") {
    return CommunityAssignmentStrategy::bin_packing;
  }
  LOG << "No valid community assignment strategy.";
  exit(0);
  return CommunityAssignmentStrategy::UNDEFINED;
}

static LouvainEdgeWeight louvainEdgeWeightFromString(const std::string& type) {
  if (type == "hybrid") {
    return LouvainEdgeWeight::hybrid;
  } else if (type == "uniform") {
    return LouvainEdgeWeight::uniform;
  } else if (type == "non_uniform") {
    return LouvainEdgeWeight::non_uniform;
  } else if (type == "degree") {
    return LouvainEdgeWeight::degree;
  }
  LOG << "No valid louvain edge weight.";
  exit(0);
  return LouvainEdgeWeight::UNDEFINED;
}

static CommunityLoadBalancingStrategy communityLoadBalancingStrategyFromString(const std::string& strategy) {
  if (strategy == "size_constraint") {
    return CommunityLoadBalancingStrategy::size_constraint;
  } else if (strategy == "none") {
    return CommunityLoadBalancingStrategy::none;
  }
  LOG << "No valid louvain edge weight.";
  exit(0);
  return CommunityLoadBalancingStrategy::none;
}

static CoarseningAlgorithm coarseningAlgorithmFromString(const std::string& type) {
  if (type == "community_coarsener") {
    return CoarseningAlgorithm::community_coarsener;
  }
  LOG << "Illegal option:" << type;
  exit(0);
  return CoarseningAlgorithm::UNDEFINED;
}

static HeavyNodePenaltyPolicy heavyNodePenaltyFromString(const std::string& penalty) {
  if (penalty == "multiplicative") {
    return HeavyNodePenaltyPolicy::multiplicative_penalty;
  } else if (penalty == "no_penalty") {
    return HeavyNodePenaltyPolicy::no_penalty;
  } else if (penalty == "edge_frequency_penalty") {
    return HeavyNodePenaltyPolicy::edge_frequency_penalty;
    // omit default case to trigger compiler warning for missing cases
  }
  LOG << "No valid edge penalty policy for rating.";
  exit(0);
  return HeavyNodePenaltyPolicy::UNDEFINED;
}

static AcceptancePolicy acceptanceCriterionFromString(const std::string& crit) {
  if (crit == "best") {
    return AcceptancePolicy::best;
  } else if (crit == "best_prefer_unmatched") {
    return AcceptancePolicy::best_prefer_unmatched;
  }
  LOG << "No valid acceptance criterion for rating.";
  exit(0);
}

static RatingFunction ratingFunctionFromString(const std::string& function) {
  if (function == "heavy_edge") {
    return RatingFunction::heavy_edge;
  }
  LOG << "No valid rating function for rating.";
  exit(0);
  return RatingFunction::UNDEFINED;
}

static InitialPartitioningMode initialPartitioningModeFromString(const std::string& mode) {
  if (mode == "direct") {
    return InitialPartitioningMode::direct;
  } else if (mode == "recursive") {
    return InitialPartitioningMode::recursive;
  }
  LOG << "Illegal option:" << mode;
  exit(0);
  return InitialPartitioningMode::UNDEFINED;
}

static LabelPropagationAlgorithm labelPropagationAlgorithmFromString(const std::string& type) {
  if (type == "label_propagation_km1") {
    return LabelPropagationAlgorithm::label_propagation_km1;
  } else if (type == "label_propagation_cut") {
    return LabelPropagationAlgorithm::label_propagation_cut;
  } else if (type == "do_nothing") {
    return LabelPropagationAlgorithm::do_nothing;
  }
  LOG << "Illegal option:" << type;
  exit(0);
  return LabelPropagationAlgorithm::do_nothing;
}

static ExecutionType executionTypeFromString(const std::string& type) {
  if (type == "exponential") {
    return ExecutionType::exponential;
  } else if (type == "multilevel") {
    return ExecutionType::multilevel;
  } else if (type == "constant") {
    return ExecutionType::constant;
  }
  LOG << "Illegal option:" << type;
  exit(0);
  return ExecutionType::UNDEFINED;
}

} // namesapce mt_kahypar
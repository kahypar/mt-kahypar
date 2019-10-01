/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "kahypar/macros.h"

namespace mt_kahypar {

enum class Type : int8_t {
  Unweighted = 0,
  EdgeWeights = 1,
  NodeWeights = 10,
  EdgeAndNodeWeights = 11,
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

static CommunityAssignmentStrategy communityAssignmentStrategyFromString(const std::string& objective) {
  if (objective == "bin_packing") {
    return CommunityAssignmentStrategy::bin_packing;
  }
  LOG << "No valid community assignment strategy.";
  exit(0);
  return CommunityAssignmentStrategy::UNDEFINED;
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
  return HeavyNodePenaltyPolicy::multiplicative_penalty;
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

} // namesapce mt_kahypar
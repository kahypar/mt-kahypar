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

enum class CoarseningVertexOrder : uint8_t {
  non_randomized,
  random_shuffle,
  pseudo_random_shuffle,
  increasing_degree_order,
  decreasing_degree_order,
  UNDEFINED
};

enum class ContractionOrder : uint8_t {
  rating_order,
  degree_order,
  heavy_node_first,
  heavy_node_second,
  UNDEFINED
};

enum class UncontractionOrder : uint8_t {
  max_subtree_size,
  min_subtree_size,
  max_vertex_degree,
  min_vertex_degree,
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
  fm_multitry,
  fm_boundary,
  do_nothing
};


std::ostream & operator<< (std::ostream& os, const CoarseningVertexOrder& order) {
  switch (order) {
    case CoarseningVertexOrder::non_randomized: return os << "non_randomized";
    case CoarseningVertexOrder::random_shuffle: return os << "random_shuffle";
    case CoarseningVertexOrder::pseudo_random_shuffle: return os << "pseudo_random_shuffle";
    case CoarseningVertexOrder::increasing_degree_order: return os << "increasing_degree_order";
    case CoarseningVertexOrder::decreasing_degree_order: return os << "decreasing_degree_order";
    case CoarseningVertexOrder::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(order);
}

std::ostream & operator<< (std::ostream& os, const ContractionOrder& order) {
  switch (order) {
    case ContractionOrder::rating_order: return os << "rating_order";
    case ContractionOrder::degree_order: return os << "degree_order";
    case ContractionOrder::heavy_node_first: return os << "heavy_node_first";
    case ContractionOrder::heavy_node_second: return os << "heavy_node_second";
    case ContractionOrder::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(order);
}

std::ostream & operator<< (std::ostream& os, const UncontractionOrder& order) {
  switch (order) {
    case UncontractionOrder::max_subtree_size: return os << "max_subtree_size";
    case UncontractionOrder::min_subtree_size: return os << "min_subtree_size";
    case UncontractionOrder::max_vertex_degree: return os << "max_vertex_degree";
    case UncontractionOrder::min_vertex_degree: return os << "min_vertex_degree";
    case UncontractionOrder::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(order);
}

static CoarseningVertexOrder coarseningVertexOrderFromString(const std::string& order) {
  if (order == "non_randomized") {
    return CoarseningVertexOrder::non_randomized;
  } else if (order == "random_shuffle") {
    return CoarseningVertexOrder::random_shuffle;
  } else if (order == "pseudo_random_shuffle") {
    return CoarseningVertexOrder::pseudo_random_shuffle;
  } else if (order == "increasing_degree_order") {
    return CoarseningVertexOrder::increasing_degree_order;
  } else if (order == "decreasing_degree_order") {
    return CoarseningVertexOrder::decreasing_degree_order;
  }
  ERROR("No valid coarsening vertex order.");
}

static ContractionOrder contractionOrderFromString(const std::string& order) {
  if (order == "rating_order") {
    return ContractionOrder::rating_order;
  } else if (order == "degree_order") {
    return ContractionOrder::degree_order;
  } else if (order == "heavy_node_first") {
    return ContractionOrder::heavy_node_first;
  } else if (order == "heavy_node_second") {
    return ContractionOrder::heavy_node_second;
  }
  ERROR("No valid contraction order.");
}

static UncontractionOrder uncontractionOrderFromString(const std::string& order) {
  if (order == "max_subtree_size") {
    return UncontractionOrder::max_subtree_size;
  } else if (order == "min_subtree_size") {
    return UncontractionOrder::min_subtree_size;
  } else if (order == "max_vertex_degree") {
    return UncontractionOrder::max_vertex_degree;
  } else if (order == "min_vertex_degree") {
    return UncontractionOrder::min_vertex_degree;
  }
  ERROR("No valid uncontraction order.");
}

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

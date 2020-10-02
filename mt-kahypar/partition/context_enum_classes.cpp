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

#include "context_enum_classes.h"

#include "mt-kahypar/macros.h"

namespace mt_kahypar {

  std::ostream & operator<< (std::ostream& os, const Type& type) {
    switch (type) {
      case Type::Unweighted: return os << "unweighted";
      case Type::EdgeWeights: return os << "edge_weights";
      case Type::NodeWeights: return os << "node_weights";
      case Type::EdgeAndNodeWeights: return os << "edge_and_node_weights";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(type);
  }

  std::ostream & operator<< (std::ostream& os, const Paradigm& paradigm) {
    switch (paradigm) {
      case Paradigm::multilevel: return os << "multilevel";
      case Paradigm::nlevel: return os << "nlevel";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(paradigm);
  }

  std::ostream & operator<< (std::ostream& os, const LouvainEdgeWeight& type) {
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

  std::ostream & operator<< (std::ostream& os, const SimiliarNetCombinerStrategy& strategy) {
    switch (strategy) {
      case SimiliarNetCombinerStrategy::union_nets: return os << "union";
      case SimiliarNetCombinerStrategy::max_size: return os << "max_size";
      case SimiliarNetCombinerStrategy::importance: return os << "importance";
      case SimiliarNetCombinerStrategy::UNDEFINED: return os << "UNDEFINED";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(strategy);
  }

  std::ostream & operator<< (std::ostream& os, const CoarseningAlgorithm& algo) {
    switch (algo) {
      case CoarseningAlgorithm::multilevel_coarsener: return os << "multilevel_coarsener";
      case CoarseningAlgorithm::nlevel_coarsener: return os << "nlevel_coarsener";
      case CoarseningAlgorithm::UNDEFINED: return os << "UNDEFINED";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(algo);
  }

  std::ostream & operator<< (std::ostream& os, const HeavyNodePenaltyPolicy& heavy_hn_policy) {
    switch (heavy_hn_policy) {
      case HeavyNodePenaltyPolicy::multiplicative_penalty: return os << "multiplicative";
      case HeavyNodePenaltyPolicy::no_penalty: return os << "no_penalty";
      case HeavyNodePenaltyPolicy::additive: return os << "additive";
      case HeavyNodePenaltyPolicy::UNDEFINED: return os << "UNDEFINED";
    }
    return os << static_cast<uint8_t>(heavy_hn_policy);
  }

  std::ostream & operator<< (std::ostream& os, const AcceptancePolicy& acceptance_policy) {
    switch (acceptance_policy) {
      case AcceptancePolicy::best: return os << "best";
      case AcceptancePolicy::best_prefer_unmatched: return os << "best_prefer_unmatched";
      case AcceptancePolicy::UNDEFINED: return os << "UNDEFINED";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(acceptance_policy);
  }

  std::ostream & operator<< (std::ostream& os, const RatingFunction& func) {
    switch (func) {
      case RatingFunction::heavy_edge: return os << "heavy_edge";
      case RatingFunction::sameness: return os << "sameness";
      case RatingFunction::UNDEFINED: return os << "UNDEFINED";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(func);
  }

  std::ostream & operator<< (std::ostream& os, const InitialPartitioningAlgorithm& algo) {
    switch (algo) {
      case InitialPartitioningAlgorithm::random: return os << "random";
      case InitialPartitioningAlgorithm::bfs: return os << "bfs";
      case InitialPartitioningAlgorithm::greedy_round_robin_fm: return os << "greedy_round_robin_fm";
      case InitialPartitioningAlgorithm::greedy_global_fm: return os << "greedy_global_fm";
      case InitialPartitioningAlgorithm::greedy_sequential_fm: return os << "greedy_sequential_fm";
      case InitialPartitioningAlgorithm::greedy_round_robin_max_net: return os << "greedy_round_robin_max_net";
      case InitialPartitioningAlgorithm::greedy_global_max_net: return os << "greedy_global_max_net";
      case InitialPartitioningAlgorithm::greedy_sequential_max_net: return os << "greedy_sequential_max_net";
      case InitialPartitioningAlgorithm::label_propagation: return os << "label_propagation";
      case InitialPartitioningAlgorithm::UNDEFINED: return os << "UNDEFINED";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(algo);
  }

  std::ostream & operator<< (std::ostream& os, const InitialPartitioningMode& mode) {
    switch (mode) {
      case InitialPartitioningMode::direct: return os << "direct";
      case InitialPartitioningMode::recursive: return os << "recursive";
      case InitialPartitioningMode::recursive_bisection: return os << "recursive_bisection";
      case InitialPartitioningMode::UNDEFINED: return os << "UNDEFINED";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(mode);
  }

  std::ostream & operator<< (std::ostream& os, const LabelPropagationAlgorithm& algo) {
    switch (algo) {
      case LabelPropagationAlgorithm::label_propagation_km1: return os << "label_propagation_km1";
      case LabelPropagationAlgorithm::label_propagation_cut: return os << "label_propagation_cut";
      case LabelPropagationAlgorithm::do_nothing: return os << "lp_do_nothing";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(algo);
  }

  std::ostream & operator<< (std::ostream& os, const FMAlgorithm& algo) {
    switch (algo) {
      case FMAlgorithm::fm_gain_cache: return os << "fm_gain_cache";
      case FMAlgorithm::fm_gain_cache_on_demand : return os << "fm_gain_cache_on_demand";
      case FMAlgorithm::fm_gain_delta: return os << "fm_gain_delta";
      case FMAlgorithm::fm_recompute_gain: return os << "fm_recompute_gain";
      case FMAlgorithm::do_nothing: return os << "fm_do_nothing";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(algo);
  }

  LouvainEdgeWeight louvainEdgeWeightFromString(const std::string& type) {
    if (type == "hybrid") {
      return LouvainEdgeWeight::hybrid;
    } else if (type == "uniform") {
      return LouvainEdgeWeight::uniform;
    } else if (type == "non_uniform") {
      return LouvainEdgeWeight::non_uniform;
    } else if (type == "degree") {
      return LouvainEdgeWeight::degree;
    }
    ERROR("No valid louvain edge weight.");
    return LouvainEdgeWeight::UNDEFINED;
  }

  SimiliarNetCombinerStrategy similiarNetCombinerStrategyFromString(const std::string& type) {
    if (type == "union") {
      return SimiliarNetCombinerStrategy::union_nets;
    } else if (type == "max_size") {
      return SimiliarNetCombinerStrategy::max_size;
    } else if (type == "importance") {
      return SimiliarNetCombinerStrategy::importance;
    }
    ERROR("No valid similiar net unifier strategy.");
    return SimiliarNetCombinerStrategy::UNDEFINED;
  }

  CoarseningAlgorithm coarseningAlgorithmFromString(const std::string& type) {
    if (type == "multilevel_coarsener") {
      return CoarseningAlgorithm::multilevel_coarsener;
    } else if (type == "nlevel_coarsener") {
      return CoarseningAlgorithm::nlevel_coarsener;
    }
    ERROR("Illegal option: " + type);
    return CoarseningAlgorithm::UNDEFINED;
  }

  HeavyNodePenaltyPolicy heavyNodePenaltyFromString(const std::string& penalty) {
    if (penalty == "multiplicative") {
      return HeavyNodePenaltyPolicy::multiplicative_penalty;
    } else if (penalty == "no_penalty") {
      return HeavyNodePenaltyPolicy::no_penalty;
    } else if (penalty == "additive") {
      return HeavyNodePenaltyPolicy::additive;
      // omit default case to trigger compiler warning for missing cases
    }
    ERROR("No valid edge penalty policy for rating.");
    return HeavyNodePenaltyPolicy::UNDEFINED;
  }

  AcceptancePolicy acceptanceCriterionFromString(const std::string& crit) {
    if (crit == "best") {
      return AcceptancePolicy::best;
    } else if (crit == "best_prefer_unmatched") {
      return AcceptancePolicy::best_prefer_unmatched;
    }
    ERROR("No valid acceptance criterion for rating.");
  }

  RatingFunction ratingFunctionFromString(const std::string& function) {
    if (function == "heavy_edge") {
      return RatingFunction::heavy_edge;
    } else  if (function == "sameness") {
      return RatingFunction::sameness;
    }
    ERROR("No valid rating function for rating.");
    return RatingFunction::UNDEFINED;
  }

  InitialPartitioningAlgorithm initialPartitioningAlgorithmFromString(const std::string& algo) {
    if (algo == "random") {
      return InitialPartitioningAlgorithm::random;
    } else if (algo == "bfs") {
      return InitialPartitioningAlgorithm::bfs;
    } else if (algo == "greedy_round_robin_fm") {
      return InitialPartitioningAlgorithm::greedy_round_robin_fm;
    } else if (algo == "greedy_global_fm") {
      return InitialPartitioningAlgorithm::greedy_global_fm;
    } else if (algo == "greedy_sequential_fm") {
      return InitialPartitioningAlgorithm::greedy_sequential_fm;
    } else if (algo == "greedy_round_robin_max_net") {
      return InitialPartitioningAlgorithm::greedy_round_robin_max_net;
    } else if (algo == "greedy_global_max_net") {
      return InitialPartitioningAlgorithm::greedy_global_max_net;
    } else if (algo == "greedy_sequential_max_net") {
      return InitialPartitioningAlgorithm::greedy_sequential_max_net;
    } else if (algo == "label_propagation") {
      return InitialPartitioningAlgorithm::label_propagation;
    }
    ERROR("Illegal option: " + algo);
    return InitialPartitioningAlgorithm::UNDEFINED;
  }

  InitialPartitioningMode initialPartitioningModeFromString(const std::string& mode) {
    if (mode == "direct") {
      return InitialPartitioningMode::direct;
    } else if (mode == "recursive") {
      return InitialPartitioningMode::recursive;
    } else if (mode == "recursive_bisection") {
      return InitialPartitioningMode::recursive_bisection;
    }
    ERROR("Illegal option: " + mode);
    return InitialPartitioningMode::UNDEFINED;
  }

  LabelPropagationAlgorithm labelPropagationAlgorithmFromString(const std::string& type) {
    if (type == "label_propagation_km1") {
      return LabelPropagationAlgorithm::label_propagation_km1;
    } else if (type == "label_propagation_cut") {
      return LabelPropagationAlgorithm::label_propagation_cut;
    } else if (type == "do_nothing") {
      return LabelPropagationAlgorithm::do_nothing;
    }
    ERROR("Illegal option: " + type);
    return LabelPropagationAlgorithm::do_nothing;
  }

  FMAlgorithm fmAlgorithmFromString(const std::string& type) {
    if (type == "fm_gain_cache") {
      return FMAlgorithm::fm_gain_cache;
    } else if (type == "fm_gain_cache_on_demand") {
      return FMAlgorithm::fm_gain_cache_on_demand;
    } else if (type == "fm_gain_delta") {
      return FMAlgorithm::fm_gain_delta;
    } else if (type == "fm_recompute_gain") {
      return FMAlgorithm::fm_recompute_gain;
    } else if (type == "do_nothing") {
      return FMAlgorithm::do_nothing;
    }
    ERROR("Illegal option: " + type);
    return FMAlgorithm::do_nothing;
  }

}
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

  std::ostream & operator<< (std::ostream& os, const FileFormat& format) {
    switch (format) {
      case FileFormat::hMetis: return os << "hMetis";
      case FileFormat::Metis: return os << "Metis";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(format);
  }

  std::ostream & operator<< (std::ostream& os, const InstanceType& type) {
    switch (type) {
      case InstanceType::graph: return os << "graph";
      case InstanceType::hypergraph: return os << "hypergraph";
      case InstanceType::UNDEFINED: return os << "UNDEFINED";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(type);
  }

  std::ostream & operator<< (std::ostream& os, const PresetType& type) {
    switch (type) {
      case PresetType::deterministic: return os << "deterministic";
      case PresetType::default_preset: return os << "default";
      case PresetType::default_flows: return os << "default_flows";
      case PresetType::quality_preset: return os << "quality";
      case PresetType::quality_flows: return os << "quality_flows";
      case PresetType::UNDEFINED: return os << "UNDEFINED";
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

  std::ostream& operator<< (std::ostream& os, const ContextType& type) {
    if (type == ContextType::main) {
      return os << "main";
    } else {
      return os << "ip";
    }
    return os << static_cast<uint8_t>(type);
  }

  std::ostream & operator<< (std::ostream& os, const Mode& mode) {
    switch (mode) {
      case Mode::recursive_bipartitioning: return os << "recursive_bipartitioning";
      case Mode::direct: return os << "direct_kway";
      case Mode::deep_multilevel: return os << "deep_multilevel";
      case Mode::UNDEFINED: return os << "UNDEFINED";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(mode);
  }

  std::ostream& operator<< (std::ostream& os, const Objective& objective) {
    switch (objective) {
      case Objective::cut: return os << "cut";
      case Objective::km1: return os << "km1";
      case Objective::UNDEFINED: return os << "UNDEFINED";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(objective);
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
      case CoarseningAlgorithm::deterministic_multilevel_coarsener: return os << "deterministic_multilevel_coarsener";
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

  std::ostream & operator<< (std::ostream& os, const LabelPropagationAlgorithm& algo) {
    switch (algo) {
      case LabelPropagationAlgorithm::label_propagation_km1: return os << "label_propagation_km1";
      case LabelPropagationAlgorithm::label_propagation_cut: return os << "label_propagation_cut";
      case LabelPropagationAlgorithm::deterministic: return os << "deterministic";
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

  std::ostream & operator<< (std::ostream& os, const FlowAlgorithm& algo) {
    switch (algo) {
      case FlowAlgorithm::flow_cutter: return os << "flow_cutter";
      case FlowAlgorithm::mock: return os << "mock";
      case FlowAlgorithm::do_nothing: return os << "do_nothing";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(algo);
  }

  Mode modeFromString(const std::string& mode) {
    if (mode == "rb") {
      return Mode::recursive_bipartitioning;
    } else if (mode == "direct") {
      return Mode::direct;
    } else if (mode == "deep") {
      return Mode::deep_multilevel;
    }
    ERROR("Illegal option: " + mode);
    return Mode::UNDEFINED;
  }

  InstanceType instanceTypeFromString(const std::string& type) {
    if (type == "graph") {
      return InstanceType::graph;
    } else if (type == "hypergraph") {
      return InstanceType::hypergraph;
    }
    ERROR("Illegal option: " + type);
    return InstanceType::UNDEFINED;
  }

  PresetType presetTypeFromString(const std::string& type) {
    if (type == "deterministic") {
      return PresetType::deterministic;
    } else if (type == "default") {
      return PresetType::default_preset;
    } else if (type == "default_flows") {
      return PresetType::default_flows;
    } else if (type == "quality") {
      return PresetType::quality_preset;
    } else if (type == "quality_flows") {
      return PresetType::quality_flows;
    }
    ERROR("Illegal option: " + type);
    return PresetType::UNDEFINED;
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
    } else if (type == "deterministic_multilevel_coarsener") {
      return CoarseningAlgorithm::deterministic_multilevel_coarsener;
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

  LabelPropagationAlgorithm labelPropagationAlgorithmFromString(const std::string& type) {
    if (type == "label_propagation_km1") {
      return LabelPropagationAlgorithm::label_propagation_km1;
    } else if (type == "label_propagation_cut") {
      return LabelPropagationAlgorithm::label_propagation_cut;
    } else if (type == "deterministic") {
      return LabelPropagationAlgorithm::deterministic;
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

  FlowAlgorithm flowAlgorithmFromString(const std::string& type) {
    if (type == "flow_cutter") {
      return FlowAlgorithm::flow_cutter;
    } else if (type == "do_nothing") {
      return FlowAlgorithm::do_nothing;
    }
    ERROR("Illegal option: " + type);
    return FlowAlgorithm::do_nothing;
  }
}
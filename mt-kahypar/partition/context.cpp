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

#include "context.h"

#include <algorithm>

#include "mt-kahypar/utils/exception.h"
#include "mt-kahypar/partition/conversion.h"

namespace mt_kahypar {

  std::ostream & operator<< (std::ostream& str, const PartitioningParameters& params) {
    str << "Partitioning Parameters:" << std::endl;
    str << "  Hypergraph:                         " << params.graph_filename << std::endl;
    if ( params.fixed_vertex_filename != "" ) {
      str << "  Fixed Vertex File:                  " << params.fixed_vertex_filename << std::endl;
    }
    if ( params.write_partition_file ) {
      str << "  Partition File:                     " << params.graph_partition_filename << std::endl;
    }
    str << "  Mode:                               " << params.mode << std::endl;
    str << "  Objective:                          " << params.objective << std::endl;
    str << "  Gain Policy:                        " << params.gain_policy << std::endl;
    str << "  Input File Format:                  " << params.file_format << std::endl;
    if ( params.instance_type != InstanceType::UNDEFINED ) {
      str << "  Instance Type:                      " << params.instance_type << std::endl;
    }
    if ( params.preset_type != PresetType::UNDEFINED ) {
      str << "  Preset Type:                        " << params.preset_type << std::endl;
    }
    str << "  Partition Type:                     " << params.partition_type << std::endl;
    str << "  k:                                  " << params.k << std::endl;
    str << "  epsilon:                            " << params.epsilon << std::endl;
    str << "  seed:                               " << params.seed << std::endl;
    str << "  Number of V-Cycles:                 " << params.num_vcycles << std::endl;
    str << "  Ignore HE Size Threshold:           " << params.ignore_hyperedge_size_threshold << std::endl;
    str << "  Large HE Size Threshold:            " << params.large_hyperedge_size_threshold << std::endl;
    if ( params.use_individual_part_weights ) {
      str << "  Individual Part Weights:            ";
      for ( const HypernodeWeight& w : params.max_part_weights ) {
        str << w << " ";
      }
      str << std::endl;
    }
    if ( params.mode == Mode::deep_multilevel ) {
      str << "  Perform Parallel Recursion:         " << std::boolalpha
          << params.perform_parallel_recursion_in_deep_multilevel << std::endl;
    }
    return str;
  }

  std::ostream & operator<< (std::ostream& str, const CommunityDetectionParameters& params) {
    str << "  Community Detection Parameters:" << std::endl;
    str << "    Edge Weight Function:                " << params.edge_weight_function << std::endl;
    str << "    Maximum Louvain-Pass Iterations:     " << params.max_pass_iterations << std::endl;
    str << "    Minimum Vertex Move Fraction:        " << params.min_vertex_move_fraction << std::endl;
    str << "    Vertex Degree Sampling Threshold:    " << params.vertex_degree_sampling_threshold << std::endl;
    str << "    Number of subrounds (deterministic): " << params.num_sub_rounds_deterministic << std::endl;
    return str;
  }

  std::ostream & operator<< (std::ostream& str, const PreprocessingParameters& params) {
    str << "Preprocessing Parameters:" << std::endl;
    str << "  Use Community Detection:            " << std::boolalpha << params.use_community_detection << std::endl;
    str << "  Disable C. D. for Mesh Graphs:      " << std::boolalpha << params.disable_community_detection_for_mesh_graphs << std::endl;
    if (params.use_community_detection) {
      str << std::endl << params.community_detection;
    }
    return str;
  }

  std::ostream & operator<< (std::ostream& str, const RatingParameters& params) {
    str << "  Rating Parameters:" << std::endl;
    str << "    Rating Function:                  " << params.rating_function << std::endl;
    str << "    Heavy Node Penalty:               " << params.heavy_node_penalty_policy << std::endl;
    str << "    Acceptance Policy:                " << params.acceptance_policy << std::endl;
    return str;
  }

  std::ostream & operator<< (std::ostream& str, const CoarseningParameters& params) {
    str << "Coarsening Parameters:" << std::endl;
    str << "  Algorithm:                          " << params.algorithm << std::endl;
    str << "  Use Adaptive Edge Size:             " << std::boolalpha << params.use_adaptive_edge_size << std::endl;
    str << "  Max Allowed Weight Multiplier:      " << params.max_allowed_weight_multiplier << std::endl;
    str << "  Maximum Allowed Hypernode Weight:   " << params.max_allowed_node_weight << std::endl;
    str << "  Contraction Limit Multiplier:       " << params.contraction_limit_multiplier << std::endl;
    str << "  Deep ML Contraction Limit Multi.:   " << params.deep_ml_contraction_limit_multiplier << std::endl;
    str << "  Contraction Limit:                  " << params.contraction_limit << std::endl;
    str << "  Minimum Shrink Factor:              " << params.minimum_shrink_factor << std::endl;
    str << "  Maximum Shrink Factor:              " << params.maximum_shrink_factor << std::endl;
    str << "  Vertex Degree Sampling Threshold:   " << params.vertex_degree_sampling_threshold << std::endl;
    if ( params.algorithm == CoarseningAlgorithm::deterministic_multilevel_coarsener ) {
      str << "  Number of Subrounds:                " << params.num_sub_rounds_deterministic << std::endl;
      str << "  Resolve Node Swaps:                 " << std::boolalpha << params.det_resolve_swaps << std::endl;
    }
    if ( params.algorithm == CoarseningAlgorithm::multilevel_coarsener || params.algorithm == CoarseningAlgorithm::nlevel_coarsener ) {
      str << std::endl << params.rating;
    }
    return str;
  }

  std::ostream & operator<< (std::ostream& str, const LabelPropagationParameters& params) {
    str << "  Label Propagation Parameters:" << std::endl;
    str << "    Algorithm:                        " << params.algorithm << std::endl;
    if ( params.algorithm != LabelPropagationAlgorithm::do_nothing ) {
      str << "    Maximum Iterations:               " << params.maximum_iterations << std::endl;
      str << "    Unconstrained:                    " << std::boolalpha << params.unconstrained << std::endl;
      str << "    Rebalancing:                      " << std::boolalpha << params.rebalancing << std::endl;
      str << "    HE Size Activation Threshold:     " << std::boolalpha << params.hyperedge_size_activation_threshold << std::endl;
      str << "    Relative Improvement Threshold:   " << params.relative_improvement_threshold << std::endl;
    }
    return str;
  }

  std::ostream & operator<< (std::ostream& str, const JetParameters& params) {
    str << "  Jet Parameters:" << std::endl;
    str << "    Algorithm:                        " << params.algorithm << std::endl;
    if ( params.algorithm != JetAlgorithm::do_nothing ) {
      str << "    Iterations without Improvement:   " << params.num_iterations << std::endl;
      str << "    Relative Improvement Threshold:   " << params.relative_improvement_threshold << std::endl;
      str << "    Dynamic Rounds:                   " << params.dynamic_rounds << std::endl;
      str << "    Initial Negative Gain Factor:     " << params.initial_negative_gain_factor << std::endl;
      str << "    Final Negative Gain Factor:       " << params.final_negative_gain_factor << std::endl;
    }
    return str;
  }

  std::ostream& operator<<(std::ostream& out, const FMParameters& params) {
    out << "  FM Parameters: \n";
    out << "    Algorithm:                        " << params.algorithm << std::endl;
    if ( params.algorithm != FMAlgorithm::do_nothing ) {
      out << "    Multitry Rounds:                  " << params.multitry_rounds << std::endl;
      out << "    Parallel Global Rollbacks:        " << std::boolalpha << params.rollback_parallel << std::endl;
      out << "    Rollback Bal. Violation Factor:   " << params.rollback_balance_violation_factor << std::endl;
      out << "    Num Seed Nodes:                   " << params.num_seed_nodes << std::endl;
      out << "    Enable Random Shuffle:            " << std::boolalpha << params.shuffle << std::endl;
      out << "    Obey Minimal Parallelism:         " << std::boolalpha << params.obey_minimal_parallelism << std::endl;
      out << "    Minimum Improvement Factor:       " << params.min_improvement << std::endl;
      out << "    Release Nodes:                    " << std::boolalpha << params.release_nodes << std::endl;
      out << "    Time Limit Factor:                " << params.time_limit_factor << std::endl;
    }
    if ( params.algorithm == FMAlgorithm::unconstrained_fm ) {
      out << "    Unconstrained Rounds:             " << params.unconstrained_rounds << std::endl;
      out << "    Threshold Border Node Inclusion:  " << params.treshold_border_node_inclusion << std::endl;
      out << "    Minimum Imbalance Penalty Factor: " << params.imbalance_penalty_min << std::endl;
      out << "    Maximum Imbalance Penalty Factor: " << params.imbalance_penalty_max << std::endl;
      out << "    Start Upper Bound for Unc.:       " << params.unconstrained_upper_bound << std::endl;
      out << "    Final Upper Bound for Unc.:       " << params.unconstrained_upper_bound_min << std::endl;
      out << "    Unc. Minimum Improvement Factor:  " << params.unconstrained_min_improvement << std::endl;
      out << "    Activate Unc. Dynamically:        " << std::boolalpha << params.activate_unconstrained_dynamically << std::endl;
      if ( params.activate_unconstrained_dynamically ) {
        out << "    Penalty for Activation Test:      " << params.penalty_for_activation_test << std::endl;
      }
    }
    out << std::flush;
    return out;
  }

  std::ostream& operator<<(std::ostream& out, const NLevelGlobalRefinementParameters& params) {
    if ( params.use_global_refinement ) {
      out << "  Global Refinement Parameters:" << std::endl;
      out << "    Refine Until No Improvement:      " << std::boolalpha << params.refine_until_no_improvement << std::endl;
      out << "    FM Algorithm:                     " << params.fm_algorithm << std::endl;
      out << "    FM Num Seed Nodes:                " << params.fm_num_seed_nodes << std::endl;
      out << "    FM Obey Minimal Parallelism:      " << std::boolalpha << params.fm_obey_minimal_parallelism << std::endl;
      out << "    LP Algorithm:                     " << params.lp_algorithm << std::endl;
      out << "    LP Unconstrained:                 " << std::boolalpha << params.lp_unconstrained << std::endl;
    }
    return out;
  }

  std::ostream& operator<<(std::ostream& out, const FlowParameters& params) {
    out << "  Flow Parameters: \n";
    out << "    Algorithm:                        " << params.algorithm << std::endl;
    if ( params.algorithm != FlowAlgorithm::do_nothing ) {
      out << "    Flow Scaling:                     " << params.alpha << std::endl;
      out << "    Maximum Number of Pins:           " << params.max_num_pins << std::endl;
      out << "    Find Most Balanced Cut:           " << std::boolalpha << params.find_most_balanced_cut << std::endl;
      out << "    Determine Distance From Cut:      " << std::boolalpha << params.determine_distance_from_cut << std::endl;
      out << "    Number of Parallel Searches:      " << params.num_parallel_searches << std::endl;
      out << "    Maximum BFS Distance:             " << params.max_bfs_distance << std::endl;
      out << "    Min Rel. Improvement Per Round:   " << params.min_relative_improvement_per_round << std::endl;
      out << "    Time Limit Factor:                " << params.time_limit_factor << std::endl;
      out << "    Skip Small Cuts:                  " << std::boolalpha << params.skip_small_cuts << std::endl;
      out << "    Skip Unpromising Blocks:          " << std::boolalpha << params.skip_unpromising_blocks << std::endl;
      out << "    Pierce in Bulk:                   " << std::boolalpha << params.pierce_in_bulk << std::endl;
      out << "    Steiner Tree Policy:              " << params.steiner_tree_policy << std::endl;
      out << std::flush;
    }
    return out;
  }

  std::ostream& operator<<(std::ostream& out, const DeterministicRefinementParameters& params) {
    out << "    Number of sub-rounds for Sync LP:  " << params.num_sub_rounds_sync_lp << std::endl;
    out << "    Use active node set:               " << std::boolalpha << params.use_active_node_set << std::endl;
    return out;
  }

  std::ostream& operator<<(std::ostream& out, const RebalancingParameters& params) {
    out << "  Rebalancing Parameters:" << std::endl;
    out << "    Algorithm:                        " << params.algorithm << std::endl;
    if ( params.algorithm == RebalancingAlgorithm::deterministic ) {
      out << "    Heavy vertex exclusion factor:    " << params.det_heavy_vertex_exclusion_factor << std::endl;
      out << "    Relative deadone size:            " << params.det_relative_deadzone_size << std::endl;
      out << "    Max rounds:                       " << params.det_max_rounds << std::endl;
    }
    return out;
  }

  std::ostream & operator<< (std::ostream& str, const RefinementParameters& params) {
    str << "Refinement Parameters:" << std::endl;
    str << "  Refine Until No Improvement:        " << std::boolalpha << params.refine_until_no_improvement << std::endl;
    str << "  Relative Improvement Threshold:     " << params.relative_improvement_threshold << std::endl;
    str << "  Maximum Batch Size:                 " << params.max_batch_size << std::endl;
    str << "  Min Border Vertices Per Thread:     " << params.min_border_vertices_per_thread << std::endl;
    str << "\n" << params.label_propagation;
    str << "\n" << params.jet;
    str << "\n" << params.fm;
    if ( params.global.use_global_refinement ) {
      str << "\n" << params.global;
    }
    str << "\n" << params.flows;
    str << "\n" << params.rebalancing;
    return str;
  }

  std::ostream & operator<< (std::ostream& str, const InitialPartitioningParameters& params) {
    str << "Initial Partitioning Parameters:" << std::endl;
    str << "  Initial Partitioning Mode:          " << params.mode << std::endl;
    str << "  Number of Runs:                     " << params.runs << std::endl;
    str << "  Use Adaptive IP Runs:               " << std::boolalpha << params.use_adaptive_ip_runs << std::endl;
    if ( params.use_adaptive_ip_runs ) {
      str << "  Min Adaptive IP Runs:               " << params.min_adaptive_ip_runs << std::endl;
    }
    str << "  Perform Refinement On Best:         " << std::boolalpha << params.perform_refinement_on_best_partitions << std::endl;
    str << "  Fm Refinement Rounds:               " << params.fm_refinment_rounds << std::endl;
    str << "  Remove Degree-Zero HNs Before IP:   " << std::boolalpha << params.remove_degree_zero_hns_before_ip << std::endl;
    str << "  Maximum Iterations of LP IP:        " << params.lp_maximum_iterations << std::endl;
    str << "  Initial Block Size of LP IP:        " << params.lp_initial_block_size << std::endl;
    str << "\nInitial Partitioning ";
    str << params.refinement << std::endl;
    return str;
  }

  std::ostream & operator<< (std::ostream& str, const MappingParameters& params) {
    str << "Mapping Parameters:                   " << std::endl;
    str << "  Target Graph File:                  " << params.target_graph_file << std::endl;
    str << "  One-To-One Mapping Strategy:        " << params.strategy << std::endl;
    str << "  Use Local Search:                   " << std::boolalpha << params.use_local_search << std::endl;
    str << "  Use Two-Phase Approach:             " << std::boolalpha << params.use_two_phase_approach << std::endl;
    str << "  Max Precomputed Steiner Tree Size:  " << params.max_steiner_tree_size << std::endl;
    str << "  Large HE Size Threshold:            " << params.large_he_threshold << std::endl;
    return str;
  }

  std::ostream & operator<< (std::ostream& str, const SharedMemoryParameters& params) {
    str << "Shared Memory Parameters:             " << std::endl;
    str << "  Number of Threads:                  " << params.num_threads << std::endl;
    if constexpr (TBBInitializer::provides_numa_information) {
      str << "  Number of used NUMA nodes:          " << TBBInitializer::instance().num_used_numa_nodes() << std::endl;
    }
    str << "  Use Localized Random Shuffle:       " << std::boolalpha << params.use_localized_random_shuffle << std::endl;
    str << "  Random Shuffle Block Size:          " << params.shuffle_block_size << std::endl;
    return str;
  }

  bool Context::isNLevelPartitioning() const {
    return partition.partition_type == N_LEVEL_GRAPH_PARTITIONING ||
      partition.partition_type == N_LEVEL_HYPERGRAPH_PARTITIONING;
  }

  bool Context::forceGainCacheUpdates() const {
    return isNLevelPartitioning() ||
      partition.mode == Mode::deep_multilevel ||
      refinement.refine_until_no_improvement;
  }

  void Context::setupPartWeights(const HypernodeWeight total_hypergraph_weight) {
    if (partition.use_individual_part_weights) {
      ASSERT(static_cast<size_t>(partition.k) == partition.max_part_weights.size());
      const HypernodeWeight max_part_weights_sum = std::accumulate(partition.max_part_weights.cbegin(),
                                                                   partition.max_part_weights.cend(), 0);
      double weight_fraction = total_hypergraph_weight / static_cast<double>(max_part_weights_sum);
      HypernodeWeight perfect_part_weights_sum = 0;
      partition.perfect_balance_part_weights.clear();
      for (const HyperedgeWeight& part_weight : partition.max_part_weights) {
        const HypernodeWeight perfect_weight = ceil(weight_fraction * part_weight);
        partition.perfect_balance_part_weights.push_back(perfect_weight);
        perfect_part_weights_sum += perfect_weight;
      }

      if (max_part_weights_sum < total_hypergraph_weight) {
        throw InvalidInputException(
          "Sum of individual part weights is less than the total hypergraph weight. "
          "Finding a valid partition is not possible.\n"
          "Total hypergraph weight: " + std::to_string(total_hypergraph_weight) + "\n"
          "Sum of part weights:     " + std::to_string(max_part_weights_sum));
      } else {
        // To avoid rounding issues, epsilon should be calculated using the sum of the perfect part weights instead of
        // the total hypergraph weight. See also recursive_bipartitioning_initial_partitioner
        partition.epsilon = std::min(0.99, max_part_weights_sum / static_cast<double>(std::max(perfect_part_weights_sum, 1)) - 1);
      }
    } else {
      partition.perfect_balance_part_weights.clear();
      partition.perfect_balance_part_weights.push_back(ceil(
              total_hypergraph_weight
              / static_cast<double>(partition.k)));
      for (PartitionID part = 1; part != partition.k; ++part) {
        partition.perfect_balance_part_weights.push_back(
                partition.perfect_balance_part_weights[0]);
      }
      partition.max_part_weights.clear();
      partition.max_part_weights.push_back((1 + partition.epsilon)
                                          * partition.perfect_balance_part_weights[0]);
      for (PartitionID part = 1; part != partition.k; ++part) {
        partition.max_part_weights.push_back(partition.max_part_weights[0]);
      }
    }
  }

  void Context::setupContractionLimit(const HypernodeWeight total_hypergraph_weight) {
    // Setup contraction limit
    if (initial_partitioning.mode == Mode::deep_multilevel) {
      coarsening.contraction_limit =
        std::max(2 * shared_memory.num_threads, static_cast<size_t>(partition.k)) *
          coarsening.contraction_limit_multiplier;
    } else {
      coarsening.contraction_limit =
              coarsening.contraction_limit_multiplier * partition.k;
    }

    // Setup maximum allowed vertex and high-degree vertex weight
    setupMaximumAllowedNodeWeight(total_hypergraph_weight);
  }

  void Context::setupMaximumAllowedNodeWeight(const HypernodeWeight total_hypergraph_weight) {
    HypernodeWeight min_block_weight = std::numeric_limits<HypernodeWeight>::max();
    for ( PartitionID part_id = 0; part_id < partition.k; ++part_id ) {
      min_block_weight = std::min(min_block_weight, partition.max_part_weights[part_id]);
    }

    double hypernode_weight_fraction =
            coarsening.max_allowed_weight_multiplier
            / coarsening.contraction_limit;
    coarsening.max_allowed_node_weight =
            std::ceil(hypernode_weight_fraction * total_hypergraph_weight);
    coarsening.max_allowed_node_weight =
            std::min(coarsening.max_allowed_node_weight, min_block_weight);
  }

  void Context::sanityCheck(const TargetGraph* target_graph) {
    if ( isNLevelPartitioning() && coarsening.algorithm == CoarseningAlgorithm::multilevel_coarsener ) {
        ALGO_SWITCH("Coarsening algorithm" << coarsening.algorithm << "is only supported in multilevel mode."
                                           << "Do you want to use the n-level version instead (Y/N)?",
                    "Partitioning with" << coarsening.algorithm
                                        << "coarsener in n-level mode is not supported!",
                    coarsening.algorithm,
                    CoarseningAlgorithm::nlevel_coarsener);
    } else if ( !isNLevelPartitioning() && coarsening.algorithm == CoarseningAlgorithm::nlevel_coarsener ) {
        ALGO_SWITCH("Coarsening algorithm" << coarsening.algorithm << "is only supported in n-Level mode."
                                           << "Do you want to use the multilevel version instead (Y/N)?",
                    "Partitioning with" << coarsening.algorithm
                                        << "coarsener in multilevel mode is not supported!",
                    coarsening.algorithm,
                    CoarseningAlgorithm::multilevel_coarsener);
    }

    ASSERT(partition.use_individual_part_weights != partition.max_part_weights.empty());
    if (partition.use_individual_part_weights && static_cast<size_t>(partition.k) != partition.max_part_weights.size()) {
      ALGO_SWITCH("Individual part weights specified, but number of parts doesn't match k."
                          << "Do you want to use k =" << partition.max_part_weights.size() << "instead (Y/N)?",
                  "Number of parts is not equal to k!",
                  partition.k,
                  partition.max_part_weights.size());
    }

    shared_memory.static_balancing_work_packages = std::clamp(shared_memory.static_balancing_work_packages, size_t(4), size_t(256));

    if ( partition.objective == Objective::steiner_tree ) {
      if ( partition.preset_type == PresetType::large_k ) {
        // steiner trees scale really badly with k (cubic with no parallelization), so we don't want to support this
        throw UnsupportedOperationException("Large k partitioning is not supported for steiner tree metric.");
      } else if ( partition.k > 64 && partition.instance_type == InstanceType::hypergraph ) {
        // larger k currently don't work correctly due to collisions in the hash table
        throw UnsupportedOperationException("Steiner tree metric on hypergraphs is currently only supported for k <= 64.");
      }
      if ( !target_graph ) {
        partition.objective = Objective::km1;
        INFO("No target graph provided for steiner tree metric. Switching to km1 metric.");
      } else {
        if ( partition.mode == Mode::deep_multilevel ) {
          ALGO_SWITCH("Partitioning mode" << partition.mode << "is not supported for steiner tree metric."
                                          << "Do you want to use the multilevel mode instead (Y/N)?",
                      "Partitioning mode" << partition.mode
                                          << "is not supported for steiner tree metric!",
                      partition.mode, Mode::direct);
        }
        if ( initial_partitioning.mode == Mode::deep_multilevel ) {
          ALGO_SWITCH("Initial partitioning mode" << partition.mode << "is not supported for steiner tree metric."
                                            << "Do you want to use the multilevel mode instead (Y/N)?",
                      "Initial partitioning mode" << partition.mode
                                            << "is not supported for steiner tree metric!",
                      partition.mode, Mode::direct);
        }
      }
      if (mapping.max_steiner_tree_size < 2) {
        mapping.max_steiner_tree_size = 2;
        INFO("For steiner tree metric, max-steiner-tree-size needs to be at least 2. Setting value to 2.");
      }
    }


    shared_memory.static_balancing_work_packages = std::clamp(shared_memory.static_balancing_work_packages, UL(4), UL(256));

    if ( partition.deterministic ) {
      // disable adaptive IP
      if ( initial_partitioning.use_adaptive_ip_runs ) {
        initial_partitioning.use_adaptive_ip_runs = false;
        WARNING("Disabling adaptive initial partitioning runs since deterministic mode is active");
      }

      // disable FM since there is no deterministic version
      if ( refinement.fm.algorithm != FMAlgorithm::do_nothing || initial_partitioning.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
        refinement.fm.algorithm = FMAlgorithm::do_nothing;
        initial_partitioning.refinement.fm.algorithm = FMAlgorithm::do_nothing;
        WARNING("Disabling FM refinement since deterministic mode is active");
      }

      // switch to deterministic algorithms
      bool switched = false;

      auto coarsening_algo = coarsening.algorithm;
      if ( coarsening_algo != CoarseningAlgorithm::do_nothing_coarsener && coarsening_algo != CoarseningAlgorithm::deterministic_multilevel_coarsener ) {
        coarsening.algorithm = CoarseningAlgorithm::deterministic_multilevel_coarsener;
        switched = true;
      }

      // refinement
      auto lp_algo = refinement.label_propagation.algorithm;
      if ( lp_algo != LabelPropagationAlgorithm::do_nothing && lp_algo != LabelPropagationAlgorithm::deterministic ) {
        refinement.label_propagation.algorithm = LabelPropagationAlgorithm::deterministic;
        switched = true;
      }
      auto jet_algo = refinement.jet.algorithm;
      if ( jet_algo != JetAlgorithm::do_nothing && jet_algo != JetAlgorithm::deterministic ) {
        refinement.jet.algorithm = JetAlgorithm::deterministic;
        switched = true;
      }
      auto rebalancing_algo = refinement.rebalancing.algorithm;
      if ( rebalancing_algo != RebalancingAlgorithm::do_nothing && rebalancing_algo != RebalancingAlgorithm::deterministic ) {
        refinement.rebalancing.algorithm = RebalancingAlgorithm::deterministic;
        switched = true;
      }

      // refinement during initial partitioning
      lp_algo = initial_partitioning.refinement.label_propagation.algorithm;
      if ( lp_algo != LabelPropagationAlgorithm::do_nothing && lp_algo != LabelPropagationAlgorithm::deterministic ) {
        initial_partitioning.refinement.label_propagation.algorithm = LabelPropagationAlgorithm::deterministic;
        switched = true;
      }
      jet_algo = initial_partitioning.refinement.jet.algorithm;
      if ( jet_algo != JetAlgorithm::do_nothing && jet_algo != JetAlgorithm::deterministic ) {
        initial_partitioning.refinement.jet.algorithm = JetAlgorithm::deterministic;
        switched = true;
      }
      rebalancing_algo = initial_partitioning.refinement.rebalancing.algorithm;
      if ( rebalancing_algo != RebalancingAlgorithm::do_nothing && rebalancing_algo != RebalancingAlgorithm::deterministic ) {
        initial_partitioning.refinement.rebalancing.algorithm = RebalancingAlgorithm::deterministic;
        switched = true;
      }

      if (switched) {
        WARNING("Switching to deterministic algorithm variants since deterministic mode is active");
      }
    }

    if ( partition.instance_type == InstanceType::UNDEFINED ) {
      partition.instance_type = to_instance_type(partition.file_format);
    }

    // Set correct gain policy type
    setupGainPolicy();

    if ( partition.preset_type == PresetType::large_k ) {
      // Silently switch to deep multilevel scheme for large k partitioning
      partition.mode = Mode::deep_multilevel;
    }
  }

  void Context::setupThreadsPerFlowSearch() {
    if ( refinement.flows.algorithm == FlowAlgorithm::flow_cutter ) {
      // = min{t, k, k * (k - 1) / 2}
      // t = number of threads
      // k * (k - 1) / 2 = maximum number of edges in the quotient graph
      refinement.flows.num_parallel_searches = std::min(shared_memory.num_threads,
        std::min(static_cast<size_t>(partition.k),
          static_cast<size_t>((partition.k * (partition.k - 1)) / 2) ));
    }
  }

  void Context::setupGainPolicy() {
    #ifndef KAHYPAR_ENABLE_SOED_METRIC
    if ( partition.objective == Objective::soed ) {
      throw InvalidParameterException(
        "SOED metric is deactivated. Add -DKAHYPAR_ENABLE_SOED_METRIC=ON to the "
        "cmake command and rebuild Mt-KaHyPar.");
    }
    #endif

    #ifndef KAHYPAR_ENABLE_STEINER_TREE_METRIC
    if ( partition.objective == Objective::steiner_tree ) {
      throw InvalidParameterException(
        "Steiner tree metric is deactivated. Add -DKAHYPAR_ENABLE_STEINER_TREE_METRIC=ON "
        "to the cmake command and rebuild Mt-KaHyPar.");
    }
    #endif

    if ( partition.instance_type == InstanceType::hypergraph ) {
      switch ( partition.objective ) {
        case Objective::km1: partition.gain_policy = GainPolicy::km1; break;
        case Objective::cut: partition.gain_policy = GainPolicy::cut; break;
        case Objective::soed: partition.gain_policy = GainPolicy::soed; break;
        case Objective::steiner_tree: partition.gain_policy = GainPolicy::steiner_tree; break;
        case Objective::UNDEFINED: partition.gain_policy = GainPolicy::none; break;
      }
    } else if ( partition.instance_type == InstanceType::graph ) {
      if ( partition.objective != Objective::cut && partition.objective != Objective::steiner_tree ) {
        partition.objective = Objective::cut;
        INFO("Current objective function is equivalent to the edge cut metric for graphs. Objective function is set to edge cut metric.");
      }
      if ( partition.objective == Objective::cut ) {
        partition.gain_policy = GainPolicy::cut_for_graphs;
      } else {
        partition.gain_policy = GainPolicy::steiner_tree_for_graphs;
      }
    }
  }

  std::ostream & operator<< (std::ostream& str, const Context& context) {
    str << "*******************************************************************************\n"
        << "*                            Partitioning Context                             *\n"
        << "*******************************************************************************\n"
        << context.partition
        << "-------------------------------------------------------------------------------\n"
        << context.preprocessing
        << "-------------------------------------------------------------------------------\n"
        << context.coarsening
        << "-------------------------------------------------------------------------------\n"
        << context.initial_partitioning
        << "-------------------------------------------------------------------------------\n"
        << context.refinement
        << "-------------------------------------------------------------------------------\n";
    if ( context.partition.objective == Objective::steiner_tree ) {
      str << context.mapping
          << "-------------------------------------------------------------------------------\n";
    }
    str << context.shared_memory
        << "-------------------------------------------------------------------------------";
    return str;
  }
}

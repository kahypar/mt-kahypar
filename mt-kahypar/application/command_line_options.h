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

#include <boost/program_options.hpp>

#if defined(_MSC_VER)
#include <process.h>
#include <Windows.h>
#else
#include <sys/ioctl.h>
#endif

#include <cctype>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/mt_kahypar.h"
#include "mt-kahypar/partition/context.h"

namespace po = boost::program_options;

namespace mt_kahypar {
namespace platform {
int getTerminalWidth() {
  int columns = 0;
#if defined(_MSC_VER)
  CONSOLE_SCREEN_BUFFER_INFO csbi;
  GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
  columns = csbi.srWindow.Right - csbi.srWindow.Left + 1;
#else
  struct winsize w = { };
  ioctl(0, TIOCGWINSZ, &w);
  columns = w.ws_col;
#endif
  return columns;
}

int getProcessID() {
#if defined(_MSC_VER)
  return _getpid();
#else
  return getpid();
#endif
}
}  // namespace platform

po::options_description createGeneralOptionsDescription(Context& context, const int num_columns) {
  po::options_description options("General Options", num_columns);
  options.add_options()
    ("seed",
    po::value<int>(&context.partition.seed)->value_name("<int>"),
    "Seed for random number generator \n"
    "(default: -1)")
    ("cmaxnet",
    po::value<HyperedgeID>(&context.partition.hyperedge_size_threshold)->value_name("<uint64_t>"),
    "Hyperedges larger than cmaxnet are ignored during partitioning process.")
    ("objective,o",
    po::value<std::string>()->value_name("<string>")->required()->notifier([&](const std::string& s) {
      if (s == "cut") {
        context.partition.objective = kahypar::Objective::cut;
      } else if (s == "km1") {
        context.partition.objective = kahypar::Objective::km1;
      }
    }),
    "Objective: \n"
    " - cut : cut-net metric \n"
    " - km1 : (lambda-1) metric")
    ("mode,m",
    po::value<std::string>()->value_name("<string>")->required()->notifier(
      [&](const std::string& mode) {
      context.partition.mode = kahypar::modeFromString(mode);
    }),
    "Partitioning mode: \n"
    " - (recursive) bisection \n"
    " - (direct) k-way");
  return options;
}

po::options_description createGenericOptionsDescription(Context& context,
                                                        const int num_columns) {
  po::options_description generic_options("Generic Options", num_columns);
  generic_options.add_options()
    ("help", "show help message")
    ("verbose,v", po::value<bool>(&context.partition.verbose_output)->value_name("<bool>"),
    "Verbose main partitioning output")
    ("quiet,q", po::value<bool>(&context.partition.quiet_mode)->value_name("<bool>"),
    "Quiet Mode: Completely suppress console output")
    ("show-detailed-timings", po::value<bool>(&context.partition.detailed_timings)->value_name("<bool>"),
    "If true, detailed timings overview is shown")
    ("enable-progress-bar", po::value<bool>(&context.partition.enable_progress_bar)->value_name("<bool>"),
    "If true, than progress bar is shown during uncoarsening")
    ("time-limit", po::value<int>(&context.partition.time_limit)->value_name("<int>"),
    "Time limit in seconds")
    ("sp-process,s", po::value<bool>(&context.partition.sp_process_output)->value_name("<bool>"),
    "Summarize partitioning results in RESULT line compatible with sqlplottools "
    "(https://github.com/bingmann/sqlplottools)");
  return generic_options;
}

po::options_description createPreprocessingOptionsDescription(Context& context, const int num_columns) {
  po::options_description options("Preprocessing Options", num_columns);
  options.add_options()
    ("p-use-community-structure-from-file",
    po::value<bool>(&context.preprocessing.use_community_structure_from_file)->value_name("<bool>"),
    "If true, than community structure is read from file <path-to-hypergraph>/<hypergraph-name>.community")
    ("p-enable-community-detection",
    po::value<bool>(&context.preprocessing.use_community_detection)->value_name("<bool>"),
    "If true, community detection is used as preprocessing step to guide contractions in coarsening phase")
    ("p-community-load-balancing-strategy",
    po::value<std::string>()->value_name("<string>")->notifier(
      [&](const std::string& strategy) {
      context.preprocessing.community_detection.load_balancing_strategy = communityLoadBalancingStrategyFromString(strategy);
    }),
    "Community load balancing strategies:\n"
    "- size_constraint\n"
    "- label_propagation\n"
    "- none")
    ("p-community-size-constraint-factor",
    po::value<size_t>(&context.preprocessing.community_detection.size_constraint_factor)->value_name("<size_t>"),
    "If load balancing strategy is 'size_constraint', than a community is not allowed to have an volume\n"
    "greater than total_volume / ( size_constraint_factor * num_threads )")
    ("p-louvain-edge-weight-function",
    po::value<std::string>()->value_name("<string>")->notifier(
      [&](const std::string& type) {
      context.preprocessing.community_detection.edge_weight_function = louvainEdgeWeightFromString(type);
    }),
    "Louvain edge weight functions:\n"
    "- hybrid\n"
    "- uniform\n"
    "- non_uniform\n"
    "- degree")
    ("p-max-louvain-pass-iterations",
    po::value<uint32_t>(&context.preprocessing.community_detection.max_pass_iterations)->value_name("<uint32_t>"),
    "Maximum number of iterations over all nodes of one louvain pass")
    ("p-louvain-min-eps-improvement",
    po::value<long double>(&context.preprocessing.community_detection.min_eps_improvement)->value_name("<long double>"),
    "Minimum improvement of quality during a louvain pass which leads to further passes")
    ("p-enable-community-redistribution",
    po::value<bool>(&context.preprocessing.use_community_redistribution)->value_name("<bool>"),
    "If true, hypergraph is redistributed based on community information to numa nodes")
    ("p-community-redistribution-objective",
    po::value<std::string>()->value_name("<string>")->notifier(
      [&](const std::string& objective) {
      context.preprocessing.community_redistribution.assignment_objective = mt_kahypar::communityAssignmentObjectiveFromString(objective);
    }),
    "Objective used during community redistribution of hypergraph: \n"
    " - vertex_objective \n"
    " - pin_objective")
    ("p-community-redistribution-strategy",
    po::value<std::string>()->value_name("<string>")->notifier(
      [&](const std::string& strategy) {
      context.preprocessing.community_redistribution.assignment_strategy = mt_kahypar::communityAssignmentStrategyFromString(strategy);
    }),
    "Strategy used during community redistribution of hypergraph: \n"
    " - bin_packing");
  return options;
}

po::options_description createCoarseningOptionsDescription(Context& context,
                                                           const int num_columns) {
  po::options_description options("Coarsening Options", num_columns);
  options.add_options()
    ("c-type",
    po::value<std::string>()->value_name("<string>")->notifier(
      [&](const std::string& ctype) {
      context.coarsening.algorithm = mt_kahypar::coarseningAlgorithmFromString(ctype);
    }),
    "Coarsening Algorithm:\n"
    " - community_coarsener\n"
    " - multilevel_coarsener")
    ("c-s",
    po::value<double>(&context.coarsening.max_allowed_weight_multiplier)->value_name("<double>"),
    "The maximum weight of a vertex in the coarsest hypergraph H is:\n"
    "(s * w(H)) / (t * k)\n")
    ("c-s-high-degree",
    po::value<double>(&context.coarsening.max_allowed_high_degree_node_weight_multiplier)->value_name("<double>"),
    "The maximum weight of a high degree vertex in the coarsest hypergraph H is:\n"
    "(s * w(H)) / (t * k)\n")
    ("c-t",
    po::value<HypernodeID>(&context.coarsening.contraction_limit_multiplier)->value_name("<int>"),
    "Coarsening stops when there are no more than t * k hypernodes left")
    ("c-shrink-factor",
    po::value<double>(&context.coarsening.multilevel_shrink_factor)->value_name("<double>"),
    "Multilevel coarsener creates a new hierarchy, if number of nodes is below |V| / shrink_factor")
    ("c-use-high-degree-vertex-threshold",
    po::value<bool>(&context.coarsening.use_high_degree_vertex_threshold)->value_name("<bool>"),
    "If true, than all hypernodes with a degree greater than mean + 5 * stdev are skipped during coarsening")
    ("c-rating-score",
    po::value<std::string>()->value_name("<string>")->notifier(
      [&](const std::string& rating_score) {
      context.coarsening.rating.rating_function =
        mt_kahypar::ratingFunctionFromString(rating_score);
    }), "Rating function used to calculate scores for vertex pairs:\n"
        "- heavy_edge")
    ("c-rating-heavy-node-penalty",
    po::value<std::string>()->value_name("<string>")->notifier(
      [&](const std::string& penalty) {
      context.coarsening.rating.heavy_node_penalty_policy =
        heavyNodePenaltyFromString(penalty);
    }),
    "Penalty function to discourage heavy vertices:\n"
    "- multiplicative\n"
    "- no_penalty\n"
    "- edge_frequency_penalty")
    ("c-rating-acceptance-criterion",
    po::value<std::string>()->value_name("<string>")->notifier(
      [&](const std::string& crit) {
      context.coarsening.rating.acceptance_policy =
        acceptanceCriterionFromString(crit);
    }),
    "Acceptance/Tiebreaking criterion for contraction partners having the same score:\n"
    "- best\n"
    "- best_prefer_unmatched");
  return options;
}

po::options_description createInitialPartitioningOptionsDescription(Context& context, const int num_columns) {
  po::options_description options("Initial Partitioning Options", num_columns);
  options.add_options()
    ("i-mode",
    po::value<std::string>()->value_name("<string>")->notifier(
      [&](const std::string& mode) {
      context.initial_partitioning.mode = initialPartitioningModeFromString(mode);
    }),
    "Mode of initial partitioning:\n"
    "- direct\n"
    "- recursive\n"
    "- recursive_bisection")
    ("i-runs",
    po::value<size_t>(&context.initial_partitioning.runs)->value_name("<size_t>"),
    "Number of runs for initial partitioner \n"
    "(default: 1)")
    ("i-use-adaptive-epsilon",
    po::value<bool>(&context.initial_partitioning.use_adaptive_epsilon)->value_name("<bool>"),
    "If true, adaptive epsilon is used during recursive initial partitioning \n"
    "(default: false)")
    ("i-lp-maximum-iterations",
    po::value<size_t>(&context.initial_partitioning.lp_maximum_iterations)->value_name("<size_t>"),
    "Maximum number of iterations of label propagation initial partitioner \n"
    "(default: 1)")
    ("i-lp-initial-block-size",
    po::value<size_t>(&context.initial_partitioning.lp_initial_block_size)->value_name("<size_t>"),
    "Initial block size used for label propagation initial partitioner \n"
    "(default: 1)");
  return options;
}

po::options_description createRefinementOptionsDescription(Context& context, const int num_columns) {
  po::options_description options("Refinement Options", num_columns);
  options.add_options()
    ("r-use-batch-uncontractions",
    po::value<bool>(&context.refinement.use_batch_uncontractions)->value_name("<bool>"),
    "If true, contractions are reversed in parallel. The batch size is determined by r-batch-size.\n"
    "(default false)")
    ("r-batch-size",
    po::value<size_t>(&context.refinement.batch_size)->value_name("<size_t>"),
    "Determines how many contractions are reversed in parallel if batch uncontractions are used.\n"
    "(default 1000)")
    ("r-lp-type",
    po::value<std::string>()->value_name("<string>")->notifier(
      [&](const std::string& type) {
      context.refinement.label_propagation.algorithm =
        labelPropagationAlgorithmFromString(type);
    }),
    "Algorithm used for label propagation:\n"
    "- label_propagation_km1\n"
    "- label_propagation_cut\n"
    "- do_nothing")
    ("r-lp-maximum-iterations",
    po::value<size_t>(&context.refinement.label_propagation.maximum_iterations)->value_name("<size_t>"),
    "Maximum number of iterations over all nodes during label propagation\n"
    "(default 1)")
    ("r-lp-part-weight-update-factor",
    po::value<double>(&context.refinement.label_propagation.part_weight_update_factor)->value_name("<double>"),
    "Determines after how many iterations the local part weights are updated\n"
    "=> steps = min(100, max(5, current_num_nodes * part_weight_update_factor))\n"
    "(default 100)")
    ("r-lp-localized",
    po::value<bool>(&context.refinement.label_propagation.localized)->value_name("<bool>"),
    "If true, label propagation is executed only on the previously uncontracted vertices)\n"
    "(default false)")
    ("r-lp-numa-aware",
    po::value<bool>(&context.refinement.label_propagation.numa_aware)->value_name("<bool>"),
    "If true, label propagation is executed numa friendly (which means that nodes are processed on its numa nodes)\n"
    "(default false)")
    ("r-lp-rebalancing",
    po::value<bool>(&context.refinement.label_propagation.rebalancing)->value_name("<bool>"),
    "If true, zero gain moves are used to rebalance solution\n"
    "(default true)")
    ("r-lp-execution-policy",
    po::value<std::string>()->value_name("<string>")->notifier(
      [&](const std::string& type) {
      context.refinement.label_propagation.execution_policy =
        executionTypeFromString(type);
    }),
    "Execution policy used for label propagation:\n"
    "- exponential\n"
    "- multilevel\n"
    "- constant")
    ("r-lp-execution-policy-alpha",
    po::value<double>(&context.refinement.label_propagation.execution_policy_alpha)->value_name("<double>"),
    "In case of execution policy 'exponential', LP is executed in each level which is a power of alpha\n"
    "In case of execution policy 'multilevel', LP is executed on each level INITIAL_NUM_NODES / alpha ^ i\n"
    "In case of execution policy 'constant', LP is executed on each level which is a multiple of alpha");
  return options;
}

po::options_description createSharedMemoryOptionsDescription(Context& context,
                                                             const int num_columns) {
  po::options_description shared_memory_options("Shared Memory Options", num_columns);
  shared_memory_options.add_options()
    ("s-num-threads",
    po::value<size_t>(&context.shared_memory.num_threads)->value_name("<size_t>"),
    "Number of threads used during shared memory hypergraph partitioning\n"
    "(default 1)")
    ("s-initial-hyperedge-distribution",
    po::value<std::string>()->value_name("<string>")->notifier(
      [&](const std::string& strategy) {
      context.shared_memory.initial_hyperedge_distribution = mt_kahypar::initialHyperedgeDistributionFromString(strategy);
    }),
    "Determines how hyperedges are distributed to numa nodes after reading hypergraph file: \n"
    " - equally\n"
    " - random\n"
    " - all_on_one");

  return shared_memory_options;
}

void processCommandLineInput(Context& context, int argc, char* argv[]) {
  const int num_columns = platform::getTerminalWidth();

  po::options_description generic_options = createGenericOptionsDescription(context, num_columns);

  po::options_description required_options("Required Options", num_columns);
  required_options.add_options()
    ("hypergraph,h",
    po::value<std::string>(&context.partition.graph_filename)->value_name("<string>")->required(),
    "Hypergraph filename")
    ("blocks,k",
    po::value<PartitionID>(&context.partition.k)->value_name("<int>")->required(),
    "Number of blocks")
    ("epsilon,e",
    po::value<double>(&context.partition.epsilon)->value_name("<double>")->required(),
    "Imbalance parameter epsilon");

  std::string context_path;
  po::options_description preset_options("Preset Options", num_columns);
  preset_options.add_options()
    ("preset,p", po::value<std::string>(&context_path)->value_name("<string>"),
    "Context Presets (see config directory):\n"
    " - <path-to-custom-ini-file>");

  po::options_description general_options = createGeneralOptionsDescription(context, num_columns);

  po::options_description preprocessing_options =
    createPreprocessingOptionsDescription(context, num_columns);
  po::options_description coarsening_options =
    createCoarseningOptionsDescription(context, num_columns);
  po::options_description initial_paritioning_options =
    createInitialPartitioningOptionsDescription(context, num_columns);
  po::options_description refinement_options =
    createRefinementOptionsDescription(context, num_columns);
  po::options_description shared_memory_options =
    createSharedMemoryOptionsDescription(context, num_columns);

  po::options_description cmd_line_options;
  cmd_line_options.add(generic_options)
  .add(required_options)
  .add(preset_options)
  .add(general_options)
  .add(preprocessing_options)
  .add(coarsening_options)
  .add(initial_paritioning_options)
  .add(refinement_options)
  .add(shared_memory_options);

  po::variables_map cmd_vm;
  po::store(po::parse_command_line(argc, argv, cmd_line_options), cmd_vm);

  // placing vm.count("help") here prevents required attributes raising an
  // error if only help was supplied
  if (cmd_vm.count("help") != 0 || argc == 1) {
    mt_kahypar::io::printBanner(context);
    LOG << cmd_line_options;
    exit(0);
  }

  po::notify(cmd_vm);

  std::ifstream file(context_path.c_str());
  if (!file) {
    ERROR("Could not load context file at: " + context_path);
  }

  po::options_description ini_line_options;
  ini_line_options.add(general_options)
  .add(preprocessing_options)
  .add(coarsening_options)
  .add(initial_paritioning_options)
  .add(refinement_options)
  .add(shared_memory_options);

  po::store(po::parse_config_file(file, ini_line_options, true), cmd_vm);
  po::notify(cmd_vm);

  std::string epsilon_str = std::to_string(context.partition.epsilon);
  epsilon_str.erase(epsilon_str.find_last_not_of('0') + 1, std::string::npos);

  context.partition.graph_partition_filename =
    context.partition.graph_filename
    + ".part"
    + std::to_string(context.partition.k)
    + ".epsilon"
    + epsilon_str
    + ".seed"
    + std::to_string(context.partition.seed)
    + ".KaHyPar";
  context.partition.graph_community_filename =
    context.partition.graph_filename + ".community";
}

void parseIniToContext(Context& context, const std::string& ini_filename) {
  std::ifstream file(ini_filename.c_str());
  if (!file) {
    ERROR("Could not load context file at: " << ini_filename);
  }
  const int num_columns = 80;

  po::options_description general_options =
    createGeneralOptionsDescription(context, num_columns);
  po::options_description preprocessing_options =
    createPreprocessingOptionsDescription(context, num_columns);
  po::options_description coarsening_options =
    createCoarseningOptionsDescription(context, num_columns);
  po::options_description initial_paritioning_options =
    createInitialPartitioningOptionsDescription(context, num_columns);
  po::options_description refinement_options =
    createRefinementOptionsDescription(context, num_columns);
  po::options_description shared_memory_options =
    createSharedMemoryOptionsDescription(context, num_columns);

  po::variables_map cmd_vm;
  po::options_description ini_line_options;
  ini_line_options.add(general_options)
  .add(preprocessing_options)
  .add(coarsening_options)
  .add(initial_paritioning_options)
  .add(refinement_options)
  .add(shared_memory_options);

  po::store(po::parse_config_file(file, ini_line_options, true), cmd_vm);
  po::notify(cmd_vm);
}
}  // namespace mt_kahypar

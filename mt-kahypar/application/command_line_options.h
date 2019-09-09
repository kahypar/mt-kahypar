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
#include <Windows.h>
#include <process.h>
#else
#include <sys/ioctl.h>
#endif

#include <cctype>
#include <limits>
#include <string>
#include <vector>
#include <fstream>

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/io/partitioning_output.h"

namespace po = boost::program_options;

namespace kahypar {
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
    po::value<HyperedgeID>(&context.partition.hyperedge_size_threshold)->value_name("<uint32_t>"),
    "Hyperedges larger than cmaxnet are ignored during partitioning process.")
    ("objective,o",
    po::value<std::string>()->value_name("<string>")->required()->notifier([&](const std::string& s) {
      if (s == "cut") {
        context.partition.objective = Objective::cut;
      } else if (s == "km1") {
        context.partition.objective = Objective::km1;
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
    ("time-limit", po::value<int>(&context.partition.time_limit)->value_name("<int>"),
    "Time limit in seconds")
    ("sp-process,s", po::value<bool>(&context.partition.sp_process_output)->value_name("<bool>"),
    "Summarize partitioning results in RESULT line compatible with sqlplottools "
    "(https://github.com/bingmann/sqlplottools)");
  return generic_options;
}

po::options_description createSharedMemoryOptionsDescription(Context& context,
                                                             const int num_columns) {
  po::options_description shared_memory_options("Shared Memory Options", num_columns);
  shared_memory_options.add_options()
    ("s-num-threads",
    po::value<size_t>(&context.shared_memory.num_threads)->value_name("<size_t>"),
    "Number of threads used during shared memory hypergraph partitioning\n"
    "(default 1)");
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

  po::options_description shared_memory_options =
    createSharedMemoryOptionsDescription(context, num_columns);

  po::options_description cmd_line_options;
  cmd_line_options.add(generic_options)
  .add(required_options)
  .add(preset_options)
  .add(general_options)
  .add(shared_memory_options);

  po::variables_map cmd_vm;
  po::store(po::parse_command_line(argc, argv, cmd_line_options), cmd_vm);

  // placing vm.count("help") here prevents required attributes raising an
  // error if only help was supplied
  if (cmd_vm.count("help") != 0 || argc == 1) {
    kahypar::io::printBanner(context);
    LOG << cmd_line_options;
    exit(0);
  }

  po::notify(cmd_vm);

  std::ifstream file(context_path.c_str());
  if (!file) {
    std::cerr << "Could not load context file at: " << context_path << std::endl;
    std::exit(-1);
  }

  po::options_description ini_line_options;
  ini_line_options.add(general_options)
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
}

} // namespace kahypar
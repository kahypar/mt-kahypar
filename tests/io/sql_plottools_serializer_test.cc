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

#include "gtest/gtest.h"

#include <vector>
#include <unordered_map>
#include <set>
#include <cstring>
#include <fstream>

#include "mt-kahypar/io/sql_plottools_serializer.h"
#include "mt-kahypar/io/csv_output.h"

using ::testing::Test;

namespace mt_kahypar {
namespace io {

std::vector<std::string> target_structs =
  { "PartitioningParameters", "CommunityDetectionParameters", "CommunityRedistributionParameters",
    "PreprocessingParameters", "RatingParameters", "CoarseningParameters", "InitialPartitioningParameters",
    "SparsificationParameters", "LabelPropagationParameters", "FMParameters", "NLevelGlobalFMParameters", "RefinementParameters",
    "SharedMemoryParameters" };

std::unordered_map<std::string, std::string> target_struct_prefix =
  { {"PartitioningParameters", ""}, {"CommunityDetectionParameters", "community_"}, {"CommunityRedistributionParameters", "community_redistribution_"},
    {"PreprocessingParameters", ""}, {"RatingParameters", "rating_"}, {"CoarseningParameters", "coarsening_"},
    {"InitialPartitioningParameters", "initial_partitioning_"}, {"SparsificationParameters", "sparsification_"},
    {"LabelPropagationParameters", "lp_"}, {"FMParameters", "fm_"}, {"NLevelGlobalFMParameters", "global_fm_"}, {"RefinementParameters", ""},
    {"SharedMemoryParameters", ""} };

std::set<std::string> excluded_members =
  { "verbose_output", "show_detailed_timings", "show_detailed_clustering_timings", "timings_output_depth", "show_memory_consumption", "show_advanced_cut_analysis", "enable_progress_bar", "sp_process_output",
    "measure_detailed_uncontraction_timings", "write_partition_file", "graph_partition_output_folder", "graph_partition_filename", "graph_community_filename", "community_detection",
    "community_redistribution", "coarsening_rating", "label_propagation", "lp_execute_sequential",
    "snapshot_interval", "initial_partitioning_refinement", "initial_partitioning_sparsification", "initial_partitioning_enabled_ip_algos",
    "stable_construction_of_incident_edges", "fm", "global_fm", "csv_output", "preset_file", "degree_of_parallelism" };

bool is_target_struct(const std::string& line) {
  for ( const std::string& target_struct : target_structs ) {
    if ( line.find("struct " + target_struct) != std::string::npos ) {
      return true;
    }
  }
  return false;
}

std::string get_target_struct_prefix(const std::string& line) {
  std::string prefix = "";
  for ( const std::string& target_struct : target_structs ) {
    if ( line.find("struct " + target_struct) != std::string::npos ) {
      return target_struct_prefix[target_struct];
    }
  }
  return "";
}

void read_all_members_of_target_struct(std::ifstream& context_file,
                                       const std::string& target_struct_line,
                                       std::vector<std::string>& members) {
  std::string prefix = get_target_struct_prefix(target_struct_line);
  std::string line;
  std::getline(context_file, line);
  while ( line != "};" ) {
    if ( line == "" || line.find("//") != std::string::npos ) {
      std::getline(context_file, line);
      continue;
    }

    char* input = new char[line.length() + 1];
    std::strcpy(input, line.c_str());
    if ( strcmp(input, "  #else") != 0 &&
         strcmp(input, "  #endif") != 0 &&
         strcmp(input, "  #ifdef KAHYPAR_USE_N_LEVEL_PARADIGM") != 0 &&
         strcmp(input, "  InitialPartitioningParameters() :") != 0 &&
         strcmp(input, "  InitialPartitioningParameters() :") != 0 &&
         strcmp(input, "    // Enable all initial partitioner per default") != 0 &&
         strcmp(input, "    enabled_ip_algos(static_cast<size_t>(InitialPartitioningAlgorithm::UNDEFINED), true) { }") != 0 ) {
      char* token = std::strtok(input, " ;");
      if ( strcmp(token, "mutable") == 0 ) {
        token = std::strtok(NULL, " ;");
      }
      // Second value is member name
      token = std::strtok(NULL, " ;");
      if ( strcmp(token, "double") == 0 ) { // long double
        token = std::strtok(NULL, " ;");
      }
      members.emplace_back(prefix + std::string(token));
    }
    std::getline(context_file, line);
  }
}

std::vector<std::string> get_all_members_in_context() {
  std::vector<std::string> members;

  std::ifstream context_file("../mt-kahypar/partition/context.h");
    if ( context_file ) {
    std::string line;
    while( std::getline(context_file, line) ) {
      if ( is_target_struct(line) ) {
        read_all_members_of_target_struct(context_file, line, members);
      }
    }
    context_file.close();
  } else {
    ERROR("Context file not found");
  }

  return members;
}

std::set<std::string> get_all_members_from_result_line(const std::string& result) {
  std::set<std::string> members;

  char* input = new char[result.length() + 1];
  std::strcpy(input, result.c_str());
  char* token = std::strtok(input, " =");
  token = std::strtok(NULL, " =");
  int idx = 0;
  while ( token != NULL ) {
    if ( idx % 2 == 0 ) {
      members.emplace(token);
    }
    token = std::strtok(NULL, " =");
    ++idx;
  }

  delete [] input;
  return members;
}

std::string map_context_member_to_result_line_member(const std::string& context_member) {
  if ( context_member == "perfect_balance_part_weights" ) {
    return "perfect_balance_part_weight";
  } else if ( context_member == "max_part_weights" ) {
    return "max_part_weight";
  } else if ( context_member == "graph_filename" ) {
    return "graph";
  } else if ( context_member == "community_redistribution_use_community_redistribution" ) {
    return "use_community_redistribution";
  } else if ( context_member == "rating_rating_function" ) {
    return "rating_function";
  }
  return context_member;
}

bool check_if_member_is_contained_in_result_line(const std::string& context_member,
                                                 const std::set<std::string>& members_result) {
  std::string mapped_context_member = map_context_member_to_result_line_member(context_member);
  if ( excluded_members.find(mapped_context_member) != excluded_members.end() ) {
    return true;
  } else if ( members_result.find(mapped_context_member) == members_result.end() ) {
    return false;
  }

  return true;
}

TEST(ASqlPlotSerializerTest, ChecksIfSomeParametersFromContextAreMissing) {
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  Hypergraph dummy_hypergraph;
  PartitionedHypergraph dummy_partitioned_hypergraph(2, dummy_hypergraph);
  Context dummy_context;
  dummy_context.partition.graph_filename = "dummy.hgr";
  dummy_context.partition.k = 0;
  dummy_context.partition.sp_process_output = true;
  dummy_context.partition.perfect_balance_part_weights.assign(2, 0);
  dummy_context.partition.max_part_weights.assign(2, 0);
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds(end - start);

  std::string result = serializer::serialize(dummy_partitioned_hypergraph, dummy_context, elapsed_seconds);
  std::set<std::string> members_result = get_all_members_from_result_line(result);
  std::vector<std::string> members_context = get_all_members_in_context();
  for ( const std::string& member : members_context ) {
    if ( !check_if_member_is_contained_in_result_line(member, members_result) ) {
      ERROR("Context member" << member << "not found in result line."
        << "Maybe it has a different name or should be excluded from this test.");
    }
  }
}

TEST(CSVTest, HeaderAndRowContainSameNumberOfColumns) {
  std::string header = csv::header();
  Hypergraph dummy_hypergraph;
  PartitionedHypergraph dummy_partitioned_hypergraph(2, dummy_hypergraph);
  Context dummy_context;
  dummy_context.partition.k = 2;
  dummy_context.partition.perfect_balance_part_weights.assign(2, 0);
  dummy_context.partition.max_part_weights.assign(2, 0);
  std::string body = csv::serialize(dummy_partitioned_hypergraph, dummy_context, std::chrono::duration<double>(0.2));
  ASSERT_EQ(std::count(body.begin(), body.end(), ','), std::count(header.begin(), header.end(), ','));

}

}  // namespace io
}  // namespace mt_kahypar

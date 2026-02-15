#pragma once

#include <charconv>
#include <mt-kahypar/io/hypergraph_factory.h>
#include <filesystem>
#include <fstream>

#include "dynamic_strategy.h"

namespace mt_kahypar::dyn {
  inline void generateFileName(Context& context) {
      context.dynamic.output_path = context.dynamic.result_folder;
      // add folder for graph name
      context.dynamic.output_path += context.partition.graph_filename.substr(context.partition.graph_filename.find_last_of("/\\") + 1);
      context.dynamic.output_path += "/";
      // add folder for change file
      context.dynamic.output_path += context.dynamic.changes_file.substr(context.dynamic.changes_file.find_last_of("/\\") + 1);
      context.dynamic.output_path += "/";
      // create folder
      std::filesystem::create_directories(context.dynamic.output_path);
      // add file name
      context.dynamic.output_path += std::to_string(context.partition.k) + "k_";
      context.dynamic.output_path += context.dynamic.getOutputFileName();
      // add
      context.dynamic.output_path += ".csv";
      // reset file
      std::ofstream file(context.dynamic.output_path);
      if (!context.dynamic.server)
      {
        std::cout << "Output file: " << context.dynamic.output_path << std::endl;
      }
    }

  inline size_t getCurrentRAMUsage() {
    std::ifstream status_file("/proc/self/status");
    std::string line;
    while (std::getline(status_file, line)) {
      if (line.rfind("VmRSS:", 0) == 0) {
        std::istringstream iss(line);
        std::string key;
        size_t value_kb;
        std::string unit;
        if (iss >> key >> value_kb >> unit) {
          return value_kb / 1024; // Convert from KB to MB
        }
      }
    }
    return 0; // Could not determine RAM usage
  }

  inline void print_progress_bar(size_t i, size_t total, size_t km1, double imbalance) {
      // clear the line
      std::cout << "\r\033[K" << std::flush;
      std::string output = "";
      output += "km1: " + std::to_string(km1) + ", imb: " + std::to_string(imbalance);
      output += "    [";
      for (size_t j = 0; j < 50; ++j) {
        if (j < i * 50 / total) {
          output += "#";
        } else {
          output += "-";
        }
      }
      output += "]";
      output += " " + std::to_string(i) + "/" + std::to_string(total);
      // output += " RAM: " + std::to_string(getCurrentRAMUsage()) + "MB";
      std::cout << output;
      std::cout.flush();
    }

inline std::vector<HypernodeID> parseIDs(const std::string& line) {
  std::vector<HypernodeID> ids;
  ids.reserve(2);

  const char* p = line.data();
  const char* end = p + line.size();

  while (p < end) {
    // skip whitespace
    while (p < end && *p <= ' ') ++p;
    if (p >= end) break;

    HypernodeID id;
    auto [np, ec] = std::from_chars(p, end, id);
    if (ec == std::errc()) {
      ids.push_back(id);
      p = np;
    } else {
      // skip malformed token
      while (p < end && *p > ' ') ++p;
    }
  }

  return ids;
}

    inline std::vector<PinChange> parsePins(const std::string& line) {
      std::vector<PinChange> pins;
      std::stringstream ss(line);
      std::string token;
      while (ss >> token) {
        size_t comma_pos = token.find(',');
        if (comma_pos != std::string::npos) {
          try {
            HypernodeID node = std::stoi(token.substr(1, comma_pos - 1)); // Skip '('
            HyperedgeID edge = std::stoi(token.substr(comma_pos + 1, token.size() - comma_pos - 2)); // Exclude ')'
            pins.push_back({node, edge});
          } catch (const std::invalid_argument&) {
            std::cerr << "Invalid pin format in token: " << token << "\n";
          }
        }
      }
      return pins;
    }

  inline std::vector<PinChange> parsePins_fast(const std::string& line) {
    std::vector<PinChange> pins;
    pins.reserve(16);

    const char* p   = line.data();
    const char* end = p + line.size();

    while (p < end) {
      // Skip whitespace
      while (p < end && *p <= ' ') ++p;
      if (p >= end) break;

      if (*p != '(') {
        // Malformed token, skip
        while (p < end && *p > ' ') ++p;
        continue;
      }
      ++p; // skip '('

      HypernodeID node;
      auto [p1, ec1] = std::from_chars(p, end, node);
      if (ec1 != std::errc()) break;
      p = p1;

      if (p >= end || *p != ',') break;
      ++p; // skip ','

      HyperedgeID edge;
      auto [p2, ec2] = std::from_chars(p, end, edge);
      if (ec2 != std::errc()) break;
      p = p2;

      if (p >= end || *p != ')') break;
      ++p; // skip ')'

      pins.push_back({node, edge});
    }

    return pins;
  }

    inline std::string getNextNonCommentLine(std::ifstream& file) {
      std::string line;
      while (std::getline(file, line)) {
        if (line[0] != '/') {
          return line;
        }
      }
      return "";
    }

    Change parseChange(std::ifstream& file) {
      return Change{
        parseIDs(getNextNonCommentLine(file)), // added nodes
        parseIDs(getNextNonCommentLine(file)), // added edges
        parsePins_fast(getNextNonCommentLine(file)), // added pins
        parseIDs(getNextNonCommentLine(file)), // removed nodes
        parseIDs(getNextNonCommentLine(file)), // removed edges
        parsePins_fast(getNextNonCommentLine(file)) // removed pins
      };
    }

    inline size_t getNumChanges(std::ifstream& file)
    {
      const std::string line = getNextNonCommentLine(file);
      return  std::stoi(line);
    }

    // parse the changes from a file
    //
    // Number of changes
    // added_nodes_id1 added_nodes_id2 ...
    // added_edges_id1 added_edges_id2 ...
    // (added_pins1.hypernode,added_pins1.hyperedge) (added_pins2.hypernode,added_pins2.hyperedge), ...
    // removed_nodes_id1 removed_nodes_id2, ...
    // removed_edges_id1 removed_edges_id2, ...
    // (removed_pins1.hypernode,removed_pins1.hyperedge) (removed_pins2.hypernode,removed_pins2.hyperedge), ...
    inline std::vector<Change> parseChanges(std::ifstream file, size_t num_changes) {
      std::vector<Change> changes;
      changes.reserve(num_changes);
      for (size_t i = 0; i < 10000000; ++i) {
        Change change = parseChange(file);
        changes.emplace_back(change);
      }
      file.close();
      return changes;
    }

    inline void log_km1_live(size_t i, size_t max_changes, Context& context, ds::PartitionedHypergraph<ds::MutableHypergraph>& hypergraph, std::chrono::duration<double> time) {

      size_t km1 = metrics::quality(hypergraph, context);
      double imbalance = metrics::imbalance(hypergraph, context);

      // TODO only open file once maybe in memory and only write at end
      std::string filename = context.dynamic.output_path;
      std::ofstream file(filename, std::ios_base::app);
      if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
      }
      //log format: change, km1, imbalance, time, km1_gain_vcycle, vcycle_duration_sum, km1_gain_rebalance_push, rebalance_push_duration_sum, km1_gain_rebalance_pull, rebalance_pull_duration_sum, km1_gain_localFM, localFM_duration_sum, processing_duration_sum, sorting_duration_sum, gain_cache_update_duration_sum
      file << i << ", " << km1 << ", " << imbalance << ", " << std::chrono::duration_cast<std::chrono::milliseconds>(time).count() << ", "
           << context.dynamic.km1_gain_vcycle << ", "
           << std::chrono::duration_cast<std::chrono::milliseconds>(context.dynamic.vcycle_duration_sum).count() << ", "
           << context.dynamic.km1_gain_rebalance_push << ", "
           << std::chrono::duration_cast<std::chrono::milliseconds>(context.dynamic.rebalance_duration_sum_push).count() << ", "
           << context.dynamic.km1_gain_rebalance_pull << ", "
           << std::chrono::duration_cast<std::chrono::milliseconds>(context.dynamic.rebalance_duration_sum_pull).count() << ", "
           << context.dynamic.km1_gain_localFM << ", "
           << std::chrono::duration_cast<std::chrono::milliseconds>(context.dynamic.localFM_duration_sum).count() << ", "
           << std::chrono::duration_cast<std::chrono::milliseconds>(context.dynamic.processing_duration_sum).count() << ", "
           << std::chrono::duration_cast<std::chrono::milliseconds>(context.dynamic.sorting_duration_sum).count() << ", "
           << std::chrono::duration_cast<std::chrono::milliseconds>(context.dynamic.gain_cache_update_duration_sum).count() << ", "
           << context.dynamic.move_count
           << std::endl;
      file.close();

      if (!context.dynamic.server) {
        print_progress_bar(i, max_changes, km1, imbalance);
      }

    }

    inline void initOutputFile(Context& context) {
      generateFileName(context);
      std::string filename = context.dynamic.output_path;
      std::ofstream file(filename);
      if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
      }
      file << "change, km1, imbalance, time, km1_gain_vcycle, vcycle_duration_sum, km1_gain_rebalance_push, rebalance_push_duration_sum, km1_gain_rebalance_pull, rebalance_pull_duration_sum, km1_gain_localFM, localFM_duration_sum, processing_duration_sum, sorting_duration_sum, gain_cache_update_duration_sum, move_count" << std::endl;
      file.close();
    }

    inline void generateErrorFile(Context& context, DynamicStrategy* strategy, std::exception& e) {
      (void) strategy;
      std::string filename = context.dynamic.output_path;
      std::ofstream file(filename + ".error");
      if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
      }
      file << "Dynamic strategy: " << context.dynamic.strategy << std::endl;
      file << "Error: " << e.what() << std::endl;

      file.close();
    }
}
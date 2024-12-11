#pragma once

#include <mt-kahypar/io/hypergraph_factory.h>
#include <filesystem>
#include "mt-kahypar/dynamic/dynamic_datastructures.h"
#include "mt-kahypar/dynamic/dynamic_strategy.h"

namespace mt_kahypar::dyn {

    void generateFileName(Context& context) {
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
      std::cout << "Output file: " << context.dynamic.output_path << std::endl;
    }

    void print_progress_bar(size_t i, size_t total, const std::vector<DynamicStrategy::PartitionResult>* history) {
      // clear the line
      std::cout << "\r\e[K" << std::flush;
      std::string output = "";
      output += "km1: " + std::to_string(history->back().km1) + ", imb: " + std::to_string(history->back().imbalance);
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
      std::cout << output;
      std::cout.flush();
    }

    std::vector<HypernodeID> parseIDs(const std::string& line) {
      std::vector<HypernodeID> ids;
      std::stringstream ss(line);
      int id;
      while (ss >> id) {
        ids.push_back(id);
      }
      return ids;
    }

    std::vector<PinChange> parsePins(const std::string& line) {
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

    // parse the changes from a file
    //
    // Number of changes
    // added_nodes_id1 added_nodes_id2 ...
    // added_edges_id1 added_edges_id2 ...
    // (added_pins1.hypernode,added_pins1.hyperedge) (added_pins2.hypernode,added_pins2.hyperedge), ...
    // removed_nodes_id1 removed_nodes_id2, ...
    // removed_edges_id1 removed_edges_id2, ...
    // (removed_pins1.hypernode,removed_pins1.hyperedge) (removed_pins2.hypernode,removed_pins2.hyperedge), ...
    std::vector<Change> parseChanges(const std::string& filename) {
      std::vector<Change> changes;
      std::ifstream file(filename);
      if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
      }

      size_t num_changes;
      file >> num_changes;
      file.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Move to the next line

      for (size_t i = 0; i < num_changes; ++i) {
        Change change;
        std::string line;

        // Parse added nodes
        std::getline(file, line);
        change.added_nodes = parseIDs(line);

        // Parse added edges
        std::getline(file, line);
        change.added_edges = parseIDs(line);

        // Parse added pins
        std::getline(file, line);
        change.added_pins = parsePins(line);

        // Parse removed nodes
        std::getline(file, line);
        change.removed_nodes = parseIDs(line);

        // Parse removed edges
        std::getline(file, line);
        change.removed_edges = parseIDs(line);

        // Parse removed pins
        std::getline(file, line);
        change.removed_pins = parsePins(line);

        // Add the parsed change to the list
        changes.push_back(change);
      }

      file.close();
      return changes;
    }


    void resetHypergraph(ds::StaticHypergraph& hypergraph_s, std::vector<Change>& changes, Context& context) {
      //iterate backwards through the changes and undo them
      for (auto it = changes.rbegin(); it != changes.rend(); ++it) {
        Change& change = *it;
        for (const PinChange& pin: change.added_pins) {
          hypergraph_s.removePin(pin.node, pin.edge);
        }
        for (const HypernodeID& hn : change.added_nodes) {
          hypergraph_s.disableHypernodeWithEdges(hn);
          if (!context.dynamic.use_final_weight) {
            hypergraph_s.decrementTotalWeight(hn);
          }
        }
        for (const HyperedgeID& he : change.added_edges) {
          hypergraph_s.removeEdge(he);
        }
        for (const PinChange& pin: change.removed_pins) {
          hypergraph_s.restorePin(pin.node, pin.edge);
        }
        for (const HypernodeID& hn : change.removed_nodes) {
          hypergraph_s.enableHypernodeWithEdges(hn);
          if (!context.dynamic.use_final_weight) {
            hypergraph_s.incrementTotalWeight(hn);
          }
        }
        for (const HyperedgeID& he : change.removed_edges) {
          hypergraph_s.restoreEdge(he);
        }
      }
    }

    std::vector<Change> generateHypernodeChanges(ds::StaticHypergraph& hypergraph_s, Context& context) {

      std::vector<HypernodeID> disabling_order(hypergraph_s.initialNumNodes());
      std::iota(disabling_order.begin(), disabling_order.end(), 0);

      // Shuffle the order of the nodes using seed using std::mt19937
      std::shuffle(disabling_order.begin(), disabling_order.end(), std::mt19937(context.partition.seed));

      size_t start_id = context.dynamic.initial_partitioning_size;

      std::vector<Change> added_nodes_changes;

      // Disable all the nodes using seed
      for ( const HypernodeID& hn : disabling_order ) {
        if ( !hypergraph_s.nodeIsEnabled(hn) ) {
          continue;
        }
        std::vector<HyperedgeID> edges = hypergraph_s.incidentEdgesCopy(hn);
        std::vector<HypernodeID> remove_prio_nodes;
        std::vector<PinChange> removed_pins;
        std::vector<HyperedgeID> removed_edges;
        for ( const HyperedgeID& he : edges ) {
          ASSERT(hypergraph_s.edgeIsEnabled(he));
          hypergraph_s.removePin(hn, he);
          removed_pins.push_back({hn, he});
          ASSERT(hypergraph_s.edgeSize(he) > 0);
          if ( hypergraph_s.edgeSize(he) == 1 ) {
            HypernodeID hn2 = *hypergraph_s.pins(he).begin();
            hypergraph_s.removePin(hn2, he);
            hypergraph_s.removeEdge(he);
            removed_pins.push_back({hn2, he});
            removed_edges.push_back(he);
            if (context.dynamic.reduce_deg_0) {
              //remove node if it is degree 0
              if (hypergraph_s.incidentEdges(hn2).empty()) {
                if (!hypergraph_s.nodeIsEnabled(hn2)) {
                  continue;
                }
                hypergraph_s.disableHypernodeWithEdges(hn2);
                if (!context.dynamic.use_final_weight) {
                  hypergraph_s.decrementTotalWeight(hn2);
                }
                remove_prio_nodes.push_back(hn2);
              }
            }
          }
        }
        hypergraph_s.disableHypernodeWithEdges(hn);
        if ( !context.dynamic.use_final_weight ) {
          hypergraph_s.decrementTotalWeight(hn);
        }
        added_nodes_changes.push_back({{hn}, removed_edges, removed_pins, {}, {}, {}});
        for ( const HypernodeID& hn2 : remove_prio_nodes ) {
          added_nodes_changes.push_back({{hn2}, {}, {}, {}, {}, {}});
        }
      }

      // reverse the order of added_nodes_changes
      std::reverse(added_nodes_changes.begin(), added_nodes_changes.end());

      // re-enable all nodes until start_id
      for ( size_t i = 0; i < start_id; ++i ) {
        Change& change = added_nodes_changes[i];
        ASSERT(change.added_nodes.size() == 1);
        HypernodeID hn = change.added_nodes[0];
        ASSERT(hn < hypergraph_s.initialNumNodes());
        hypergraph_s.enableHypernodeWithEdges(hn);
        for ( const HyperedgeID& he : change.added_edges ) {
          hypergraph_s.restoreEdge(he);
        }
        for ( const PinChange& pin: change.added_pins ) {
          hypergraph_s.restorePin(pin.node, pin.edge);
        }
        if ( !context.dynamic.use_final_weight ) {
          hypergraph_s.incrementTotalWeight(hn);
        }
      }

      //remove all re-enabled nodes from added_nodes
      added_nodes_changes.erase(added_nodes_changes.begin(), added_nodes_changes.begin() + start_id);

      return added_nodes_changes;
    }

    //currently disables and re-enables all edges
    std::vector<Change> generateHyperedgeChanges(ds::StaticHypergraph& hypergraph_s, Context& context) {
      std::vector<HyperedgeID> disabling_order(hypergraph_s.initialNumEdges());
      std::iota(disabling_order.begin(), disabling_order.end(), 0);

      // Shuffle the order of the edges using seed
      std::shuffle(disabling_order.begin(), disabling_order.end(), std::mt19937(context.partition.seed));

      std::vector<Change> added_edges_changes;

      // Disable all the edges using seed
      for ( const HyperedgeID& he : disabling_order ) {
        if ( !hypergraph_s.edgeIsEnabled(he) ) {
          continue;
        }
        std::vector<HypernodeID> pins = hypergraph_s.pinsCopy(he);
        std::vector<PinChange> pins_changes;
        for ( const HypernodeID& hn : pins ) {
          ASSERT(hypergraph_s.nodeIsEnabled(hn));
          hypergraph_s.removePin(hn, he);
          pins_changes.push_back({hn, he});
        }
        hypergraph_s.removeEdge(he);
        added_edges_changes.push_back({{}, {he}, pins_changes, {}, {}, {}});
      }

      // reverse the order of added_edges_changes
      std::reverse(added_edges_changes.begin(), added_edges_changes.end());

      // re-enable all edges until start_id
      for (size_t i = 0; i < added_edges_changes.size(); ++i ) {
        Change& change = added_edges_changes[i];
        HyperedgeID he = change.added_edges[0];
        ASSERT(he < hypergraph_s.initialNumEdges());
        hypergraph_s.restoreEdge(he);
        for ( const PinChange& pin: change.added_pins ) {
          hypergraph_s.restorePin(pin.node, pin.edge);
        }
      }

      //remove all re-enabled edges from added_edges
      added_edges_changes.clear();

      return added_edges_changes;
    }

    //currently only disables and re-enables all pins
    std::vector<Change> generatePinChanges(ds::StaticHypergraph& hypergraph_s, Context& context) {
      std::vector<HyperedgeID> disabling_order(hypergraph_s.initialNumEdges());
      std::iota(disabling_order.begin(), disabling_order.end(), 0);

      // Shuffle the order of the pins using seed
      std::shuffle(disabling_order.begin(), disabling_order.end(), std::mt19937(context.partition.seed));

      std::vector<Change> added_pins_changes;
      std::vector<PinChange> removed_pins;
      std::vector<HyperedgeID> removed_edges;

      // Disable all the pins using seed
      for ( const auto& he: disabling_order ) {
        removed_pins.clear();
        removed_edges.clear();
        if ( !hypergraph_s.edgeIsEnabled(he) ) {
          continue;
        }
        std::vector<HypernodeID> pins = hypergraph_s.pinsCopy(he);
        for ( const auto& hn: pins ) {
          hypergraph_s.removePin(hn, he);
          removed_pins.push_back({hn, he});
          ASSERT(hypergraph_s.edgeSize(he) > 0);
          if ( hypergraph_s.edgeSize(he) == 1 ) {
            HypernodeID hn2 = *hypergraph_s.pins(he).begin();
            hypergraph_s.removePin(hn2, he);
            hypergraph_s.removeEdge(he);
            removed_pins.push_back({hn2, he});
            removed_edges.push_back(he);
            break;
          }
        }
        added_pins_changes.push_back({{}, removed_edges, removed_pins, {}, {}, {}});
      }

      // reverse the order of added_pins_changes
      std::reverse(added_pins_changes.begin(), added_pins_changes.end());

      // re-enable all pins
      for ( Change& change : added_pins_changes ) {
        for ( const HyperedgeID he : change.added_edges ) {
          hypergraph_s.restoreEdge(he);
        }
        for ( const PinChange& pin: change.added_pins ) {
          hypergraph_s.restorePin(pin.node, pin.edge);
        }
      }

      //remove all re-enabled pins from added_pins
      added_pins_changes.clear();

      return added_pins_changes;
    }


    std::vector<Change> generateChanges(ds::StaticHypergraph& hypergraph_s, Context& context) {

      std::cout << "Adding nodes" << std::endl;
      std::vector<Change> added_nodes = generateHypernodeChanges(hypergraph_s, context);
      std::cout << "Adding edges" << std::endl;
      std::vector<Change> added_edges = generateHyperedgeChanges(hypergraph_s, context);
      std::cout << "Adding pins" << std::endl;
      std::vector<Change> added_pins = generatePinChanges(hypergraph_s, context);

      //merge all changes
      std::vector<Change> changes;
      for (size_t i = 0; i < std::max({added_nodes.size(), added_edges.size(), added_pins.size()}); ++i) {
        Change change;
        if (i < added_nodes.size()) {
          change.append(added_nodes[i]);
        }
        if (i < added_edges.size()) {
          change.append(added_edges[i]);
        }
        if (i < added_pins.size()) {
          change.append(added_pins[i]);
        }
        changes.push_back(change);
      }

      return changes;
    }

    void log_km1(Context& context, const std::vector<DynamicStrategy::PartitionResult>* history) {
      std::string filename = context.dynamic.output_path;
      std::ofstream file(filename);
      if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
      }
      for (size_t i = 0; i < history->size(); ++i) {
        if (!history->at(i).valid) {
          continue;
        }
        file << i << ", " << history->at(i).km1 << ", " << history->at(i).imbalance << std::endl;
      }
      file.close();
    }

    void log_km1_live(size_t i, Context& context, DynamicStrategy::PartitionResult result, std::chrono::duration<double> time) {
      if (!result.valid) {
        return;
      }
      std::string filename = context.dynamic.output_path;
      std::ofstream file(filename, std::ios_base::app);
      if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
      }
      file << i << ", " << result.km1 << ", " << result.imbalance << ", " << time.count() << std::endl;
      file.close();
    }

    void initOutputFile(Context& context) {
      generateFileName(context);
      std::string filename = context.dynamic.output_path;
      std::ofstream file(filename);
      if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
      }
      file << "change, km1, imbalance, time" << std::endl;
      file.close();
    }

    void generateErrorFile(Context& context, DynamicStrategy* strategy, std::exception& e) {
      std::string filename = context.dynamic.output_path;
      std::ofstream file(filename + ".error");
      if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
      }
      file << "Error: " << e.what() << std::endl;

      //print the history length
      file << "History length: " << strategy->history.size() << std::endl;

      // print the last result
      if (!strategy->history.empty()) {
        file << "Last result: " << strategy->history.back().km1 << ", " << strategy->history.back().imbalance << std::endl;
      }
      file.close();
    }
}
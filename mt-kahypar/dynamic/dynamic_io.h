#pragma once

#include <mt-kahypar/io/hypergraph_factory.h>
#include "mt-kahypar/dynamic/dynamic_datastructures.h"

namespace mt_kahypar::dyn {

    std::vector<HypernodeID> generateHypernodeChanges(ds::StaticHypergraph& hypergraph_s, Context& context) {

      std::vector<HypernodeID> disabling_order(hypergraph_s.initialNumNodes());
      std::iota(disabling_order.begin(), disabling_order.end(), 0);

      // Shuffle the order of the nodes using seed using std::random
      std::shuffle(disabling_order.begin(), disabling_order.end(), std::default_random_engine(context.partition.seed));

      size_t start_id = context.dynamic.initial_partitioning_size;

      std::vector<HypernodeID> added_nodes;

      // Disable all the nodes using seed
      for ( const HypernodeID& hn : disabling_order ) {
        hypergraph_s.disableHypernodeWithEdges(hn);
        if ( !context.dynamic.use_final_weight ) {
          hypergraph_s.decrementTotalWeight(hn);
        }
        added_nodes.push_back(hn);
      }

      // re-enable all nodes until start_id
      for ( size_t i = 0; i < start_id; ++i ) {
        HypernodeID hn = added_nodes[i];
        ASSERT(hn < hypergraph_s.initialNumNodes());
        hypergraph_s.enableHypernodeWithEdges(hn);
        if ( !context.dynamic.use_final_weight ) {
          hypergraph_s.incrementTotalWeight(hn);
        }
      }

      //remove all re-enabled nodes from added_nodes
      added_nodes.erase(added_nodes.begin(), added_nodes.begin() + start_id);

      return added_nodes;
    }

    //currently disables and re-enables all edges
    std::vector<HyperedgeID> generateHyperedgeChanges(ds::StaticHypergraph& hypergraph_s, Context& context) {
      std::vector<HyperedgeID> disabling_order(hypergraph_s.initialNumEdges());
      std::iota(disabling_order.begin(), disabling_order.end(), 0);

      // Shuffle the order of the edges using seed
      std::shuffle(disabling_order.begin(), disabling_order.end(), std::default_random_engine(context.partition.seed));

      //TODO: Add edge parameter
      size_t start_id = hypergraph_s.initialNumEdges();

      std::vector<HyperedgeID> added_edges;

      // Disable all the edges using seed
      for ( const HyperedgeID& he : disabling_order ) {
        hypergraph_s.removeEdge(he);
        added_edges.push_back(he);
      }

      // re-enable all edges until start_id
      for ( size_t i = 0; i < start_id; ++i ) {
        HyperedgeID& he = added_edges[i];
        ASSERT(he < hypergraph_s.initialNumEdges());
        hypergraph_s.restoreEdge(he);
      }

      //remove all re-enabled edges from added_edges
      added_edges.erase(added_edges.begin(), added_edges.begin() + start_id);

      return added_edges;
    }

    //currently only disables and re-enables all pins
    std::vector<PinChange> generatePinChanges(ds::StaticHypergraph& hypergraph_s, Context& context) {
      std::vector<HyperedgeID> disabling_order(hypergraph_s.initialNumEdges());
      std::iota(disabling_order.begin(), disabling_order.end(), 0);

      // Shuffle the order of the pins using seed
      std::shuffle(disabling_order.begin(), disabling_order.end(), std::default_random_engine(context.partition.seed));

      //TODO: Add pin parameter
      size_t start_id = hypergraph_s.initialNumPins();
      (void) start_id; // suppress warning (unused variable

      std::vector<PinChange> added_pins;

      // Disable all the pins using seed
      for ( const auto& he: disabling_order ) {
        if ( !hypergraph_s.edgeIsEnabled(he) ) {
          continue;
        }
        for ( const auto& hn: hypergraph_s.pins(he) ) {
          if ( !hypergraph_s.nodeIsEnabled(hn) ) {
            continue;
          }
          hypergraph_s.removePin(hn, he);
          added_pins.push_back({hn, he});
        }
      }

      // re-enable all pins until start_id
      for ( auto [hn, he]: added_pins ) {
        hypergraph_s.restorePin(hn, he);
        ASSERT(std::find(hypergraph_s.pins(he).begin(), hypergraph_s.pins(he).end(), hn) != hypergraph_s.pins(he).end());
      }

      //remove all re-enabled pins from added_pins
      added_pins.clear();

      return added_pins;
    }


    std::tuple<std::vector<Change>, mt_kahypar_hypergraph_t> generateChanges(Context& context) {

      // Read Hypergraph
      mt_kahypar_hypergraph_t hypergraph = io::readInputFile(
              context.partition.graph_filename, context.partition.preset_type,
              context.partition.instance_type, context.partition.file_format,
              context.preprocessing.stable_construction_of_incident_edges);

      auto& hypergraph_s = utils::cast<ds::StaticHypergraph>(hypergraph);

      std::vector<HypernodeID> added_nodes = generateHypernodeChanges(hypergraph_s, context);
      std::vector<HyperedgeID> added_edges = generateHyperedgeChanges(hypergraph_s, context);
      std::vector<PinChange> added_pins = generatePinChanges(hypergraph_s, context);

      //merge all changes
      std::vector<Change> changes;
      for (size_t i = 0; i < std::max({added_nodes.size(), added_edges.size(), added_pins.size()}); ++i) {
        Change change;
        if (i < added_nodes.size()) {
          change.added_nodes.push_back(added_nodes[i]);
        }
        if (i < added_edges.size()) {
          change.added_edges.push_back(added_edges[i]);
        }
        if (i < added_pins.size()) {
          change.added_pins.push_back(added_pins[i]);
        }
        changes.push_back(change);
      }

      return {changes, hypergraph};
    }

    void log_km1(Context& context, const std::vector<DynamicStrategy::PartitionResult>* history) {
      //TODO change initial_partitioning_size to useful value
      std::string filename = context.dynamic.result_folder + context.dynamic.strategy + "_" + std::to_string(context.dynamic.initial_partitioning_size) + "_" + std::to_string(context.partition.k) + "k" + (context.dynamic.use_final_weight ? "_final_weight" : "");
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

    void log_live_km1(Context& context, const std::vector<DynamicStrategy::PartitionResult>* history) {
      (void) context;
      ASSERT(!history->empty());
      if (history->back().valid) {
        std::cout << history->size() << ", " << history->back().km1 << ", " << history->back().imbalance << std::endl;
      }
    }
}
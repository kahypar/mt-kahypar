#pragma once

#include "mt-kahypar/dynamic/dynamic_datastructures.h"

namespace mt_kahypar::dyn {

    std::tuple<std::vector<Change>, mt_kahypar_hypergraph_t> generateChanges(Context& context) {

      // Read Hypergraph
      mt_kahypar_hypergraph_t hypergraph = io::readInputFile(
              context.partition.graph_filename, context.partition.preset_type,
              context.partition.instance_type, context.partition.file_format,
              context.preprocessing.stable_construction_of_incident_edges);

      auto& hypergraph_s = utils::cast<ds::StaticHypergraph>(hypergraph);

      std::vector<HypernodeID> disabling_order(hypergraph_s.initialNumNodes());
      std::iota(disabling_order.begin(), disabling_order.end(), 0);

      // Shuffle the order of the nodes using seed
      utils::Randomize::instance().shuffleVector(disabling_order);

      std::vector<Change> changes;
      size_t start_id = context.dynamic.initial_partitioning_size;

      // Disable all the nodes using seed
      for ( const HypernodeID& hn : disabling_order ) {
        hypergraph_s.disableHypernodeWithEdges(hn);
        if ( !context.dynamic.use_final_weight ) {
          hypergraph_s.decrementTotalWeight(hn);
        }
        changes.emplace_back(Change{{hn}});
      }

      // re-enable all nodes until start_id
      for ( size_t i = 0; i < start_id; ++i ) {
        HypernodeID hn = changes[i].added_nodes[0];
        ASSERT(hn < hypergraph_s.initialNumNodes());
        hypergraph_s.enableHypernodeWithEdges(hn);
        if ( !context.dynamic.use_final_weight ) {
          hypergraph_s.incrementTotalWeight(hn);
        }
      }
      return {changes, hypergraph};
    }
}
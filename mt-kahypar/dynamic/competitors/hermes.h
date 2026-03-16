#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>

namespace mt_kahypar::dyn {

    class Hermes : public DynamicStrategy {

    private:

        PartitionID assign_node_first_free_partition(const HypernodeID& hn) {
          for (PartitionID p = 0; p < context.partition.k; ++p) {
            if (partitioned_hypergraph_m.partWeight(p) + hypergraph_m.nodeWeight(hn) <=
                context.partition.max_part_weights[p]) {
              partitioned_hypergraph_m.addNode(hn, p);
              return p;
            }
          }
          // if no partition could accommodate the node put in the first
          partitioned_hypergraph_m.addNode(hn, 0);
          return 0;
        }

        PartitionID getTargetPart(const HypernodeID& hn, const size_t stage) {
          const PartitionID current_part = partitioned_hypergraph_m.partID(hn);
          ASSERT(current_part != kInvalidPartition);
          size_t max_gain = 0;
          PartitionID target_part = current_part;
          for ( PartitionID p = 0; p < context.partition.k; ++p ) {
            if (p != current_part) {
              if (stage == 1 && p > current_part || stage == 2 && p < current_part) {
                size_t gain = 0;
                for ( const HyperedgeID& he : hypergraph_m.incidentEdges(hn) ) {
                  if (partitioned_hypergraph_m.pinCountInPart(he, p) == 0) {
                    gain -= hypergraph_m.edgeWeight(he);
                  }
                  if (partitioned_hypergraph_m.pinCountInPart(he, current_part) == 1) {
                    gain += hypergraph_m.edgeWeight(he);
                  }
                }
                if (gain > max_gain && partitioned_hypergraph_m.partWeight(p) + hypergraph_m.nodeWeight(hn) <=
                    context.partition.max_part_weights[p]) {
                  max_gain = gain;
                  target_part = p;
                }
              }
            }
          }
          return target_part;
        }

    public:

      Hermes(ds::MutableHypergraph& hypergraph_m, Context& context)
          : DynamicStrategy(hypergraph_m, context) {}
      
      MutablePartitionedHypergraph& init() override {
          partitioned_hypergraph_m = partition_hypergraph_km1(hypergraph_m, context);
        ASSERT(partitioned_hypergraph_m.checkAllConnectivitySets());
        return partitioned_hypergraph_m;
      }

      void partition(Change& change, size_t changes_size) override {
        (void) changes_size;
        ASSERT(metrics::isBalanced(partitioned_hypergraph_m, context));
        // ASSERT(partitioned_hypergraph_m.checkAllConnectivitySets());

        // Removals are not supported in streaming mode

        ASSERT(change.removed_nodes.empty());
        ASSERT(change.removed_edges.empty());
        ASSERT(change.removed_pins.empty());

        for (const HypernodeID& hn : change.added_nodes) {
          const HypernodeID new_hn = hypergraph_m.addHypernode({}, 1);
          (void) new_hn;
          ASSERT(hn == new_hn);
          updateMaxPartWeight(context, hypergraph_m);
          const PartitionID assigned_part = assign_node_first_free_partition(hn);
          (void) assigned_part;
          ASSERT(assigned_part != kInvalidPartition);
        }

        for (const HyperedgeID& he : change.added_edges) {
          hypergraph_m.addHyperedge({}, 1);
          partitioned_hypergraph_m.addEdge(he);
        }

        for (const auto& [node, edge] : change.added_pins)
        {
          hypergraph_m.addPin(edge, node);
          partitioned_hypergraph_m.incrementPinCountOfBlockWrapper(edge, partitioned_hypergraph_m.partID(node));
        }

        // for (const HypernodeID& hn : change.added_nodes)
        for (const HypernodeID& hn : hypergraph_m.nodes())
        {
          for (size_t stage = 1; stage <= 2; ++stage)
          {
            const PartitionID target_part = getTargetPart(hn, stage);
            if (target_part != partitioned_hypergraph_m.partID(hn)) {
              context.dynamic.move_count++;
              partitioned_hypergraph_m.changeNodePart(hn, partitioned_hypergraph_m.partID(hn), target_part);
            }
          }
        }
      }
    };
}

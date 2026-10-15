#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>

namespace mt_kahypar::dyn {

    class Hermes : public DynamicStrategy {

    private:
        size_t change_count = 0;
        HyperedgeWeight prior_total_weight = 0;
        double changed_weight = 0;

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


      PartitionID move_node_to_best_partition(const HypernodeID& hn) {
          const PartitionID current_part = partitioned_hypergraph_m.partID(hn);
          ASSERT(current_part != kInvalidPartition);

          // compute for each block the number of nodes the new node connected to
          std::vector block_connectivities(context.partition.k, std::make_tuple(0,0));
          for ( PartitionID p = 0; p < context.partition.k; ++p ) {
            block_connectivities[p] = std::make_tuple(0, p);
          }
          for ( const HyperedgeID& he : hypergraph_m.incidentEdges(hn) ) {
            ASSERT(partitioned_hypergraph_m.checkConnectivitySet(he, context.partition.k));
            for ( const PartitionID& p : partitioned_hypergraph_m.connectivitySet(he) ) {
              block_connectivities[p] = std::make_tuple(std::get<0>(block_connectivities[p]) + 1, p);
            }
          }

          // sort block_connectivities in descending order
          std::sort(block_connectivities.begin(), block_connectivities.end(), [](const auto& a, const auto& b) {
            return std::get<0>(a) > std::get<0>(b);
          });

          //Add node to block with the highest connectivity if it doesn't violate max_part_weights (imbalance)
          for (const auto& block_connectivity : block_connectivities) {
            // TODO Why not stay if current is best?
            if (std::get<1>(block_connectivity) != current_part)
            {
              if (partitioned_hypergraph_m.partWeight(std::get<1>(block_connectivity)) + hypergraph_m.nodeWeight(hn) <=
                  context.partition.max_part_weights[std::get<1>(block_connectivity)]) {
                partitioned_hypergraph_m.changeNodePart(hn, current_part,  std::get<1>(block_connectivity));
                return std::get<1>(block_connectivity);
              }
            }
          }
          // if no other partition could accommodate the node, stay in current partition
          return current_part;
        }

        PartitionID getTargetPart(const HypernodeID& hn, const size_t stage) {
          const PartitionID current_part = partitioned_hypergraph_m.partID(hn);
          ASSERT(current_part != kInvalidPartition);
          int max_gain = 0;
          PartitionID target_part = current_part;
          for ( PartitionID p = 0; p < context.partition.k; ++p ) {
            if (p != current_part) {
              if ((stage == 1 && p > current_part) || (stage == 2 && p < current_part) || (stage == 0)) {
                int gain = 0;
                for ( const HyperedgeID& he : hypergraph_m.incidentEdges(hn) ) {
                  if (partitioned_hypergraph_m.pinCountInPart(he, p) == 0 && partitioned_hypergraph_m.pinCountInPart(he, current_part) > 1) {
                    gain -= hypergraph_m.edgeWeight(he);
                  }
                  if (partitioned_hypergraph_m.pinCountInPart(he, p) > 0 && partitioned_hypergraph_m.pinCountInPart(he, current_part) == 1) {
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
        prior_total_weight = hypergraph_m.totalWeight();
        changed_weight = 0;
        return partitioned_hypergraph_m;
      }

      void partition(Change& change, size_t changes_size) override {
        (void) changes_size;
        change_count++;
        ASSERT(metrics::isBalanced(partitioned_hypergraph_m, context));
        // ASSERT(partitioned_hypergraph_m.checkAllConnectivitySets());

        // Removals are not supported by hermes

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
          // Move the node to the best partition
          // PartitionID target_part = getTargetPart(node, 0);
          // if (target_part != partitioned_hypergraph_m.partID(node))
          // {
          //   context.dynamic.move_count++;
          //   partitioned_hypergraph_m.changeNodePart(node, partitioned_hypergraph_m.partID(node), target_part);
          // }

        }

        for (const HypernodeID& hn : change.added_nodes)
        {
          move_node_to_best_partition(hn);
        }

        // only trigger the "lightweight" refiner every batch_size changes since e.g. orkut would take > 24h otherwise
        changed_weight += change.added_nodes.size() + change.removed_nodes.size();
        context.dynamic.simulate_opt_vcycle = true;
        if (changed_weight > context.dynamic.vcycle_step_size_pct * prior_total_weight || (context.dynamic.simulate_opt_vcycle && change_count == changes_size))
        {
          prior_total_weight += changed_weight;
          changed_weight = 0;
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
      }
    };
}

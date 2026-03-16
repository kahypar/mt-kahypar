#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>

namespace mt_kahypar::dyn {

    class Leopard : public DynamicStrategy {

    private:

        vec<size_t> skippedComps; // number of consecutive times a node was not examined for reassignment
        double skipping_threshold = 0.1; // threshold for examining a node for reassignment as experimentally determined in the leopard paper

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


        PartitionID getTargetPartFennel(const HypernodeID& hn) {
          double fennel_alpha = sqrt(sqrt(context.partition.k) * hypergraph_m.initialNumEdges() / pow(hypergraph_m.initialNumNodes(),1.5));
          const PartitionID current_part = partitioned_hypergraph_m.partID(hn);
          ASSERT(current_part != kInvalidPartition);
          std::vector block_connectivities(context.partition.k, std::make_tuple(0,0));
          for ( PartitionID p = 0; p < context.partition.k; ++p ) {
            block_connectivities[p] = std::make_tuple(0, p);
          }
          // compute neighbor count in each block
          for ( const HyperedgeID& he : hypergraph_m.incidentEdges(hn) ) {
            for ( const HypernodeID& hn2 : hypergraph_m.pins(he) ) {
              block_connectivities[partitioned_hypergraph_m.partID(hn2)] = std::make_tuple(std::get<0>(block_connectivities[partitioned_hypergraph_m.partID(hn2)]) + 1, partitioned_hypergraph_m.partID(hn2));
            }
          }
          std::vector<std::tuple<double, PartitionID>> block_scores(context.partition.k, std::make_tuple(0,0));
          for ( PartitionID p = 0; p < context.partition.k; ++p )
          {
            constexpr double fennel_gamma = 1.5;
            double score = std::get<0>(block_connectivities[p]) - fennel_alpha * (fennel_gamma / 2) * pow(partitioned_hypergraph_m.partWeight(p), fennel_gamma - 1);
            block_scores[p] = std::make_tuple(score, p);
          }

          // sort block_scores in descending order
          std::sort(block_scores.begin(), block_scores.end(), [](const auto& a, const auto& b) {
            return std::get<0>(a) > std::get<0>(b);
          });

          // return the block with the highest score that doesn't violate max_part_weights (imbalance)
          for (const auto& block_score : block_scores)
          {
            if (std::get<1>(block_score) != current_part)
            {
              if (partitioned_hypergraph_m.partWeight(std::get<1>(block_score)) + hypergraph_m.nodeWeight(hn) <=
                  context.partition.max_part_weights[std::get<1>(block_score)])
              {
                return std::get<1>(block_score);
              }
            }
          }
          return current_part;
        }

        bool toExamineOrNot(const HypernodeID& hn, size_t threshold) {
          size_t neighbors = 0;
          for (const HyperedgeID& he : hypergraph_m.incidentEdges(hn)) {
            for (const HypernodeID& hn2 : hypergraph_m.pins(he)) {
              if (hn2 != hn) {
                neighbors++;
              }
            }
          }
          if (neighbors == 0) {
            return false;
          }
          if ((skippedComps[hn] + 1) / neighbors >= threshold) {
            skippedComps[hn] = 0;
            return true;
          } else {
            skippedComps[hn]++;
            return false;
          }
        }

    public:

      Leopard(ds::MutableHypergraph& hypergraph_m, Context& context)
          : DynamicStrategy(hypergraph_m, context) {}

      MutablePartitionedHypergraph& init() override {
          partitioned_hypergraph_m = partition_hypergraph_km1(hypergraph_m, context);
        ASSERT(partitioned_hypergraph_m.checkAllConnectivitySets());
        for (HypernodeID hn = 0; hn < hypergraph_m.initialNumNodes(); ++hn) {
          skippedComps.push_back(0);
        }
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
          skippedComps.push_back(0);
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

        std::vector<HypernodeID> nodes_to_examine;
        for (const HyperedgeID& he : change.added_edges)
        {
          for (const HypernodeID& hn : hypergraph_m.pins(he))
          {
            nodes_to_examine.push_back(hn);
          }
        }

        size_t examined_nodes  = 1;

        while (!nodes_to_examine.empty() && examined_nodes < 10)
        {
          const HypernodeID hn = nodes_to_examine.back();
          nodes_to_examine.pop_back();
          examined_nodes++;
          for (size_t stage = 1; stage <= 2; ++stage)
          {
            const PartitionID target_part = getTargetPartFennel(hn);
            if (target_part != partitioned_hypergraph_m.partID(hn)) {
              context.dynamic.move_count++;
              partitioned_hypergraph_m.changeNodePart(hn, partitioned_hypergraph_m.partID(hn), target_part);
              for (const HyperedgeID& he : hypergraph_m.incidentEdges(hn))
              {
                for (const HypernodeID& hn2 : hypergraph_m.pins(he))
                {
                  if (hn2 != hn && toExamineOrNot(hn2, skipping_threshold)) {
                    nodes_to_examine.push_back(hn2);
                  }
                }
              }
            }
          }
        }
      }
    };
}
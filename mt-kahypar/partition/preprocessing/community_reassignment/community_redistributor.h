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

#include "tbb/task_group.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include "mt-kahypar/definitions.h"

namespace mt_kahypar {
namespace preprocessing {

template< typename TypeTraits >
class CommunityRedistributorT {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;

 public:
  CommunityRedistributorT(const CommunityRedistributorT&) = delete;
  CommunityRedistributorT& operator= (const CommunityRedistributorT&) = delete;
  CommunityRedistributorT(CommunityRedistributorT&&) = delete;
  CommunityRedistributorT& operator= (CommunityRedistributorT&&) = delete;

  ~CommunityRedistributorT() = default;

  static HyperGraph redistribute(HyperGraph& hg,
                                 const PartitionID k,
                                 const std::vector<PartitionID>& community_assignment) {
    int used_numa_nodes = TBB::instance().num_used_numa_nodes();

    // Compute Node Mapping
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    std::vector<HypernodeID> node_mapping(hg.initialNumNodes(), -1);
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(0UL, hg.initialNumNodes()),
      [&](const tbb::blocked_range<HypernodeID>& range) {
        for ( HypernodeID hn = range.begin(); hn < range.end(); ++hn ) {
          node_mapping[hn] = community_assignment[hg.communityID(hg.globalNodeID(hn))];
        }
    });
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("compute_node_mapping", "Compute Node Mapping",
      "redistribution", mt_kahypar::utils::Timer::Type::PREPROCESSING, 1, std::chrono::duration<double>(end - start).count());

    // Compute Hyperedge Mapping
    start = std::chrono::high_resolution_clock::now();
    std::vector<parallel::scalable_vector<PartitionID>> hyperedge_mapping(used_numa_nodes);
    for ( int node = 0; node < used_numa_nodes; ++node ) {
      TBB::instance().numa_task_arena(node).execute([&, node] {
        TBB::instance().numa_task_group(node).run([&, node] {
          hyperedge_mapping[node].assign(hg.initialNumEdges(node), -1);
          tbb::parallel_for(tbb::blocked_range<HyperedgeID>(0UL, hg.initialNumEdges(node)),
            [&](const tbb::blocked_range<HyperedgeID>& range) {
              for ( HyperedgeID local_he = range.begin(); local_he < range.end(); ++local_he ) {
                const HyperedgeID he = StreamingHyperGraph::get_global_edge_id(node, local_he);

                // Compute for each hyperedge the pin count on each node for the new assignment
                parallel::scalable_vector<size_t> pin_count(used_numa_nodes, 0);
                for ( const HypernodeID& pin : hg.pins(he) ) {
                  ASSERT(hg.communityID(pin) < (PartitionID) community_assignment.size());
                  ++pin_count[community_assignment[hg.communityID(pin)]];
                }

                // Compute assignment based maximum pin count
                size_t max_pin_count = 0;
                PartitionID assigned_node = -1;
                for ( PartitionID current_node = 0; current_node < used_numa_nodes; ++current_node ) {
                  if ( pin_count[current_node] >= max_pin_count ) {
                    max_pin_count = pin_count[current_node];
                    assigned_node = current_node;
                  }
                }
                ASSERT(assigned_node != -1);
                hyperedge_mapping[node][local_he] = assigned_node;
              }
            });
        });
      });
    }
    TBB::instance().wait();
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("compute_hyperedge_mapping", "Compute Hyperedge Mapping",
      "redistribution", mt_kahypar::utils::Timer::Type::PREPROCESSING, 2, std::chrono::duration<double>(end - start).count());

    // Reset Pins to original node ids
    start = std::chrono::high_resolution_clock::now();
    hg.resetPinsToOriginalNodeIds();
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("reset_pins_to_original_ids", "Reset Pins to original IDs",
      "redistribution", mt_kahypar::utils::Timer::Type::PREPROCESSING, 3, std::chrono::duration<double>(end - start).count());

    start = std::chrono::high_resolution_clock::now();
    // Initialize Streaming Hypergraphs
    std::vector<StreamingHyperGraph> numa_hypergraphs;
    tbb::task_group group;
    for ( int node = 0; node < used_numa_nodes; ++node ) {
      TBB::instance().numa_task_arena(node).execute([&] {
        group.run([&] {
          numa_hypergraphs.emplace_back(node, k);
        });
      });
      TBB::instance().wait(node, group);
    }

    // Stream hyperedges into hypergraphs
    for ( int node = 0; node < used_numa_nodes; ++node ) {
      for ( int streaming_node = 0; streaming_node < used_numa_nodes; ++streaming_node ) {
        TBB::instance().numa_task_arena(streaming_node).execute([&, node, streaming_node] {
          TBB::instance().numa_task_group(node).run([&, node, streaming_node] {
            tbb::parallel_for(tbb::blocked_range<HyperedgeID>(0UL, hg.initialNumEdges(node)),
            [&, node, streaming_node](const tbb::blocked_range<HyperedgeID>& range) {
              for ( HyperedgeID local_he = range.begin(); local_he < range.end(); ++local_he ) {
                ASSERT(streaming_node == HwTopology::instance().numa_node_of_cpu(sched_getcpu()));
                if ( hyperedge_mapping[node][local_he] == streaming_node ) {
                  const HyperedgeID he = StreamingHyperGraph::get_global_edge_id(node, local_he);
                  parallel::scalable_vector<HypernodeID> hyperedge;
                  for ( const HypernodeID& pin : hg.pins(he) ) {
                    hyperedge.emplace_back(pin);
                  }
                  numa_hypergraphs[streaming_node].streamHyperedge(hyperedge, hg.originalEdgeID(he), hg.edgeWeight(he));
                }
              }
            });
          });
        });
      }
      TBB::instance().wait();
    }

    // Initialize hyperedges in numa hypergraphs
    for ( int node = 0; node < used_numa_nodes; ++node ) {
      TBB::instance().numa_task_arena(node).execute([&] {
        TBB::instance().numa_task_group(node).run([&, node] {
          numa_hypergraphs[node].initializeHyperedges(hg.initialNumNodes());
        });
      });
    }
    TBB::instance().wait();
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("stream_hyperedges", "Stream Hyperedges",
      "redistribution", mt_kahypar::utils::Timer::Type::PREPROCESSING, 4, std::chrono::duration<double>(end - start).count());

    // Initialize hypergraph
    start = std::chrono::high_resolution_clock::now();
    HyperGraph hypergraph(hg.initialNumNodes(), std::move(numa_hypergraphs), std::move(node_mapping), k);
    ASSERT(hypergraph.initialNumNodes() == hg.initialNumNodes());
    ASSERT(hypergraph.initialNumEdges() == hg.initialNumEdges());
    ASSERT(hypergraph.initialNumPins() == hg.initialNumPins());
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("initialize hypergraph", "Initialize Hypergraph",
      "redistribution", mt_kahypar::utils::Timer::Type::PREPROCESSING, 5, std::chrono::duration<double>(end - start).count());

    // Initialize Communities
    start = std::chrono::high_resolution_clock::now();
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(0UL, hypergraph.initialNumNodes()),
      [&](const tbb::blocked_range<HypernodeID>& range) {
      for ( HypernodeID hn = range.begin(); hn < range.end(); ++hn ) {
        HypernodeID old_global_id = hg.globalNodeID(hn);
        HypernodeID new_global_id = hypergraph.globalNodeID(hn);
        hypergraph.setCommunityID(new_global_id, hg.communityID(old_global_id));
      }
    });
    hypergraph.initializeCommunities();
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("initialize_communities", "Initialize Communities",
      "redistribution", mt_kahypar::utils::Timer::Type::PREPROCESSING, 6, std::chrono::duration<double>(end - start).count());

    HEAVY_PREPROCESSING_ASSERT([&] {
      for ( const HypernodeID& hn : hypergraph.nodes() ) {
        int node = StreamingHyperGraph::get_numa_node_of_vertex(hn);
        PartitionID community_id = hypergraph.communityID(hn);
        if ( community_assignment[community_id] != node ) {
          LOG << "Hypernode" << hn << "should be on numa node" << community_assignment[community_id]
              << "but is on node" << node;
          return false;
        }
      }
      return true;
    }(), "There are verticies assigned to wrong numa node");

    return hypergraph;
  }

 protected:
  CommunityRedistributorT() = default;
};

using CommunityRedistributor = CommunityRedistributorT<GlobalTypeTraits>;

} // namespace preprocessing
} // namespace mt_kahypar
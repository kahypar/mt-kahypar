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

#include <algorithm>
#include <type_traits>
#include <chrono>

#include "kahypar/macros.h"
#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/streaming_hypergraph.h"

namespace mt_kahypar {
namespace ds {

template <typename HypernodeType_ = Mandatory,
          typename HyperedgeType_ = Mandatory,
          typename HypernodeWeightType_ = Mandatory,
          typename HyperedgeWeightType_ = Mandatory,
          typename PartitionIDType_ = Mandatory,
          typename HardwareTopology = Mandatory,
          typename TBBNumaArena = Mandatory>
class Hypergraph {

  static constexpr bool debug = false;

  using HypernodeID = HypernodeType_;
  using HyperedgeID = HyperedgeType_;
  using HypernodeWeight = HypernodeWeightType_;
  using HyperedgeWeight = HyperedgeWeightType_;
  using PartitionID = PartitionIDType_;

  using StreamingHypergraph = mt_kahypar::ds::StreamingHypergraph<HypernodeID, 
                                                                  HyperedgeID, 
                                                                  HypernodeWeight, 
                                                                  HyperedgeWeight, 
                                                                  PartitionID,
                                                                  HardwareTopology,
                                                                  TBBNumaArena>;

 public:

  explicit Hypergraph() :
    _num_hypernodes(0),
    _hypergraphs(),
    _node_mapping() { }

  explicit Hypergraph(const HypernodeID num_hypernodes,
                      std::vector<StreamingHypergraph>&& hypergraphs) :
    _num_hypernodes(num_hypernodes),
    _hypergraphs(std::move(hypergraphs)),
    _node_mapping(num_hypernodes, 0) { 
    initializeHypernodes();
  }

  Hypergraph(const Hypergraph&) = delete;
  Hypergraph& operator= (const Hypergraph&) = delete;

  Hypergraph(Hypergraph&& other) = default;
  Hypergraph& operator= (Hypergraph&&) = default;

  ~Hypergraph() = default;


 private:

  void initializeHypernodes() {
    size_t num_streaming_hypergraphs = _hypergraphs.size();
    // Computes mapping for each node to a streaming hypergraph
    // A node is assigned to the streaming hypergraph where it occurs
    // most as pin.
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(0UL, _num_hypernodes), 
      [&](const tbb::blocked_range<HypernodeID>& range) {
      for ( HypernodeID hn = range.begin(); hn < range.end(); ++hn ) {
        size_t max_pins = 0;
        HypernodeID max_node_id = 0;
        for ( HypernodeID node = 1; node < num_streaming_hypergraphs; ++node ) {
          size_t num_pins = _hypergraphs[node].vertexPinCount(hn);
          if ( num_pins > max_pins ) {
            max_pins = num_pins;
            max_node_id = node;
          }
        }
        _node_mapping[hn] = max_node_id;
      }
    });

    // Stream hypernodes into corresponding streaming hypergraph, where it
    // is assigned to
    tbb::task_group group;
    for ( HypernodeID node = 0; node < num_streaming_hypergraphs; ++node ) {
      TBBNumaArena::instance().numa_task_arena(node).execute([&] {
        group.run([&, node] {
          tbb::parallel_for(tbb::blocked_range<HypernodeID>(0UL, _num_hypernodes), 
            [&](const tbb::blocked_range<HypernodeID>& range) {
            for ( HypernodeID hn = range.begin(); hn < range.end(); ++hn ) {
              if ( _node_mapping[hn] == node ) {
                _node_mapping[hn] = _hypergraphs[node].streamHypernode(hn, 1);
              }
            }
          });
        });
      });
    }
    group.wait();

    // Initialize hypernodes on each streaming hypergraph
    // NOTE, that also involves streaming local incident nets to other
    // streaming hypergraphs
    for ( size_t node = 0; node < num_streaming_hypergraphs; ++node ) {
      group.run([&, node] {
        TBBNumaArena::instance().numa_task_arena(node).execute([&] {
          _hypergraphs[node].initializeHypernodes(_hypergraphs, _node_mapping);
        });
      });
    }
    group.wait();

    // Initialize incident nets of hypernodes
    for ( size_t node = 0; node < num_streaming_hypergraphs; ++node ) {
      group.run([&, node] {
        TBBNumaArena::instance().numa_task_arena(node).execute([&] {
          _hypergraphs[node].initializeIncidentNets();
        });
      });
    }
    group.wait();

    ASSERT([&] {
      // Internally verify that incident nets are constructed correctly
      for ( size_t node = 0; node < num_streaming_hypergraphs; ++node ) {
        if ( !_hypergraphs[node].verify_incident_nets_of_hypergraph(_hypergraphs) ) {
          return false;
        }
      }
      return true;
    }(), "Initialization of incident nets failed");
  }

  HypernodeID _num_hypernodes; 

  std::vector<StreamingHypergraph> _hypergraphs;
  std::vector<HypernodeID> _node_mapping;

 
};

} // namespace ds
} // namespace mt_kahypar
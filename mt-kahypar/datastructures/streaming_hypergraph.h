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
#include <thread>
#include <type_traits>

#include "tbb/task_scheduler_observer.h"
#include "tbb/task_arena.h"
#include "tbb/task_group.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_sort.h"
#include "tbb/blocked_range.h"

#include "kahypar/macros.h"
#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/datastructures/streaming_map.h"


namespace mt_kahypar {
namespace ds {

template <typename HypernodeType_ = Mandatory,
          typename HyperedgeType_ = Mandatory,
          typename HypernodeWeightType_ = Mandatory,
          typename HyperedgeWeightType_ = Mandatory,
          typename PartitionIDType_ = Mandatory,
          typename HardwareTopology = Mandatory,
          typename TBBNumaArena = Mandatory>
class StreamingHypergraph {

  static constexpr bool debug = false;
  static constexpr size_t NUMA_NODE_INDENTIFIER = 48;

  using HypernodeID = HypernodeType_;
  using HyperedgeID = HyperedgeType_;
  using HypernodeWeight = HypernodeWeightType_;
  using HyperedgeWeight = HyperedgeWeightType_;
  using PartitionID = PartitionIDType_;

  using Self = StreamingHypergraph<HypernodeID, HyperedgeID, HypernodeWeight,
                                   HyperedgeWeight, PartitionID, HardwareTopology,
                                   TBBNumaArena>;

  using IncidentNets = std::vector<std::vector<HyperedgeID>>;

  static_assert( sizeof(HypernodeID) == 8 );
  static_assert( std::is_unsigned<HypernodeID>::value );
  static_assert( sizeof(HyperedgeID) == 8 );
  static_assert( std::is_unsigned<HyperedgeID>::value );

  class Hypernode {

    public:
      Hypernode() :
        _id(0),
        _original_id(0),
        _weight(1),
        _valid(false) { }

      Hypernode(const HypernodeID id,
                const HypernodeID original_id,
                const HypernodeWeight weight) :
        _id(id),
        _original_id(original_id),
        _weight(weight),
        _valid(false) { }

      HypernodeID nodeId() const {
        return _id;
      }

    private:
      // ! Hypernode id
      HypernodeID _id;
      // ! Original hypernode id
      HypernodeID _original_id;
      // ! Hypernode weight
      HyperedgeWeight _weight;
      // ! Flag indicating whether or not the element is active.
      bool _valid;
  };

  static_assert( std::is_trivially_copyable<Hypernode>::value );

  class Hyperedge {

    public:
      Hyperedge() :
        _begin(0),
        _size(0),
        _weight(1),
        _valid(false) { }

      Hyperedge(const size_t begin, 
                const size_t size,
                const HyperedgeWeight weight) :
        _begin(begin),
        _size(size),
        _weight(weight),
        _valid(true) { }

      // ! Disables the hypernode/hyperedge. Disable hypernodes/hyperedges will be skipped
      // ! when iterating over the set of all nodes/edges.
      void disable() {
        ASSERT(!isDisabled());
        _valid = false;
      }

      bool isDisabled() const {
        return _valid == false;
      }

      void enable() {
        ASSERT(isDisabled());
        _valid = true;
      }

      // ! Returns the index of the first element in _incidence_array
      size_t firstEntry() const {
        return _begin;
      }

      // ! Sets the index of the first element in _incidence_array to begin
      void setFirstEntry(size_t begin) {
        ASSERT(!isDisabled());
        _begin = begin;
      }

      // ! Returns the index of the first element in _incidence_array
      size_t firstInvalidEntry() const {
        return _begin + _size;
      }

      size_t size() const {
        ASSERT(!isDisabled());
        return _size;
      }

      void setSize(size_t size) {
        ASSERT(!isDisabled());
        _size = size;
      }

      void incrementSize() {
        ASSERT(!isDisabled());
        ++_size;
      }

      void decrementSize() {
        ASSERT(!isDisabled());
        ASSERT(_size > 0);
        --_size;
      }

      HyperedgeWeight weight() const {
        ASSERT(!isDisabled());
        return _weight;
      }

      void setWeight(HyperedgeWeight weight) {
        ASSERT(!isDisabled());
        _weight = weight;
      }

      bool operator== (const Hyperedge& rhs) const {
        return _begin == rhs._begin && _size == rhs._size && _weight == rhs._weight;
      }

      bool operator!= (const Hyperedge& rhs) const {
        return !operator== (this, rhs);
      }

    private:
      // ! Index of the first element in _incidence_array
      size_t _begin;
      // ! Number of _incidence_array elements
      size_t _size;
      // ! yperedge weight
      HyperedgeWeight _weight;
      // ! Flag indicating whether or not the element is active.
      bool _valid;
  };

  using CPUHypernodeBuffer = std::vector<std::vector<Hypernode>>;

 public:
  explicit StreamingHypergraph(const int node) :
    _node(node),
    _arena(TBBNumaArena::instance().numa_task_arena(node)),
    _hypernodes(),
    _incident_nets(),
    _hyperedges(),
    _incidence_array(),
    _vertex_pin_count(),
    _pin_stream(),
    _hyperedge_stream(),
    _next_node_id(0),
    _hypernode_stream(),
    _incident_net_stream() { 
    // Make sure constructor is called on corresponding numa node
    ASSERT(HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node, 
      "Only allowed to allocate numa hypergraph on node" << _node << ", but it is"
      << HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()));
  }

  StreamingHypergraph(const StreamingHypergraph&) = delete;
  StreamingHypergraph& operator= (const StreamingHypergraph&) = delete;

  StreamingHypergraph(StreamingHypergraph&& other) :
    _node(other._node),
    _arena(other._arena),
    _hypernodes(std::move(other._hypernodes)),
    _incident_nets(std::move(other._incident_nets)),
    _hyperedges(std::move(other._hyperedges)),
    _incidence_array(std::move(other._incidence_array)),
    _vertex_pin_count(std::move(other._vertex_pin_count)),
    _pin_stream(std::move(other._pin_stream)),
    _hyperedge_stream(std::move(other._hyperedge_stream)),
    _next_node_id(other._next_node_id.load()),
    _hypernode_stream(std::move(other._hypernode_stream)),
    _incident_net_stream(std::move(other._incident_net_stream)) { }

  StreamingHypergraph& operator= (StreamingHypergraph&&) = default;

  ~StreamingHypergraph() = default;

  size_t vertexPinCount(const HypernodeID hn) const {
    ASSERT(hn < _vertex_pin_count.size());
    return _vertex_pin_count[hn];
  }

  HypernodeID streamHypernode(HypernodeID original_id, HypernodeWeight weight) {
    // Make sure calling process is part of correct numa node
    ASSERT(HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node);
    HypernodeID node_id = get_global_node_id();
    _hypernode_stream.stream(node_id, original_id, weight);
    return node_id;
  }

  void streamHyperedge(const std::vector<HypernodeID>& hyperedge, const HyperedgeWeight& weight) {
    int cpu_id = sched_getcpu();
    // Make sure calling process is part of correct numa node
    ASSERT(HardwareTopology::instance().numa_node_of_cpu(cpu_id) == _node);
    _hyperedge_stream.stream(_pin_stream.size(cpu_id), hyperedge.size(), weight);
    for ( const HypernodeID& pin : hyperedge  ) {
      _pin_stream.stream(pin);
    }
  }

  void streamIncidentNet(const HypernodeID hn, const HyperedgeID he) {
    ASSERT(get_numa_node_of_vertex(hn) == _node);
    _incident_net_stream.stream(hn, he);
  }

  void initializeHyperedges(const HypernodeID num_hypernodes) {
    // Make sure calling process is part of correct numa node
    ASSERT(HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node,
      "Expected that assigned cpu is on numa node" << _node << ", but was on node"
        << HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()));

    // Copy streamed data into global vectors
    _incidence_array = _pin_stream.copy(_arena);
    _hyperedges = _hyperedge_stream.copy(_arena);

    // Update start position of each hyperedge to correct one in global incidence array
    // Note, start positions are stored relative to the local buffer there are streamed into.
    // However, memcpy does hyperedges invalidates those local positions.
    tbb::task_group group;
    _arena.execute([&] {
      for ( size_t cpu_id = 0; cpu_id < std::thread::hardware_concurrency(); ++cpu_id ) {
        group.run([&, cpu_id] {
          size_t start = _hyperedge_stream.prefix_sum(cpu_id);
          size_t end = start + _hyperedge_stream.size(cpu_id);
          size_t delta = _pin_stream.prefix_sum(cpu_id);
          ASSERT(end <= _hyperedges.size());
          for ( size_t pos = start; pos < end; ++pos ) {
            _hyperedges[pos].setFirstEntry( _hyperedges[pos].firstEntry() + delta );
          }
        });
      }
    });
    TBBNumaArena::instance().wait(_node, group);

    ASSERT([&] {
      for ( size_t cpu_id = 0; cpu_id < std::thread::hardware_concurrency(); ++cpu_id ) {
        size_t start = _hyperedge_stream.prefix_sum(cpu_id);
        for ( size_t idx = 0; idx < _hyperedge_stream.size(cpu_id); ++idx ) {
          const Hyperedge& stream_he = _hyperedge_stream.value(cpu_id, idx);
          const Hyperedge& he = _hyperedges[idx + start];
          if ( stream_he.size() != he.size() ) {
            LOG << "Unequal size";
            return false;
          }
          size_t global_pos = he.firstEntry();
          for ( size_t pos = stream_he.firstEntry(); pos < stream_he.firstEntry() + stream_he.size(); ++pos ) {
            if ( _incidence_array[global_pos] != _pin_stream.value(cpu_id, pos) ) {
              LOG << "Pins in stream and global incidence array are not equal"
                  << V(_incidence_array[global_pos]) << V(_pin_stream.value(cpu_id, pos));
              return false;
            }
            ++global_pos;
          }
        }
      }
      return true;
    }(), "Failed to copy buffer to hypergraph");

    // Compute how many times a hypernode occurs on this node
    // as pin. Will be later important to compute node assignment.
    // TODO(heuer): Think how to parallelize this
    _vertex_pin_count.resize(num_hypernodes);
    for ( size_t i = 0; i < _incidence_array.size(); ++i ) {
      const HypernodeID& pin = _incidence_array[i];
      ASSERT(pin < _vertex_pin_count.size());
      _vertex_pin_count[pin]++;
    }

    _pin_stream.clear();
    _hyperedge_stream.clear();
  }

  void initializeHypernodes(std::vector<Self>& hypergraphs,
                            const std::vector<HypernodeID>& node_mapping) {
    // Make sure calling process is part of correct numa node
    ASSERT(HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node,
      "Expected that assigned cpu is on numa node" << _node << ", but was on node"
        << HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()));

    _hypernodes = _hypernode_stream.copy(_arena);

    // Sort hypernodes in increasing order of their node id and ...
    tbb::task_group group;
    _arena.execute([&] {
      group.run([&] {
        tbb::parallel_sort(_hypernodes.begin(), _hypernodes.end(), 
          [&](const Hypernode& lhs, const Hypernode& rhs) {
          return lhs.nodeId() < rhs.nodeId();
        });
      });
    });

    // ... parallel to that remap node ids in incidence array to new node ids
    _arena.execute([&] {
      group.run([&] {
        tbb::parallel_for(tbb::blocked_range<size_t>(0UL, _incidence_array.size()),
          [&](const tbb::blocked_range<size_t>& range) {
          for ( size_t pos = range.begin(); pos < range.end(); ++pos ) {
            HypernodeID pin = _incidence_array[pos];
            ASSERT(pin < node_mapping.size());
            _incidence_array[pos] = node_mapping[pin];
          }
        });
      });
    });
    TBBNumaArena::instance().wait(_node, group);

    ASSERT([&]{
      for ( size_t i = 0; i < _hypernodes.size(); ++i ) {
        const Hypernode& hn = _hypernodes[i];
        const int node = get_numa_node_of_vertex(hn.nodeId()); 
        const HypernodeID local_id = get_local_node_id_of_vertex(hn.nodeId());
        if ( node != _node ) {
          LOG << "Hypernode" << hn.nodeId() << "should be on numa node" << node
              << ", but is on numa node" << _node;
          return false;
        } else if ( i != local_id ) {
          LOG << "Hypernode" << local_id << "should have local id" << i;
          return false;
        }

        if ( i > 0 ) {
          const Hypernode& previous_hn = _hypernodes[i-1];
          const HypernodeID previous_local_id = get_local_node_id_of_vertex(previous_hn.nodeId());
          if ( local_id <= previous_local_id ) {
            LOG << "Local node ids not sorted" << V(previous_local_id) << V(local_id);
            return false;
          } else if ( previous_local_id + 1 != local_id ) {
            LOG << "Local node ids not consecutive" << V(previous_local_id) << V(local_id);
            return false;
          }
        }
      }
      return true;
    }(), "Initialization of hypernodes failed");

    // Stream incident nets
    _arena.execute([&] {
      group.run([&] {
        tbb::parallel_for(tbb::blocked_range<size_t>(0UL, _hyperedges.size()),
          [&](const tbb::blocked_range<size_t>& range) {
          for ( size_t pos = range.begin(); pos < range.end(); ++pos ) {
            const Hyperedge& he = _hyperedges[pos];
            const HyperedgeID he_id = get_global_edge_id(pos);
            for ( size_t incidence_array_pos = he.firstEntry();
                  incidence_array_pos < he.firstEntry() + he.size();
                  ++incidence_array_pos ) {
              const HypernodeID hn = _incidence_array[incidence_array_pos];
              const int node = get_numa_node_of_vertex(hn);
              ASSERT(node < (int) hypergraphs.size());
              hypergraphs[node].streamIncidentNet(hn, he_id);
            }
          }
        });
      });
    });
    TBBNumaArena::instance().wait(_node, group);

    _hypernode_stream.clear();
  }

  void initializeIncidentNets() {
    // Make sure calling process is part of correct numa node
    ASSERT(HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node,
      "Expected that assigned cpu is on numa node" << _node << ", but was on node"
        << HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()));

    _incident_nets.resize(_hypernodes.size());
    _incident_net_stream.copy(_arena, _incident_nets, [&](const HypernodeID& u) {
      ASSERT(get_numa_node_of_vertex(u) == _node);
      return get_local_node_id_of_vertex(u);
    });

    _incident_net_stream.clear();
  }

  bool verify_incident_nets_of_hypergraph(const std::vector<Self>& hypergraphs) const {
    for ( size_t pos = 0; pos < _incident_nets.size(); ++pos ) {
      const HypernodeID& hn = _hypernodes[pos].nodeId();
      for ( const HyperedgeID& he : _incident_nets[pos] ) {
        const int node = get_numa_node_of_hyperedge(he);
        const HyperedgeID local_edge_id = get_local_edge_id_of_hyperedge(he);
        ASSERT(node < (int) hypergraphs.size());
        ASSERT(local_edge_id < hypergraphs[node]._hyperedges.size());
        const Hyperedge& e = hypergraphs[node]._hyperedges[local_edge_id];
        const auto first = hypergraphs[node]._incidence_array.begin() + e.firstEntry();
        const auto last = first + e.size();
        if ( std::find(first, last, hn) == last ) {
          LOG << "Hypernode" << hn << "not part of hyperedge" << he << "on numa node" << node;
          return false;
        }
      }
    }
    return true;
  }

 private:

  HypernodeID get_global_node_id() {
    HypernodeID local_node_id = _next_node_id++;
    return ( ( (HypernodeID) _node ) << NUMA_NODE_INDENTIFIER ) | local_node_id;
  }

  HyperedgeID get_global_edge_id(const size_t edge_pos) {
    return ( ( (HyperedgeID) _node ) << NUMA_NODE_INDENTIFIER ) | edge_pos;
  }

  static int get_numa_node_of_vertex(const HypernodeID u) {
    return (int) (u >> NUMA_NODE_INDENTIFIER); 
  }

  static HypernodeID get_local_node_id_of_vertex(const HypernodeID u) {
    return ( (1UL << NUMA_NODE_INDENTIFIER) - 1 ) & u; 
  }

  static int get_numa_node_of_hyperedge(const HyperedgeID e) {
    return (int) (e >> NUMA_NODE_INDENTIFIER); 
  }

  static HypernodeID get_local_edge_id_of_hyperedge(const HyperedgeID e) {
    return ( (1UL << NUMA_NODE_INDENTIFIER) - 1 ) & e; 
  }

  const int _node;
  // ! task arena for numa node
  tbb::task_arena& _arena;

  // ! Hypernodes
  std::vector<Hypernode> _hypernodes;
  // ! Incident nets of hypernodes
  IncidentNets _incident_nets;
  // ! Hyperedges
  std::vector<Hyperedge> _hyperedges; 
  // ! Pins of hyperedges
  std::vector<HypernodeID> _incidence_array;

  // ! Contains for each hypernode, how many time it
  // ! occurs as pin in incidence array
  std::vector<size_t> _vertex_pin_count; 
  
  // ! Stream for pins of hyperedges
  StreamingVector<HypernodeID> _pin_stream;
  // ! Stream for hyperedges
  StreamingVector<Hyperedge> _hyperedge_stream;

  // ! Atomic counter for unique local node ids
  std::atomic<HypernodeID> _next_node_id;
  // ! Stream for hypernodes
  StreamingVector<Hypernode> _hypernode_stream;

  // ! Stream for incident nets
  StreamingMap<HypernodeID, HyperedgeID> _incident_net_stream;
};

} // namespace ds
} // namespace mt_kahypar
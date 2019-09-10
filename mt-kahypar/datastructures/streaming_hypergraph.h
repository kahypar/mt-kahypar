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
#include "tbb/blocked_range.h"

#include "kahypar/macros.h"
#include "kahypar/meta/mandatory.h"

namespace kahypar {
namespace ds {

template <typename HypernodeType_ = Mandatory,
          typename HyperedgeType_ = Mandatory,
          typename HypernodeWeightType_ = Mandatory,
          typename HyperedgeWeightType_ = Mandatory,
          typename PartitionIDType_ = Mandatory,
          typename HardwareTopology = Mandatory,
          typename TBBNumaArena = Mandatory>
class StreamingHypergraph {

  static constexpr bool debug = true;

  using HypernodeID = HypernodeType_;
  using HyperedgeID = HyperedgeType_;
  using HypernodeWeight = HypernodeWeightType_;
  using HyperedgeWeight = HyperedgeWeightType_;
  using PartitionID = PartitionIDType_;

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

  static_assert( std::is_trivially_copyable<Hyperedge>::value );

  using HyperedgeStream = std::vector<HypernodeID>;
  using CPUHyperedgeBuffer = std::vector<HyperedgeStream>;
  using CPUHyperedgeIndexBuffer = std::vector<std::vector<Hyperedge>>;

 public:
  explicit StreamingHypergraph(const int node) :
    _node(node),
    _arena(TBBNumaArena::instance().numa_task_arena(node)),
    _hyperedges(),
    _incidence_array(),
    _vertex_pin_count(),
    _buffer(std::thread::hardware_concurrency()),
    _index_buffer(std::thread::hardware_concurrency()) { 
    // Make sure constructor is called on corresponding numa node
    ASSERT(HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node, 
      "Only allowed to allocate numa hypergraph on node" << _node << ", but it is"
      << HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()));
  }

  StreamingHypergraph(const StreamingHypergraph&) = delete;
  StreamingHypergraph& operator= (const StreamingHypergraph&) = delete;

  StreamingHypergraph(StreamingHypergraph&& other) = default;
  StreamingHypergraph& operator= (StreamingHypergraph&&) = default;

  ~StreamingHypergraph() = default;

  void stream(const HyperedgeStream& hyperedge, const HyperedgeWeight& weight) {
    int cpu_id = sched_getcpu();
    // Make sure calling process is part of correct numa node
    ASSERT(HardwareTopology::instance().numa_node_of_cpu(cpu_id) == _node);
    _index_buffer[cpu_id].emplace_back(_buffer[cpu_id].size(), hyperedge.size(), weight);
    _buffer[cpu_id].insert(_buffer[cpu_id].end(), hyperedge.begin(), hyperedge.end());
  }

  void initialize(const HypernodeID num_hypernodes) {
    // Make sure calling process is part of correct numa node
    ASSERT(HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node,
      "Expected that assigned cpu is on numa node" << _node << ", but was on node"
        << HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()));

    size_t num_pins = 0;
    size_t num_hyperedges = 0;
    // Compute positions for copying data of buffers to incidence
    // and hyperedge array
    std::vector<size_t> prefix_pin_sum(_buffer.size());
    std::vector<size_t> prefix_hyperedge_sum(_index_buffer.size());
    for ( size_t i = 0; i < _buffer.size(); ++i ) {
      prefix_pin_sum[i] = num_pins;
      prefix_hyperedge_sum[i] += num_hyperedges;
      num_pins += _buffer[i].size();
      num_hyperedges += _index_buffer[i].size();
    }
    _hyperedges.resize(num_hyperedges);
    _incidence_array.resize(num_pins);

    // Copy buffers to incidence and hyperedge array
    tbb::task_group group;
    _arena.execute([&] {
      for ( size_t cpu_id = 0; cpu_id < std::thread::hardware_concurrency(); ++cpu_id ) {
        group.run([&, cpu_id] {
          copyHyperedgesFromCpuBufferToHypergraph(cpu_id, prefix_hyperedge_sum[cpu_id], prefix_pin_sum[cpu_id]);
        });
      }
    });
    TBBNumaArena::instance().wait(_node, group);

    // Compute how many times a hypernode occurs on this node
    // as pin. Will be later important to compute node assignment.
    // TODO(heuer): Think how to parallelize this
    _vertex_pin_count.resize(num_hypernodes);
    for ( size_t i = 0; i < _incidence_array.size(); ++i ) {
      const HypernodeID& pin = _incidence_array[i];
      ASSERT(pin < _incidence_array.size());
      _vertex_pin_count[pin]++;
    }
  }

 private:

  void copyHyperedgesFromCpuBufferToHypergraph(const int cpu_id,
                                               const size_t hyperedge_pos,
                                               const size_t incidence_array_pos) {
    DBG << "Copy buffer of cpu" << cpu_id << "of size" << _buffer[cpu_id].size()
        << "on numa node" << _node << "(CPU =" << sched_getcpu() 
        << ") to incidence array to position" << incidence_array_pos;
    memcpy(_incidence_array.data() + incidence_array_pos,
            _buffer[cpu_id].data(), _buffer[cpu_id].size() * sizeof(HypernodeID));
    memcpy(_hyperedges.data() + hyperedge_pos,
            _index_buffer[cpu_id].data(), _index_buffer[cpu_id].size() * sizeof(Hyperedge));
    size_t start = hyperedge_pos;
    size_t end = start + _index_buffer[cpu_id].size();
    for ( size_t i = start; i < end; ++i ) {
      _hyperedges[i].setFirstEntry( _hyperedges[i].firstEntry() + incidence_array_pos );
    }

    ASSERT([&] {
      size_t num_buffer_hyperedges = _index_buffer[cpu_id].size();
      size_t global_pos = start;
      for ( size_t i = 0; i < num_buffer_hyperedges; ++i ) {
        const Hyperedge& buffer_he = _index_buffer[cpu_id][i];
        const Hyperedge& he = _hyperedges[global_pos];
        size_t global_j = he.firstEntry();
        for ( size_t j = buffer_he.firstEntry(); j < buffer_he.firstInvalidEntry(); ++j ) {
          if ( _buffer[cpu_id][j] != _incidence_array[global_j] ) {
            LOG << "Expected pin" << _buffer[cpu_id][j] << "on position" << global_j
                << "in incidence array, but was" << _incidence_array[global_j];
            return false;
          }
          ++global_j;
        }
        ++global_pos;
      }
      return true;
    }(), "Failed to copy buffer to hypergraph");

    _buffer[cpu_id].clear();
    _index_buffer[cpu_id].clear();
  }

  const int _node;
  tbb::task_arena& _arena;

  std::vector<Hyperedge> _hyperedges; 
  std::vector<HypernodeID> _incidence_array;
  std::vector<size_t> _vertex_pin_count; 
  
  CPUHyperedgeBuffer _buffer;
  CPUHyperedgeIndexBuffer _index_buffer;
 
};

} // namespace ds
} // namespace kahypar
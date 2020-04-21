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

#include <mutex>
#include <shared_mutex>
#include <memory>
#include <unordered_map>

#include "tbb/parallel_for.h"
#include "tbb/scalable_allocator.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/stl/scalable_unique_ptr.h"
#include "mt-kahypar/utils/memory_tree.h"

namespace mt_kahypar {
namespace parallel {

/*!
 * Singleton that handles huge memory allocations.
 * Memory chunks can be registered with a key and all memory
 * chunks can be collectively allocated in parallel.
 */
class MemoryPool {
  static constexpr bool debug = false;

  // ! Represents a memory chunk.
  struct MemoryChunk {

    explicit MemoryChunk(const size_t num_elements,
                         const size_t size) :
      _chunk_mutex(),
      _num_elements(num_elements),
      _size(size),
      _data(nullptr),
      _is_assigned(false) { }

    MemoryChunk(MemoryChunk&& other) :
      _chunk_mutex(),
      _num_elements(other._num_elements),
      _size(other._size),
      _data(std::move(other._data)),
      _is_assigned(other._is_assigned) {
      other._data = nullptr;
      other._is_assigned = false;
    }

    // ! Requests the memory chunk.
    // ! Note, successive calls to this method will return
    // ! nullptr until release_chunk() is called.
    char* request_chunk() {
      std::lock_guard<std::mutex> lock(_chunk_mutex);
      if ( !_data ) {
        allocate();
      }
      if ( !_is_assigned ) {
        _is_assigned = true;
        return _data;
      } else {
        return nullptr;
      }
    }

    // ! Releases the memory chunks
    void release_chunk() {
      std::lock_guard<std::mutex> lock(_chunk_mutex);
      _is_assigned = false;
    }

    // ! Allocates the memory chunk
    // ! Note, the memory chunk is zero initialized.
    void allocate() {
      if ( !_data ) {
        _data = (char*) scalable_calloc(_num_elements, _size);
      }
    }

    // ! Frees the memory chunk
    void free() {
      if ( _data ) {
        scalable_free(_data);
        _data = nullptr;
      }
    }

    // ! Returns the size in bytes of the memory chunk
    size_t size_in_bytes() const {
      size_t size = 0;
      if ( _data ) {
        size = _num_elements * _size;
      }
      return size;
    }

    std::mutex _chunk_mutex;
    const size_t _num_elements; // Number of elements to allocate
    const size_t _size; // Data type size in bytes
    char* _data;
    bool _is_assigned; // True, if already assigned to a vector
  };


 public:
  MemoryPool(const MemoryPool&) = delete;
  MemoryPool & operator= (const MemoryPool &) = delete;

  MemoryPool(MemoryPool&&) = delete;
  MemoryPool & operator= (MemoryPool &&) = delete;

  ~MemoryPool() {
    free_memory_chunks();
  }

  static MemoryPool& instance() {
    static MemoryPool instance;
    return instance;
  }

  // ! Registers a memory chunk in the memory pool. The memory chunk is
  // ! associated with a memory group and a unique key within that group.
  // ! Note, that the memory chunk is not immediatly allocated. One has to call
  // ! allocate_memory_chunks() to collectively allocate all memory chunks.
  void register_memory_chunk(const std::string& group,
                             const std::string& key,
                             const size_t num_elements,
                             const size_t size) {
    std::unique_lock<std::shared_timed_mutex> lock(_memory_mutex);
    const size_t memory_id = _memory_chunks.size();
    if ( _memory_id_map.find(group) == _memory_id_map.end() ) {
      _memory_id_map.insert(std::make_pair(
        group, std::unordered_map<std::string, size_t>()));
    }
    std::unordered_map<std::string, size_t>& key_to_memory_id = _memory_id_map[group];
    ASSERT(key_to_memory_id.find(key) == key_to_memory_id.end());
    key_to_memory_id.insert(std::make_pair(key, memory_id));
    _memory_chunks.emplace_back(num_elements, size);
  }

  // ! Allocates all registered memory chunks in parallel
  void allocate_memory_chunks() {
    std::shared_lock<std::shared_timed_mutex> lock(_memory_mutex);
    const size_t num_memory_segments = _memory_chunks.size();
    tbb::parallel_for(0UL, num_memory_segments, [&](const size_t i) {
      _memory_chunks[i].allocate();
    });
  }

  // ! Returns the memory chunk registered under the corresponding
  // ! group with the specified key. If the memory chunk is already
  // ! requested, the size of the memory chunk is smaller than the
  // ! requested size or the requested memory chunk does not exist,
  // ! than nullptr is returned.
  char* request_mem_chunk(const std::string& group,
                          const std::string& key,
                          const size_t num_elements,
                          const size_t size) {
    std::shared_lock<std::shared_timed_mutex> lock(_memory_mutex);
    MemoryChunk* chunk = find_memory_chunk(group, key);
    if ( chunk && num_elements * size <= chunk->size_in_bytes() ) {
      return chunk->request_chunk();
    }
    return nullptr;
  }

  // ! Releases the memory chunk under the corresponding group with
  // ! the specified key. Afterwards, memory chunk is available for
  // ! further requests.
  void release_mem_chunk(const std::string& group,
                         const std::string& key) {
    std::shared_lock<std::shared_timed_mutex> lock(_memory_mutex);
    MemoryChunk* chunk = find_memory_chunk(group, key);
    if ( chunk ) {
      chunk->release_chunk();
    }
  }

  // ! Returns the memory chunk under the corresponding group with
  // ! the specified key. In contrast to assign_mem_chunk, no explicit
  // ! checks are performed, if chunk is already assigned.
  char* mem_chunk(const std::string& group,
                  const std::string& key) {
    std::shared_lock<std::shared_timed_mutex> lock(_memory_mutex);
    MemoryChunk* chunk = find_memory_chunk(group, key);
    if ( chunk )   {
      return chunk->_data;
    } else {
      return nullptr;
    }
  }

  // ! Frees all memory chunks in parallel
  void free_memory_chunks() {
    std::unique_lock<std::shared_timed_mutex> lock(_memory_mutex);
    const size_t num_memory_segments = _memory_chunks.size();
    tbb::parallel_for(0UL, num_memory_segments, [&](const size_t i) {
      _memory_chunks[i].free();
    });
    _memory_chunks.clear();
    _memory_id_map.clear();
  }

  // ! Builds a memory tree that reflects the memory
  // ! consumption of the memory pool
  void memory_consumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);
    std::shared_lock<std::shared_timed_mutex> lock(_memory_mutex);
    for ( const auto& group_element : _memory_id_map ) {
      const std::string& group = group_element.first;
      const auto& key_to_memory_id = group_element.second;
      utils::MemoryTreeNode* group_node = parent->addChild(group);
      for ( const auto& element : key_to_memory_id ) {
        const std::string& key = element.first;
        const size_t memory_id = element.second;
        ASSERT(memory_id < _memory_chunks.size());
        group_node->addChild(key, _memory_chunks[memory_id].size_in_bytes());
      }
    }
  }

 private:
  explicit MemoryPool() :
    _memory_mutex(),
    _memory_id_map(),
    _memory_chunks() { }

  // ! Returns a pointer to memory chunk under the corresponding group with
  // ! the specified key.
  MemoryChunk* find_memory_chunk(const std::string& group,
                                 const std::string& key) {

    if ( _memory_id_map.find(group) != _memory_id_map.end() &&
         _memory_id_map[group].find(key) != _memory_id_map[group].end() ) {
      const size_t memory_id = _memory_id_map[group][key];
      ASSERT(memory_id < _memory_chunks.size());
      return &_memory_chunks[memory_id];
    }
    return nullptr;
  }

  // ! Read-Write Lock for memory pool
  mutable std::shared_timed_mutex _memory_mutex;
  // ! Mapping from group-key to a memory chunk id
  // ! The memory chunk id maps points to the memory chunk vector
  std::unordered_map<std::string, std::unordered_map<std::string, size_t>> _memory_id_map;
  // ! Memory chunks
  std::vector<MemoryChunk> _memory_chunks;
};

}  // namespace parallel
}  // namespace mt_kahypar

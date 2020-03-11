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

#include <cmath>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace ds {

/*!
 * Stores for each block of the partition its size and weight.
 * Internally, we store for each block ceil(log2(p)/log2(k)) local block
 * info objects, where p is the number of processors and k is the number
 * of blocks we want to partition the hypergraph into. A block info
 * objects stores the block size and weight as an atomic. To determine
 * the block weight or size, we have to sum up over the local block info
 * objects for the corresponding block.
 * The idea of using log2(p)/log2(k) local block info objects is to prevent
 * congestion when performing atomic updates of block weights and sizes concurrently
 * compared to an variant where we only use one atomic for the weights and sizes.
 * Furthermore, by using only a logarithmic number of local block infos, we
 * still can guarantee * reasonable fast queries to get the global part
 * weights and sizes.
 * If k is small compared to the number of threads, than the workload
 * is scattered onto a logarithmic number of local block infos. If k
 * is equal or even greater than number of threads, we use only one
 * local part info object for each block, which also represents than
 * the global block weight and size.
 */
class PartitionInfo {

  using AtomicWeight = parallel::IntegralAtomicWrapper<HypernodeWeight>;
  using AtomicSize = parallel::IntegralAtomicWrapper<HypernodeID>;

 public:
  /*!
   * For each block \f$V_i\f$ of the \f$k\f$-way partition
   * \f$\mathrm{\Pi} = \{V_1, \dots, V_k\}\f$, a BlockInfo object stores
   * the number of hypernodes currently assigned to block \f$V_i\f$
   * as well as the sum of their weights.
   */
  class BlockInfo {
   public:
    explicit BlockInfo() :
      weight(0),
      size(0) { }

    bool operator== (const BlockInfo& other) const {
      return weight == other.weight && size == other.size;
    }

    AtomicWeight weight;
    AtomicSize size;
  };

 private:
  // For one block, we store several BlockInfo object in order
  // to minimize congestion on the atomics.
  using LocalBlockInfo = parallel::scalable_vector<BlockInfo>;

 public:
  PartitionInfo() :
    _k(0),
    _num_local_block_infos(0),
    _part_info() { }

  PartitionInfo(const size_t num_threads, const PartitionID k) :
    _k(k),
    _num_local_block_infos(computeNumLocalBlockInfos(num_threads, k)),
    _part_info(k, LocalBlockInfo(computeNumLocalBlockInfos(num_threads, k))) {
    ASSERT(_num_local_block_infos > 0);
  }

  PartitionInfo(const PartitionInfo&) = delete;
  PartitionInfo & operator= (const PartitionInfo &) = delete;

  PartitionInfo(PartitionInfo&& other) :
    _k(other._k),
    _num_local_block_infos(other._num_local_block_infos),
    _part_info(std::move(other._part_info)) { }

  PartitionInfo & operator= (PartitionInfo&& other) {
    _k = other._k;
    _num_local_block_infos = other._num_local_block_infos;
    _part_info = std::move(other._part_info);
    return *this;
  }

  HypernodeID numLocalBlockInfos() const {
    return _num_local_block_infos;
  }

  inline void setBlockSize(const PartitionID id,
                           const HypernodeID size) {
    ASSERT(id != kInvalidPartition && id < _k);
    _part_info[id][0].size = size;
    ASSERT(partSize(id) == size);
  }

  inline void setBlockWeight(const PartitionID id,
                             const HypernodeWeight weight) {
    ASSERT(id != kInvalidPartition && id < _k);
    _part_info[id][0].weight = weight;
    ASSERT(partWeight(id) == weight);
  }

  inline void updateNodeWeight(const HypernodeID u,
                               const PartitionID id,
                               const HypernodeWeight delta) {
    ASSERT(id != kInvalidPartition && id < _k);
    const HypernodeID local_bucket = u % _num_local_block_infos;
    _part_info[id][local_bucket].weight += delta;
  }


  inline void setNodePart(const HypernodeID u,
                          const PartitionID id,
                          const HypernodeWeight weight) {
    ASSERT(id != kInvalidPartition && id < _k);
    const HypernodeID local_bucket = u % _num_local_block_infos;
    ++_part_info[id][local_bucket].size;
    _part_info[id][local_bucket].weight += weight;
  }

  inline void changeNodePart(const HypernodeID u,
                             const PartitionID from,
                             const PartitionID to,
                             const HypernodeWeight weight) {
    ASSERT(from != kInvalidPartition && from < _k);
    ASSERT(to != kInvalidPartition && to < _k);
    ASSERT(from != to);
    const HypernodeID local_bucket = u % _num_local_block_infos;
    --_part_info[from][local_bucket].size;
    _part_info[from][local_bucket].weight -= weight;
    ++_part_info[to][local_bucket].size;
    _part_info[to][local_bucket].weight += weight;
  }

  inline HypernodeWeight partWeight(const PartitionID id) const {
    ASSERT(id != kInvalidPartition && id < _k);
    HypernodeWeight weight = 0;
    for ( size_t i = 0; i < _num_local_block_infos; ++i ) {
      weight += _part_info[id][i].weight;
    }
    return weight;
  }

  inline HypernodeID partSize(const PartitionID id) const {
    ASSERT(id != kInvalidPartition && id < _k);
    HypernodeID size = 0;
    for ( size_t i = 0; i < _num_local_block_infos; ++i ) {
      size += _part_info[id][i].size;
    }
    return size;
  }

  void reset() {
    for ( PartitionID block = 0; block < _k; ++block ) {
      for ( HypernodeID local_bucket = 0; local_bucket < _num_local_block_infos; ++local_bucket ) {
        _part_info[block][local_bucket].size = ID(0);
        _part_info[block][local_bucket].weight = 0;
      }
    }
  }

 private:
  HypernodeID computeNumLocalBlockInfos(const size_t num_threads, const PartitionID k) {
    const double log2_num_threads = std::log2(static_cast<double>(num_threads));
    const double log2_k = std::log2(static_cast<double>(k));
    return std::max(std::ceil(log2_num_threads / log2_k), 1.0);
  }

  PartitionID _k;
  HypernodeID _num_local_block_infos;
  parallel::scalable_vector<LocalBlockInfo> _part_info;
};
}  // namespace ds
}  // namespace mt_kahypar

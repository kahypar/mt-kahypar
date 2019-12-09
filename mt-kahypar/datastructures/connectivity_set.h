/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include <atomic>
#include <type_traits>

#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/bit_ops.h"
#include "mt-kahypar/utils/range.h"

#include "mt-kahypar/macros.h"

namespace mt_kahypar {
namespace ds {


class ConnectivitySets {
public:

  static constexpr bool debug = true;
  using PartitionID = uint32_t;
  using HyperedgeID = uint64_t;	// TODO how to keep synced with definitions.h, other than templates?


  ConnectivitySets(const HyperedgeID numEdges, const PartitionID k) : k(k), numEdges(numEdges),
                                                                      numBlocksPerHyperedge(k / bits_per_block + (k % bits_per_block != 0)),
                                                                      bits(numEdges * numBlocksPerHyperedge),
                                                                      connectivity_cache(numEdges)
  {
  }

  void add(const HyperedgeID he, const PartitionID p) {
    toggle(he, p);
    connectivity_cache[he].fetch_add(1, std::memory_order_relaxed);
  }

  void remove(const HyperedgeID he, const PartitionID p) {
    toggle(he, p);
    connectivity_cache[he].fetch_sub(1, std::memory_order_relaxed);
  }

  bool contains(const HyperedgeID he, const PartitionID p) const {
    const size_t div = p / bits_per_block;
    const size_t rem = p % bits_per_block;
    return bits[he * numBlocksPerHyperedge + div].load(std::memory_order_relaxed) & (UnsafeBlock(1) << rem);
  }

  // not threadsafe
  void clear(const HyperedgeID he) {
    for (size_t i = he * numBlocksPerHyperedge; i < (he + 1) * numBlocksPerHyperedge; ++i) {
      bits[i].store(0, std::memory_order_relaxed);
    }
    connectivity_cache[he].store(0, std::memory_order_relaxed);
  }

  // Note: this might differ from the lookup
  PartitionID computeConnectivity(const HyperedgeID he) const {
    PartitionID conn = 0;
    for (size_t i = he * numBlocksPerHyperedge; i < (he + 1) * numBlocksPerHyperedge; ++i) {
      conn += utils::popcount_64(bits[i].load(std::memory_order_relaxed));
    }
    return conn;
  }

  // Note: this might differ from the number of actually set bits
  PartitionID connectivity(const HyperedgeID he) const {
    return connectivity_cache[he].load(std::memory_order_relaxed);
  }


private:
  using UnsafeBlock = uint64_t;
  static constexpr int bits_per_block = std::numeric_limits<UnsafeBlock>::digits;
  using Block = parallel::IntegralAtomicWrapper<UnsafeBlock>;
  using BlockIterator = std::vector<Block>::const_iterator;

	
	PartitionID k;
	HyperedgeID numEdges;
	PartitionID numBlocksPerHyperedge;
	std::vector<Block> bits;
	std::vector< parallel::IntegralAtomicWrapper<PartitionID> > connectivity_cache;

	void toggle(const HyperedgeID he, const PartitionID p) {
	  assert(p < k);
	  assert(he < numEdges);
    const size_t div = p / bits_per_block, rem = p % bits_per_block;
    const size_t idx = he * numBlocksPerHyperedge + div;
    assert(idx < bits.size());
	  bits[idx].fetch_xor(UnsafeBlock(1) << rem, std::memory_order_relaxed);
	}

public:

  class Iterator : public std::iterator<std::forward_iterator_tag, PartitionID, std::ptrdiff_t, const PartitionID*, PartitionID> {
  public:
    Iterator(BlockIterator first, PartitionID part, PartitionID k) : currentPartition(part), k(k), firstBlockIt(first) {
      findNextBit();
    }

    PartitionID operator*() const {
      return currentPartition;
    }

    Iterator& operator++() {
      findNextBit();
      return *this;
    }

    Iterator operator++(int ) {
      const Iterator res = *this;
      findNextBit();
      return res;
    }

    bool operator==(const Iterator& o) const {
      return currentPartition == o.currentPartition && firstBlockIt == o.firstBlockIt;
    }

    bool operator!=(const Iterator& o) const {
      return !operator==(o);
    }

  private:
    PartitionID currentPartition;
    PartitionID k;
    BlockIterator firstBlockIt;

    void findNextBit() {
      ++currentPartition;
      while (currentBlock()->load(std::memory_order_relaxed) >> (currentPartition % bits_per_block) == 0 && currentPartition < k) {
        currentPartition += (bits_per_block - (currentPartition % bits_per_block));   // skip rest of block
      }
      if (currentPartition < k) {
        UnsafeBlock b = currentBlock()->load(std::memory_order_relaxed) >> (currentPartition % bits_per_block);
        if (b != 0) {
          currentPartition += utils::lowest_set_bit_64(b);
        } else {
          // this should only happen if another thread has removed the incident parts set in the currentBlock() since exiting the while loop
          currentPartition += (bits_per_block - (currentPartition % bits_per_block));   // skip rest of block
        }
      } else {
        currentPartition = k;
      }
    }

    BlockIterator currentBlock() const {
      return firstBlockIt + currentPartition / bits_per_block;
    }

  };

	Iterator hyperedgeBegin(const HyperedgeID he) const {
	  return Iterator(bits.begin() + he * numBlocksPerHyperedge, -1, k);
	}

	Iterator hyperedgeEnd(const HyperedgeID he) const {
	  return Iterator(bits.begin() + he * numBlocksPerHyperedge, k-1, k);
	}


  IteratorRange<Iterator> connectivitySet(const HyperedgeID he) const {
    return IteratorRange<Iterator>(hyperedgeBegin(he), hyperedgeEnd(he));
  }

};
	


}  // namespace ds
}  // namespace mt_kahypar

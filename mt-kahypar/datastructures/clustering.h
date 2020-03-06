/*******************************************************************************
 * This file is part of KaHyPar.
 *
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

#include <cassert>
#include <functional>

#include <tbb/parallel_for.h>

#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include <mt-kahypar/definitions.h>

namespace mt_kahypar {
namespace ds {
class Clustering : public std::vector<PartitionID> {
 public:
  using Base = std::vector<PartitionID>;

  explicit Clustering(size_t n) :
    Base(n) { }

  // make Clustering a callable, so we don't need to wrap it in other callables.
  const PartitionID & operator() (const size_t x) const {
    return operator[] (x);
  }
  PartitionID & operator() (const size_t x) {
    return operator[] (x);
  }

  void assignSingleton() {
    tbb::parallel_for(PartitionID(0), static_cast<PartitionID>(size()), [&](PartitionID i) {
          (*this)[i] = i;
        });
  }

  size_t compactify(PartitionID upperIDBound = -1, size_t numTasks = 1) {
    if (upperIDBound < 0)
      upperIDBound = static_cast<PartitionID>(size()) - 1;
    const PartitionID res = numTasks > 1 ?
      parallelCompactify(upperIDBound, numTasks) :
      sequentialCompactify(upperIDBound);
    return static_cast<size_t>(res);
  }

 private:
  PartitionID sequentialCompactify(PartitionID upperIDBound) {
    std::vector<PartitionID> mapping(upperIDBound + 1, -1);
    PartitionID i = 0;
    for (PartitionID& c : *this) {
      if (mapping[c] == -1)
        mapping[c] = i++;
      c = mapping[c];
    }
    return i;
  }

  PartitionID parallelCompactify(PartitionID upperIDBound, size_t numTasks) {
    // TODO implement RoutingKit style rank bitvector and parallelize

#ifdef KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS
    Clustering seq = *this;
    PartitionID numClustersFromSeq = seq.sequentialCompactify(upperIDBound);
#endif

    std::vector<PartitionID> mapping(upperIDBound + 1, 0);
    tbb::parallel_for_each(*this, [&](const PartitionID& c) {
          mapping[c] = 1;
        });

    parallel::PrefixSum::parallelTwoPhase(mapping.begin(), mapping.end(), mapping.begin(), std::plus<PartitionID>(), PartitionID(0), numTasks);
    // PrefixSum::parallelTBBNative(mapping.begin(), mapping.end(), mapping.begin(), std::plus<PartitionID>(), PartitionID(0), numTasks);
    // NOTE Benchmark!

    tbb::parallel_for_each(*this, [&](PartitionID& c) {
          c = mapping[c];
        });

#ifdef KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS
    assert(numClustersFromSeq == mapping.back());
    for (size_t i = 0; i < size(); ++i)
      assert((*this)[i] == seq[i]);
#endif
    return mapping.back();
  }
};
}  // namespace ds
}  // namespace mt_kahypar

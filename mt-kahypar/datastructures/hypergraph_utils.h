/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2026 Nikolai Maas <nikolai.maas@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#pragma once

#include "mt-kahypar/datastructures/hypergraph_common.h"

#include <atomic>
#include <limits>

#include <tbb/parallel_reduce.h>
#include <tbb/enumerable_thread_specific.h>


namespace mt_kahypar {

namespace impl {
  inline HNWeightScalar add(HNWeightScalar lhs, HNWeightScalar rhs, std::atomic_bool& error_flag) {
    if (lhs <= std::numeric_limits<HNWeightScalar>::max() - rhs) {
      return lhs + rhs;
    } else {
      error_flag.store(true, std::memory_order_relaxed);
      return lhs;  // prevent UB from signed overflow
    }
  }

  struct safe_addition {
    std::atomic_bool error_flag;

    HNWeightScalar operator()(const HNWeightScalar& lhs, const HNWeightScalar& rhs) {
      return add(lhs, rhs, error_flag);
    }
  };
}

template<typename Hypergraph>
void computeTotalNodeWeightParallel(const Hypergraph& hypergraph, AllocatedHNWeight& total_weight, AllocatedHNWeight& max_weight) {
  // For some reason, TBB has difficulty handling the exception if it is thrown from within the parallel loop
  // (crashes in debug mode, usually works in release mode but hangs in rare cases). Therefore, we instead set
  // a flag and throw the exception after the calculation is finished.

  impl::safe_addition adder {false};
  tbb::enumerable_thread_specific<AllocatedHNWeight> local_sum(hypergraph.dimension(), 0);
  tbb::enumerable_thread_specific<AllocatedHNWeight> local_max(hypergraph.dimension(), 0);
  hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
    auto hn_weight = hypergraph.nodeWeight(hn);
    local_sum.local() = weight::mapBinary(local_sum.local(), hn_weight,
      [&](HNWeightScalar lhs, HNWeightScalar rhs) {
        return adder(lhs, rhs);
      });
    local_max.local() = weight::max(local_max.local(), hn_weight);
  });

  total_weight = weight::broadcast(0, hypergraph.dimension());
  for (const auto& weight: local_sum) {
    total_weight += weight;
  }
  max_weight = weight::broadcast(0, hypergraph.dimension());
  for (const auto& weight: local_max) {
    max_weight = weight::max(max_weight, weight);
  }

  // TODO: multi-constraint overflow and epsilon
  if (adder.error_flag) {
    throw InvalidInputException("total node weight overflows weight data type");
  }
}

} // namespace mt_kahypar

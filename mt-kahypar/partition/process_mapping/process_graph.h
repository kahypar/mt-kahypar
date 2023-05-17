/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/static_graph.h"
#include "mt-kahypar/datastructures/static_bitset.h"

namespace mt_kahypar {

class ProcessGraph {

  static constexpr size_t MEMORY_LIMIT = 100000000;

 public:
  explicit ProcessGraph(ds::StaticGraph&& graph) :
    _is_initialized(true),
    _k(kInvalidPartition),
    _graph(std::move(graph)),
    _max_precomputed_connectitivty(0),
    _distances() {
    _k = _graph.initialNumNodes();
  }

  ProcessGraph(const ProcessGraph&) = delete;
  ProcessGraph & operator= (const ProcessGraph &) = delete;

  ProcessGraph(ProcessGraph&&) = default;
  ProcessGraph & operator= (ProcessGraph &&) = default;

  PartitionID numBlocks() const {
    return _k;
  }

  void precomputeDistances(const size_t max_conectivity);

  HyperedgeWeight distance(const ds::StaticBitset& connectivity_set);

  HyperedgeWeight distance(const PartitionID i, const PartitionID j) {
    ASSERT(_is_initialized);
    return _distances[index(i, j)];
  }

 private:
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE size_t index(const PartitionID i,
                                                  const PartitionID j) {
    ASSERT(i < _k && j < _k);
    return i + j * _k;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE size_t index(const ds::StaticBitset& connectivity_set) {
    size_t index = 0;
    PartitionID multiplier = 1;
    PartitionID last_block = kInvalidPartition;
    for ( const PartitionID block : connectivity_set ) {
      index += multiplier * block;
      multiplier *= _k;
      last_block = block;
    }
    return index + (multiplier == _k ? last_block * _k : 0);
  }

  bool _is_initialized;
  PartitionID _k;
  ds::StaticGraph _graph;
  PartitionID _max_precomputed_connectitivty;
  vec<HyperedgeWeight> _distances;
};

}  // namespace kahypar

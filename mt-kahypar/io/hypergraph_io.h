/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <string>
#include <cstdint>

#include "utils/hash/murmur2_hash.hpp"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace io {
  using Hyperedge = vec<HypernodeID>;
  using HyperedgeVector = vec<Hyperedge>;

  void readHypergraphFile(const std::string& filename,
                          HyperedgeID& num_hyperedges,
                          HypernodeID& num_hypernodes,
                          HyperedgeID& num_removed_single_pin_hyperedges,
                          HyperedgeVector& hyperedges,
                          vec<HyperedgeWeight>& hyperedges_weight,
                          vec<HypernodeWeight>& hypernodes_weight,
                          const bool remove_single_pin_hes = true);

  void readGraphFile(const std::string& filename,
                     HyperedgeID& num_hyperedges,
                     HypernodeID& num_hypernodes,
                     HyperedgeVector& hyperedges,
                     vec<HyperedgeWeight>& hyperedges_weight,
                     vec<HypernodeWeight>& hypernodes_weight);

  void readPartitionFile(const std::string& filename, HypernodeID num_nodes, std::vector<PartitionID>& partition);
  void readPartitionFile(const std::string& filename, HypernodeID num_nodes, PartitionID* partition);

  void readFrequencyFile(const std::string& filename,
                         ds::DynamicSparseMap<__uint128_t, float>& frequencies,
                         double binary_threshold,
                         double ternary_threshold,
                         double ternary_value);

  inline __uint128_t hashedEdgeKey(const utils_tm::hash_tm::murmur2_hash& hasher, HypernodeID u, HypernodeID v) {
    // this is a bit ugly, but we need it to avoid collisions (since we don't have a good hashmap with a non-trivial hash function)
    uint32_t u_val = u, v_val = v;
    uint64_t key_left = (static_cast<uint64_t>(u_val) << 32) | static_cast<uint64_t>(v_val);
    return (static_cast<__uint128_t>(key_left) << 64) | static_cast<__uint128_t>(hasher(key_left));
  }

  template<typename PartitionedHypergraph>
  void writePartitionFile(const PartitionedHypergraph& phg, const std::string& filename);

}  // namespace io
}  // namespace mt_kahypar

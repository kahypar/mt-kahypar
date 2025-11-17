/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2025 Nikolai Maas <nikolai.maas@kit.edu>
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
#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {
namespace rebalancer {

struct PotentialMove {
  HypernodeID node;
  PartitionID to;
  float rating;
};

template <typename GraphAndGainTypes>
class Fallback {
 private:
  using PartitionedHypergraph = typename GraphAndGainTypes::PartitionedHypergraph;
  using GainCache = typename GraphAndGainTypes::GainCache;
  using AttributedGains = typename GraphAndGainTypes::AttributedGains;

  static constexpr bool debug = false;

 public:
  static std::pair<int64_t, size_t> runDeadlockFallback(PartitionedHypergraph& phg,
                                                        GainCache& gain_cache,
                                                        const Context& context,
                                                        vec<ds::StreamingVector<rebalancer::PotentialMove>>& tmp_potential_moves,
                                                        ds::Array<Move>& moves,
                                                        ds::Array<MoveID>& move_id_of_node,
                                                        ds::Array<uint8_t>& node_is_locked,
                                                        const vec<double>& weight_normalizer,
                                                        size_t& global_move_id);
};

} // namespace rebalancer
}  // namespace mt_kahypar

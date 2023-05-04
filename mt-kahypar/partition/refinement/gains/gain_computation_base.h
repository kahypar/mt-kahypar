/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
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

#include <vector>

#include "kahypar/meta/mandatory.h"

#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

template <class Derived = Mandatory>
class GainComputationBase {
  using DeltaGain = tbb::enumerable_thread_specific<Gain>;

 public:
  using RatingMap = ds::SparseMap<PartitionID, Gain>;
  using TmpScores = tbb::enumerable_thread_specific<RatingMap>;

  GainComputationBase(const Context& context) :
    _context(context),
    _deltas(0),
    _tmp_scores(context.partition.k) { }

  template<typename PartitionedHypergraph>
  Move computeMaxGainMove(const PartitionedHypergraph& phg,
                          const HypernodeID hn,
                          const bool rebalance = false,
                          const bool consider_non_adjacent_blocks = false) {
    return static_cast<Derived*>(this)->computeMaxGainMoveImpl(
      phg, hn, rebalance, consider_non_adjacent_blocks);
  }

  inline void computeDeltaForHyperedge(const HyperedgeID he,
                                       const HyperedgeWeight edge_weight,
                                       const HypernodeID edge_size,
                                       const HypernodeID pin_count_in_from_part_after,
                                       const HypernodeID pin_count_in_to_part_after) {
    static_cast<Derived*>(this)->computeDeltaForHyperedgeImpl(
      he, edge_weight, edge_size,
      pin_count_in_from_part_after, pin_count_in_to_part_after);
  }

  // ! Returns the delta in the objective function for all moves
  // ! performed by the calling thread relative to the last call
  // ! reset()
  Gain localDelta() {
    return _deltas.local();
  }

  // ! Returns the overall delta of all moves performed by
  // ! all threads relative to the last call of reset()
  Gain delta() const {
    Gain overall_delta = 0;
    for (const Gain& delta : _deltas) {
      overall_delta += delta;
    }
    return overall_delta;
  }

  void reset() {
    for (Gain& delta : _deltas) {
      delta = 0;
    }
  }

  void changeNumberOfBlocks(const PartitionID new_k) {
    for ( auto& tmp_score : _tmp_scores ) {
      if ( static_cast<size_t>(new_k) > tmp_score.size() ) {
        tmp_score = RatingMap(new_k);
      }
    }
  }

 protected:
  const Context& _context;
  DeltaGain _deltas;
  TmpScores _tmp_scores;
};

}  // namespace mt_kahypar

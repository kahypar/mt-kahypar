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

#include "mt-kahypar/partition/refinement/gains/gain_computation_base.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {

class Km1GainComputation : public GainComputationBase<Km1GainComputation> {
  using Base = GainComputationBase<Km1GainComputation>;
  using RatingMap = typename Base::RatingMap;

  static constexpr bool enable_heavy_assert = false;

 public:
  Km1GainComputation(const Context& context,
                     bool disable_randomization = false) :
    Base(context),
    _disable_randomization(disable_randomization) { }

  template<typename PartitionedHypergraph>
  Move computeMaxGainMoveImpl(const PartitionedHypergraph& phg,
                              const HypernodeID hn,
                              const bool rebalance,
                              const bool consider_non_adjacent_blocks) {
    RatingMap& tmp_scores = _tmp_scores.local();
    ASSERT(tmp_scores.size() == 0, "Rating map not empty");

    PartitionID from = phg.partID(hn);
    Gain internal_weight = 0;
    for (const HyperedgeID& he : phg.incidentEdges(hn)) {
      HypernodeID pin_count_in_from_part = phg.pinCountInPart(he, from);
      HyperedgeWeight he_weight = phg.edgeWeight(he);

      // In case, there is more one than one pin left in from part, we would
      // increase the connectivity, if we would move the pin to one block
      // no contained in the connectivity set. In such cases, we can only
      // increase the connectivity of a hyperedge and therefore gather
      // the edge weight of all those edges and add it later to move gain
      // to all other blocks.
      if ( pin_count_in_from_part > 1 ) {
        internal_weight += he_weight;
      }

      // Substract edge weight of all incident blocks.
      // Note, in case the pin count in from part is greater than one
      // we will later add that edge weight to the gain (see internal_weight).
      for (const PartitionID& to : phg.connectivitySet(he)) {
        if (from != to) {
          tmp_scores[to] -= he_weight;
        }
      }
    }

    Move best_move { from, from, hn, rebalance ? std::numeric_limits<Gain>::max() : 0 };
    HypernodeWeight hn_weight = phg.nodeWeight(hn);
    int cpu_id = SCHED_GETCPU;
    utils::Randomize& rand = utils::Randomize::instance();
    auto test_and_apply = [&](const PartitionID to,
                              const Gain score,
                              const bool no_tie_breaking = false) {
      bool new_best_gain = (score < best_move.gain) ||
                            (score == best_move.gain &&
                            !_disable_randomization &&
                            (no_tie_breaking || rand.flipCoin(cpu_id)));
      if (new_best_gain && phg.partWeight(to) + hn_weight <=
          _context.partition.max_part_weights[to]) {
        best_move.to = to;
        best_move.gain = score;
        return true;
      } else {
        return false;
      }
    };

    for ( const auto& entry : tmp_scores ) {
      const PartitionID to = entry.key;
      if (from != to) {
        const Gain score = entry.value + internal_weight;
        test_and_apply(to, score);
      }
    }

    if ( consider_non_adjacent_blocks && best_move.from == from ) {
      // This is important for our rebalancer as the last fallback strategy
      vec<PartitionID> non_adjacent_block;
      for ( PartitionID to = 0; to < _context.partition.k; ++to ) {
        if ( from != to && !tmp_scores.contains(to) ) {
          // This block is not adjacent to the current node
          if ( test_and_apply(to, internal_weight, true /* no tie breaking */ ) ) {
            non_adjacent_block.push_back(to);
          }
        }
      }

      if ( non_adjacent_block.size() > 0 ) {
        // Choose one at random
        const PartitionID to = non_adjacent_block[
          rand.getRandomInt(0, static_cast<int>(non_adjacent_block.size() - 1), cpu_id)];
        best_move.to = to;
        best_move.gain = internal_weight;
      }
    }

    tmp_scores.clear();
    return best_move;
  }

  inline void computeDeltaForHyperedgeImpl(const HyperedgeID he,
                                           const HyperedgeWeight edge_weight,
                                           const HypernodeID edge_size,
                                           const HypernodeID pin_count_in_from_part_after,
                                           const HypernodeID pin_count_in_to_part_after) {
    _deltas.local() += km1Delta(he, edge_weight, edge_size,
      pin_count_in_from_part_after, pin_count_in_to_part_after);
  }

  using Base::_context;
  using Base::_deltas;
  using Base::_tmp_scores;
  bool _disable_randomization;
};

}  // namespace mt_kahypar

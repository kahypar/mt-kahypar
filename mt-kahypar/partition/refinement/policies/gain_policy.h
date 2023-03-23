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

#include "tbb/enumerable_thread_specific.h"

#include "kahypar/meta/policy_registry.h"
#include "kahypar/meta/typelist.h"

#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

template <class Derived = Mandatory,
          class HyperGraph = Mandatory>
class GainPolicy : public kahypar::meta::PolicyBase {
  using DeltaGain = tbb::enumerable_thread_specific<Gain>;

 public:
  using RatingMap = ds::SparseMap<PartitionID, Gain>;
  using TmpScores = tbb::enumerable_thread_specific<RatingMap>;

  GainPolicy(const Context& context) :
    _context(context),
    _deltas(0),
    _tmp_scores(context.partition.k) { }

  Move computeMaxGainMove(const HyperGraph& hypergraph,
                          const HypernodeID hn,
                          const bool rebalance = false,
                          const bool consider_non_adjacent_blocks = false) {
    return static_cast<Derived*>(this)->computeMaxGainMoveImpl(
      hypergraph, hn, rebalance, consider_non_adjacent_blocks);
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

template <class HyperGraph = Mandatory>
class Km1Policy : public GainPolicy<Km1Policy<HyperGraph>, HyperGraph> {
  using Base = GainPolicy<Km1Policy<HyperGraph>, HyperGraph>;
  using RatingMap = typename Base::RatingMap;

  static constexpr bool enable_heavy_assert = false;

 public:
  Km1Policy(const Context& context,
            bool disable_randomization = false) :
    Base(context),
    _disable_randomization(disable_randomization) { }

  Move computeMaxGainMoveImpl(const HyperGraph& hypergraph,
                              const HypernodeID hn,
                              const bool rebalance,
                              const bool consider_non_adjacent_blocks) {
    RatingMap& tmp_scores = _tmp_scores.local();
    ASSERT(tmp_scores.size() == 0, "Rating map not empty");

    PartitionID from = hypergraph.partID(hn);
    Gain internal_weight = 0;
    for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
      HypernodeID pin_count_in_from_part = hypergraph.pinCountInPart(he, from);
      HyperedgeWeight he_weight = hypergraph.edgeWeight(he);

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
      for (const PartitionID& to : hypergraph.connectivitySet(he)) {
        if (from != to) {
          tmp_scores[to] -= he_weight;
        }
      }
    }

    Move best_move { from, from, hn, rebalance ? std::numeric_limits<Gain>::max() : 0 };
    HypernodeWeight hn_weight = hypergraph.nodeWeight(hn);
    int cpu_id = SCHED_GETCPU;
    utils::Randomize& rand = utils::Randomize::instance();
    auto test_and_apply = [&](const PartitionID to,
                              const Gain score,
                              const bool no_tie_breaking = false) {
      bool new_best_gain = (score < best_move.gain) ||
                            (score == best_move.gain &&
                            !_disable_randomization &&
                            (no_tie_breaking || rand.flipCoin(cpu_id)));
      if (new_best_gain && hypergraph.partWeight(to) + hn_weight <=
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
                                            pin_count_in_from_part_after,
                                            pin_count_in_to_part_after);
  }

  using Base::_context;
  using Base::_deltas;
  using Base::_tmp_scores;
  bool _disable_randomization;
};

template <class HyperGraph = Mandatory>
class CutPolicy : public GainPolicy<CutPolicy<HyperGraph>, HyperGraph> {
  using Base = GainPolicy<CutPolicy<HyperGraph>, HyperGraph>;
  using RatingMap = typename Base::RatingMap;

  static constexpr bool enable_heavy_assert = false;

 public:
  CutPolicy(const Context& context,
            bool disable_randomization = false) :
    Base(context),
    _disable_randomization(disable_randomization) { }

  Move computeMaxGainMoveImpl(const HyperGraph& hypergraph,
                              const HypernodeID hn,
                              const bool rebalance,
                              const bool consider_non_adjacent_blocks) {
    RatingMap& tmp_scores = _tmp_scores.local();
    ASSERT(tmp_scores.size() == 0, "Rating map not empty");

    PartitionID from = hypergraph.partID(hn);
    Gain internal_weight = 0;
    for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
      PartitionID connectivity = hypergraph.connectivity(he);
      HypernodeID pin_count_in_from_part = hypergraph.pinCountInPart(he, from);
      HyperedgeWeight weight = hypergraph.edgeWeight(he);
      if (connectivity == 1) {
        ASSERT(hypergraph.edgeSize(he) > 1);
        // In case, the hyperedge is a non-cut hyperedge, we would increase
        // the cut, if we move vertex hn to an other block.
        internal_weight += weight;
      } else if (connectivity == 2 && pin_count_in_from_part == 1) {
        for (const PartitionID& to : hypergraph.connectivitySet(he)) {
          // In case there are only two blocks contained in the current
          // hyperedge and only one pin left in the from part of the hyperedge,
          // we would make the current hyperedge a non-cut hyperedge when moving
          // vertex hn to the other block.
          if (from != to) {
            tmp_scores[to] -= weight;
          }
        }
      }
    }

    Move best_move { from, from, hn, rebalance ? std::numeric_limits<Gain>::max() : 0 };
    HypernodeWeight hn_weight = hypergraph.nodeWeight(hn);
    int cpu_id = SCHED_GETCPU;
    utils::Randomize& rand = utils::Randomize::instance();
    auto test_and_apply = [&](const PartitionID to,
                              const Gain score,
                              const bool no_tie_breaking = false) {
      bool new_best_gain = (score < best_move.gain) ||
                            (score == best_move.gain &&
                            !_disable_randomization &&
                            (no_tie_breaking || rand.flipCoin(cpu_id)));
      if (new_best_gain && hypergraph.partWeight(to) + hn_weight <=
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
    _deltas.local() += cutDelta(he, edge_weight, edge_size,
                                            pin_count_in_from_part_after,
                                            pin_count_in_to_part_after);
  }

  using Base::_context;
  using Base::_deltas;
  using Base::_tmp_scores;
  bool _disable_randomization;
};
}  // namespace mt_kahypar

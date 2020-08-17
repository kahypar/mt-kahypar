/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include <boost/dynamic_bitset.hpp>

#include "tbb/parallel_invoke.h"
#include "tbb/parallel_scan.h"

#include "mt-kahypar/partition/refinement/fm/fm_commons.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

struct BalanceAndBestIndexScan {
  const PartitionedHypergraph& phg;
  const vec<Move>& moves;

  struct Prefix {
    Gain gain = 0;                           /** gain when using valid moves up to best_index */
    MoveID best_index = 0;                   /** local ID of first move to revert */
    HypernodeWeight heaviest_weight =
            std::numeric_limits<HypernodeWeight>::max();   /** weight of the heaviest part */

    bool operator<(const Prefix& o) const {
      return gain > o.gain ||
             (gain == o.gain && std::tie(heaviest_weight, best_index) < std::tie(o.heaviest_weight, o.best_index));
    }
  };
  std::shared_ptr< tbb::enumerable_thread_specific<Prefix> > local_best;

  Gain gain_sum = 0;

  vec<HypernodeWeight> part_weights;
  const std::vector<HypernodeWeight>& max_part_weights;

  BalanceAndBestIndexScan(BalanceAndBestIndexScan& b, tbb::split) :
          phg(b.phg),
          moves(b.moves),
          local_best(b.local_best),
          gain_sum(0),
          part_weights(b.part_weights.size(), 0),
          max_part_weights(b.max_part_weights) { }


  BalanceAndBestIndexScan(const PartitionedHypergraph& phg,
                          const vec<Move>& moves,
                          const vec<HypernodeWeight>& part_weights,
                          const std::vector<HypernodeWeight>& max_part_weights) :
          phg(phg),
          moves(moves),
          local_best(std::make_shared< tbb::enumerable_thread_specific<Prefix> >()),
          part_weights(part_weights),
          max_part_weights(max_part_weights)
  {
  }


  void operator()(const tbb::blocked_range<MoveID>& r, tbb::pre_scan_tag ) {
    for (MoveID i = r.begin(); i < r.end(); ++i) {
      const Move& m = moves[i];
      if (m.gain != invalidGain) {  // skip locally reverted moves
        gain_sum += m.gain;
        part_weights[m.from] -= phg.nodeWeight(m.node);
        part_weights[m.to] += phg.nodeWeight(m.node);
      }
    }
  }

  // subranges a | b | c | d . assuming this ran pre_scan on c,
  // then lhs ran pre_scan on b and final_scan of this will be on d
  void reverse_join(BalanceAndBestIndexScan& lhs) {
    for (size_t i = 0; i < part_weights.size(); ++i) {
      part_weights[i] += lhs.part_weights[i];
    }
    gain_sum += lhs.gain_sum;
  }

  void operator()(const tbb::blocked_range<MoveID>& r, tbb::final_scan_tag ) {
    size_t overloaded = 0;
    for (size_t i = 0; i < part_weights.size(); ++i) {
      if (part_weights[i] > max_part_weights[i]) {
        overloaded++;
      }
    }

    Prefix current;
    for (MoveID i = r.begin(); i < r.end(); ++i) {
      const Move& m = moves[i];

      if (m.gain != invalidGain) {  // skip locally reverted moves
        gain_sum += m.gain;

        const bool from_overloaded = part_weights[m.from] > max_part_weights[m.from];
        part_weights[m.from] -= phg.nodeWeight(m.node);
        if (from_overloaded && part_weights[m.from] <= max_part_weights[m.from]) {
          overloaded--;
        }
        const bool to_overloaded = part_weights[m.to] > max_part_weights[m.to];
        part_weights[m.to] += phg.nodeWeight(m.node);
        if (!to_overloaded && part_weights[m.to] > max_part_weights[m.to]) {
          overloaded++;
        }

        if (overloaded == 0 && gain_sum >= current.gain) {
          Prefix new_prefix = { gain_sum, i + 1, *std::max_element(part_weights.begin(), part_weights.end()) };
          current = std::min(current, new_prefix);
        }
      }
    }

    if (current.best_index != 0) {
      Prefix& lb = local_best->local();
      lb = std::min(lb, current);
    }
  }

  void assign(BalanceAndBestIndexScan& b) {
    gain_sum = b.gain_sum;
  }

  Prefix finalize(const vec<HypernodeWeight>& initial_part_weights) {
    Prefix res { 0, 0, *std::max_element(initial_part_weights.begin(), initial_part_weights.end()) };
    for (const Prefix& x : *local_best) {
      res = std::min(res, x);
    }
    return res;
  }
};


class GlobalRollback {
  static constexpr bool enable_heavy_assert = false;
public:
  explicit GlobalRollback(const Hypergraph& hg, const Context& context, PartitionID numParts) :
    context(context),
    maxPartWeightScaling(context.refinement.fm.rollback_balance_violation_factor),
    numParts(numParts),
    remaining_original_pins(),
    first_move_in(),
    last_move_out() {

    tbb::parallel_invoke([&] {
      remaining_original_pins.resize("Refinement", "remaining_original_pins",
        static_cast<size_t>(hg.numNonGraphEdges()) * numParts);
    }, [&] {
      first_move_in.resize("Refinement", "first_move_in",
        static_cast<size_t>(hg.numNonGraphEdges()) * numParts);
    }, [&] {
      last_move_out.resize("Refinement", "last_move_out",
        static_cast<size_t>(hg.numNonGraphEdges()) * numParts);
    });
    resetStoredMoveIDs();
  }

  HyperedgeWeight revertToBestPrefix(PartitionedHypergraph& phg,
                                     FMSharedData& sharedData,
                                     const vec<HypernodeWeight>& partWeights);

  HyperedgeWeight revertToBestPrefixParallel(PartitionedHypergraph& phg,
                                             FMSharedData& sharedData,
                                             const vec<HypernodeWeight>& partWeights,
                                             const std::vector<HypernodeWeight>& maxPartWeights);


  void recalculateGains(PartitionedHypergraph& phg, FMSharedData& sharedData);

  MoveID lastMoveOut(HyperedgeID he, PartitionID block) const {
    return last_move_out[size_t(he) * numParts + block].load(std::memory_order_relaxed);
  }

  MoveID firstMoveIn(HyperedgeID he, PartitionID block) const {
    return first_move_in[size_t(he) * numParts + block].load(std::memory_order_relaxed);
  }

  HypernodeID remainingPinsFromBeginningOfMovePhase(HyperedgeID he, PartitionID block) const {
    return remaining_original_pins[size_t(he) * numParts + block].load(std::memory_order_relaxed);
  }

  bool verifyPinCountMatches(const PartitionedHypergraph& phg) const ;

  bool verifyGains(PartitionedHypergraph& phg, FMSharedData& sharedData);

  void resetStoredMoveIDs();

  void setRemainingOriginalPins(PartitionedHypergraph& phg);

  void memoryConsumption(utils::MemoryTreeNode* parent) const;

  const Context& context;

  // ! Factor to multiply max part weight with, in order to relax or disable the balance criterion. Set to zero for disabling
  double maxPartWeightScaling;

  PartitionID numParts;

  // ! For each hyperedge and block, the number of pins_in_part
  // ! at the beginning of a move phase minus the number of moved out pins
  // ! A value of zero means at some point a gain improvement might have occured, which then has to be checked
  // ! with the move ID flagging
  ds::Array<CAtomic<HypernodeID>> remaining_original_pins;
  // ! For each hyperedge and each block, the ID of the first move to place a pin in that block
  ds::Array<CAtomic<MoveID>> first_move_in;
  // ! For each hyperedge and each block, the ID of the last move to remove a pin from that block
  ds::Array<CAtomic<MoveID>> last_move_out;

};

}
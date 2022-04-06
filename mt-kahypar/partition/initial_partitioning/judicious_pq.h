/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2022 Noah Wahl <noah.wahl@student.kit.edu>
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
#include <mt-kahypar/datastructures/priority_queue.h>
#include <mt-kahypar/definitions.h>
#include <mt-kahypar/partition/context.h>

namespace mt_kahypar {
class JudiciousPQ final {
public:
  // using positive penalties as keys to make secondary criteria of lighter
  // blocks easier
  using PriorityQueue = ds::ExclusiveHandleHeap<
      ds::Heap<HypernodeWeight, PartitionID, std::greater<>>>;
  using BlockPQ = ds::ExclusiveHandleHeap<
      ds::Heap<std::pair<HypernodeWeight, HypernodeWeight>, PartitionID,
               std::greater<>>>;

  explicit JudiciousPQ(const Context &context, HypernodeID num_nodes)
      : _context(context), _toPQs(static_cast<size_t>(context.partition.k),
                                  PriorityQueue(num_nodes)),
        _blockPQ(static_cast<size_t>(context.partition.k)) {}

  void insert(const PartitionedHypergraph &phg, const HypernodeID v) {
    const PartitionID pv = phg.partID(v);
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      if (i == pv)
        continue;
      Gain gain = 0;
      if (pv == kInvalidPartition) {
        gain = computeGainForInvalidFrom(phg, v, i);
      } else {
        // const Gain gain = computeGain(phg, v, pv, i);
      }
      _toPQs[i].insert(v, gain);
    }
  }

  void initBlockPQ(PartitionedHypergraph &phg) {
    ASSERT(_blockPQ.empty());
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      if (!_toPQs[i].empty()) {
        _blockPQ.insert(i, std::make_pair(phg.partLoad(i), _toPQs[i].topKey()));
      }
    }
  }

  void updateGain(const PartitionedHypergraph &phg, const HypernodeID v,
                  const HyperedgeID he, const PartitionID to) {
    if (phg.partID(v) == kInvalidPartition) {
      _toPQs[to].increaseKey(v, _toPQs[to].keyOf(v) - phg.edgeWeight(he));
    }
  }

  bool findNextMove(const PartitionedHypergraph &phg, Move &m) {
    if (!updatePQs(phg)) {
      return false;
    }
    ASSERT(!_blockPQ.empty());
    const PartitionID to = _blockPQ.top();
    ASSERT(!_toPQs[to].empty());
    const HypernodeID u = _toPQs[to].top();
    const Gain gain = -_toPQs[to].topKey();
    m.node = u;
    m.from = phg.partID(u);
    m.to = to;
    m.gain = gain;
    for (auto &pq : _toPQs) {
      pq.remove(u);
    }
    return true;
  }

private:
  bool updatePQs(const PartitionedHypergraph &phg) {
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      updateOrRemoveToPQFromBlocks(i, phg);
    }
    return !_blockPQ.empty();
  }

  void updateOrRemoveToPQFromBlocks(const PartitionID i,
                                    const PartitionedHypergraph &phg) {
    if (!_toPQs[i].empty()) {
      _blockPQ.insertOrAdjustKey(
          i, std::make_pair(phg.partLoad(i), _toPQs[i].topKey()));
    } else if (_blockPQ.contains(i)) {
      _blockPQ.remove(i);
    }
  }

  Gain computeGainForInvalidFrom(const PartitionedHypergraph &phg,
                                 const HypernodeID u, const PartitionID to) {
    Gain penalty = 0;
    for (const auto &he : phg.incidentEdges(u)) {
      const HyperedgeWeight pin_count_in_to = phg.pinCountInPart(he, to);
      if (pin_count_in_to == 0) {
        penalty += phg.edgeWeight(he);
      }
    }
    penalty += phg.weightOfDisabledEdges(u);
    return penalty;
  }

  // std::pair<PartitionID, HyperedgeWeight>
  // computeBestTargetBlock(const PartitionedHypergraph &phg, const HypernodeID
  // u,
  //                        const PartitionID from) {
  //   const HyperedgeWeight from_load =
  //       from != kInvalidPartition ? phg.partLoad(from) : 0;
  //   PartitionID to = kInvalidPartition;
  //   HyperedgeWeight to_load = std::numeric_limits<HyperedgeWeight>::max();
  //   HyperedgeWeight to_load_after =
  //   std::numeric_limits<HyperedgeWeight>::max(); HyperedgeWeight to_penalty =
  //   std::numeric_limits<HyperedgeWeight>::max(); Gain benefit = (from !=
  //   kInvalidPartition ? phg.moveFromBenefit(u) : 0) +
  //                  phg.weightOfDisabledEdges(u);
  //   HyperedgeWeight from_load_after = from_load - benefit;
  //   bool found_positive_target = false;
  //   for (PartitionID i = 0; i < _context.partition.k; ++i) {
  //     if (i != from) {
  //       const HyperedgeWeight load = phg.partLoad(i);
  //       HyperedgeWeight penalty = 0;
  //       if (from == kInvalidPartition) {
  //         for (const auto &he : phg.incidentEdges(u)) {
  //           if (phg.pinCountInPart(he, i) == 0) {
  //             penalty += phg.edgeWeight(he);
  //           }
  //         }
  //       } else {
  //         penalty += phg.moveToPenalty(u, i);
  //       }
  //       HyperedgeWeight load_after = load + penalty;
  //       if (load_after < from_load_after && penalty < to_penalty) {
  //         to = i;
  //         to_load = load;
  //         to_load_after = load_after;
  //         to_penalty = penalty;
  //         found_positive_target = true;
  //       } else if (!found_positive_target &&
  //                  (load_after < to_load_after ||
  //                   (load_after == to_load_after && load < to_load))) {
  //         to = i;
  //         to_load = load;
  //         to_load_after = load_after;
  //         to_penalty = penalty;
  //       }
  //     }
  //   }
  //   ASSERT(to != kInvalidPartition);
  //   to_load_after += phg.weightOfDisabledEdges(u);
  //   const Gain gain =
  //       found_positive_target ? benefit : from_load - to_load_after;
  //
  //   return std::make_pair(to, gain);
  // }

  // Gain computeGain(const PartitionedHypergraph &phg, const HypernodeID u,
  //                         const PartitionID from, const PartitionID to, const
  //                         PartitionID heaviestPart) {
  //   HyperedgeWeight highest_load = phg.partLoad(heaviestPart);
  //   Gain gain = 0;
  //   Gain penalty = 0;
  //   if (from == kInvalidPartition) {
  //     for (const auto &he : phg.incidentEdges(u)) {
  //       const HyperedgeWeight pin_count_in_to_after =
  //           phg.pinCountInPart(he, to);
  //       if (pin_count_in_to_after == 0) {
  //         penalty += phg.edgeWeight(he);
  //       }
  //     }
  //   } else {
  //     penalty = phg.moveToPenalty(u, he);
  //   }
  //   penalty += phg.weightOfDisabledEdges(u);
  //   HyperedgeWeight to_load_after = partLoad(to) + penalty;
  //   if (from != kInvalidPartition && to_load_after > highest_load) {
  //     gain = phg.moveFromBenefit(u) + phg.weightOfDisabledEdges(u);
  //   } else {
  //     gain = -penalty;
  //   }
  //   return gain;
  // }

  const Context &_context;
  vec<PriorityQueue> _toPQs;
  BlockPQ _blockPQ;
};
} // namespace mt_kahypar

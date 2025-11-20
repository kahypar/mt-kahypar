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

#include "mt-kahypar/partition/refinement/rebalancing/fallback.h"

#include <set>

#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/parallel/scalable_sort.h"
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer_common.h"
#include "mt-kahypar/weight/hypernode_weight_common.h"

namespace mt_kahypar {

namespace impl {
  // factoring out this function seems to significantly improve compile time,
  // probably because the calling function becomes overly complex if it is inlined
  MT_KAHYPAR_ATTRIBUTE_NO_INLINE void scalableSortMoves(vec<rebalancer::PotentialMove>& moves) {
    parallel::scalable_sort(moves, [](const rebalancer::PotentialMove& a, const rebalancer::PotentialMove& b) {
      return a.rating > b.rating;
    });
  }
}

namespace rebalancer {

static constexpr MoveID kInvalidMove = std::numeric_limits<MoveID>::max();

template <typename GraphAndGainTypes>
std::pair<int64_t, size_t> Fallback<GraphAndGainTypes>::runDeadlockFallback(PartitionedHypergraph& phg,
                                                                            GainCache& gain_cache,
                                                                            const Context& context,
                                                                            vec<ds::StreamingVector<rebalancer::PotentialMove>>& tmp_potential_moves,
                                                                            ds::Array<Move>& moves,
                                                                            ds::Array<MoveID>& move_id_of_node,
                                                                            ds::Array<uint8_t>& node_is_locked,
                                                                            const vec<double>& weight_normalizer,
                                                                            size_t& global_move_id) {
  tmp_potential_moves.resize(context.partition.k);
  for (auto& moves: tmp_potential_moves) {
    moves.clear_sequential();
  }

  HypernodeWeightArray max_weight_per_block(context.partition.k, phg.dimension(), 0, false);
  for (PartitionID block = 0; block < context.partition.k; ++block) {
    max_weight_per_block[block] = context.refinement.rebalancing.fallback_weight_threshold * context.partition.max_part_weights[block];
  }

  // compute the extremal dimensions and weight normalizers of all blocks, to speed up later calculations
  vec<uint8_t> is_max_block_dimension(context.partition.k * phg.dimension(), static_cast<bool>(false));
  vec<uint8_t> is_min_block_dimension(context.partition.k * phg.dimension(), static_cast<bool>(false));
  vec<vec<double>> block_weight_normalizers(context.partition.k, vec<double>{});
  for (PartitionID block = 0; block < context.partition.k; ++block) {
    vec<double> normalizer(phg.dimension(), 0);
    for (Dimension d = 0; d < phg.dimension(); ++d) {
      normalizer[d] = 1 / static_cast<double>(context.partition.max_part_weights[block].at(d));
    }
    const HNWeightConstRef block_weight = weight::toNonAtomic(phg.partWeight(block));
    impl::getExtremalDimensions(block_weight, normalizer,
                                &is_max_block_dimension[block * phg.dimension()], true);
    impl::getExtremalDimensions(block_weight, normalizer,
                                &is_min_block_dimension[block * phg.dimension()], false);
    block_weight_normalizers[block] = std::move(normalizer);
  }
  // it can happen (due to degree 0 nodes or individual block weights) that some min dimension is not present for any block
  // ---> compute block with most fitting imbalance instead
  for (Dimension d = 0; d < phg.dimension(); ++d) {
    bool has_block_target = false;
    for (PartitionID block = 0; block < context.partition.k; ++block) {
      if (is_min_block_dimension[block * phg.dimension() + d]) {
        has_block_target = true;
        break;
      }
    }
    if (!has_block_target) {
      float smallest_ratio = 2;
      PartitionID best_block = kInvalidPartition;
      for (PartitionID block = 0; block < context.partition.k; ++block) {
        const HNWeightConstRef block_weight = weight::toNonAtomic(phg.partWeight(block));
        float summed_weight = impl::normalizedSum(block_weight, block_weight_normalizers[block]);
        float weight_of_dim = block_weight_normalizers[block].at(d) * block_weight.at(d);
        float ratio = weight_of_dim / summed_weight;
        if (ratio < smallest_ratio) {
          best_block = block;
          smallest_ratio = ratio;
        }
      }
      ASSERT(best_block != kInvalidPartition);
      is_min_block_dimension[best_block * phg.dimension() + d] = static_cast<uint8_t>(true);
    }
  }

  auto compute_best_target = [&](const HypernodeID hn, const PartitionID from, const HNWeightConstRef weight, const vec<uint8_t>& max_dimensions, bool skip_blocks) {
    PartitionID best_to_part = kInvalidPartition;
    float best_rating = std::numeric_limits<float>::min();
    HyperedgeWeight best_gain = std::numeric_limits<HyperedgeWeight>::min();
    const HyperedgeWeight penalty_term = gain_cache.penaltyTerm(hn, phg.partID(hn));
    HyperedgeWeight adjacent_weight = 0;
    if (context.refinement.rebalancing.fallback_relative_block_priority) {
      for (HyperedgeID he: phg.incidentEdges(hn)) {
        adjacent_weight += phg.edgeWeight(he);
      }
    }

    for (PartitionID to = 0; to < context.partition.k; ++to) {
      if (from == to) continue;
      if (skip_blocks && !impl::hasMatchingDimension(weight, max_dimensions.data(), &is_min_block_dimension[to * phg.dimension()])) continue;

      float rating = 1.0;
      const HNWeightConstRef to_weight = weight::toNonAtomic(phg.partWeight(to));
      switch (context.refinement.rebalancing.fallback_block_selection) {
        case RbFallbackBlockSelectionPolicy::any_fitting_min_dimension: break;
        case RbFallbackBlockSelectionPolicy::by_internal_imbalance: {
          float matchingBlockWeight = impl::weightOfMatchingDimension(to_weight, max_dimensions.data(), block_weight_normalizers[to]);
          float weight_sum = impl::normalizedSum(to_weight, block_weight_normalizers[to]);
          rating = (weight_sum - matchingBlockWeight) / (0.01 * weight_sum + matchingBlockWeight);
          break;
        }
        case RbFallbackBlockSelectionPolicy::by_dot_product: {
          rating = std::max<float>(impl::dotProduct(weight, context.partition.max_part_weights[to] - to_weight, block_weight_normalizers[to]), 0);
          break;
        }
        case RbFallbackBlockSelectionPolicy::by_progress: {
          rating = impl::computeBalanceProgress(weight, weight::toNonAtomic(phg.partWeight(from)), context.partition.max_part_weights[from],
                                                to_weight.load(std::memory_order_relaxed), context.partition.max_part_weights[to], weight_normalizer).first;
          float sum = impl::normalizedSum(weight, weight_normalizer);
          rating = (rating + sum) / sum;
          break;
        }
      };

      HyperedgeWeight gain = 0;
      if (context.refinement.rebalancing.fallback_relative_block_priority || rating == best_rating) {
        HyperedgeWeight benefit;
        if (gain_cache.blockIsAdjacent(hn, to)) {
          benefit = gain_cache.benefitTerm(hn, to);
        } else {
          benefit = gain_cache.recomputeBenefitTerm(phg, hn, to);
        }
        gain = benefit - penalty_term;
      }
      if (context.refinement.rebalancing.fallback_relative_block_priority) {
        ASSERT(adjacent_weight + gain >= 0);
        rating *= (2 * adjacent_weight + gain);
      }
      ASSERT(rating >= 0);
      if (rating > best_rating || (rating == best_rating && gain > best_gain)) {
        best_to_part = to;
        best_rating = rating;
        best_gain = gain;
      }
    }
    ASSERT(best_to_part != kInvalidPartition);
    return std::make_pair(best_to_part, best_rating);
  };

  // compute moves for nodes
  tbb::enumerable_thread_specific<vec<uint8_t>> local_max_dimensions(phg.dimension(), static_cast<bool>(false));
  phg.doParallelForAllNodes([&](const HypernodeID hn) {
    const PartitionID from = phg.partID(hn);
    const HNWeightConstRef from_weight = weight::toNonAtomic(phg.partWeight(from));
    const HNWeightConstRef weight = phg.nodeWeight(hn);
    if (from_weight <= context.partition.max_part_weights[from]) return;
    if (!(weight <= max_weight_per_block[from])) return;

    auto& max_dimensions = local_max_dimensions.local();
    max_dimensions.assign(phg.dimension(), static_cast<bool>(false));
    impl::getExtremalDimensions(weight, block_weight_normalizers[from], max_dimensions.data(), true);
    if (!impl::hasMatchingDimension(weight, max_dimensions.data(), &is_max_block_dimension[from * phg.dimension()])) return;

    auto [to_part, rating] = compute_best_target(hn, from, weight, max_dimensions, true);
    ASSERT(to_part != kInvalidPartition);

    float node_rating = 1.0;
    switch (context.refinement.rebalancing.fallback_node_selection) {
      case RbFallbackNodeSelectionPolicy::any_fitting_max_dimension: break;
      case RbFallbackNodeSelectionPolicy::by_internal_imbalance: {
        float matchingNodeWeight = impl::weightOfMatchingDimension(weight, max_dimensions.data(), block_weight_normalizers[from]);
        node_rating = matchingNodeWeight / (1.01 * impl::normalizedSum(weight, block_weight_normalizers[from]) - matchingNodeWeight);
        break;
      }
      case RbFallbackNodeSelectionPolicy::by_dot_product: {
        node_rating = impl::dotProduct(weight, from_weight, block_weight_normalizers[from]);
        break;
      }
    };
    if (context.refinement.rebalancing.fallback_node_priority_by_weight) {
      const auto diff_to_block = weight::max(from_weight - context.partition.max_part_weights[from], weight::broadcast(0, phg.dimension()));
      float penalty = std::max(impl::normalizedSum(weight, block_weight_normalizers[from]), impl::normalizedSum(diff_to_block, block_weight_normalizers[from]));
      rating /= penalty;
    }
    if (context.refinement.rebalancing.fallback_relative_node_priority) {
      rating *= node_rating;
    } else {
      rating = node_rating;
    }
    tmp_potential_moves[from].stream(rebalancer::PotentialMove{hn, to_part, rating});
  });

  ds::StreamingVector<rebalancer::PotentialMove> all_moves;
  tbb::parallel_for(0, context.partition.k, [&](const PartitionID from) {
    if (tmp_potential_moves[from].size() > 0) {
      vec<rebalancer::PotentialMove> moves_of_block = tmp_potential_moves[from].copy_parallel();

      // sort the moves from each overweight part by rating
      impl::scalableSortMoves(moves_of_block);

      // select prefix of moves that should be applied
      AllocatedHNWeight from_weight;
      from_weight = phg.partWeight(from);
      auto should_stop = [&]() {
        switch (context.refinement.rebalancing.fallback_node_count) {
          case RbFallbackNodeCountPolicy::only_one:
            return true;
          case RbFallbackNodeCountPolicy::until_balanced:
            return from_weight <= context.partition.max_part_weights[from];
          case RbFallbackNodeCountPolicy::below_threshold:
            return from_weight <= context.refinement.rebalancing.fallback_node_count_threshold * context.partition.max_part_weights[from];
        };
        return true;
      };
      for (const auto& m: moves_of_block) {
        ASSERT(phg.partID(m.node) == from);
        from_weight -= phg.nodeWeight(m.node);
        all_moves.stream(m);

        if (should_stop()) break;
      }
    }
  });
  vec<rebalancer::PotentialMove> move_list = all_moves.copy_parallel();
  impl::scalableSortMoves(move_list);
  // ASSERT(!move_list.empty());  --> weight threshold
  DBG << "Applying" << move_list.size() << "moves to break deadlock";

  ASSERT([&]{
    std::set<HypernodeID> move_set;
    for (const auto& m: move_list) {
      if (move_set.find(m.node) != move_set.end()) {
        LOG << "Node" << m.node << "appears twice in move set!";
        return false;
      }
      move_set.insert(m.node);
    }
    return true;
  }());

  if (context.refinement.rebalancing.fallback_use_locking) {
    node_is_locked.assign(phg.initialNumNodes(), static_cast<uint8_t>(false), true);
  }

  const size_t old_id = global_move_id;
  int64_t attributed_gain = 0;
  for (const auto& m: move_list) {
    const PartitionID from = phg.partID(m.node);
    if (phg.partWeight(from) <= context.partition.max_part_weights[from]) continue;

    // recompute target
    const HNWeightConstRef weight = phg.nodeWeight(m.node);
    auto& max_dimensions = local_max_dimensions.local();
    max_dimensions.assign(phg.dimension(), static_cast<bool>(false));
    impl::getExtremalDimensions(weight, block_weight_normalizers[from], max_dimensions.data(), true);
    auto [best_to_part, _] = compute_best_target(m.node, from, weight, max_dimensions, false);
    ASSERT(best_to_part != kInvalidPartition);

    DBG << V(weight) << V(from) << V(phg.partWeight(from)) << V(best_to_part) << V(phg.partWeight(best_to_part));

    int64_t gain = 0;
    phg.changeNodePart(gain_cache, m.node, from, best_to_part,
      [&](const SynchronizedEdgeUpdate& sync_update) {
        gain = AttributedGains::gain(sync_update);
        attributed_gain += gain;
      });
    if (move_id_of_node[m.node] == kInvalidMove) {
      move_id_of_node[m.node] = global_move_id;
    }
    moves[global_move_id++] = Move{from, best_to_part, m.node, static_cast<Gain>(gain)};
    if (context.refinement.rebalancing.fallback_use_locking) {
      node_is_locked[m.node] = static_cast<uint8_t>(true);
    }
  }
  return {attributed_gain, global_move_id - old_id};
}

// explicitly instantiate so the compiler can generate them when compiling this cpp file
namespace {
  #define FALLBACK(X) Fallback<X>
}

// explicitly instantiate so the compiler can generate them when compiling this cpp file
INSTANTIATE_CLASS_WITH_VALID_TRAITS(FALLBACK)

} // namespace rebalancer
}  // namespace mt_kahypar

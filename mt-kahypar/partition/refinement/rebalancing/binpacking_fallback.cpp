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

#include "mt-kahypar/partition/refinement/rebalancing/binpacking_fallback.h"

#include <algorithm>
#include <atomic>
#include <iostream>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/buffered_vector.h"
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer_common.h"
#include "mt-kahypar/utils/exception.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
namespace bp {

static constexpr bool debug = false;

enum class BPAlgorithm: uint8_t {
  first_fit = 0,
  best_fit = 1,
  worst_fit = 2,
  random_fit = 3,
  old_block_or_best = 4,
  best_by_dot_product = 5,
  best_by_internal_imbalance = 6,
  UNDEFINED = 7
};

std::ostream& operator<< (std::ostream& os, const BPAlgorithm& algo) {
  switch (algo) {
    case BPAlgorithm::first_fit: return os << "first_fit";
    case BPAlgorithm::best_fit: return os << "best_fit";
    case BPAlgorithm::worst_fit: return os << "worst_fit";
    case BPAlgorithm::random_fit: return os << "random_fit";
    case BPAlgorithm::old_block_or_best: return os << "old_block_or_best";
    case BPAlgorithm::best_by_dot_product: return os << "best_by_dot_product";
    case BPAlgorithm::best_by_internal_imbalance: return os << "best_by_internal_imbalance";
    case BPAlgorithm::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(algo);
}

struct BPResult {
  BPResult() = default;

  BPResult(const BPResult&) = delete;
  BPResult& operator=(const BPResult&) = delete;

  BPResult(BPResult&&) = default;
  BPResult& operator=(BPResult&&) = default;

  bool applyBetterResult(const Context& context,
                         const vec<RebalancingNode>& other_assignment,
                         double other_rating,
                         HyperedgeWeight other_loss,
                         BPAlgorithm other_algo,
                         bool other_permuted) {
    unused(context);
    unused(other_loss);
    if (assignment.empty()
        || other_rating < balance_rating
        || (other_rating == balance_rating && utils::Randomize::instance().flipCoin(THREAD_ID))) {
      assignment.resize(other_assignment.size());
      std::copy(other_assignment.cbegin(), other_assignment.cend(), assignment.begin());
      balance_rating = other_rating;
      estimated_quality_loss = other_loss;
      algo = other_algo;
      permuted = other_permuted;
      return true;
    } else {
      return false;
    }
  }

  bool applyBetterResult(const Context& context,
                         const BPResult& other) {
    return applyBetterResult(context, other.assignment, other.balance_rating,
                             other.estimated_quality_loss, other.algo, other.permuted);
  }

  bool isValid() const {
    return !assignment.empty();
  }

  void reset() {
    assignment.clear();
    balance_rating = 0;
    estimated_quality_loss = 0;
    algo = BPAlgorithm::UNDEFINED;
    permuted = false;
  }

  vec<RebalancingNode> assignment;
  double balance_rating = 0;
  HyperedgeWeight estimated_quality_loss = 0;
  BPAlgorithm algo = BPAlgorithm::UNDEFINED;
  bool permuted = false;
};

std::ostream& operator<< (std::ostream& os, const BPResult& result) {
  os << "[" << result.algo;
  if (result.permuted) {
    os << " (permute)";
  }
  return os << ": " << result.balance_rating << "]";
}

struct BlockInfo {
  PartitionID id;
  HNWeightConstRef current_weight;
  HNWeightConstRef max_weight;
  const vec<double>* weight_normalizer = nullptr;
};

struct BinPacker {
  BinPacker(const Context& context, const vec<vec<double>>& block_weight_normalizers):
    target_bound(0),
    max_fill_grade(context.dimension(), 0),
    block_order(),
    block_weights(context.partition.k, context.dimension(), 0, false),
    max_part_weights(context.partition.max_part_weights),
    block_weight_normalizers(block_weight_normalizers) {
    for (PartitionID block = 0; block < context.partition.k; ++block) {
      block_order.emplace_back(block);
    }
  }

  BinPacker(const BinPacker&) = delete;
  BinPacker& operator=(const BinPacker&) = delete;

  Dimension dimension() const {
    return block_weights.dimension();
  }

  void initialize(double bound) {
    ASSERT(bound <= 1 && bound > 0);
    target_bound = bound;
    max_fill_grade.assign(dimension(), 0);
    block_weights.assign(block_weights.size(), 0, false);
    utils::Randomize::instance().shuffleVector(block_order);
  }

  bool isEmpty(PartitionID block) const {
    ASSERT(block != kInvalidPartition && UL(block) < block_weights.size());
    return weight::isZero(block_weights[block]);
  }

  bool fits(const RebalancingNode& node, PartitionID to, bool is_empty) const {
    ASSERT(target_bound > 0);
    ASSERT(to != kInvalidPartition && UL(to) < max_part_weights.size());
    if (is_empty) {
      // TODO: individual part weights...
      return true;
    }

    const auto new_weight = block_weights[to] + node.weight;
    const auto limit = target_bound * max_part_weights[to];
    for (Dimension d = 0; d < dimension(); ++d) {
      if (node.weight.at(d) > 0 && new_weight.at(d) > limit.at(d)) {
        return false;
      }
    }
    return true;
  }

  bool fits(const RebalancingNode& node, PartitionID to) const {
    return fits(node, to, isEmpty(to));
  }

  bool assignNode(RebalancingNode& node, PartitionID to) {
    ASSERT(target_bound > 0);
    if (!fits(node, to)) return false;

    assignFitting(node, to);
    return true;
  }

  void assignFitting(RebalancingNode& node, PartitionID to) {
    ASSERT(fits(node, to));

    const bool is_empty = isEmpty(to);
    block_weights[to] += node.weight;
    node.to = to;
    if (!is_empty) {
      // update max_fill_grade
      for (Dimension d = 0; d < dimension(); ++d) {
        if (node.weight.at(d) > 0) {
          double fill_grade = block_weight_normalizers[to][d] * block_weights[to].at(d);
          max_fill_grade[d] = std::max(max_fill_grade[d], fill_grade);
        }
      }
    }
  }

  template<typename F>
  bool assignByRating(RebalancingNode& node, const F& rating_fn) {
    PartitionID best_block = kInvalidPartition;
    float best_rating = 0;
    bool has_seen_empty = false;
    for (PartitionID block: block_order) {
      const bool is_empty = isEmpty(block);
      if (is_empty && has_seen_empty && block != node.from) continue;
      has_seen_empty = is_empty;

      if (fits(node, block, is_empty)) {
        BlockInfo current_info{block, block_weights[block], max_part_weights[block], &block_weight_normalizers[block]};
        float rating = rating_fn(current_info);
        if (best_block == kInvalidPartition || rating < best_rating || (rating == best_rating && block == node.from)) {
          best_block = block;
          best_rating = rating;
        }
      }
    }
    if (best_block != kInvalidPartition) {
      assignFitting(node, best_block);
      return true;
    } else {
      return false;
    }
  }

  double target_bound;
  vec<double> max_fill_grade;
  vec<PartitionID> block_order;
  HypernodeWeightArray block_weights;
  const HypernodeWeightArray& max_part_weights;
  const vec<vec<double>>& block_weight_normalizers;
};

double computeBalanceRating(const Context& context, const BinPacker& result) {
  unused(context);
  double squares = 0;
  for (double val: result.max_fill_grade) {
    squares += val * val;
  }
  return squares;
}

bool assignNextNode(BinPacker& packer, RebalancingNode& node, BPAlgorithm algo, utils::Randomize& rand) {
  if (algo == BPAlgorithm::random_fit) {
    const int cpu_id = THREAD_ID;
    const size_t max_attempts = packer.block_order.size() + 1;
    for (size_t i = 0; i < max_attempts; ++i) {
      PartitionID to = rand.getRandomInt(0, packer.block_order.size() - 1, cpu_id);
      if (packer.assignNode(node, to)) return true;
    }
    return false;
  }
  if (algo == BPAlgorithm::old_block_or_best) {
    if (packer.assignNode(node, node.from)) return true;
  }

  const auto compute_max = [](const RebalancingNode& node, const BlockInfo& info) {
    float max = 0;
    for (Dimension d = 0; d < node.weight.dimension(); ++d) {
      max = std::max(max, static_cast<float>(info.weight_normalizer->at(d)) * static_cast<float>(info.current_weight.at(d) + node.weight.at(d)));
    }
    return max;
  };

  switch (algo) {
    // Note: small rating is good
    case BPAlgorithm::first_fit:
      return packer.assignByRating(node, [&](const BlockInfo&) {
        return 1.0;
      });
    case BPAlgorithm::old_block_or_best:  // fallthrough
    case BPAlgorithm::best_fit:
      return packer.assignByRating(node, [&](const BlockInfo& current) {
        return compute_max(node, current);
      });
    case BPAlgorithm::worst_fit:
      return packer.assignByRating(node, [&](const BlockInfo& current) {
        return -compute_max(node, current);
      });
    case BPAlgorithm::best_by_dot_product:
      return packer.assignByRating(node, [&](const BlockInfo& current) {
        // TODO: better normalization?
        return -impl::dotProduct(node.weight, current.max_weight - current.current_weight, *current.weight_normalizer);
      });
    case BPAlgorithm::best_by_internal_imbalance:
      return packer.assignByRating(node, [&](const BlockInfo& current) {
        return impl::internalImbalance(node.weight + current.current_weight, *current.weight_normalizer);
      });
    case BPAlgorithm::random_fit: throw InvalidParameterException("invalid bin packing algorithm");
    case BPAlgorithm::UNDEFINED: throw InvalidParameterException("invalid bin packing algorithm");
      // omit default case to trigger compiler warning for missing cases
  }
  throw InvalidParameterException("invalid bin packing algorithm");
}

bool computePackingFromNodeOrder(BinPacker& packer, vec<RebalancingNode>& nodes, const vec<double>& weight_normalizer, BPAlgorithm algo, bool permute) {
  auto& randomize = utils::Randomize::instance();
  std::sort(nodes.begin(), nodes.end(), [&](const RebalancingNode& lhs, const RebalancingNode& rhs) {
    return lhs.normalized_weight > rhs.normalized_weight;
  });

  if (permute) {
    const int cpu_id = THREAD_ID;
    for (size_t i = randomize.flipCoin(cpu_id) ? 0 : 1; i + 1 < nodes.size(); ++i) {
      if (randomize.flipCoin(cpu_id)) {
        std::swap(nodes[i], nodes[i + 1]);
      }
    }
  }

  for (RebalancingNode& node: nodes) {
    bool success = assignNextNode(packer, node, algo, randomize);
    if (!success) {
      // DBG << "Failed to compute packing: " << V(algo) << V(permute);
      return false;
    };
  }
  // if (debug) {
  //   std::cout << "Successfully computed packing with " << V(algo) << ", " << V(permute) << " - resulting fill grades: ";
  //   for (Dimension d = 0; d < packer.dimension(); ++d) {
  //     std::cout << packer.max_fill_grade[d] << " ";
  //   }
  //   std::cout << std::endl;
  // }
  return true;
}

bool computeSinglePacking(BinPacker& packer, vec<RebalancingNode>& nodes, const vec<double>& weight_normalizer, BPAlgorithm algo, double bound) {
  packer.initialize(bound);
  return computePackingFromNodeOrder(packer, nodes, weight_normalizer, algo, false);
}

bool computeManyPackings(const Context& context,
                         tbb::enumerable_thread_specific<BinPacker>& local_packer,
                         tbb::enumerable_thread_specific<BPResult>& local_best,
                         tbb::enumerable_thread_specific<vec<RebalancingNode>>& local_nodes,
                         BPResult& best_result,
                         const vec<double>& weight_normalizer,
                         double bound) {
  for (BPResult& result: local_best) {
    result.reset();
  }
  const size_t total_num_tasks = context.refinement.rebalancing.bin_packing_tasks * static_cast<size_t>(BPAlgorithm::UNDEFINED);
  vec<std::pair<BPAlgorithm, bool>> tasks;
  tasks.reserve(total_num_tasks);
  for (uint8_t algo = 0; algo < static_cast<uint8_t>(BPAlgorithm::UNDEFINED); ++algo) {
    for (size_t id = 0; id < context.refinement.rebalancing.bin_packing_tasks; ++id) {
      tasks.emplace_back(static_cast<BPAlgorithm>(algo), id != 0);
    }
  }

  std::atomic_bool overall_success = false;
  tbb::parallel_for(UL(0), tasks.size(), [&](const size_t task_id) {
    auto [algo, permute] = tasks[task_id];
    auto& packer = local_packer.local();
    auto& nodes = local_nodes.local();
    packer.initialize(bound);
    bool success = computePackingFromNodeOrder(packer, nodes, weight_normalizer, algo, permute);

    if (success) {
      const double balance_rating = computeBalanceRating(context, packer);
      auto& local_best_result = local_best.local();
      local_best_result.applyBetterResult(context, nodes, balance_rating, 0, algo, permute);
      overall_success.store(true, std::memory_order_relaxed);
    }
  });
  if (overall_success) {
    for (const BPResult& result: local_best) {
      best_result.applyBetterResult(context, result);
    }
  }
  DBG << CYAN << "Best result: " << best_result << END;
  return overall_success;
}

template<typename PartitionedHypergraph>
vec<RebalancingNode> determineNodesForRebalancing(PartitionedHypergraph& phg, const Context& context) {
  // TODO: perhaps this should depend on epsilon
  const double selection_factor = context.refinement.rebalancing.bin_packing_selection_threshold;
  const vec<vec<double>> block_weight_normalizers = impl::computeBlockWeightNormalizers(context);
  const size_t max_count = phg.dimension() * std::ceil((1.0 + context.partition.epsilon) / selection_factor * context.partition.k);

  ds::BufferedVector<RebalancingNode> result(std::min(max_count, static_cast<size_t>(phg.initialNumNodes())));
  phg.doParallelForAllNodes([&](const HypernodeID hn) {
    const HNWeightConstRef weight = phg.nodeWeight(hn);
    const PartitionID from = phg.partID(hn);
    const vec<double>& normalizer = block_weight_normalizers[from];
    for (Dimension d = 0; d < phg.dimension(); ++d) {
      if (normalizer[d] * weight.at(d) > selection_factor) {
        result.push_back_buffered(RebalancingNode(hn, weight, from));
        break;
      }
    }
  });
  result.finalize();
  auto data = result.getData();
  data.resize(result.size());
  // TODO: fallback if result is empty / too small?
  return data;
}

bool computeBinPacking(const Context& context, vec<RebalancingNode>& nodes, const vec<double>& weight_normalizer) {
  constexpr size_t max_search_steps = 8;
  ASSERT(!nodes.empty());

  // precompute the normalized weights
  tbb::parallel_for(UL(0), nodes.size(), [&](const size_t i) {
    nodes[i].normalized_weight = impl::normalizedSum(nodes[i].weight, weight_normalizer);
  });

  const vec<vec<double>> block_weight_normalizers = impl::computeBlockWeightNormalizers(context);
  tbb::enumerable_thread_specific<BinPacker> local_packer(context, block_weight_normalizers);
  tbb::enumerable_thread_specific<BPResult> local_best;
  tbb::enumerable_thread_specific<vec<RebalancingNode>> local_nodes(nodes);
  BPResult global_best_result;

  double upper_bound = 2.0;  // make first step at 1.0, reducing special casing
  double lower_bound = 0;
  double best_bound = 2.0;
  bool best_bound_was_multi_packed = false;
  for (size_t i = 0; i < max_search_steps; ++i) {
    const double current_bound = (upper_bound + lower_bound) / 2;
    const bool skip_single_packing = current_bound < best_bound && best_bound_was_multi_packed;
    DBG << "Binary search step: " << V(current_bound);

    bool success = false;
    if (!skip_single_packing) {
      auto& current_packer = local_packer.local();
      auto& current_nodes = local_nodes.local();
      const BPAlgorithm algo = BPAlgorithm::best_fit;
      success = computeSinglePacking(current_packer, current_nodes, weight_normalizer, algo, current_bound);
      if (success) {
        const double balance_rating = computeBalanceRating(context, current_packer);
        global_best_result.applyBetterResult(context, current_nodes, balance_rating, 0, algo, false);
        best_bound_was_multi_packed = false;
      }
    }
    if (!success) {
      success = computeManyPackings(context, local_packer, local_best, local_nodes,
                                    global_best_result, weight_normalizer, current_bound);
      if (success) {
        best_bound_was_multi_packed = true;
      }
    }
    if (success) {
      ASSERT(current_bound < best_bound);
      best_bound = current_bound;
      upper_bound = current_bound;
    } else if (i == 0) {
      // we couldn't find a packing within the max part weight
      return false;
    } else {
      lower_bound = current_bound;
    }
  }

  if (!best_bound_was_multi_packed) {
    DBG << "Global best result before multi packing: " << global_best_result;
    computeManyPackings(context, local_packer, local_best, local_nodes,
                        global_best_result, weight_normalizer, best_bound);
  }
  DBG << CYAN << "Global best result: " << global_best_result << END;

  nodes = std::move(global_best_result.assignment);
  return true;
}

namespace {
#define DETERMINE_NODES_FOR_REBALANCING(X) vec<RebalancingNode> determineNodesForRebalancing(X& phg, const Context& context)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(DETERMINE_NODES_FOR_REBALANCING)

}
}  // namespace mt_kahypar

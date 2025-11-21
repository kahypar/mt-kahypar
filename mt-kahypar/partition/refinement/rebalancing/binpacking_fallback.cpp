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
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
namespace bp {

static constexpr bool debug = true;

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
  }

  vec<RebalancingNode> assignment;
  double balance_rating;
  HyperedgeWeight estimated_quality_loss;
  BPAlgorithm algo = BPAlgorithm::UNDEFINED;
  bool permuted;
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
  const vec<double>* weight_normalizer;
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

  bool fits(const RebalancingNode& node, PartitionID to) const {
    ASSERT(target_bound > 0);
    ASSERT(to != kInvalidPartition && UL(to) < max_part_weights.size());
    if (isEmpty(to)) {
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

  bool assignNode(RebalancingNode& node, PartitionID to) {
    ASSERT(target_bound > 0);
    if (!fits(node, to)) return false;

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
    return true;
  }

  template<typename F>
  bool assignByComparator(RebalancingNode& node, F compare_fn) {
    PartitionID best_block = kInvalidPartition;
    BlockInfo best_info;
    for (PartitionID block: block_order) {
      if (fits(node, block)) {
        BlockInfo current_info{block, block_weights[block], max_part_weights[block], &block_weight_normalizers[block]};
        if (best_block == kInvalidPartition || compare_fn(best_info, current_info)) {
          best_block = block;
          best_info = current_info;
        }
      }
    }
    if (best_block != kInvalidPartition) {
      return assignNode(node, best_block);
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

bool assignNextNode(BinPacker& packer, RebalancingNode& node, BPAlgorithm algo) {
  return packer.assignByComparator(node, [&](const BlockInfo& best, const BlockInfo& current) {
    const auto best_new_weight = best.current_weight + node.weight;
    const auto current_new_weight = current.current_weight + node.weight;
    double best_max = 0;
    double current_max = 0;
    for (Dimension d = 0; d < node.weight.dimension(); ++d) {
      best_max = std::max(best_max, best.weight_normalizer->at(d) * best_new_weight.at(d));
      current_max = std::max(current_max, current.weight_normalizer->at(d) * current_new_weight.at(d));
    }
    ASSERT(best_max <= 1 && current_max <= 1);
    return current_max < best_max || (current_max == best_max && current.id == node.from);
  });
}

bool computePackingFromNodeOrder(BinPacker& packer, vec<RebalancingNode>& nodes, const vec<double>& weight_normalizer, BPAlgorithm algo, bool permute) {
  std::sort(nodes.begin(), nodes.end(), [&](const RebalancingNode& lhs, const RebalancingNode& rhs) {
    return impl::normalizedSum(lhs.weight, weight_normalizer) > impl::normalizedSum(rhs.weight, weight_normalizer);
  });

  if (permute) {
    auto& randomize = utils::Randomize::instance();
    const int cpu_id = THREAD_ID;
    for (size_t i = randomize.flipCoin(cpu_id) ? 0 : 1; i + 1 < nodes.size(); ++i) {
      if (randomize.flipCoin(cpu_id)) {
        std::swap(nodes[i], nodes[i + 1]);
      }
    }
  }

  for (RebalancingNode& node: nodes) {
    bool success = assignNextNode(packer, node, algo);
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
  for (const BPResult& result: local_best) {
    best_result.applyBetterResult(context, result);
  }
  DBG << CYAN << "Best result: " << best_result << END;
  return overall_success;
}

template<typename PartitionedHypergraph>
vec<RebalancingNode> determineNodesForRebalancing(PartitionedHypergraph& phg, const Context& context) {
  constexpr double selection_factor = 0.02;
  const vec<vec<double>> block_weight_normalizers = impl::computeBlockWeightNormalizers(context);

  ds::BufferedVector<RebalancingNode> result(std::ceil(1.0 / selection_factor * context.partition.k));
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

  const vec<vec<double>> block_weight_normalizers = impl::computeBlockWeightNormalizers(context);
  tbb::enumerable_thread_specific<BinPacker> local_packer(context, block_weight_normalizers);
  tbb::enumerable_thread_specific<BPResult> local_best;
  tbb::enumerable_thread_specific<vec<RebalancingNode>> local_nodes(nodes);
  BPResult global_best_result;
  BPResult last_best_result;
  BPResult current_best_result;

  double upper_bound = 2.0;  // make first step at 1.0, reducing special casing
  double lower_bound = 0;
  double best_bound = 2.0;
  bool best_bound_was_multi_packed = false;
  for (size_t i = 0; i < max_search_steps; ++i) {
    const double current_bound = (upper_bound + lower_bound) / 2;
    const bool skip_single_packing = current_bound < best_bound && best_bound_was_multi_packed;
    current_best_result.reset();
    DBG << "Binary search step: " << V(current_bound);

    bool success = false;
    if (!skip_single_packing) {
      auto& current_packer = local_packer.local();
      auto& current_nodes = local_nodes.local();
      const BPAlgorithm algo = BPAlgorithm::best_fit;
      success = computeSinglePacking(current_packer, current_nodes, weight_normalizer, algo, current_bound);
      if (success) {
        const double balance_rating = computeBalanceRating(context, current_packer);
        bool check = current_best_result.applyBetterResult(context, current_nodes, balance_rating, 0, algo, false);
        ASSERT(check);
        best_bound_was_multi_packed = false;
      }
    }
    if (!success) {
      success = computeManyPackings(context, local_packer, local_best, local_nodes,
                                    current_best_result, weight_normalizer, current_bound);
      if (success) {
        best_bound_was_multi_packed = true;
      }
    }
    if (success) {
      ASSERT(current_bound < best_bound);
      best_bound = current_bound;
      upper_bound = current_bound;
      global_best_result.applyBetterResult(context, current_best_result);
      last_best_result = std::move(current_best_result);
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
  DBG << "Last best result: " << last_best_result;
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

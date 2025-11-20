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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/buffered_vector.h"
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer_common.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
namespace bp {

static constexpr bool debug = true;

struct BlockInfo {
  PartitionID id;
  HNWeightConstRef current_weight;
  HNWeightConstRef max_weight;
  const vec<double>* weight_normalizer;
};

struct BinPacker {
  BinPacker(const Context& context, const vec<vec<double>>& block_weight_normalizers):
    block_order(),
    block_weights(context.partition.k, context.dimension(), 0, false),
    max_part_weights(context.partition.max_part_weights),
    block_weight_normalizers(block_weight_normalizers) {
    for (PartitionID block = 0; block < context.partition.k; ++block) {
      block_order.emplace_back(block);
    }
    utils::Randomize::instance().shuffleVector(block_order);
  }

  bool fits(const RebalancingNode& node, PartitionID to) const {
    ASSERT(to != kInvalidPartition && UL(to) < max_part_weights.size());
    return block_weights[to] + node.weight <= max_part_weights[to];
  }

  bool assignNode(RebalancingNode& node, PartitionID to) {
    if (!fits(node, to)) return false;

    block_weights[to] += node.weight;
    node.to = to;
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

  vec<PartitionID> block_order;
  HypernodeWeightArray block_weights;
  const HypernodeWeightArray& max_part_weights;
  const vec<vec<double>>& block_weight_normalizers;
};

bool computeSinglePacking(BinPacker& packer, vec<RebalancingNode>& nodes, const vec<double>& weight_normalizer) {
  std::sort(nodes.begin(), nodes.end(), [&](const RebalancingNode& lhs, const RebalancingNode& rhs) {
    return impl::normalizedSum(lhs.weight, weight_normalizer) > impl::normalizedSum(rhs.weight, weight_normalizer);
  });
  for (RebalancingNode& node: nodes) {
    DBG << "Trying to pack node" << node.id << "with weight" << node.weight;
    bool success = packer.assignByComparator(node, [&](const BlockInfo& best, const BlockInfo& current) {
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
    if (!success) return false;
  }
  return true;
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
  return data;
}

bool computeBinPacking(const Context& context, vec<RebalancingNode>& nodes, const vec<double>& weight_normalizer) {
  const vec<vec<double>> block_weight_normalizers = impl::computeBlockWeightNormalizers(context);
  BinPacker packer(context, block_weight_normalizers);
  return computeSinglePacking(packer, nodes, weight_normalizer);
}

namespace {
#define DETERMINE_NODES_FOR_REBALANCING(X) vec<RebalancingNode> determineNodesForRebalancing(X& phg, const Context& context)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(DETERMINE_NODES_FOR_REBALANCING)

}
}  // namespace mt_kahypar

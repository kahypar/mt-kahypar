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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/buffered_vector.h"
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer_common.h"

namespace mt_kahypar {
namespace bp {

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

bool computeBinPacking(const Context& context, vec<RebalancingNode>& nodes) {
  return false;
}

namespace {
#define DETERMINE_NODES_FOR_REBALANCING(X) vec<RebalancingNode> determineNodesForRebalancing(X& phg, const Context& context)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(DETERMINE_NODES_FOR_REBALANCING)

}
}  // namespace mt_kahypar

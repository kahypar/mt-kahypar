/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2026 Nikolai Maas <nikolai.maas@student.kit.edu>
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

#include "mt-kahypar/macros.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/weight/hypernode_weight_common.h"

namespace mt_kahypar {

namespace bipartitioning {

struct OriginalHypergraphInfo {

  // The initial allowed imbalance cannot be used for each bipartition as this could result in an
  // imbalanced k-way partition when performing recursive bipartitioning. We therefore adaptively
  // adjust the allowed imbalance for each bipartition individually based on the adaptive imbalance
  // definition described in our papers.
  double computeAdaptiveEpsilon(const HNWeightScalar current_hypergraph_weight,
                                const PartitionID current_k) const {
    // TODO: epsilon could be multi-dimensional ??
    if ( current_hypergraph_weight == 0 ) {
      // In recursive bipartitioning, it can happen that a block becomes too light that
      // all nodes of the block fit into one block in a subsequent bipartitioning step.
      // This will create an empty block, which we fix later in a rebalancing step.
      return 0.0;
    } else {
      double base = ceil(static_cast<double>(original_hypergraph_weight) / original_k)
        / ceil(static_cast<double>(current_hypergraph_weight) / current_k)
        * (1.0 + original_epsilon);
      double adaptive_epsilon = std::min(0.99, std::max(std::pow(base, 1.0 /
        ceil(log2(static_cast<double>(current_k)))) - 1.0, 0.0));
      return adaptive_epsilon;
    }
  }

  const HNWeightScalar original_hypergraph_weight;
  const PartitionID original_k;
  const double original_epsilon;
};

MT_KAHYPAR_ATTRIBUTE_NO_INLINE
inline void setupBipartitioningContext(Context& b_context,
                                       const Context& context,
                                       const OriginalHypergraphInfo& info,
                                       const PartitionID start_k,
                                       const PartitionID end_k,
                                       const HNWeightConstRef total_weight,
                                       const bool is_graph) {
  ASSERT(end_k - start_k >= 2);
  const Dimension dimension = context.dimension();

  b_context.partition.k = 2;
  b_context.partition.objective = Objective::cut;
  b_context.partition.gain_policy = is_graph ? GainPolicy::cut_for_graphs : GainPolicy::cut;
  b_context.partition.enable_logging = false;
  b_context.initial_partitioning.mode = Mode::direct;
  b_context.type = ContextType::initial_partitioning;

  // Setup Part Weights
  const PartitionID k = end_k - start_k;
  const PartitionID k0 = k / 2 + (k % 2 != 0 ? 1 : 0);
  const PartitionID k1 = k / 2;
  ASSERT(k0 + k1 == k);
  if ( context.partition.use_individual_part_weights ) {
    AllocatedHNWeight max_part_weights_sum(dimension, 0);
    std::vector<double> weight_fraction(dimension, 0);

    for (const auto& part_weight: context.partition.max_part_weights) {
      max_part_weights_sum += part_weight;
    }
    for (size_t d = 0; d < dimension; ++d) {
      weight_fraction[d] = total_weight.at(d) / static_cast<double>(max_part_weights_sum.at(d));
      // assertion doesn't hold if previous bipartitions are imbalanced...
      // ASSERT(weight_fraction[d] <= 1.0);
    }

    AllocatedHNWeight perfect_weight_p0(dimension, 0);
    for ( PartitionID i = 0; i < k0; ++i ) {
      for (size_t d = 0; d < dimension; ++d) {
        perfect_weight_p0.at(d) += ceil(weight_fraction[d] * context.partition.max_part_weights[i].at(d));
      }
    }
    AllocatedHNWeight perfect_weight_p1(dimension, 0);
    for ( PartitionID i = k0; i < k; ++i ) {
      for (size_t d = 0; d < dimension; ++d) {
        perfect_weight_p1.at(d) += ceil(weight_fraction[d] * context.partition.max_part_weights[i].at(d));
      }
    }
    // In the case of individual part weights, the usual adaptive epsilon formula is not applicable because it
    // assumes equal part weights. However, by observing that ceil(current_weight / current_k) is the current
    // perfect part weight and (1 + epsilon)ceil(original_weight / original_k) is the maximum part weight,
    // we can derive an equivalent formula using the sum of the perfect part weights and the sum of the
    // maximum part weights.
    // Note that the sum of the perfect part weights might be unequal to the hypergraph weight due to rounding.
    // Thus, we need to use the former instead of using the hypergraph weight directly, as otherwise it could
    // happen that (1 + epsilon)perfect_part_weight > max_part_weight because of rounding issues.
    std::vector<double> epsilon_per_dimension(dimension, 0);
    double summed_epsilon = 0;
    for (size_t d = 0; d < dimension; ++d) {
      const double base = max_part_weights_sum.at(d) / static_cast<double>((perfect_weight_p0 + perfect_weight_p1).at(d));
      epsilon_per_dimension[d] = total_weight.at(d) == 0 ? 0 : std::clamp(
                                    std::pow(base, 1.0 / ceil(log2(static_cast<double>(k)))) - 1.0, 0.0, 0.99);
      summed_epsilon += epsilon_per_dimension[d];
    }

    b_context.partition.epsilon = summed_epsilon / dimension;  // TODO
    b_context.partition.perfect_balance_part_weights.replaceWith(2, dimension, 0, false);
    b_context.partition.max_part_weights.replaceWith(2, dimension, 0, false);
    b_context.partition.perfect_balance_part_weights[0] = perfect_weight_p0;
    b_context.partition.perfect_balance_part_weights[1] = perfect_weight_p1;
    for (size_t d = 0; d < dimension; ++d) {
      b_context.partition.max_part_weights[0].set(d,
              round((1 + epsilon_per_dimension[d]) * perfect_weight_p0.at(d)));
      b_context.partition.max_part_weights[1].set(d,
              round((1 + epsilon_per_dimension[d]) * perfect_weight_p1.at(d)));
    }
  } else {
    b_context.partition.epsilon = info.computeAdaptiveEpsilon(weight::sum(total_weight), k);

    b_context.partition.perfect_balance_part_weights.replaceWith(2, dimension, 0, false);
    b_context.partition.max_part_weights.replaceWith(2, dimension, 0, false);
    for (size_t d = 0; d < dimension; ++d) {
      b_context.partition.perfect_balance_part_weights[0].set(d,
              std::ceil(k0 / static_cast<double>(k) * static_cast<double>(total_weight.at(d))));
      b_context.partition.perfect_balance_part_weights[1].set(d,
              std::ceil(k1 / static_cast<double>(k) * static_cast<double>(total_weight.at(d))));
      b_context.partition.max_part_weights[0].set(d,
              (1 + b_context.partition.epsilon) * b_context.partition.perfect_balance_part_weights[0].at(d));
      b_context.partition.max_part_weights[1].set(d,
              (1 + b_context.partition.epsilon) * b_context.partition.perfect_balance_part_weights[1].at(d));
    }
  }
  b_context.setupContractionLimit(total_weight);
  b_context.setupThreadsPerFlowSearch();
}


}  // namespace bipartitioning
}  // namespace mt_kahypar

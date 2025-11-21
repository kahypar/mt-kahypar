/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/weight/hypernode_weight_common.h"

namespace mt_kahypar {
namespace impl {

inline void getExtremalDimensions(HNWeightConstRef weight, const vec<double>& weight_normalizer, uint8_t* out_ptr, bool maximimize) {
  float max_weight = maximimize ? 0 : 1;
  for (Dimension d = 0; d < weight.dimension(); ++d) {
    float normalized_weight = weight_normalizer[d] * weight.at(d);
    max_weight = maximimize ? std::max(max_weight, normalized_weight) : std::min(max_weight, normalized_weight);
  }
  for (Dimension d = 0; d < weight.dimension(); ++d) {
    float normalized_weight = weight_normalizer[d] * weight.at(d);
    if (maximimize ? (normalized_weight >= max_weight) : (normalized_weight <= max_weight)) {
    out_ptr[d] = static_cast<uint8_t>(true);
    }
  }
}

inline bool hasMatchingDimension(HNWeightConstRef weight, const uint8_t* lhs, const uint8_t* rhs) {
  for (Dimension d = 0; d < weight.dimension(); ++d) {
    if (lhs[d] && rhs[d]) {
    return true;
    }
  }
  return false;
}

inline float weightOfMatchingDimension(HNWeightConstRef weight, const uint8_t* valid_dims, const vec<double>& weight_normalizer) {
  double sum = 0;
  for (Dimension d = 0; d < weight.dimension(); ++d) {
    if (valid_dims[d] && sum == 0) {
    sum += weight_normalizer[d] * static_cast<double>(weight.at(d));
    }
  }
  return sum;
}

template<typename HNWeightExpression>
float dotProduct(HNWeightConstRef weight, const HNWeightExpression& expr, const vec<double>& weight_normalizer) {
  double sum = 0;
  for (Dimension d = 0; d < weight.dimension(); ++d) {
    sum += weight_normalizer[d] * weight_normalizer[d]
            * static_cast<double>(weight.at(d))* static_cast<double>(expr.at(d));
  }
  return sum;
}

template<typename HNWeightExpression>
float normalizedSum(const HNWeightExpression& weight, const vec<double>& weight_normalizer) {
  double sum = 0;
  for (Dimension d = 0; d < weight.dimension(); ++d) {
    sum += weight_normalizer[d] * static_cast<double>(weight.at(d));
  }
  return sum;
}

inline float internalImbalanceByDim(HNWeightConstRef weight, const uint8_t* valid_dims, const vec<double>& weight_normalizer) {
  float matching_weight = impl::weightOfMatchingDimension(weight, valid_dims, weight_normalizer);
  float weight_sum = impl::normalizedSum(weight, weight_normalizer);
  return (0.01 * weight_sum + matching_weight) / (1.01 * weight_sum - matching_weight);
}

template<typename HNWeightExpression>
inline float internalImbalance(const HNWeightExpression& weight, const vec<double>& weight_normalizer) {
  float max_weight = 0;
  for (Dimension d = 0; d < weight.dimension(); ++d) {
    max_weight = std::max(max_weight, static_cast<float>(weight_normalizer[d] * weight.at(d)));
  }

  float weight_sum = impl::normalizedSum(weight, weight_normalizer);
  return (0.01 * weight_sum + max_weight) / (1.01 * weight_sum - max_weight);
}

inline double imbalance(const HypernodeWeightArray& part_weights, const Context& context) {
  double max_balance = 0;
  for (PartitionID i = 0; i < context.partition.k; ++i) {
    for (Dimension d = 0; d < part_weights.dimension(); ++d) {
      const double curr_balance =
              (part_weights[i].at(d) /
                  static_cast<double>(context.partition.perfect_balance_part_weights[i].at(d)));
      max_balance = std::max(max_balance, curr_balance);
    }
  }

  return max_balance - 1.0;
}

inline double imbalanceSum(const HypernodeWeightArray& part_weights, const Context& context, const vec<double>& weight_normalizer) {
double sum = 0;
  for (PartitionID i = 0; i < context.partition.k; ++i) {
    for (Dimension d = 0; d < part_weights.dimension(); ++d) {
    HNWeightScalar diff = part_weights[i].at(d) - context.partition.max_part_weights[i].at(d);
      if (diff > 0) {
        sum += weight_normalizer[d] * static_cast<double>(diff);
      }
    }
  }
  return sum;
}

inline float transformGainFromProgress(Gain gain_, float progress, bool has_negative_progress, double negative_progress_penalty) {
  // here: positive gain means improvement
  float gain = gain_;
  if (has_negative_progress) {
    progress /= negative_progress_penalty;
  }
  if (gain > 0) {
    gain *= progress;
  } else if (gain < 0) {
    gain /= progress;
  }
  return gain;
}

inline float transformGain(Gain gain_, HNWeightConstRef wu, HNWeightAtomicCRef from_weight, HNWeightConstRef max_part_weight) {
  // here: positive gain means improvement
  if (wu.dimension() == 1) {
    float gain = gain_;
    if (gain > 0) {
    gain *= wu.at(0);
    } else if (gain < 0) {
    gain /= wu.at(0);
    }
    return gain;
  } else {
    float relevant_weight_fraction = 0;
    for (Dimension d = 0; d < wu.dimension(); ++d) {
      if (from_weight.at(d) > max_part_weight.at(d)) {
        relevant_weight_fraction += wu.at(d) / static_cast<float>(max_part_weight.at(d));
      }
    }
    return transformGainFromProgress(gain_, relevant_weight_fraction, false, 1.0);
  }
}

inline std::pair<float, bool> computeBalanceProgress(HNWeightConstRef wu, HNWeightConstRef from_weight, HNWeightConstRef max_part_weight_from,
                                                     HNWeightAtomicCRef to_weight, HNWeightConstRef max_part_weight_to, const vec<double>& weight_normalizer) {
  ASSERT(wu.dimension() > 1);
  float relative_progress = 0;
  const auto new_to_weight = wu + to_weight;
  bool has_negative_progress = false;
  for (Dimension d = 0; d < wu.dimension(); ++d) {
    const HNWeightScalar reduced_weight = std::min(from_weight.at(d) - max_part_weight_from.at(d), wu.at(d));
    if (reduced_weight > 0) {
      relative_progress += weight_normalizer[d] * reduced_weight;
    }
    const HNWeightScalar overweight = std::min(new_to_weight.at(d) - max_part_weight_to.at(d), wu.at(d));
    if (overweight > 0) {
      relative_progress -= weight_normalizer[d] * overweight;
      has_negative_progress = true;
    }
  }
  return {relative_progress, has_negative_progress};
}

inline vec<vec<double>> computeBlockWeightNormalizers(const Context& context) {
  vec<vec<double>> block_weight_normalizers(context.partition.k, vec<double>{});
  for (PartitionID block = 0; block < context.partition.k; ++block) {
    vec<double> normalizer(context.dimension(), 0);
    for (Dimension d = 0; d < context.dimension(); ++d) {
      normalizer[d] = 1 / static_cast<double>(context.partition.max_part_weights[block].at(d));
    }
    block_weight_normalizers[block] = std::move(normalizer);
  }
  return block_weight_normalizers;
}

} // namespace impl
}  // namespace mt_kahypar

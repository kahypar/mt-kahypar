/*******************************************************************************
 * MIT License
 *
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2017 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include "kahypar-resources/meta/policy_registry.h"
#include "kahypar-resources/meta/typelist.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/macros.h"

namespace mt_kahypar {

class NoWeightPenalty final : public kahypar::meta::PolicyBase {
 public:
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HypernodeWeight penalty(const HypernodeWeight, const HypernodeWeight, const Context) {
    return 1;
  }
};

class MultiplicativePenalty final : public kahypar::meta::PolicyBase {
 public:
   MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HypernodeWeight penalty(const HypernodeWeight weight_u, const HypernodeWeight weight_v, const Context) {
    return weight_u * weight_v;
  }
};

class AdditivePenalty final : public kahypar::meta::PolicyBase {
 public:
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HypernodeWeight penalty(const HypernodeWeight weight_u, const HypernodeWeight weight_v,const Context) {
    return weight_u + weight_v;
  }
};

class ScalingAdditivePenalty final : public kahypar::meta::PolicyBase {
 public:
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HypernodeWeight penalty(const HypernodeWeight weight_u, const HypernodeWeight weight_v, const Context &context) {
    HypernodeWeight threshold = context.coarsening.max_allowed_node_weight *
                                context.coarsening.penalty_max_weight_factor;
    return std::pow(std::max(1, static_cast<int>(weight_u + weight_v - threshold)), context.coarsening.gamma);
  }
};

using HeavyNodePenaltyPolicies = kahypar::meta::Typelist<MultiplicativePenalty,
                                                         NoWeightPenalty,
                                                         AdditivePenalty,
                                                         ScalingAdditivePenalty>;
}  // namespace mt_kahypar

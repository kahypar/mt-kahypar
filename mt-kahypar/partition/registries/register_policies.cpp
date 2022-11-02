/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "kahypar/meta/policy_registry.h"
#include "kahypar/meta/registrar.h"

#include "mt-kahypar/partition/coarsening/policies/rating_acceptance_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_heavy_node_penalty_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_score_policy.h"
#include "mt-kahypar/partition/context_enum_classes.h"

#define REGISTER_POLICY(policy, id, policy_class)                                                    \
  static kahypar::meta::Registrar<kahypar::meta::PolicyRegistry<policy> > register_ ## policy_class( \
    id, new policy_class())

namespace mt_kahypar {
// //////////////////////////////////////////////////////////////////////////////
//                       Coarsening / Rating Policies
// //////////////////////////////////////////////////////////////////////////////
REGISTER_POLICY(RatingFunction, RatingFunction::heavy_edge,
                HeavyEdgeScore);
REGISTER_POLICY(RatingFunction, RatingFunction::sameness,
                SamenessScore);

REGISTER_POLICY(HeavyNodePenaltyPolicy, HeavyNodePenaltyPolicy::no_penalty,
                NoWeightPenalty);
REGISTER_POLICY(HeavyNodePenaltyPolicy, HeavyNodePenaltyPolicy::multiplicative_penalty,
                MultiplicativePenalty);
REGISTER_POLICY(HeavyNodePenaltyPolicy, HeavyNodePenaltyPolicy::additive,
                AdditivePenalty);

REGISTER_POLICY(AcceptancePolicy, AcceptancePolicy::best,
                BestRatingWithTieBreaking);
REGISTER_POLICY(AcceptancePolicy, AcceptancePolicy::best_prefer_unmatched,
                BestRatingPreferringUnmatched);
}  // namespace mt_kahypar

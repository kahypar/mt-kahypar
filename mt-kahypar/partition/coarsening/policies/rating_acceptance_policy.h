/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2017 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "kahypar/meta/policy_registry.h"
#include "kahypar/meta/typelist.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/coarsening/policies/rating_tie_breaking_policy.h"

namespace mt_kahypar {

class BestRatingWithoutTieBreaking final : public kahypar::meta::PolicyBase {
 public:
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static bool acceptRating(const RatingType tmp,
                                                              const RatingType max_rating,
                                                              const HypernodeID u,
                                                              const HypernodeID v,
                                                              const int,
                                                              const kahypar::ds::FastResetFlagArray<> &) {
    return max_rating < tmp || ( max_rating == tmp && u < v );
  }
};

class BestRatingWithTieBreaking final : public kahypar::meta::PolicyBase {
 public:
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static bool acceptRating(const RatingType tmp,
                                                              const RatingType max_rating,
                                                              const HypernodeID,
                                                              const HypernodeID,
                                                              const int cpu_id,
                                                              const kahypar::ds::FastResetFlagArray<> &) {
    return max_rating < tmp || (max_rating == tmp && RandomRatingWins::acceptEqual(cpu_id));
  }
};

class BestRatingPreferringUnmatched final : public kahypar::meta::PolicyBase {
 public:
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static bool acceptRating(const RatingType tmp,
                                                              const RatingType max_rating,
                                                              const HypernodeID old_target,
                                                              const HypernodeID new_target,
                                                              const int cpu_id,
                                                              const kahypar::ds::FastResetFlagArray<>& already_matched) {
    return max_rating < tmp ||
           ((max_rating == tmp) &&
            ((already_matched[old_target] && !already_matched[new_target]) ||
             (already_matched[old_target] && already_matched[new_target] &&
              RandomRatingWins::acceptEqual(cpu_id)) ||
             (!already_matched[old_target] && !already_matched[new_target] &&
              RandomRatingWins::acceptEqual(cpu_id))));
  }
};

using AcceptancePolicies = kahypar::meta::Typelist<BestRatingWithTieBreaking, BestRatingPreferringUnmatched>;
}  // namespace mt_kahypar

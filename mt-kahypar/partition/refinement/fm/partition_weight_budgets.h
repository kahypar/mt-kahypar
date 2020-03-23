/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/definitions.h"


namespace mt_kahypar {
namespace refinement {
class PartitionWeightBudgets {
public:

  // Initially distribute budget evenly. It is allowed to steal everything from another processor.
  // To ensure we use all of the available budget, even if we use less than max_num_threads,
  // we steal from the highest search IDs first.

  PartitionWeightBudgets(size_t k, size_t max_num_threads) : search_local_budgets(max_num_threads, vec<CAtomic<HypernodeWeight>>(k, CAtomic<HypernodeWeight>(0))) { }


  void initialize(vec<HypernodeWeight>& global_budgets) {
    for (PartitionID p = 0; p < global_budgets.size(); ++p) {
      HypernodeWeight local_budget = static_cast<HypernodeWeight>(global_budgets[p] / maxNumThreads());
      size_t num_threads_with_one_additional = global_budgets[p] % maxNumThreads();

      for (size_t thread_id = 0; thread_id < num_threads_with_one_additional; ++thread_id) {
        search_local_budgets[thread_id][p].store(local_budget + 1);
      }
      for (size_t thread_id = num_threads_with_one_additional; thread_id < maxNumThreads(); ++thread_id) {
        search_local_budgets[thread_id][p].store(local_budget);
      }
    }
  }

  size_t maxNumThreads() const {
    return search_local_budgets.size();
  }

  vec<vec<CAtomic<HypernodeWeight>>> search_local_budgets;

  vec<CAtomic<HypernodeWeight>>& localBudgets(uint32_t my_search) {
    return search_local_budgets[my_search];
  }

  // don't try stealing every time a move is infeasible
  HypernodeWeight steal(uint32_t my_search, PartitionID to, HypernodeWeight least_desired_amount) {
    for (uint32_t i = static_cast<SearchID>(search_local_budgets.size()); i > 0; --i) {
      if (i == my_search) continue;
    }
    return 0;
  }

};

}
}
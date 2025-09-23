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

#pragma once

#include <cstdint>

#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/coarsening/coarsening_commons.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {

template<typename T>
struct Statistic {
  float avg = 0.0;
  float sd = 0.0;
  float skew = 0.0;
  float entropy = 0.0;
  T min = 0;
  T q1 = 0;
  T med = 0;
  T q3 = 0;
  T max = 0;
};

struct GlobalFeatures {
  uint32_t n = 0;
  uint32_t m = 0;
  double irregularity = 0.0;
  uint32_t exp_median_degree = 0;
  Statistic<uint32_t> degree_stats;
  uint32_t n_communities_0 = 0;
  uint32_t n_communities_1 = 0;
  uint32_t n_communities_2 = 0;
  double modularity_0 = 0.0;
  double modularity_1 = 0.0;
  double modularity_2 = 0.0;
};

struct N1Features {
  uint32_t degree = 0;
  float inverse_degree = 0;
  float degree_quantile = 0;
  Statistic<uint32_t> degree_stats;
  uint32_t to_n1_edges = 0;
  uint32_t to_n2_edges = 0;
  uint32_t d1_nodes = 0;
  float modularity = 0;
  float max_modularity = 0;
  uint32_t max_modularity_size = 0;
  float min_contracted_degree = 0;
  uint32_t min_contracted_degree_size = 0;
  float clustering_coefficient = 0;
  float chi_squared_degree_deviation = 0;
};

struct EdgeFeatures {
  float strawman_similarity = 0;
  bool comm_0_equal = false;
  bool comm_1_equal = false;
  bool comm_2_equal = false;
};

std::tuple<GlobalFeatures, ds::Array<N1Features>, bool> computeFeatures(const ds::StaticGraph& graph, const Context& context);

inline EdgeFeatures computeEdgeFeatures(const ds::StaticGraph& /*graph*/, const Context& context, HypernodeID u, const N1Features& u_f, HypernodeID v, const N1Features& v_f, bool skip_comm_1) {
  EdgeFeatures result;
  result.strawman_similarity = (u_f.inverse_degree == 1.0 || v_f.inverse_degree == 1.0) ? 1.0 : u_f.inverse_degree * v_f.inverse_degree;

  auto comm_ptr = context.preprocessing.community_stack;
  if (comm_ptr != nullptr) {
    const auto& community_stack = *comm_ptr;

    auto equal_communities = [&](size_t i) {
      ASSERT(community_stack.size() - i - 1 < community_stack.size());
      const auto& [clustering, _, _m] = community_stack[community_stack.size() - i - 1];
      return clustering[u] == clustering[v];
    };

    result.comm_0_equal = equal_communities(0);
    if (skip_comm_1) {
      result.comm_1_equal = equal_communities(2);
      result.comm_2_equal = equal_communities(3);
    } else {
      result.comm_1_equal = equal_communities(1);
      result.comm_2_equal = equal_communities(2);
    }
  }
  return result;
}

}  // namespace mt_kahypar

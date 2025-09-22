/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Nikolai Maas <nikolai.maas@kit.edu>
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
  double avg = 0.0;
  double sd = 0.0;
  double skew = 0.0;
  double entropy = 0.0;
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
  double degree_quantile = 0;
  Statistic<uint32_t> degree_stats;
  uint32_t to_n1_edges = 0;
  uint32_t to_n2_edges = 0;
  uint32_t d1_nodes = 0;
  double modularity = 0;
  double max_modularity = 0;
  uint32_t max_modularity_size = 0;
  double min_contracted_degree = 0;
  uint32_t min_contracted_degree_size = 0;
  double clustering_coefficient = 0;
  double chi_squared_degree_deviation = 0;
};

std::pair<GlobalFeatures, ds::Array<N1Features>> computeFeatures(const ds::StaticGraph& graph, const Context& context);

}  // namespace mt_kahypar

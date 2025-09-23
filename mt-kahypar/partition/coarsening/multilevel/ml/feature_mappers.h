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

#include "mt-kahypar/partition/coarsening/multilevel/ml/compute_features.h"

namespace mt_kahypar {
namespace features {

// macros for simplifying the definitions
#define  JOIN(a, b)  a ## b

#define GLOBAL_MAPPER(FEATURE_NAME) \
  struct FEATURE_NAME { \
    static float map(const GlobalFeatures& global, const N1Features&, const N1Features&, const EdgeFeatures&) { \
      return global.FEATURE_NAME; \
    } \
  }

#define DEGREE_MAPPER(FEATURE_NAME) \
  struct JOIN(FEATURE_NAME, _degree) { \
    static float map(const GlobalFeatures& global, const N1Features&, const N1Features&, const EdgeFeatures&) { \
      return global.degree_stats.FEATURE_NAME; \
    } \
    \
    static float map_n1(const N1Features& features) { \
      return features.degree_stats.FEATURE_NAME; \
    } \
  }

#define EDGE_MAPPER(FEATURE_NAME) \
  struct FEATURE_NAME { \
    static float map(const GlobalFeatures&, const N1Features&, const N1Features&, const EdgeFeatures& edge) { \
      return edge.FEATURE_NAME; \
    } \
  }

#define N1_MAPPER(FEATURE_NAME) \
  struct FEATURE_NAME { \
    static float map_n1(const N1Features& features) { \
      return features.FEATURE_NAME; \
    } \
  }

#define N1_DEGREE_MAPPER(FEATURE_NAME) \
  struct FEATURE_NAME { \
    static float map_n1(const N1Features& features) { \
      return features.FEATURE_NAME; \
    } \
  }


// define the mappings for all features
template<typename N1_MAPPER>
struct N0 {
  static float map(const GlobalFeatures&, const N1Features& features, const N1Features&, const EdgeFeatures&) {
    return N1_MAPPER::map_n1(features);
  }
};

template<typename N1_MAPPER>
struct N1 {
  static float map(const GlobalFeatures&, const N1Features&, const N1Features& features, const EdgeFeatures&) {
    return N1_MAPPER::map_n1(features);
  }
};

GLOBAL_MAPPER(n);
GLOBAL_MAPPER(m);
GLOBAL_MAPPER(irregularity);
GLOBAL_MAPPER(exp_median_degree);
GLOBAL_MAPPER(n_communities_0);
GLOBAL_MAPPER(n_communities_1);
GLOBAL_MAPPER(n_communities_2);
GLOBAL_MAPPER(modularity_0);
GLOBAL_MAPPER(modularity_1);
GLOBAL_MAPPER(modularity_2);

DEGREE_MAPPER(avg);
DEGREE_MAPPER(sd);
DEGREE_MAPPER(skew);
DEGREE_MAPPER(entropy);
DEGREE_MAPPER(min);
DEGREE_MAPPER(q1);
DEGREE_MAPPER(med);
DEGREE_MAPPER(q3);
DEGREE_MAPPER(max);

EDGE_MAPPER(strawman_similarity);
EDGE_MAPPER(comm_0_equal);
EDGE_MAPPER(comm_1_equal);
EDGE_MAPPER(comm_2_equal);

N1_MAPPER(degree);
N1_MAPPER(degree_quantile);
N1_MAPPER(to_n1_edges);
N1_MAPPER(to_n2_edges);
N1_MAPPER(d1_nodes);
N1_MAPPER(modularity);
N1_MAPPER(max_modularity);
N1_MAPPER(max_modularity_size);
N1_MAPPER(min_contracted_degree);
N1_MAPPER(min_contracted_degree_size);
N1_MAPPER(clustering_coefficient);
N1_MAPPER(chi_squared_degree_deviation);

}  // namespace features
}  // namespace mt_kahypar

/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <boost/program_options.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <cstdint>
#include <iomanip>
#include <type_traits>
#include <limits>
#include <chrono>

#include "tbb/parallel_sort.h"
#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_reduce.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/static_graph.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/preprocessing/community_detection/parallel_louvain.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/delete.h"
#include "mt-kahypar/utils/hypergraph_statistics.h"

#include "kahypar-resources/utils/math.h"

#include "neighborhood_computation.h"

using namespace mt_kahypar;
namespace po = boost::program_options;

using StaticGraph = ds::StaticGraph;
using LouvainGraph = ds::Graph<ds::StaticGraph>;

enum class FeatureType {
  floatingpoint,
  integer,
  boolean,
};

struct Feature {
  double value;
  FeatureType type;

  Feature(double value, FeatureType type): value(value), type(type) {
    ALWAYS_ASSERT(type == FeatureType::floatingpoint);
  }
  Feature(uint64_t value, FeatureType type): value(static_cast<double>(value)), type(type) {
    ALWAYS_ASSERT(type != FeatureType::floatingpoint);
    if (type == FeatureType::boolean) {
      ALWAYS_ASSERT(value == 0 || value == 1);
    }
  }
  Feature(uint32_t value, FeatureType type): Feature(static_cast<uint64_t>(value), type) { }
  Feature(uint8_t value, FeatureType type): Feature(static_cast<uint64_t>(value), type) { }
};


void printHeader(std::ostream& os, const std::vector<std::string>& header) {
  for (size_t i = 0; i < header.size(); ++i) {
    if (i > 0) {
      os << ",";
    }
    os << header[i];
  }
  os << std::endl;
}

void printFeatures(std::ostream& os, const std::vector<Feature>& features) {
  for (size_t i = 0; i < features.size(); ++i) {
    if (i > 0) {
      os << ",";
    }
    if (features[i].type != FeatureType::floatingpoint) {
      os << static_cast<uint64_t>(features[i].value);
    } else {
      os << std::setprecision(4) << features[i].value;
    }
  }
  os << std::endl;
}


// ################ Feature Definitions ################

struct Invalid {};

template<typename T = Invalid>
struct Statistic {
  static constexpr uint64_t num_entries = 9;

  double avg = 0.0;
  double sd = 0.0;
  double skew = 0.0;
  double entropy = 0.0;
  T min = 0;
  T q1 = 0;
  T med = 0;
  T q3 = 0;
  T max = 0;
  // TODO: chi^2 ?

  std::vector<Feature> featureList() const {
    FeatureType type = std::is_floating_point_v<T> ? FeatureType::floatingpoint : FeatureType::integer;
    std::vector<Feature> result {
      {avg, FeatureType::floatingpoint},
      {sd, FeatureType::floatingpoint},
      {skew, FeatureType::floatingpoint},
      {entropy, FeatureType::floatingpoint},
      {min, type},
      {q1, type},
      {med, type},
      {q3, type},
      {max, type},
    };
    ALWAYS_ASSERT(result.size() == num_entries, "header info was not properly updated");
    return result;
  }

  static std::vector<std::string> header(const char* suffix) {
    std::vector<std::string> result {
      std::string("avg_") + suffix,
      std::string("sd_") + suffix,
      std::string("skew_") + suffix,
      std::string("entropy_") + suffix,
      std::string("min_") + suffix,
      std::string("q1_") + suffix,
      std::string("med_") + suffix,
      std::string("q3_") + suffix,
      std::string("max_") + suffix,
    };
    ALWAYS_ASSERT(result.size() == num_entries, "header info was not properly updated");
    return result;
  }
};

struct GlobalFeatures {
  static constexpr uint64_t num_entries = 2 * Statistic<>::num_entries + 10;

  uint64_t n = 0;
  uint64_t m = 0;
  double irregularity = 0.0;
  uint64_t exp_median_degree = 0;
  Statistic<uint64_t> degree_stats;  // over nodes
  Statistic<double> locality_stats;  // over edges
  uint64_t n_communities_0 = 0;
  uint64_t n_communities_1 = 0;
  uint64_t n_communities_2 = 0;
  double modularity_0 = 0.0;
  double modularity_1 = 0.0;
  double modularity_2 = 0.0;

  std::vector<Feature> featureList() const {
    std::vector<Feature> result_1 {
      {n, FeatureType::integer},
      {m, FeatureType::integer},
      {irregularity, FeatureType::floatingpoint},
      {exp_median_degree, FeatureType::integer},
    };
    std::vector<Feature> degree_features = degree_stats.featureList();
    std::vector<Feature> locality_features = locality_stats.featureList();
    result_1.insert(result_1.end(), degree_features.begin(), degree_features.end());
    result_1.insert(result_1.end(), locality_features.begin(), locality_features.end());
    std::vector<Feature> result_2 {
      {n_communities_0, FeatureType::integer},
      {n_communities_1, FeatureType::integer},
      {n_communities_2, FeatureType::integer},
      {modularity_0, FeatureType::floatingpoint},
      {modularity_1, FeatureType::floatingpoint},
      {modularity_2, FeatureType::floatingpoint},
    };
    result_1.insert(result_1.end(), result_2.begin(), result_2.end());
    ALWAYS_ASSERT(result_1.size() == num_entries, "header info was not properly updated");
    return result_1;
  }

  static std::vector<std::string> header() {
    std::vector<std::string> result_1 {"n", "m", "irregularity", "exp_median_degree"};
    std::vector<std::string> degree_header = Statistic<>::header("degree");
    std::vector<std::string> locality_header = Statistic<>::header("locality");
    result_1.insert(result_1.end(), degree_header.begin(), degree_header.end());
    result_1.insert(result_1.end(), locality_header.begin(), locality_header.end());
    std::vector<std::string> result_2 {
      "n_communities_0", "n_communities_1", "n_communities_2",
      "modularity_0", "modularity_1", "modularity_2",
    };
    result_1.insert(result_1.end(), result_2.begin(), result_2.end());
    ALWAYS_ASSERT(result_1.size() == num_entries, "header info was not properly updated");
    return result_1;
  }
};

struct N1Features {
  static constexpr uint64_t num_entries = 2 * Statistic<>::num_entries + 13;

  uint64_t degree = 0;
  double degree_quantile = 0;
  Statistic<uint64_t> degree_stats;  // over nodes
  Statistic<double> locality_stats;  // over nodes
  uint64_t to_n1_edges = 0;
  uint64_t to_n2_edges = 0;
  uint64_t d1_nodes = 0;
  double modularity = 0;
  double max_modularity = 0;
  uint64_t max_modularity_size = 0;
  double min_contracted_degree = 0;
  uint64_t min_contracted_degree_size = 0;
  double clustering_coefficient = 0;
  double chi_squared_degree_deviation = 0;
  uint64_t max_clique = 0;
  // TODO: community overlap?

  std::vector<Feature> featureList() const {
    std::vector<Feature> result_1 {
      {degree, FeatureType::integer},
      {degree_quantile, FeatureType::floatingpoint},
    };
    std::vector<Feature> degree_features = degree_stats.featureList();
    std::vector<Feature> locality_features = locality_stats.featureList();
    result_1.insert(result_1.end(), degree_features.begin(), degree_features.end());
    result_1.insert(result_1.end(), locality_features.begin(), locality_features.end());
    std::vector<Feature> result_2 {
      {to_n1_edges, FeatureType::integer},
      {to_n2_edges, FeatureType::integer},
      {d1_nodes, FeatureType::integer},
      {modularity, FeatureType::floatingpoint},
      {max_modularity, FeatureType::floatingpoint},
      {max_modularity_size, FeatureType::integer},
      {min_contracted_degree, FeatureType::floatingpoint},
      {min_contracted_degree_size, FeatureType::integer},
      {clustering_coefficient, FeatureType::floatingpoint},
      {chi_squared_degree_deviation, FeatureType::floatingpoint},
      {max_clique, FeatureType::integer},
    };
    result_1.insert(result_1.end(), result_2.begin(), result_2.end());
    ALWAYS_ASSERT(result_1.size() == num_entries, "header info was not properly updated");
    return result_1;
  }

  static std::vector<std::string> header() {
    std::vector<std::string> result_1 {"degree", "degree_quantile"};
    std::vector<std::string> degree_header = Statistic<>::header("n1_degree");
    std::vector<std::string> locality_header = Statistic<>::header("n1_locality");
    result_1.insert(result_1.end(), degree_header.begin(), degree_header.end());
    result_1.insert(result_1.end(), locality_header.begin(), locality_header.end());
    std::vector<std::string> result_2 {
      "n1_to_n1_edges", "n1_to_n2_edges", "n1_d1_nodes", "n1_modularity", "n1_max_modularity",
      "n1_max_modularity_size", "n1_min_contracted_degree", "n1_min_contracted_degree_size", "n1_clustering_coefficient",
      "chi_squared_degree_deviation", "n1_max_clique",
    };
    result_1.insert(result_1.end(), result_2.begin(), result_2.end());
    ALWAYS_ASSERT(result_1.size() == num_entries, "header info was not properly updated");
    return result_1;
  }

  static std::vector<std::string> header(const char* prefix) {
    std::string p(prefix);
    std::vector<std::string> result = header();
    for (std::string& val: result) {
      val = p + val;
    }
    return result;
  }
};

struct N2Features {
  static constexpr uint64_t num_entries = 2 * Statistic<>::num_entries + 6;

  uint64_t size = 0;
  Statistic<uint64_t> degree_stats;  // over nodes
  Statistic<double> locality_stats;  // over nodes
  uint64_t to_n1n2_edges = 0;
  uint64_t out_edges = 0;
  double modularity = 0;
  double max_modularity = 0;
  uint64_t max_modularity_size = 0;
  // TODO: community overlap?

  std::vector<Feature> featureList() const {
    std::vector<Feature> result_1 { {size, FeatureType::integer} };
    std::vector<Feature> degree_features = degree_stats.featureList();
    std::vector<Feature> locality_features = locality_stats.featureList();
    result_1.insert(result_1.end(), degree_features.begin(), degree_features.end());
    result_1.insert(result_1.end(), locality_features.begin(), locality_features.end());
    std::vector<Feature> result_2 {
      {to_n1n2_edges, FeatureType::integer},
      {out_edges, FeatureType::integer},
      {modularity, FeatureType::floatingpoint},
      {max_modularity, FeatureType::floatingpoint},
      {max_modularity_size, FeatureType::integer},
    };
    result_1.insert(result_1.end(), result_2.begin(), result_2.end());
    ALWAYS_ASSERT(result_1.size() == num_entries, "header info was not properly updated");
    return result_1;
  }

  static std::vector<std::string> header() {
    std::vector<std::string> result_1 {"n2_size"};
    std::vector<std::string> degree_header = Statistic<>::header("n2_degree");
    std::vector<std::string> locality_header = Statistic<>::header("n2_locality");
    result_1.insert(result_1.end(), degree_header.begin(), degree_header.end());
    result_1.insert(result_1.end(), locality_header.begin(), locality_header.end());
    std::vector<std::string> result_2 {"n2_to_n1n2_edges", "n2_out_edges", "n2_modularity", "n2_max_modularity", "n2_max_modularity_size"};
    result_1.insert(result_1.end(), result_2.begin(), result_2.end());
    ALWAYS_ASSERT(result_1.size() == num_entries, "header info was not properly updated");
    return result_1;
  }
};

struct EdgeFeatures {
  static constexpr uint64_t num_entries = 2 * N1Features::num_entries + 8;

  double locality = 0;
  double jaccard_index = 0;
  double cosine_similarity = 0;
  double dice_similarity = 0;
  double strawman_similarity = 0;
  N1Features intersect_features;
  N1Features union_features;
  // TODO: community overlap?
  uint8_t comm_0_equal = 0;
  uint8_t comm_1_equal = 0;
  uint8_t comm_2_equal = 0;

  std::vector<Feature> featureList() const {
    std::vector<Feature> result_1 {
      {locality, FeatureType::floatingpoint},
      {jaccard_index, FeatureType::floatingpoint},
      {cosine_similarity, FeatureType::floatingpoint},
      {dice_similarity, FeatureType::floatingpoint},
      {strawman_similarity, FeatureType::floatingpoint},
    };
    std::vector<Feature> intersect_list = intersect_features.featureList();
    std::vector<Feature> union_list = union_features.featureList();
    result_1.insert(result_1.end(), intersect_list.begin(), intersect_list.end());
    result_1.insert(result_1.end(), union_list.begin(), union_list.end());
    std::vector<Feature> result_2 {
      {comm_0_equal, FeatureType::boolean},
      {comm_1_equal, FeatureType::boolean},
      {comm_2_equal, FeatureType::boolean},
    };
    result_1.insert(result_1.end(), result_2.begin(), result_2.end());
    ALWAYS_ASSERT(result_1.size() == num_entries, "header info was not properly updated");
    return result_1;
  }

  static std::vector<std::string> header() {
    std::vector<std::string> result_1 {
      "locality", "jaccard_index", "cosine_similarity", "dice_similarity", "strawman_similarity",
    };
    std::vector<std::string> intersect_header = N1Features::header("intersect_");
    std::vector<std::string> union_header = N1Features::header("union_");
    result_1.insert(result_1.end(), intersect_header.begin(), intersect_header.end());
    result_1.insert(result_1.end(), union_header.begin(), union_header.end());
    std::vector<std::string> result_2 {"comm_0_equal", "comm_1_equal", "comm_2_equal"};
    result_1.insert(result_1.end(), result_2.begin(), result_2.end());
    ALWAYS_ASSERT(result_1.size() == num_entries, "header info was not properly updated");
    return result_1;
  }
};



// ################ Feature Computation ################
double scaledEntropyFromOccurenceCounts(const ds::DynamicSparseMap<int64_t, int64_t>& occurence, size_t total) {
  // collect and sort summands
  std::vector<double> summands;
  for (auto& element : occurence) {
      double p_x = (double)element.value / (double)total;
      double summand = p_x * log2(p_x);
      // double summand = (pair.second * log2(pair.second) - pair.second * log2(total)) / total;
      summands.push_back(summand);
  }
  std::sort(summands.begin(), summands.end(), [] (double a, double b) { return abs(a) < abs(b); });
  // calculate entropy
  double entropy = 0;
  for (double summand : summands) {
      entropy -= summand;
  }
  // scale by log of number of categories    
  return log2(summands.size()) == 0 ? 0 : (double)entropy / log2(summands.size());
}

double scaledEntropy(const std::vector<double>& distribution) {
  ds::DynamicSparseMap<int64_t, int64_t> occurence;
  for (double value : distribution) {
    // snap to 2 digits after decimal point
    int64_t snap = static_cast<int64_t>(std::round(100 * value));
    occurence[snap]++;
  }
  return scaledEntropyFromOccurenceCounts(occurence, distribution.size());
}

double scaledEntropy(const std::vector<uint64_t>& distribution) {
  ds::DynamicSparseMap<int64_t, int64_t> occurence;
  for (auto value : distribution) {
    occurence[value]++;
  }
  return scaledEntropyFromOccurenceCounts(occurence, distribution.size());
}



template <typename T>
Statistic<T> createStats(std::vector<T>& vec, bool parallel) {
  Statistic<T> stats;
  if (!vec.empty()) {
    double avg = 0;
    double stdev = 0;
    double skew = 0;
    if (parallel) {
      tbb::parallel_sort(vec.begin(), vec.end());
      avg = utils::parallel_avg(vec, vec.size());
      stdev = utils::parallel_stdev(vec, avg, vec.size());
      skew = utils::parallel_skew(vec, avg, stdev, vec.size());
    } else {
      std::sort(vec.begin(), vec.end());
      for (auto val: vec) {
        avg += val;
      }
      avg = avg / static_cast<double>(vec.size());
      for (auto val: vec) {
        stdev += (val - avg) * (val - avg);
        skew += (val - avg) * (val - avg) * (val - avg);
      }
      stdev = std::sqrt(stdev / static_cast<double>(vec.size()));
      skew = stdev == 0 ? 0.0 : skew / (static_cast<double>(vec.size()) * std::pow(stdev, 3));
    }
    const auto quartiles = kahypar::math::firstAndThirdQuartile(vec);
    stats.min = vec.front();
    stats.q1 = quartiles.first;
    stats.med = kahypar::math::median(vec);
    stats.q3 = quartiles.second;
    stats.max = vec.back();
    stats.avg = avg;
    stats.sd = stdev;
    stats.skew = skew;
    stats.entropy = scaledEntropy(vec);
  }
  return stats;
}


bool float_eq(double left, double right) {
  return left <= 1.001 * right && right <= 1.001 * left;
}


std::tuple<GlobalFeatures, std::vector<uint64_t>, bool> computeGlobalFeatures(const StaticGraph& graph,
    std::vector<std::pair<ds::Clustering, double>>& community_stack) {
  GlobalFeatures features;

  std::vector<uint64_t> hn_degrees;
  hn_degrees.resize(graph.initialNumNodes());
  graph.doParallelForAllNodes([&](const HypernodeID& node) {
    hn_degrees[node] = graph.nodeDegree(node);
    ALWAYS_ASSERT(graph.nodeWeight(node) == 1);
  });

  HypernodeID num_nodes = graph.initialNumNodes();
  HyperedgeID num_edges = StaticGraph::is_graph ? graph.initialNumEdges() / 2 : graph.initialNumEdges();
  Statistic degree_stats = createStats(hn_degrees, true);
  features.n = num_nodes;
  features.m = num_edges;
  features.degree_stats = degree_stats;
  features.irregularity = degree_stats.sd / degree_stats.avg;

  // compute exp_median_degree via suffix sum
  uint64_t pins = graph.initialNumPins();
  uint64_t count = 0;
  for (size_t i = num_nodes; i > 0; --i) {
    count += hn_degrees[i-1];
    if (count >= pins / 2) {
      features.exp_median_degree = hn_degrees[i-1];
      break;
    }
  }

  // modularity features
  ds::DynamicSparseMap<PartitionID, uint32_t> comm_set;
  auto modularity_features = [&](size_t i) {
    const auto& [clustering, modularity] = community_stack.at(community_stack.size() - i - 1);
    comm_set.clear();
    for (PartitionID c: clustering) {
      comm_set[c] = 0;
    }
    uint64_t n_comms = 0;
    for (auto _: comm_set) {
      n_comms++;
    }
    return std::make_pair(n_comms, modularity);
  };
  bool skip_comm_1 = false;
  std::tie(features.n_communities_0, features.modularity_0) = modularity_features(0);
  std::tie(features.n_communities_1, features.modularity_1) = modularity_features(1);
  std::tie(features.n_communities_2, features.modularity_2) = modularity_features(2);
  if (community_stack.size() > 3 && features.n_communities_1 < 2 * features.n_communities_0) {
    // small hack to get more meaningful features
    std::tie(features.n_communities_1, features.modularity_1) = modularity_features(2);
    std::tie(features.n_communities_2, features.modularity_2) = modularity_features(3);
    skip_comm_1 = true;
  }

  return {features, hn_degrees, skip_comm_1};
}

// modularity, max_modularity, n_removed
std::tuple<double, double, uint64_t> maxModularitySubset(const StaticGraph& graph, uint64_t internal_edges, uint64_t all_edges,
                                                         FastResetArray& membership, const std::vector<HypernodeID>& node_list) {
  if (node_list.empty()) {
    return {0, 0, 0};
  }

  double factor = 1.0 / static_cast<double>(graph.initialNumEdges());
  uint64_t curr_internal_edges = internal_edges;
  uint64_t curr_all_edges = all_edges;
  double best = static_cast<double>(curr_internal_edges) - factor * all_edges * all_edges;
  int64_t n_removed = 0;
  const double modularity = best / static_cast<double>(all_edges);
  bool changed = true;
  for (size_t round = 0; changed && round < 5; ++round) {
    changed = false;
    for (HypernodeID node: node_list) {
      uint64_t local_edges = 0;
      for (HyperedgeID edge: graph.incidentEdges(node)) {
        if (membership[graph.edgeTarget(edge)]) {
          local_edges++;
        }
      }
      uint64_t updated_internal;
      uint64_t updated_all;
      if (membership[node]) {
        updated_internal = curr_internal_edges - 2 * local_edges;
        updated_all = curr_all_edges - graph.nodeDegree(node);
      } else {
        updated_internal = curr_internal_edges + 2 * local_edges;
        updated_all = curr_all_edges + graph.nodeDegree(node);
      }
      double updated_best = static_cast<double>(updated_internal) - factor * updated_all * updated_all;
      if (updated_best > best) {
        best = updated_best;
        curr_internal_edges = updated_internal;
        curr_all_edges = updated_all;
        n_removed += membership[node] ? 1 : -1;
        membership.set(node, !membership[node]);
        changed = true;
      }
    }
  }
  ALWAYS_ASSERT(n_removed >= 0);
  return {modularity, best / static_cast<double>(all_edges), n_removed};
}

N1Features n1FeaturesFromNeighborhood(const StaticGraph& graph, const GlobalFeatures& global, const std::vector<uint64_t>& global_degrees, const NeighborhoodResult& data,
                                      CliqueComputation* c_comp, HyperedgeID duplicated_degree, FastResetArray& membership) {
  membership.reset();
  for (HypernodeID root: data.roots) {
    membership.set(root);
  }

  N1Features result;
  HypernodeID num_nodes = data.n1_list.size();
  result.degree = num_nodes;
  size_t lower = std::lower_bound(global_degrees.begin(), global_degrees.end(), result.degree) - global_degrees.begin();  // always < n
  size_t upper = std::upper_bound(global_degrees.begin(), global_degrees.end(), result.degree) - global_degrees.begin();  // always >= 1
  result.degree_quantile = (static_cast<double>(lower + upper)) / (2 * static_cast<double>(global_degrees.size()) - 1);

  // compute degree stats
  std::vector<uint64_t> degrees;
  degrees.reserve(num_nodes);
  for (HypernodeID node: data.n1_list) {
    degrees.push_back(graph.nodeDegree(node));
    membership.set(node);
  }
  result.degree_stats = createStats(degrees, degrees.size() >= 20000);
  double chi_sq = 0.0;
  for (uint64_t d: degrees) {
    chi_sq += (static_cast<double>(d) - global.degree_stats.avg) * (static_cast<double>(d) - global.degree_stats.avg) / global.degree_stats.avg;
  }
  result.chi_squared_degree_deviation = chi_sq;

  // TODO: edges between neighbors are ignored in the same way as in
  // the coarsening heuristic. Is this what we want?
  HyperedgeWeight out_edges = result.degree;
  HypernodeWeight node_weight = data.roots[0] == data.roots[1] ? 1 : 2;
  for (HyperedgeID d: degrees) {
    if ((static_cast<HyperedgeWeight>(d) - 2) * node_weight < out_edges) {
      out_edges += d - 2;
      node_weight++;
    }
  }
  result.min_contracted_degree = static_cast<double>(out_edges) / static_cast<double>(node_weight);
  result.min_contracted_degree_size = node_weight;

  // compute locality stats and related values
  std::vector<double> locality_values;
  locality_values.reserve(num_nodes);
  for (HypernodeID node: data.n1_list) {
    uint64_t local_n1_edges = 0;
    for (HyperedgeID edge: graph.incidentEdges(node)) {
      HypernodeID neighbor = graph.edgeTarget(edge);
      if (data.isInN1Exactly(neighbor)) {
        result.to_n1_edges++;
        local_n1_edges++;
      } else if (!data.isRoot(neighbor)) {
        result.to_n2_edges++;
      }
    }
    HypernodeID node_degree = graph.nodeDegree(node);
    uint64_t divisor = std::min(num_nodes, node_degree) - 1;
    if (divisor != 0) {
      locality_values.push_back(static_cast<double>(local_n1_edges) / static_cast<double>(divisor));
    } else if (node_degree == 1) {
      result.d1_nodes++;
    }
  }
  result.to_n1_edges /= 2;  // doubly counted
  result.clustering_coefficient = (result.degree <= 1) ? 0 : static_cast<double>(2 * result.to_n1_edges) / static_cast<double>(result.degree * (result.degree -1));
  result.locality_stats = createStats(locality_values, locality_values.size() >= 20000);

  if (c_comp != nullptr) {
    result.max_clique = c_comp->computeMaxCliqueSize(graph, data.n1_list);
  }

  // modularity computations
  uint64_t internal_edges = 2 * (result.degree + duplicated_degree + result.to_n1_edges);
  uint64_t all_edges = internal_edges + result.to_n2_edges;
  auto [modularity, max_modularity, n_removed] = maxModularitySubset(graph, internal_edges, all_edges, membership, data.n1_list);
  result.modularity = modularity;
  result.max_modularity = max_modularity;
  result.max_modularity_size = result.degree - n_removed;

  // TODO: community overlap?
  return result;
}

N2Features n2FeaturesFromNeighborhood(const StaticGraph& graph, const NeighborhoodResult& data, const N1Features& n1_data,
                                     HyperedgeID duplicated_degree, FastResetArray& membership) {
  ALWAYS_ASSERT(data.includes_two_hop);
  membership.reset();
  for (HypernodeID root: data.roots) {
    membership.set(root);
  }
  for (HypernodeID node: data.n1_list) {
    membership.set(node);
  }

  N2Features result;
  HypernodeID num_nodes = data.n2_list.size();
  result.size = num_nodes;

  // compute degree stats
  std::vector<uint64_t> degrees;
  degrees.reserve(num_nodes);
  for (HypernodeID node: data.n2_list) {
    degrees.push_back(graph.nodeDegree(node));
    membership.set(node);
  }
  result.degree_stats = createStats(degrees, degrees.size() >= 20000);

  // compute locality stats and related values
  std::vector<double> locality_values;
  locality_values.reserve(num_nodes);
  for (HypernodeID node: data.n2_list) {
    uint64_t local_edges = 0;
    for (HyperedgeID edge: graph.incidentEdges(node)) {
      HypernodeID neighbor = graph.edgeTarget(edge);
      if (data.isInN2(neighbor)) {
        result.to_n1n2_edges++;
        local_edges++;
        if (!data.isInN2Exactly(neighbor)) {
          result.to_n1n2_edges++;
        }
      } else {
        result.out_edges++;
      }
    }
    HypernodeID node_degree = graph.nodeDegree(node);
    uint64_t divisor = std::min(num_nodes + static_cast<HypernodeID>(data.n1_list.size()), node_degree);
    ALWAYS_ASSERT(divisor != 0);
    locality_values.push_back(static_cast<double>(local_edges) / static_cast<double>(divisor));
  }
  result.to_n1n2_edges /= 2;  // doubly counted
  result.locality_stats = createStats(locality_values, locality_values.size() >= 20000);

  // modularity computations
  uint64_t internal_edges = 2 * (n1_data.degree + duplicated_degree + n1_data.to_n1_edges + result.to_n1n2_edges);
  uint64_t all_edges = internal_edges + result.out_edges;
  std::vector<HypernodeID> combined_neighborhood(data.n2_list);
  combined_neighborhood.insert(combined_neighborhood.end(), data.n1_list.cbegin(), data.n1_list.cend());
  auto [modularity, max_modularity, n_removed] = maxModularitySubset(graph, internal_edges, all_edges, membership, combined_neighborhood);
  result.modularity = modularity;
  result.max_modularity = max_modularity;
  result.max_modularity_size = data.n1_list.size() + result.size - n_removed;

  // TODO: community overlap?
  return result;
}

std::vector<std::tuple<HypernodeID, N1Features, N2Features>> computeNodeFeatures(const StaticGraph& graph, const GlobalFeatures& global, const std::vector<uint64_t>& global_degrees) {
  std::vector<std::tuple<HypernodeID, N1Features, N2Features>> result;
  result.resize(graph.initialNumNodes());

  tbb::enumerable_thread_specific<NeighborhoodComputation> n_comps(graph.initialNumNodes());
  tbb::enumerable_thread_specific<CliqueComputation> clique_comps(graph.initialNumNodes());
  tbb::enumerable_thread_specific<FastResetArray> membership_buffer(graph.initialNumNodes());
  graph.doParallelForAllNodes([&](HypernodeID node) {
    NeighborhoodComputation& local_compute = n_comps.local();
    local_compute.reset();
    NeighborhoodResult neighborhood = local_compute.computeNeighborhood(graph, std::array{node}, true);
    N1Features n1_features = n1FeaturesFromNeighborhood(graph, global, global_degrees, neighborhood, &clique_comps.local(), 0, membership_buffer.local());
    N2Features n2_features = n2FeaturesFromNeighborhood(graph, neighborhood, n1_features, 0, membership_buffer.local());
    result[node] = {node, n1_features, n2_features};
  });

  return result;
}

std::vector<std::tuple<HypernodeID, HypernodeID, EdgeFeatures>> computeEdgeFeatures(const StaticGraph& graph, const GlobalFeatures& global, const std::vector<uint64_t>& global_degrees,
                                                                                    const std::vector<std::pair<ds::Clustering, double>>& community_stack, bool skip_comm_1) {
  tbb::enumerable_thread_specific<std::vector<std::tuple<HypernodeID, HypernodeID, EdgeFeatures>>> result_list;
  tbb::enumerable_thread_specific<NeighborhoodComputation> base_neighborhood(graph.initialNumNodes());
  tbb::enumerable_thread_specific<NeighborhoodComputation> result_neighborhood(graph.initialNumNodes());
  tbb::enumerable_thread_specific<CliqueComputation> clique_comps(graph.initialNumNodes());
  tbb::enumerable_thread_specific<FastResetArray> membership_buffer(graph.initialNumNodes());
  graph.doParallelForAllNodes([&](HypernodeID u) {
    NeighborhoodComputation& base_compute = base_neighborhood.local();
    base_compute.reset();
    NeighborhoodResult u_neighborhood = base_compute.computeNeighborhood(graph, std::array{u}, false);

    for (HyperedgeID edge: graph.incidentEdges(u)) {
      ALWAYS_ASSERT(graph.edgeSize(edge) == 2);
      ALWAYS_ASSERT(graph.edgeWeight(edge) == 1);
      HypernodeID v = graph.edgeTarget(edge);
      HypernodeID deg_u = graph.nodeDegree(u);
      HypernodeID deg_v = graph.nodeDegree(v);
      if (deg_u < deg_v || (deg_u == deg_v && u < v)) {
        continue;
      }

      EdgeFeatures result;
      NeighborhoodComputation& result_compute = result_neighborhood.local();
      result_compute.reset();
      NeighborhoodResult union_result = result_compute.computeNeighborhood(graph, std::array{u, v}, false);
      result.union_features = n1FeaturesFromNeighborhood(graph, global, global_degrees, union_result, nullptr,
                                                         graph.nodeDegree(u) + graph.nodeDegree(v) - union_result.n1_list.size(), membership_buffer.local());
      result_compute.reset();
      // our slightly hacky intersect computation
      NeighborhoodResult intersect_result = result_compute.computeNeighborhood(graph, std::array{v}, false,
        [&](HypernodeID node){ return u_neighborhood.isInN1Exactly(node); });
      intersect_result.roots[0] = u;
      result.intersect_features = n1FeaturesFromNeighborhood(graph, global, global_degrees, intersect_result, &clique_comps.local(), 0, membership_buffer.local());

      double intersect_size = result.intersect_features.degree;
      if (deg_u == 1 || deg_v == 1) {
        result.locality = 1;
        result.strawman_similarity = 1;
      } else {
        result.locality = intersect_size / static_cast<double>(std::min(deg_u, deg_v) - 1);
        result.strawman_similarity = 1.0 / static_cast<double>((deg_u - 1) * (deg_v - 1));
      }
      result.jaccard_index = (result.union_features.degree == 0) ? 1 : intersect_size / static_cast<double>(result.union_features.degree);
      result.cosine_similarity = intersect_size / std::sqrt(deg_u * deg_v);
      if (result.intersect_features.degree == 0) {
        result.dice_similarity = 0;
      } else {
        HypernodeID dice_divisor = result.intersect_features.degree + result.intersect_features.to_n1_edges + result.intersect_features.to_n2_edges;
        result.dice_similarity = intersect_size / static_cast<double>(dice_divisor);
      }

      // community detection
      auto equal_communities = [&](size_t i) {
        const auto& [clustering, _] = community_stack.at(community_stack.size() - i - 1);
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

      result_list.local().emplace_back(u, v, result);
    }
  });

  std::vector<std::tuple<HypernodeID, HypernodeID, EdgeFeatures>> result;
  for (const auto& list: result_list) {
    result.insert(result.end(), list.begin(), list.end());
  }
  return result;
}



// ################    Main   ################


int main(int argc, char* argv[]) {
  auto start = std::chrono::high_resolution_clock::now();

  Context context;
  std::string global_out;
  std::string nodes_out;
  std::string edges_out;

  po::options_description options("Options");
  options.add_options()
          ("help", "show help message")
          ("hypergraph,h",
           po::value<std::string>(&context.partition.graph_filename)->value_name("<string>")->required(),
           "Graph Filename")
          ("global,g",
           po::value<std::string>(&global_out)->value_name("<string>")->required(),
           "Output file for global features")
          ("nodes,n",
           po::value<std::string>(&nodes_out)->value_name("<string>")->required(),
           "Output file for node features")
          ("edges,e",
           po::value<std::string>(&edges_out)->value_name("<string>")->required(),
           "Output file for edge features")
          ("input-file-format",
            po::value<std::string>()->value_name("<string>")->notifier([&](const std::string& s) {
              if (s == "hmetis") {
                context.partition.file_format = FileFormat::hMetis;
              } else if (s == "metis") {
                context.partition.file_format = FileFormat::Metis;
              }
            })->default_value("metis"),
            "Input file format: \n"
            " - hmetis : hMETIS hypergraph file format \n"
            " - metis : METIS graph file format")
            ("p-louvain-min-vertex-move-fraction",
             po::value<long double>(&context.preprocessing.community_detection.min_vertex_move_fraction)->value_name(
                     "<long double>")->default_value(0.01),
             "Louvain pass terminates if less than that fraction of nodes moves during a pass")
            ("p-max-louvain-pass-iterations",
             po::value<uint32_t>(&context.preprocessing.community_detection.max_pass_iterations)->value_name(
                     "<uint32_t>")->default_value(5),
             "Maximum number of iterations over all nodes of one louvain pass");

  po::variables_map cmd_vm;
  po::store(po::parse_command_line(argc, argv, options), cmd_vm);
  if (cmd_vm.count("help") != 0 || argc == 1) {
    LOG << options;
    exit(0);
  }
  po::notify(cmd_vm);

  // context.shared_memory.num_threads = 1;
  // TBBInitializer::instance(context.shared_memory.num_threads);

  std::ofstream global(global_out);
  std::ofstream nodes(nodes_out);
  std::ofstream edges(edges_out);

  // Read Hypergraph
  mt_kahypar_hypergraph_t hypergraph =
    mt_kahypar::io::readInputFile(
      context.partition.graph_filename, PresetType::default_preset,
      InstanceType::graph, context.partition.file_format, true);
  StaticGraph& graph = utils::cast<StaticGraph>(hypergraph);

  double time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
  std::cout << "Starting global feature computation [" << time << "s]" << std::endl;
  LouvainGraph louvain_graph(graph, LouvainEdgeWeight::uniform, StaticGraph::is_graph);
  auto community_stack = community_detection::run_parallel_louvain(louvain_graph, context);
  ALWAYS_ASSERT(community_stack.size() > 0);
  while (community_stack.size() < 3) {
    community_stack.insert(community_stack.begin(), community_stack.front());
  }
  auto [global_features, degrees, skip_comm_1] = computeGlobalFeatures(graph, community_stack);  // does not contain locality
  time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
  std::cout << "Starting node feature computation [" << time << "s]" << std::endl;
  auto node_features = computeNodeFeatures(graph, global_features, degrees);
  time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
  std::cout << "Starting Edge feature computation [" << time << "s]" << std::endl;
  auto edge_features = computeEdgeFeatures(graph, global_features, degrees, community_stack, skip_comm_1);
  time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
  std::cout << "Feature computation complete [" << time << "s]" << std::endl;
  

  tbb::parallel_invoke([&]{
    // compute locality stats
    std::vector<double> locality_values;
    for (auto [u, v, features]: edge_features) {
      if (graph.nodeDegree(u) > 1 && graph.nodeDegree(v) > 1) {
        locality_values.push_back(features.locality);
      }
    }
    global_features.locality_stats = createStats(locality_values, true);

    // output global features
    auto header = global_features.header();
    auto features = global_features.featureList();
    ALWAYS_ASSERT(header.size() == features.size());
    printHeader(global, header);
    printFeatures(global, features);
  }, [&]{
    // output node features
    std::vector<std::string> header {"node_id"};
    std::vector<std::string> n1_header = N1Features::header();
    header.insert(header.end(), n1_header.begin(), n1_header.end());
    std::vector<std::string> n2_header = N2Features::header();
    header.insert(header.end(), n2_header.begin(), n2_header.end());
    printHeader(nodes, header);

    std::vector<Feature> features;
    for (const auto& [id, n1_features, n2_features]: node_features) {
      features.clear();
      features.emplace_back(id, FeatureType::integer);  // (not actually a feature)
      std::vector<Feature> n1 = n1_features.featureList();
      features.insert(features.end(), n1.begin(), n1.end());
      std::vector<Feature> n2 = n2_features.featureList();
      features.insert(features.end(), n2.begin(), n2.end());
      ALWAYS_ASSERT(header.size() == features.size());
      printFeatures(nodes, features);
    }
  }, [&]{
    // output edge features
    std::vector<std::string> header {"id_high_degree", "id_low_degree"};
    std::vector<std::string> edge_header = EdgeFeatures::header();
    header.insert(header.end(), edge_header.begin(), edge_header.end());
    printHeader(edges, header);

    std::vector<Feature> features;
    for (const auto& [u, v, e_features]: edge_features) {
      features.clear();
      features.emplace_back(u, FeatureType::integer);  // (not actually a feature)
      features.emplace_back(v, FeatureType::integer);  // (not actually a feature)
      std::vector<Feature> results = e_features.featureList();
      features.insert(features.end(), results.begin(), results.end());
      ALWAYS_ASSERT(header.size() == features.size());
      printFeatures(edges, features);
    }
  });

  time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
  std::cout << "Done! [" << time << "s]" << std::endl;
  // std::string graph_name = context.partition.graph_filename.substr(
  //   context.partition.graph_filename.find_last_of("/") + 1);

  utils::delete_hypergraph(hypergraph);

  return 0;
}

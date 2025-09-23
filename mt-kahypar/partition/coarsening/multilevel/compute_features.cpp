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

#include "compute_features.h"

#include <array>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_invoke.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_sort.h>

#include "kahypar-resources/datastructure/fast_reset_flag_array.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/static_graph.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/parallel/scalable_sort.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/hypergraph_statistics.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
using FastResetArray = kahypar::ds::FastResetFlagArray<>;

constexpr size_t CACHE_SIZE = 200;

template<typename Container>
double parallel_skew(const Container& data, const double avg, const double stdev, const size_t n) {
    if (stdev == 0) {
        return 0.0;
    }
    return tbb::parallel_reduce(
            tbb::blocked_range<size_t>(UL(0), data.size()), 0.0,
            [&](tbb::blocked_range<size_t>& range, double init) -> double {
            double tmp_skew = init;
            for ( size_t i = range.begin(); i < range.end(); ++i ) {
                tmp_skew += (data[i] - avg) * (data[i] - avg) * (data[i] - avg);
            }
            return tmp_skew;
            }, std::plus<double>()) / (static_cast<double>(n) * std::pow(stdev, 3));
}

template<typename Container>
double parallel_stdev(const Container& data, const double avg, const size_t n) {
    return std::sqrt(tbb::parallel_reduce(
            tbb::blocked_range<size_t>(UL(0), data.size()), 0.0,
            [&](const tbb::blocked_range<size_t>& range, double init) -> double {
            double tmp_stdev = init;
            for ( size_t i = range.begin(); i < range.end(); ++i ) {
                tmp_stdev += (data[i] - avg) * (data[i] - avg);
            }
            return tmp_stdev;
            }, std::plus<double>()) / static_cast<double>(n));
}

template<typename Container>
double parallel_avg(const Container& data, const size_t n) {
    return tbb::parallel_reduce(
            tbb::blocked_range<size_t>(UL(0), data.size()), 0.0,
            [&](const tbb::blocked_range<size_t>& range, double init) -> double {
            double tmp_avg = init;
            for ( size_t i = range.begin(); i < range.end(); ++i ) {
                tmp_avg += static_cast<double>(data[i]);
            }
            return tmp_avg;
            }, std::plus<double>()) / static_cast<double>(n);
}

template <typename Container>
double median(const Container& vec) {
  double median = 0.0;
  if ((vec.size() % 2) == 0) {
    median = static_cast<double>((vec[vec.size() / 2] + vec[(vec.size() / 2) - 1])) / 2.0;
  } else {
    median = vec[vec.size() / 2];
  }
  return median;
}

template <typename Container>
std::pair<double, double> firstAndThirdQuartile(const Container& vec) {
  if (vec.size() > 1) {
    const size_t size_mod_4 = vec.size() % 4;
    const size_t M = vec.size() / 2;
    const size_t ML = M / 2;
    const size_t MU = M + ML;
    double first_quartile = 0.0;
    double third_quartile = 0.0;
    if (size_mod_4 == 0 || size_mod_4 == 1) {
      first_quartile = (vec[ML] + vec[ML - 1]) / 2;
      third_quartile = (vec[MU] + vec[MU - 1]) / 2;
    } else if (size_mod_4 == 2 || size_mod_4 == 3) {
      first_quartile = vec[ML];
      third_quartile = vec[MU];
    }
    return std::make_pair(first_quartile, third_quartile);
  } else {
    return std::make_pair(0.0, 0.0);
  }
}

double ilog2(size_t value, const std::array<double, CACHE_SIZE>& logCache) {
  return value < logCache.size() ? logCache[value] : log2(value);
}

double scaledEntropyFromOccurenceCounts(const ds::DynamicSparseMap<int32_t, int32_t>& occurence, size_t total,
                                        const std::array<double, CACHE_SIZE>& logCache) {
  double entropy = 0;
  size_t count = 0;
  const double factor = 1.0 / static_cast<double>(total);
  const double sub = ilog2(total, logCache);
  for (auto& element : occurence) {
    double log_p_x = ilog2(element.value, logCache) - sub;
    double summand = factor * static_cast<double>(element.value) * log_p_x;
    entropy -= summand;
    count++;
  }

  // scale by log of number of categories
  double scaling = ilog2(count, logCache);
  return scaling == 0 ? 0 : (double)entropy / scaling;
}

// double scaledEntropy(const std::vector<double>& distribution) {
//   ds::DynamicSparseMap<int64_t, int64_t> occurence;
//   for (double value : distribution) {
//     // snap to 2 digits after decimal point
//     int64_t snap = static_cast<int64_t>(std::round(100 * value));
//     occurence[snap]++;
//   }
//   return scaledEntropyFromOccurenceCounts(occurence, distribution.size());
// }

template <typename Container>
double scaledEntropy(const Container& distribution, ds::DynamicSparseMap<int32_t, int32_t>& occurenceBuffer, bool parallel,
                     const std::array<double, CACHE_SIZE>& logCache) {
  // TODO: optimize this by exploiting sortedness? ==> doesn't seem to result in a speedup
  if (parallel) {
    tbb::enumerable_thread_specific<ds::DynamicSparseMap<int32_t, int32_t>> localOccurenceBuffer;
    tbb::parallel_for(UL(0), distribution.size(), [&](size_t i) {
      auto value = distribution[i];
      localOccurenceBuffer.local()[value]++;
    });

    occurenceBuffer.clear();
    for (const ds::DynamicSparseMap<int32_t, int32_t>& buffer: localOccurenceBuffer) {
      for (auto& element : buffer) {
        occurenceBuffer[element.key] += element.value;
      }
    }
  } else {
    occurenceBuffer.clear();
    for (auto value : IteratorRange(distribution.cbegin(), distribution.cend())) {
      occurenceBuffer[value]++;
    }
  }
  return scaledEntropyFromOccurenceCounts(occurenceBuffer, distribution.size(), logCache);
}

template <typename Container, typename T>
Statistic<T> createStats(Container& data, bool parallel, ds::DynamicSparseMap<int32_t, int32_t>& buffer, const std::array<double, CACHE_SIZE>& logCache) {
  Statistic<T> stats;
  if (!data.empty()) {
    double avg = 0;
    double stdev = 0;
    double skew = 0;
    if (parallel) {
      parallel::scalable_sort(data, std::less<>());
      avg = parallel_avg(data, data.size());
      stdev = parallel_stdev(data, avg, data.size());
      skew = parallel_skew(data, avg, stdev, data.size());
    } else {
      std::sort(data.begin(), data.end());
      for (auto val: data) {
        avg += val;
      }
      avg = avg / static_cast<double>(data.size());
      for (auto val: data) {
        stdev += (val - avg) * (val - avg);
        skew += (val - avg) * (val - avg) * (val - avg);
      }
      stdev = std::sqrt(stdev / static_cast<double>(data.size()));
      skew = stdev == 0 ? 0.0 : skew / (static_cast<double>(data.size()) * std::pow(stdev, 3));
    }

    const auto quartiles = firstAndThirdQuartile(data);
    stats.min = data[0];
    stats.q1 = quartiles.first;
    stats.med = median(data);
    stats.q3 = quartiles.second;
    stats.max = data.back();
    stats.avg = avg;
    stats.sd = stdev;
    stats.skew = skew;
    stats.entropy = scaledEntropy(data, buffer, parallel, logCache);
  }
  return stats;
}

double degreeQuantile(const ds::Array<uint32_t>& global_degrees, HypernodeID degree) {
  size_t lower = std::lower_bound(global_degrees.cbegin(), global_degrees.cend(), degree) - global_degrees.cbegin();  // always < n
  size_t upper = std::upper_bound(global_degrees.cbegin(), global_degrees.cend(), degree) - global_degrees.cbegin();  // always >= 1
  return (static_cast<double>(lower + upper)) / (2 * static_cast<double>(global_degrees.size()) - 1);
}

std::array<double, CACHE_SIZE> precomputeLogs() {
  std::array<double, CACHE_SIZE> result;
  for (size_t i = 0; i < result.size(); ++i) {
    result[i] = log2(static_cast<double>(i));
  }
  return result;
}

std::array<float, CACHE_SIZE> precomputeQuantiles(const ds::Array<uint32_t>& global_degrees) {
  std::array<float, CACHE_SIZE> result;
  for (size_t i = 0; i < result.size(); ++i) {
    result[i] = degreeQuantile(global_degrees, i);
  }
  return result;
}


// modularity, max_modularity, n_removed
std::tuple<double, double, uint64_t> maxModularitySubset(const ds::StaticGraph& graph, uint64_t internal_edges, uint64_t all_edges,
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

void cheapN1Features(const ds::StaticGraph& graph, N1Features& result, HypernodeID node,
                     const GlobalFeatures& global, const ds::Array<uint32_t>& global_degrees,
                     ds::DynamicSparseMap<int32_t, int32_t>& occurenceBuffer, vec<uint32_t>& degreeBuffer,
                     const std::array<float, CACHE_SIZE>& quantileCache, const std::array<double, CACHE_SIZE>& logCache) {
  HypernodeID num_nodes = graph.nodeDegree(node);
  result.degree = num_nodes;
  result.inverse_degree = num_nodes <= 1 ? 1.0 : 1.0 / static_cast<double>(num_nodes - 1);
  result.degree_quantile = result.degree < quantileCache.size() ? quantileCache[result.degree] : degreeQuantile(global_degrees, result.degree);

  // compute degree stats
  degreeBuffer.clear();
  degreeBuffer.reserve(num_nodes);
  for (HyperedgeID edge: graph.incidentEdges(node)) {
    HypernodeID curr_node = graph.edgeTarget(edge);
    degreeBuffer.push_back(graph.nodeDegree(curr_node));
  }
  result.degree_stats = createStats<vec<uint32_t>, uint32_t>(degreeBuffer, degreeBuffer.size() >= 20000, occurenceBuffer, logCache);
  double chi_sq = 0.0;
  for (uint32_t d: degreeBuffer) {
    chi_sq += (static_cast<double>(d) - global.degree_stats.avg) * (static_cast<double>(d) - global.degree_stats.avg) / global.degree_stats.avg;
  }
  result.chi_squared_degree_deviation = chi_sq;

  // min contracted degree stats
  HyperedgeWeight out_edges = result.degree;
  HypernodeWeight node_weight = 1;
  for (HyperedgeID d: degreeBuffer) {
    if ((static_cast<HyperedgeWeight>(d) - 2) * node_weight < out_edges) {
      out_edges += d - 2;
      node_weight++;
    }
  }
  result.min_contracted_degree = static_cast<double>(out_edges) / static_cast<double>(node_weight);
  result.min_contracted_degree_size = node_weight;
}

void computeCheapN1Features(const ds::StaticGraph& graph, ds::Array<N1Features>& n1_features, const GlobalFeatures& global,
                            const ds::Array<uint32_t>& global_degrees, const std::array<float, CACHE_SIZE>& quantileCache, const std::array<double, CACHE_SIZE>& logCache) {
  tbb::enumerable_thread_specific<ds::DynamicSparseMap<int32_t, int32_t>> localOccurenceBuffer;
  tbb::enumerable_thread_specific<vec<uint32_t>> localDegreeBuffer;
  graph.doParallelForAllNodes([&](const HypernodeID hn) {
    ASSERT(hn < n1_features.size());
    n1_features[hn] = N1Features{};
    cheapN1Features(graph, n1_features[hn], hn, global, global_degrees, localOccurenceBuffer.local(), localDegreeBuffer.local(), quantileCache, logCache);
  });
}

void expensiveN1Features(const ds::StaticGraph& graph, N1Features& result, HypernodeID node, FastResetArray& membership) {
  membership.reset();

  std::vector<HypernodeID> node_list;
  for (HyperedgeID edge: graph.incidentEdges(node)) {
    HypernodeID curr_node = graph.edgeTarget(edge);
    node_list.push_back(curr_node);
    membership.set(curr_node);
  }

  // compute locality stats and related values
  for (HypernodeID node: node_list) {
    uint64_t local_n1_edges = 0;
    for (HyperedgeID edge: graph.incidentEdges(node)) {
      HypernodeID neighbor = graph.edgeTarget(edge);
      if (membership[neighbor]) {
        result.to_n1_edges++;
        local_n1_edges++;
      } else if (neighbor != node) {
        result.to_n2_edges++;
      }
    }
    HypernodeID node_degree = graph.nodeDegree(node);
    if (node_degree == 1) {
      result.d1_nodes++;
    }
  }
  result.to_n1_edges /= 2;  // doubly counted
  result.clustering_coefficient = (result.degree <= 1) ? 0 : static_cast<double>(2 * result.to_n1_edges) / static_cast<double>(result.degree * (result.degree - 1));

  membership.set(node);

  // modularity computations
  uint64_t internal_edges = 2 * (result.degree + result.to_n1_edges);
  uint64_t all_edges = internal_edges + result.to_n2_edges;
  auto [modularity, max_modularity, n_removed] = maxModularitySubset(graph, internal_edges, all_edges, membership, node_list);
  result.modularity = modularity;
  result.max_modularity = max_modularity;
  result.max_modularity_size = result.degree - n_removed;
}

void computeExpensiveN1Features(const ds::StaticGraph& graph, ds::Array<N1Features>& n1_features) {
  tbb::enumerable_thread_specific<FastResetArray> localMembership(graph.initialNumNodes());
  graph.doParallelForAllNodes([&](const HypernodeID hn) {
    ASSERT(hn < n1_features.size());
    expensiveN1Features(graph, n1_features[hn], hn, localMembership.local());
  });
}

std::tuple<GlobalFeatures, bool> computeGlobalFeatures(const ds::StaticGraph& graph, ds::Array<uint32_t>& node_degrees,
                                                       const std::vector<std::tuple<ds::Clustering, HypernodeID, double>>* comm_ptr,
                                                       const std::array<double, CACHE_SIZE>& logCache) {
  GlobalFeatures features;
  ds::DynamicSparseMap<int32_t, int32_t> occurence_buffer;

  HypernodeID num_nodes = graph.initialNumNodes();
  HyperedgeID num_edges = ds::StaticGraph::is_graph ? graph.initialNumEdges() / 2 : graph.initialNumEdges();
  ALWAYS_ASSERT(num_nodes == node_degrees.size());
  Statistic<uint32_t> degree_stats = createStats<ds::Array<uint32_t>, uint32_t>(node_degrees, true, occurence_buffer, logCache);
  features.n = num_nodes;
  features.m = num_edges;
  features.degree_stats = degree_stats;
  features.irregularity = degree_stats.sd / degree_stats.avg;

  // compute exp_median_degree via suffix sum
  uint64_t pins = graph.initialNumPins();
  uint64_t count = 0;
  for (size_t i = num_nodes; i > 0; --i) {
    count += node_degrees[i-1];
    if (count >= pins / 2) {
      features.exp_median_degree = node_degrees[i-1];
      break;
    }
  }

  // modularity features
  bool skip_comm_1 = false;
  if (comm_ptr != nullptr) {
    const auto& community_stack = *comm_ptr;
    ALWAYS_ASSERT(community_stack.size() >= 3);

    auto modularity_features = [&](size_t i) {
      const auto& [_c, comm_count, modularity] = community_stack.at(community_stack.size() - i - 1);
      return std::make_pair(comm_count, modularity);
    };

    std::tie(features.n_communities_0, features.modularity_0) = modularity_features(0);
    std::tie(features.n_communities_1, features.modularity_1) = modularity_features(1);
    std::tie(features.n_communities_2, features.modularity_2) = modularity_features(2);
    if (community_stack.size() > 3 && features.n_communities_1 < 2 * features.n_communities_0) {
      // small hack to get more meaningful features
      std::tie(features.n_communities_1, features.modularity_1) = std::tie(features.n_communities_2, features.modularity_2);
      std::tie(features.n_communities_2, features.modularity_2) = modularity_features(3);
      skip_comm_1 = true;
    }
  }

  return {features, skip_comm_1};
}

std::tuple<GlobalFeatures, ds::Array<N1Features>, bool> computeFeatures(const ds::StaticGraph& graph, const Context& context) {
  ds::Array<N1Features> n1_features;
  ds::Array<uint32_t> node_degrees;

  utils::Timer& timer = utils::Utilities::instance().getTimer(context.utility_id);

  timer.start_timer("features", "Compute Features");

  timer.start_timer("features_initialize", "Features Initialize Data");
  n1_features.resizeNoAssign(graph.initialNumNodes());  // cheap since no writes are performed
  node_degrees.resize(graph.initialNumNodes());
  graph.doParallelForAllNodes([&](const HypernodeID& node) {
    node_degrees[node] = graph.nodeDegree(node);
    ASSERT(graph.nodeWeight(node) == 1);
  });
  timer.stop_timer("features_initialize");

  timer.start_timer("features_global", "Compute Global Features");
  auto logCache = precomputeLogs();
  auto [global_features, skip_comm_1] = computeGlobalFeatures(graph, node_degrees, context.preprocessing.community_stack, logCache);
  timer.stop_timer("features_global");

  timer.start_timer("features_cheap", "Compute Cheap Features");
  auto quantileCache = precomputeQuantiles(node_degrees);
  computeCheapN1Features(graph, n1_features, global_features, node_degrees, quantileCache, logCache);
  timer.stop_timer("features_cheap");

  timer.start_timer("features_expensive", "Compute Expensive Features");
  // note: this may only be called after `computeCheapN1Features`
  computeExpensiveN1Features(graph, n1_features);
  timer.stop_timer("features_expensive");

  timer.stop_timer("features");

  return {global_features, std::move(n1_features), skip_comm_1};
}

}  // namespace mt_kahypar

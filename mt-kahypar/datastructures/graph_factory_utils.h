#pragma once

#include <algorithm>
#include <cstddef>
#include <span>
#include <utility>

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/enumerable_thread_specific.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_scan.h>

#include "kahypar-resources/macros.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/exception.h"

namespace mt_kahypar {
namespace ds {

inline void apply_inclusive_prefix_sum(
    std::span<HyperedgeID> node_degrees, const HypernodeID num_nodes,
    [[maybe_unused]] const HypernodeID num_edges) {
  parallel::TBBPrefixSum<HyperedgeID, std::span> inclusive_prefix_sum(
      node_degrees);
  tbb::parallel_scan(tbb::blocked_range<size_t>(UL(0), UI64(num_nodes)),
                     inclusive_prefix_sum);
  ASSERT(inclusive_prefix_sum.total_sum() == 2 * num_edges);
}

inline vec<HyperedgeID>
compute_degree_prefix_sum_from_degrees(const HyperedgeID *node_degrees,
                                       const HypernodeID num_nodes,
                                       const HypernodeID num_edges) {
  // Copy degrees shifted by one element
  vec<HyperedgeID> degree_prefix_sum(num_nodes + 1);
  tbb::parallel_for(tbb::blocked_range<size_t>(UL(0), num_nodes),
                    [&](const tbb::blocked_range<size_t> &range) {
                      auto node_degrees_begin = node_degrees + range.begin();
                      auto node_degrees_end = node_degrees + range.end();
                      auto degree_prefix_sum_begin =
                          degree_prefix_sum.begin() + range.begin() + 1;
                      std::copy(node_degrees_begin, node_degrees_end,
                                degree_prefix_sum_begin);
                    });

  // Apply inclusive prefix sum without the first element
  std::span<HyperedgeID> num_degree_per_vertex =
      std::span(degree_prefix_sum).subspan(1);
  apply_inclusive_prefix_sum(num_degree_per_vertex, num_nodes,
                                    num_edges);
  return degree_prefix_sum;
}

inline vec<HyperedgeID> compute_degree_prefix_sum_from_edge_vector(
    const vec<std::pair<HypernodeID, HypernodeID>> &edge_vector,
    const HypernodeID num_nodes, const HypernodeID num_edges) {
  vec<HyperedgeID> degree_prefix_sum(num_nodes + 1);
  tbb::enumerable_thread_specific<vec<HyperedgeID>> local_degree_per_vertex(
      num_nodes);
  tbb::parallel_for(ID(0), num_edges, [&](const size_t pos) {
    vec<HyperedgeID> &num_degree_per_vertex = local_degree_per_vertex.local();
    const auto &[source, target] = edge_vector[pos];
    if (source == target) {
      throw InvalidInputException("Edge " + STR(pos) +
                                  " has identical source and target: " +
                                  STR(source) + " " + STR(target));
    }

    if (source >= num_nodes) {
      throw InvalidInputException("Edge " + STR(pos) +
                                  " points to invalid node: " + STR(source));
    }
    ++num_degree_per_vertex[source];
    if (target >= num_nodes) {
      throw InvalidInputException("Edge " + STR(pos) +
                                  " points to invalid node: " + STR(target));
    }
    ++num_degree_per_vertex[target];
  });

  // Attribute the degrees to the succeeding vertex
  for (vec<HyperedgeID> &local_degrees : local_degree_per_vertex) {
    tbb::parallel_for(ID(0), num_nodes, [&](const size_t node_id) {
      degree_prefix_sum[node_id + 1] += local_degrees[node_id];
    });
  }

  // Compute the exclusive prefix sum over the degrees. At this point,
  // degree_prefix_sum constains the vertex degrees shifted by one element.
  std::span<HyperedgeID> num_degree_per_vertex =
      std::span(degree_prefix_sum).subspan(1);
  apply_inclusive_prefix_sum(num_degree_per_vertex, num_nodes,
                                    num_edges);
  return degree_prefix_sum;
}

} // namespace ds
} // namespace mt_kahypar

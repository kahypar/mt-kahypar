/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2022 Nikolai Maas <nikolai.maas@student.kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "detect_low_degree_nodes.h"

#include "mt-kahypar/utils/hypergraph_statistics.h"
#include "mt-kahypar/utils/stats.h"

#include <tbb/parallel_sort.h>

#include <algorithm>

namespace mt_kahypar {
namespace star_partitioning {

void detectLowDegreeNodes(const Hypergraph& hg, const Context& context, ds::Clustering& communities) {
    const HypernodeWeight total_weight = hg.totalWeight();
    tbb::enumerable_thread_specific<HyperedgeWeight> incident_weight_local(0);
    std::vector<std::pair<double, HypernodeWeight>> weight_to_gain_rates;
    weight_to_gain_rates.assign(hg.initialNumNodes(), std::make_pair(0, 0));
    tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](HypernodeID u) {
      weight_to_gain_rates[u] = std::make_pair(static_cast<double>(hg.nodeWeight(u))
                                / std::max(static_cast<double>(hg.incidentWeight(u)), 1.0),
                                hg.nodeWeight(u));
      incident_weight_local.local() += hg.incidentWeight(u);
    });

    const HyperedgeWeight total_incident_weight = incident_weight_local.combine(std::plus<>());
    auto map_to_first = std::function([](const std::pair<double, HypernodeWeight>& pair) { return pair.first; });

    const double avg_rate = static_cast<double>(total_weight) / static_cast<double>(total_incident_weight);
    const double stdev_rate = utils::parallel_weighted_stdev(weight_to_gain_rates, avg_rate, total_incident_weight,
      [&](size_t i) {
        return hg.incidentWeight(i);
      }, map_to_first);

    const double stdev_factor = context.preprocessing.community_detection.isolated_nodes_threshold_stdev_factor;
    const double stdev_factor_min = context.preprocessing.community_detection.isolated_nodes_threshold_stdev_factor_min;

    double min_gain = 1 / (avg_rate + stdev_factor * stdev_rate);
    const double min_gain_limit = 1 / (avg_rate + stdev_factor_min * stdev_rate);
    if (context.preprocessing.community_detection.adjust_sd_factor) {
      tbb::parallel_sort(weight_to_gain_rates.begin(), weight_to_gain_rates.end(), [](auto l, auto r) { return l.first > r.first; });
      HypernodeWeight sum = 0;
      for (auto pair: weight_to_gain_rates) {
        sum += pair.second;
        if (sum >= (1 - context.preprocessing.community_detection.sd_factor_core_size_target) * hg.totalWeight()) {
          const double max_inv_gain = pair.first;
          min_gain = std::max(min_gain, 1 / max_inv_gain);
          min_gain = std::min(min_gain, min_gain_limit);
          break;
        }
      }
    }

    tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](HypernodeID u) {
      const double incident_weight = hg.incidentWeight(u);
      if (incident_weight <= 0 || incident_weight <= static_cast<double>(hg.nodeWeight(u)) * min_gain) {
        communities[u] = u;
      }
    });
}

void detectViaObjectiveFunction(const Hypergraph& hg, const Context& context, ds::Clustering& communities,
                                std::function<double(double, double)> objective) {
  const HypernodeWeight total_weight = hg.totalWeight();
  tbb::enumerable_thread_specific<HyperedgeWeight> incident_weight_local(0);
  std::vector<std::pair<double, HypernodeID>> weight_to_gain_rates;
  std::vector<unsigned char> removed_nodes;
  tbb::parallel_invoke([&] {
    weight_to_gain_rates.assign(hg.initialNumNodes(), std::make_pair(0, 0));
  }, [&] {
    removed_nodes.assign(hg.initialNumNodes(), 0);
  });
  tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](HypernodeID u) {
    weight_to_gain_rates[u] = std::make_pair(static_cast<double>(hg.nodeWeight(u))
                              / std::max(static_cast<double>(hg.incidentWeight(u)), 1.0), u);
    incident_weight_local.local() += hg.incidentWeight(u);
  });

  const HyperedgeWeight total_incident_weight = incident_weight_local.combine(std::plus<>());
  auto map_to_first = std::function([](const std::pair<double, HypernodeID>& pair) { return pair.first; });

  const double avg_rate = static_cast<double>(total_weight) / static_cast<double>(total_incident_weight);
  const double stdev_rate = utils::parallel_weighted_stdev(weight_to_gain_rates, avg_rate, total_incident_weight,
    [&](size_t i) {
      return hg.incidentWeight(i);
    }, map_to_first);

  const double stdev_factor = context.preprocessing.community_detection.isolated_nodes_threshold_stdev_factor;
  const double max_inv_gain = avg_rate + stdev_factor * stdev_rate;
  tbb::enumerable_thread_specific<HypernodeWeight> removed_node_weight_local(0);
  tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](HypernodeID u) {
    const double incident_weight = hg.incidentWeight(u);
    if (incident_weight <= 0 || incident_weight <= static_cast<double>(hg.nodeWeight(u)) / max_inv_gain) {
      communities[u] = u;
      removed_nodes[u] = 1;
      removed_node_weight_local.local() += hg.nodeWeight(u);
    }
  });

  tbb::enumerable_thread_specific<HyperedgeWeight> internal_weight_local(0);
  tbb::enumerable_thread_specific<HyperedgeWeight> cut_weight_local(0);
  tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](HypernodeID u) {
    if (removed_nodes[u]) {
      for (const HyperedgeID& e: hg.incidentEdges(u)) {
        if (removed_nodes[hg.edgeTarget(e)]) {
          internal_weight_local.local() += hg.edgeWeight(e);
        } else {
          cut_weight_local.local() += 2 * hg.edgeWeight(e);
        }
      }
    }
  });

  HyperedgeWeight removed_node_weight = removed_node_weight_local.combine(std::plus<>());
  HyperedgeWeight internal_weight = internal_weight_local.combine(std::plus<>());
  HyperedgeWeight cut_weight = cut_weight_local.combine(std::plus<>());
  auto get_objective = [&] {
    const HyperedgeWeight core_edge_weight = total_incident_weight - cut_weight - internal_weight;
    return objective(core_edge_weight / static_cast<double>(hg.totalWeight() - removed_node_weight),
                     internal_weight / static_cast<double>(removed_node_weight));
  };

  tbb::parallel_sort(weight_to_gain_rates.begin(), weight_to_gain_rates.end(),
                      [](auto l, auto r) { return l.first > r.first; });

  auto iter = std::upper_bound(weight_to_gain_rates.begin(), weight_to_gain_rates.end(),
                                std::make_pair(max_inv_gain, 0), [](auto l, auto r) { return l.first > r.first; });
  size_t num_iter = 0;
  for (; iter != weight_to_gain_rates.end(); ++iter) {
    const double old_internal_weight = internal_weight;
    const double old_cut_weight = cut_weight;
    const double old_obj = get_objective();
    const HypernodeID u = iter->second;
    if (!removed_nodes[u]) {
      for (const HyperedgeID& e: hg.incidentEdges(u)) {
        if (removed_nodes[hg.edgeTarget(e)]) {
          internal_weight += 2 * hg.edgeWeight(e);
          cut_weight -= 2 * hg.edgeWeight(e);
        } else {
          cut_weight += 2 * hg.edgeWeight(e);
        }
      }
      removed_node_weight += hg.nodeWeight(u);
    }
    if (get_objective() >= old_obj) {
      internal_weight = old_internal_weight;
      cut_weight = old_cut_weight;
      removed_node_weight -= hg.nodeWeight(u);
    } else {
      communities[u] = u;
      ++num_iter;
    }
  }
}

} // namepace star_partitioning
} // namepace mt_kahypar

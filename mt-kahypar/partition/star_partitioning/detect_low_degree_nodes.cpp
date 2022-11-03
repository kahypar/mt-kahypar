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

namespace mt_kahypar {
namespace star_partitioning {

void detectLowDegreeNodes(const Hypergraph& hg, const Context& context, ds::Clustering& communities) {
    const HypernodeWeight total_weight = hg.totalWeight();
    tbb::enumerable_thread_specific<HyperedgeWeight> incident_weight_local(0);
    std::vector<double> weight_to_gain_rates;
    weight_to_gain_rates.assign(hg.initialNumNodes(), 0);
    tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](HypernodeID u) {
      weight_to_gain_rates[u] = static_cast<double>(hg.nodeWeight(u))
                                / std::max(static_cast<double>(hg.incidentWeight(u)), 1.0);
      incident_weight_local.local() += hg.incidentWeight(u);
    });

    const HyperedgeWeight total_incident_weight = incident_weight_local.combine(std::plus<>());

    const double avg_rate = static_cast<double>(total_weight) / static_cast<double>(total_incident_weight);
    const double stdev_rate = utils::parallel_weighted_stdev(weight_to_gain_rates, avg_rate, total_incident_weight,
      [&](size_t i) {
        return hg.incidentWeight(i);
      });

    const double stdev_factor = context.preprocessing.community_detection.isolated_nodes_threshold_stdev_factor;

    tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](HypernodeID u) {
      const double incident_weight = hg.incidentWeight(u);
      if (incident_weight <= 0 || incident_weight < static_cast<double>(hg.nodeWeight(u)) / (avg_rate + stdev_factor * stdev_rate)) {
        communities[u] = u;
      }
    });
}

} // namepace star_partitioning
} // namepace mt_kahypar

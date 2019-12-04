/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include <tbb/tick_count.h>

#include "mt-kahypar/datastructures/clustering.h"
#include "mt-kahypar/datastructures/graph.h"

namespace mt_kahypar {
class ClusteringStatistics {
 private:
  static constexpr bool debug = false;

  template <typename T>
  static T percentile(double fraction, std::vector<T>& elements) {
    fraction = std::max(fraction, 0.0);
    long ind = std::lround(elements.size() * fraction);
    if (ind < 0) {
      ind = 0;
    }
    if (static_cast<size_t>(ind) >= elements.size()) {
      ind = elements.size() - 1;
    }
    return elements[ind];
  }

  template <typename T>
  static double avg(std::vector<T>& elements) {
    T sum = std::accumulate(elements.begin(), elements.end(), T());
    return static_cast<double>(sum) / static_cast<double>(elements.size());
  }

 public:
  static void printLocalMovingStats(const ds::AdjListGraph& G, ds::Clustering& C) {
    std::vector<size_t> sizeOfCluster(G.numNodes());
    PartitionID maxClusterID = 0;
    for (const PartitionID& c : C) {
      sizeOfCluster[c]++;
      maxClusterID = std::max(maxClusterID, c);
    }
    std::vector<size_t> clusterSizes;
    PartitionID numClusters = maxClusterID + 1;
    for (PartitionID c = 0; c < numClusters; ++c) {
      clusterSizes.push_back(sizeOfCluster[c]);
    }
    std::sort(clusterSizes.begin(), clusterSizes.end());

    std::vector<double> percentiles = { 0.0, 0.05, 0.2, 0.5, 0.8, 0.95, 1.0 };

    DBG << "Local Moving Done";
    DBG << V(G.numNodes()) << V(G.numArcs()) << V(numClusters);
    DBG << "Avg Cluster Size " << avg(clusterSizes);
    std::stringstream ss;
    ss << "Percentile cluster sizes : ";
    for (double p : percentiles) {
      ss << "[p = " << p << " size = " << percentile(p, clusterSizes) << "]" << "  ";
    }
    ss << std::endl;
    DBG << ss.str();
  }
};

class HypergraphCommunityAssignmentStatistics {
 public:
  static void print(const Hypergraph& hg) {
    (void)hg;
    // Quantile and Avg num nodes per comm, num inter/intra hyperedges and pins
    std::cout << "Implement some community analysis numbers" << std::endl;
  }
};
}

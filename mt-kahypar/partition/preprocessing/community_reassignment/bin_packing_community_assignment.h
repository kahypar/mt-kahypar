/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <cmath>

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/preprocessing/community_reassignment/i_community_assignment.h"

namespace mt_kahypar {
namespace preprocessing {
template <typename TypeTraits = Mandatory,
          typename Objective = Mandatory>
class BinPackingCommunityAssignmentT : public ICommunityAssignment {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using TBB = typename TypeTraits::TBB;

  static constexpr bool debug = false;

  struct Community {
    PartitionID community_id;
    HypernodeID objective;
  };

 public:
  BinPackingCommunityAssignmentT(HyperGraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context) { }

  BinPackingCommunityAssignmentT(const BinPackingCommunityAssignmentT&) = delete;
  BinPackingCommunityAssignmentT & operator= (const BinPackingCommunityAssignmentT &) = delete;

  BinPackingCommunityAssignmentT(BinPackingCommunityAssignmentT&&) = delete;
  BinPackingCommunityAssignmentT & operator= (BinPackingCommunityAssignmentT &&) = delete;

 private:
  parallel::scalable_vector<PartitionID> computeAssignmentImpl() {
    // Compute Bin Capacities
    utils::Timer::instance().start_timer("compute_community_mapping", "Compute Community Mapping");
    int used_numa_nodes = TBB::instance().num_used_numa_nodes();
    std::vector<HypernodeID> bin_capacities(used_numa_nodes, 0);
    for (int node = 0; node < used_numa_nodes; ++node) {
      bin_capacities[node] = bin_capacity(node);
      DBG << "Bin Capacity for node" << node << "is" << bin_capacities[node];
    }

    // Sort communities in increasing order of their objective (either vertex or pin count)
    std::vector<Community> communities;
    for (PartitionID community = 0; community < _hg.numCommunities(); ++community) {
      communities.emplace_back(Community { community, Objective::objective(_hg, community) });
    }
    std::sort(communities.begin(), communities.end(),
              [](const Community& lhs, const Community& rhs) {
          return lhs.objective < rhs.objective;
        });

    // Assign communities in round-robin fashion to bins
    std::vector<HypernodeID> current_capacity(used_numa_nodes, 0);
    parallel::scalable_vector<PartitionID> community_assignment(_hg.numCommunities(), -1);
    int current_bin = 0;
    while (!communities.empty()) {
      const Community community = communities.back();
      communities.pop_back();

      // Search for bin with enough space to hold current community
      int start_bin = current_bin;
      bool found = false;
      do {
        ASSERT(current_bin < (int)current_capacity.size());
        if (current_capacity[current_bin] + community.objective <= bin_capacities[current_bin]) {
          found = true;
          break;
        }
        current_bin = (current_bin + 1) % used_numa_nodes;
      } while (current_bin != start_bin);

      if (found) {
        ASSERT(current_capacity[current_bin] + community.objective <= bin_capacities[current_bin]);
        current_capacity[current_bin] += community.objective;
        community_assignment[community.community_id] = current_bin;
        DBG << "Assign community" << community.community_id << "of size" << community.objective << "to node" << current_bin
            << "(Remaining Capacity =" << (bin_capacities[current_bin] - current_capacity[current_bin]) << ")";
      } else {
        // If there is no remaining space to assign community to a node,
        // we assign the community to the node where the overflow is minimal
        HypernodeID overflow = std::numeric_limits<HypernodeID>::max();
        int assigned_bin = -1;
        for (int bin = 0; bin < used_numa_nodes; ++bin) {
          ASSERT(current_capacity[bin] + community.objective > bin_capacities[bin]);
          HypernodeID current_overflow = (current_capacity[bin] + community.objective) - bin_capacities[bin];
          if (current_overflow < overflow) {
            overflow = current_overflow;
            assigned_bin = bin;
          }
        }
        ASSERT(assigned_bin != -1);
        current_capacity[assigned_bin] += community.objective;
        community_assignment[community.community_id] = assigned_bin;
        DBG << "Assign community" << community.community_id << "of size" << community.objective << "to node" << assigned_bin
            << "(Overflow Capacity =" << (current_capacity[assigned_bin] - bin_capacities[assigned_bin]) << ")";
      }

      current_bin = (current_bin + 1) % used_numa_nodes;
    }
    utils::Timer::instance().stop_timer("compute_community_mapping");

    ASSERT(std::count(community_assignment.begin(), community_assignment.end(), -1) == 0, "There are unassigned communities");
    return community_assignment;
  }

  HypernodeID bin_capacity(const int node) const {
    ASSERT(node < TBB::instance().num_used_numa_nodes());
    int total_threads = TBB::instance().total_number_of_threads();
    int threads_of_node = TBB::instance().number_of_threads_on_numa_node(node);
    return std::ceil((((double)threads_of_node) / total_threads) * Objective::total(_hg));
  }

  HyperGraph& _hg;
  const Context& _context;
};

template <typename Objective = Mandatory>
using BinPackingCommunityAssignment = BinPackingCommunityAssignmentT<GlobalTypeTraits, Objective>;
}  // namespace preprocessing
}  // namespace mt_kahypar

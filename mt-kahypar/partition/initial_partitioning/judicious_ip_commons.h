/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2022 Noah Wahl <noah.wahl@student.kit.edu>
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
#include <algorithm>
#include <mt-kahypar/definitions.h>

namespace mt_kahypar {
struct GreedyJudiciousInitialPartitionerStats {

  static constexpr bool debug = false;
  static constexpr bool track_gain_sequence = true;

  size_t num_moved_nodes = 0;
  vec<Gain> gain_sequence;
  vec<size_t> update_hist;

  GreedyJudiciousInitialPartitionerStats(const HypernodeID num_nodes)
      : update_hist(num_nodes, 0) {}

  void print() {
    if (!track_gain_sequence) {
      return;
    }
    ASSERT(num_moved_nodes == gain_sequence.size());
    double zero_gain_ratio = static_cast<double>(std::count_if(gain_sequence.begin(), gain_sequence.end(), [](auto gain) { return gain == 0; })) / gain_sequence.size();
    LOG << zero_gain_ratio;
  }
};

struct GreedyJudiciousInitialPartitionerConfig {

  GreedyJudiciousInitialPartitionerConfig();

  GreedyJudiciousInitialPartitionerConfig(const bool preassign_nodes,
                                          const bool random_selection,
                                          const bool use_judicious_increase,
                                          const bool use_block_load_only)
      : preassign_nodes(preassign_nodes), random_selection(random_selection),
        use_judicious_increase(use_judicious_increase),
        use_block_load_only(use_block_load_only) {
    ASSERT(!(use_judicious_increase && use_block_load_only));
  }

  GreedyJudiciousInitialPartitionerConfig(const Context &context)
      : preassign_nodes(context.initial_partitioning.preassign_nodes),
        random_selection(context.initial_partitioning.random_selection),
        use_judicious_increase(
            context.initial_partitioning.use_judicious_increase),
        use_block_load_only(context.initial_partitioning.use_block_load_only) {
    ASSERT(!(use_judicious_increase && use_block_load_only));
  }

  bool preassign_nodes = false;
  bool random_selection = false;
  bool use_judicious_increase = false;
  bool use_block_load_only = false;
};
inline std::ostream &
operator<<(std::ostream &str,
           const GreedyJudiciousInitialPartitionerConfig &config) {
  str << "Config [ ";
  str << "Preassign Nodes: " << std::boolalpha << config.preassign_nodes;
  str << ", Random Selection: " << std::boolalpha << config.random_selection;
  str << ", Strategy: "
      << (config.use_judicious_increase
              ? "Judicious Increase"
              : (config.use_block_load_only ? "Block Load" : "Penalty"));
  str << " ]" << std::endl;
  return str;
}
} // namespace mt_kahypar

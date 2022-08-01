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

#include "gmock/gmock.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/star_partitioning/approximate.h"
#include "mt-kahypar/datastructures/array.h"

#include <tbb/parallel_sort.h>

using ::testing::Test;

namespace mt_kahypar {
using ds::Array;

class AnApproximate : public Test {
 public:
  void initialize(std::vector<HypernodeWeight>&& weights, std::vector<HyperedgeWeight>&& gains) {
    ASSERT(weights.size() == gains.size());
    _weights = std::move(weights);
    _gains = std::move(gains);

    _sorted_nodes.clear();
    _sorted_nodes.assign(_weights.size(), 0);
    tbb::parallel_for(0UL, _sorted_nodes.size(), [&](const size_t i) {
        _sorted_nodes[i] = i;
    });
    tbb::parallel_sort(_sorted_nodes.begin(), _sorted_nodes.end(),
      [&](const HypernodeID& left, const HypernodeID& right) {
        const HypernodeWeight weight_left = _weights[left];
        const HypernodeWeight weight_right = _weights[right];
        if (weight_left == 0) {
          return false;
        } else if (weight_right == 0) {
          return true;
        }
        return static_cast<double>(_gains[left]) / weight_left
               < static_cast<double>(_gains[right]) / weight_right;
      }
    );
  }

  void applyMinKnapsack(const HypernodeWeight& capacity, const std::vector<HypernodeID>& expected) {
    auto get_gain = [&](HypernodeID node) { return _gains[node]; };
    auto get_weight = [&](HypernodeID node) { return _weights[node]; };
    auto result = star_partitioning::minKnapsack(_sorted_nodes, capacity, get_gain, get_weight);

    ASSERT_EQ(result.size(), expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
        ASSERT_EQ(result[i], expected[i]);
    }
  }

 private:
  parallel::scalable_vector<HypernodeID> _sorted_nodes;
  std::vector<HypernodeWeight> _weights;
  std::vector<HyperedgeWeight> _gains;
};

TEST_F(AnApproximate, EdgeCase) {
    initialize({}, {});
    applyMinKnapsack(0, {});
    applyMinKnapsack(10, {});

    initialize({1, 2}, {2, 1});
    applyMinKnapsack(0, {1, 0});
    applyMinKnapsack(3, {});
}

TEST_F(AnApproximate, SimpleCases) {
    initialize({4, 3, 2}, {5, 4, 3});
    applyMinKnapsack(2, {0, 1});
    applyMinKnapsack(3, {0, 2});
    applyMinKnapsack(4, {0, 2});
    applyMinKnapsack(5, {0});
    applyMinKnapsack(6, {1});
    applyMinKnapsack(7, {2});
    applyMinKnapsack(8, {2});
    applyMinKnapsack(9, {});
    applyMinKnapsack(10, {});
}

TEST_F(AnApproximate, WeightZero) {
    initialize({1, 1, 0, 0}, {1, 2, 1, 0});
    applyMinKnapsack(0, {0, 1});
    applyMinKnapsack(1, {0});
    applyMinKnapsack(2, {});
}

TEST_F(AnApproximate, TestMultipleParts) {
    Array<HypernodeWeight> part_weights;
    part_weights.assign(2, 0);
    std::vector<HypernodeWeight> max_part_weights{5, 5};
    std::vector<HypernodeWeight> weights{3, 1, 3, 2, 1};
    std::vector<HyperedgeWeight> gains{10, 9, 2, 1, 4, 0, 3, 0, 2, 0};

    star_partitioning::Approximate ap(2);
    ap.partition(5, part_weights, max_part_weights,
      [&](HyperedgeWeight* w, const HypernodeID node) { w[0] = gains[2 * node]; w[1] = gains[2 * node + 1]; },
      [&](const HypernodeID node) { return weights[node]; },
      [&](const HypernodeID node, const PartitionID part) {
        if (node == 0 || node == 1 || node == 4) {
            ASSERT_EQ(part, 1);
        } else {
            ASSERT_EQ(part, 0);
        }
      }
    );
}

TEST_F(AnApproximate, TestMultiplePartsParallel) {
    Array<HypernodeWeight> part_weights;
    part_weights.assign(2, 0);
    std::vector<HypernodeWeight> max_part_weights{5, 5};
    std::vector<HypernodeWeight> weights{3, 1, 3, 2, 1};
    std::vector<HyperedgeWeight> gains{10, 9, 2, 1, 4, 0, 3, 0, 2, 0};

    star_partitioning::Approximate ap(2);
    ap.partition(5, part_weights, max_part_weights,
      [&](HyperedgeWeight* w, const HypernodeID node) { w[0] = gains[2 * node]; w[1] = gains[2 * node + 1]; },
      [&](const HypernodeID node) { return weights[node]; },
      [&](const HypernodeID node, const PartitionID part) {
        if (node == 0 || node == 1 || node == 4) {
            ASSERT_EQ(part, 1);
        } else {
            ASSERT_EQ(part, 0);
        }
      }, parallel_tag_t()
    );
}

}  // namespace mt_kahypar

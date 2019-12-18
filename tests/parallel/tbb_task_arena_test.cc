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

#include "gmock/gmock.h"

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/task_group.h"

#include "mt-kahypar/parallel/hardware_topology.h"
#include "mt-kahypar/parallel/tbb_numa_arena.h"
#include "tests/parallel/topology_mock.h"

using ::testing::Test;

namespace mt_kahypar {
namespace parallel {
template <int NUM_NUMA_NODES>
struct Numa {
  constexpr static int NUMA_NODES = NUM_NUMA_NODES;
};

template <typename Numa>
class ATBBNumaArenaTest : public Test {
 private:
  using TopoMock = mt_kahypar::parallel::TopologyMock<Numa::NUMA_NODES>;
  using HwTopology = mt_kahypar::parallel::HardwareTopology<TopoMock, topology_t, node_t>;
  using TBBArena = mt_kahypar::parallel::TBBNumaArena<HwTopology>;
  using SplittedTBBNumaArena = std::pair<std::unique_ptr<TBBArena>, std::unique_ptr<TBBArena>>;

 public:
  ATBBNumaArenaTest() :
    num_threads(HwTopology::instance().num_cpus()) { }

  int expected_num_numa_nodes() const {
    return Numa::NUMA_NODES;
  }

  int num_numa_nodes() const {
    return HwTopology::instance().num_numa_nodes();
  }

  int num_cpus_on_numa_node(int node) const {
    return HwTopology::instance().num_cpus_on_numa_node(node);
  }

  std::vector<int> get_cpus_of_numa_node(int node) const {
    return HwTopology::instance().get_cpus_of_numa_node(node);
  }

  int expected_number_of_threads() const {
    return num_threads;
  }

  int total_number_of_threads() const {
    return TBBArena::instance(num_threads).total_number_of_threads();
  }

  int number_of_threads_on_numa_node(int node) const {
    return TBBArena::instance(num_threads).number_of_threads_on_numa_node(node);
  }

  tbb::task_arena& numa_task_arena(int node) {
    return TBBArena::instance(num_threads).numa_task_arena(node);
  }

  SplittedTBBNumaArena split_tbb_numa_arena(const size_t num_threads_0, const size_t num_threads_1) {
    return TBBArena::instance(num_threads).split_tbb_numa_arena(num_threads_0, num_threads_1);
  }

  void wait(const int node, tbb::task_group& group) {
    TBBArena::instance(num_threads).wait(node, group);
  }

  void terminate() {
    TBBArena::instance(num_threads).terminate();
  }

 private:
  int num_threads;
};

#define SYSTEM_HAS_MORE_THAN_FOUR_CORES true
typedef ::testing::Types<Numa<1>, Numa<2>
                         #if SYSTEM_HAS_MORE_THAN_FOUR_CORES
                         , Numa<4>
                         #endif
                         > NumaNodesTemplate;

TYPED_TEST_CASE(ATBBNumaArenaTest, NumaNodesTemplate);

TYPED_TEST(ATBBNumaArenaTest, ChecksTBBArenaInitialization) {
  ASSERT_EQ(this->expected_num_numa_nodes(), this->num_numa_nodes());
  ASSERT_EQ(this->expected_number_of_threads(), this->total_number_of_threads());
  int total_threads = 0;
  for (int node = 0; node < this->expected_num_numa_nodes(); ++node) {
    ASSERT_EQ(this->num_cpus_on_numa_node(node), this->number_of_threads_on_numa_node(node));
    total_threads += this->number_of_threads_on_numa_node(node);
  }
  ASSERT_EQ(this->expected_number_of_threads(), total_threads);
}

TYPED_TEST(ATBBNumaArenaTest, ChecksThreadsToNumaNodeAssignment) {
  for (int node = 0; node < this->expected_num_numa_nodes(); ++node) {
    std::vector<bool> cpus(this->expected_number_of_threads(), false);
    tbb::task_arena& arena = this->numa_task_arena(node);
    tbb::task_group group;
    arena.execute([&group, &cpus]() {
          group.run([&cpus]() {
            tbb::parallel_for(0, 1000000, [&cpus](const size_t&) {
              cpus[sched_getcpu()] = true;
            });
          });
        });

    arena.execute([&] {
          group.wait();
        });

    std::vector<int> expected_cpus = this->get_cpus_of_numa_node(node);
    std::sort(expected_cpus.begin(), expected_cpus.end());
    std::reverse(expected_cpus.begin(), expected_cpus.end());
    int num_threads = 0;
    for (int cpu_id = 0; cpu_id < this->expected_number_of_threads(); ++cpu_id) {
      if (cpu_id > expected_cpus.back()) {
        expected_cpus.pop_back();
      }
      if (cpus[cpu_id]) {
        ASSERT_FALSE(expected_cpus.empty()) << V(cpu_id) << V(node);
        ASSERT_EQ(cpu_id, expected_cpus.back()) << V(cpu_id) << V(node);
        num_threads++;
      }
    }
    ASSERT_GE(num_threads, 1);
  }
  this->terminate();
}

TYPED_TEST(ATBBNumaArenaTest, VerifiesNumaArenaSplitWithEqualNumberOfThreads) {
  size_t num_threads = this->total_number_of_threads();
  ASSERT_EQ(0, num_threads % 2); // if it fails, your system is weird

  auto splitted_arena = this->split_tbb_numa_arena(num_threads / 2, num_threads / 2);
  auto tbb_arena_0 = std::move(splitted_arena.first);
  auto tbb_arena_1 = std::move(splitted_arena.second);

  ASSERT_EQ(num_threads / 2, tbb_arena_0->total_number_of_threads());
  ASSERT_EQ(num_threads / 2, tbb_arena_1->total_number_of_threads());
  ASSERT_EQ(this->num_numa_nodes(), tbb_arena_0->num_used_numa_nodes());
  ASSERT_EQ(this->num_numa_nodes(), tbb_arena_1->num_used_numa_nodes());
  for ( int node = 0; node < this->num_numa_nodes(); ++node ) {
    ASSERT_EQ(this->num_cpus_on_numa_node(node) / 2, tbb_arena_0->number_of_threads_on_numa_node(node));
    ASSERT_EQ(this->num_cpus_on_numa_node(node) / 2, tbb_arena_1->number_of_threads_on_numa_node(node));
  }
}

TYPED_TEST(ATBBNumaArenaTest, VerifiesNumaArenaSplitWithUnequalNumberOfThreads) {
  size_t num_threads = this->total_number_of_threads();
  int num_numa_nodes = this->num_numa_nodes();
  ASSERT_EQ(0, num_threads % 2); // if it fails, your system is weird

  auto splitted_arena = this->split_tbb_numa_arena(num_threads - 1, 1);
  auto tbb_arena_0 = std::move(splitted_arena.first);
  auto tbb_arena_1 = std::move(splitted_arena.second);

  ASSERT_EQ(num_threads - 1, tbb_arena_0->total_number_of_threads());
  ASSERT_EQ(num_numa_nodes, tbb_arena_1->total_number_of_threads());
  ASSERT_EQ(num_numa_nodes, tbb_arena_0->num_used_numa_nodes());
  ASSERT_EQ(num_numa_nodes, tbb_arena_1->num_used_numa_nodes());
  for ( int node = 0; node < num_numa_nodes - 1; ++node ) {
    ASSERT_EQ(this->num_cpus_on_numa_node(node), tbb_arena_0->number_of_threads_on_numa_node(node));
    ASSERT_EQ(1, tbb_arena_1->number_of_threads_on_numa_node(node));
  }
  ASSERT_EQ(this->num_cpus_on_numa_node(num_numa_nodes - 1) - 1,
            tbb_arena_0->number_of_threads_on_numa_node(num_numa_nodes - 1));
  ASSERT_EQ(1,  tbb_arena_1->number_of_threads_on_numa_node(num_numa_nodes - 1));
}

TYPED_TEST(ATBBNumaArenaTest, VerifiesRecursiveNumaArenaSplit) {
  size_t num_threads = this->total_number_of_threads();
  size_t num_numa_nodes = this->num_numa_nodes();
  ASSERT_EQ(0, num_threads % 2); // if it fails, your system is weird

  auto splitted_arena = this->split_tbb_numa_arena(num_threads / 2, num_threads / 2);
  ASSERT_EQ(num_threads / 2, splitted_arena.first->total_number_of_threads());
  num_threads = splitted_arena.first->total_number_of_threads();
  auto first_splitted_arena = splitted_arena.first->split_tbb_numa_arena(num_threads / 2, num_threads / 2);
  auto tbb_arena_0 = std::move(splitted_arena.first);
  auto tbb_arena_00 = std::move(first_splitted_arena.first);
  auto tbb_arena_01 = std::move(first_splitted_arena.second);

  ASSERT_EQ(std::max(num_threads / 2, num_numa_nodes), tbb_arena_00->total_number_of_threads());
  ASSERT_EQ(std::max(num_threads / 2, num_numa_nodes), tbb_arena_01->total_number_of_threads());
  ASSERT_EQ(this->num_numa_nodes(), tbb_arena_00->num_used_numa_nodes());
  ASSERT_EQ(this->num_numa_nodes(), tbb_arena_01->num_used_numa_nodes());
  for ( int node = 0; node < this->num_numa_nodes(); ++node ) {
    ASSERT_EQ(std::max(tbb_arena_0->number_of_threads_on_numa_node(node) / 2, 1),
              tbb_arena_00->number_of_threads_on_numa_node(node));
    ASSERT_EQ(std::max(tbb_arena_0->number_of_threads_on_numa_node(node) / 2, 1),
              tbb_arena_01->number_of_threads_on_numa_node(node));
  }
}
}  // namespace parallel
}  // namespace mt_kahypar

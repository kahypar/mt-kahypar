//
// Created by mlaupichler on 05.07.21.
//
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "gtest/gtest-death-test.h"

#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/datastructures/async/depth_priority_queue.h"

namespace mt_kahypar::ds {

TEST(ADepthPriorityQueue, EmptyAfterConstruction) {

  uint32_t num_depths = 5;
  auto num_per_depth = std::vector<ContractionGroupID>(num_depths, 1);

  auto dpq = DepthPriorityQueue(num_depths, std::move(num_per_depth));

  ASSERT_TRUE(dpq.unsafe_empty());
  ASSERT_EQ(dpq.unsafe_size(), 0);

}

TEST(ADepthPriorityQueue, SequentialPushPop) {

      uint32_t num_depths = 5;
      auto num_per_depth = std::vector<ContractionGroupID>(num_depths, 10);

      auto dpq = DepthPriorityQueue(num_depths, std::move(num_per_depth));

      dpq.push(1,1);
      dpq.push(2,2);
      dpq.push(0,0);

      ASSERT_EQ(dpq.unsafe_size(), 3);
      ASSERT_FALSE(dpq.unsafe_empty());

      ContractionGroupID popped_id;
      bool popped_success = dpq.try_pop(popped_id);
      ASSERT_TRUE(popped_success);
      ASSERT_EQ(popped_id, 0);

      popped_success = dpq.try_pop(popped_id);
      ASSERT_TRUE(popped_success);
      ASSERT_EQ(popped_id, 1);

      dpq.push(3, 0);
      popped_success = dpq.try_pop(popped_id);
      ASSERT_TRUE(popped_success);
      ASSERT_EQ(popped_id, 3);
      ASSERT_EQ(dpq.unsafe_size(), 1);

      popped_success = dpq.try_pop(popped_id);
      ASSERT_TRUE(popped_success);
      ASSERT_EQ(popped_id, 2);

      ASSERT_TRUE(dpq.unsafe_empty());

      // No success popping if everything empty
      popped_success = dpq.try_pop(popped_id);
      ASSERT_FALSE(popped_success);

      // Success popping again after a push
      dpq.push(4, 0);
      popped_success = dpq.try_pop(popped_id);
      ASSERT_TRUE(popped_success);
      ASSERT_EQ(popped_id, 4);

}

    TEST(ADepthPriorityQueue, ParallelPushPop) {

      uint32_t num_depths = 10;
      ContractionGroupID max_size_per_depth = 10000;
      auto num_per_depth = std::vector<ContractionGroupID>(num_depths, max_size_per_depth);

      auto dpq = DepthPriorityQueue(num_depths, std::move(num_per_depth));

      auto num_repetitions = 20;

      for (int i = 0; i < num_repetitions; ++i) {

        tbb::parallel_invoke([&] {
            tbb::parallel_for(ID(0), max_size_per_depth, [&](const ContractionGroupID id) {
                uint32_t depth = utils::Randomize::instance().getRandomInt(0,num_depths-1, sched_getcpu());
                dpq.push(id, depth);
            });
        }, [&] {
            tbb::parallel_for(ID(0), max_size_per_depth, [&](const ContractionGroupID id) {
                ContractionGroupID dest;
                dpq.try_pop(dest);
            });
        });

        // Get current number of elements in the PQ and spawn as many parallel pops. Make sure they are all successful.
        auto size = dpq.unsafe_size();
        ASSERT(size >= 0);
        ContractionGroupID pos_size = static_cast<ContractionGroupID>(size);
        tbb::parallel_for(ID(0), pos_size, [&](const ContractionGroupID) {
            ContractionGroupID dest;
            bool popped_success = dpq.try_pop(dest);
            ASSERT_TRUE(popped_success);
        });

        ASSERT_TRUE(dpq.unsafe_empty());

        dpq.reset();
      }

    }

    TEST(ADepthPriorityQueue, PushToCompletedDepthDeathTest) {
      testing::FLAGS_gtest_death_test_style="threadsafe";
      uint32_t num_depths = 1;
      auto num_per_depth = std::vector<ContractionGroupID>(num_depths, 1);

      auto dpq = DepthPriorityQueue(num_depths, std::move(num_per_depth));
      dpq.increment_finished(0);

      ASSERT_DEATH(dpq.push(0,0), "");
    }

} // end namespace
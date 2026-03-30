#include <vector>

#include "gmock/gmock.h"

#include "mt-kahypar/partition/evolutionary/evo_logs.h"
#include "mt-kahypar/partition/evolutionary/stopping_criteria.h"

namespace mt_kahypar {

using evolutionary::ImprovementLogEntry;
using evolutionary::IterationLogEntry;

TEST(EvolutionStoppingCriteriaTest, EmptyIterationLogDoesNotStop) {
  std::vector<ImprovementLogEntry> improvements;
  std::vector<IterationLogEntry> iterations;
  double early_rate = -1.0;

  const bool should_stop = evolutionary::stopping::sliding_window_improvement_rate_stop(
    -1,
    improvements,
    iterations,
    2,
    2,
    0.2,
    10,
    early_rate);

  EXPECT_FALSE(should_stop);
}

TEST(EvolutionStoppingCriteriaTest, PlateauFallbackStopsWhenNoImprovementForTooLong) {
  const std::vector<ImprovementLogEntry> improvements = {
    {1000, 10, 100.0, "Insert"}
  };
  const std::vector<IterationLogEntry> iterations = {
    {5, 900, 120.0},
    {25, 1900, 95.0}
  };
  double early_rate = -1.0;

  const bool should_stop = evolutionary::stopping::sliding_window_improvement_rate_stop(
    10,
    improvements,
    iterations,
    2,
    2,
    0.2,
    10,
    early_rate);

  EXPECT_TRUE(should_stop);
}

TEST(EvolutionStoppingCriteriaTest, InsufficientImprovementsDoesNotStop) {
  const std::vector<ImprovementLogEntry> improvements = {
    {1000, 1, 100.0, "Insert"}
  };
  const std::vector<IterationLogEntry> iterations = {
    {1, 1000, 100.0},
    {2, 2000, 99.0}
  };
  double early_rate = -1.0;

  const bool should_stop = evolutionary::stopping::sliding_window_improvement_rate_stop(
    1,
    improvements,
    iterations,
    2,
    2,
    0.2,
    100,
    early_rate);

  EXPECT_FALSE(should_stop);
}

TEST(EvolutionStoppingCriteriaTest, InitialOnlyImprovementsDoNotTriggerRateBasedStop) {
  const std::vector<ImprovementLogEntry> improvements = {
    {1000, 1, 100.0, "Initial"},
    {2000, 2, 90.0, "Initial"},
    {3000, 3, 89.0, "Initial"}
  };
  const std::vector<IterationLogEntry> iterations = {
    {1, 1000, 100.0},
    {2, 2000, 90.0},
    {3, 3000, 89.0}
  };
  double early_rate = -1.0;

  const bool should_stop = evolutionary::stopping::sliding_window_improvement_rate_stop(
    3,
    improvements,
    iterations,
    2,
    2,
    0.2,
    100,
    early_rate);

  EXPECT_FALSE(should_stop);
  EXPECT_LT(early_rate, 0.0);
}

TEST(EvolutionStoppingCriteriaTest, LowRecentRateComparedToEarlyRateStops) {
  const std::vector<ImprovementLogEntry> improvements = {
    {1000, 1, 100.0, "Insert"},
    {2000, 2, 90.0, "Insert"},
    {10000, 10, 89.0, "Insert"}
  };
  const std::vector<IterationLogEntry> iterations = {
    {1, 1000, 100.0},
    {2, 2000, 90.0},
    {10, 10000, 89.0}
  };
  double early_rate = -1.0;

  const bool should_stop = evolutionary::stopping::sliding_window_improvement_rate_stop(
    10,
    improvements,
    iterations,
    2,
    2,
    0.2,
    100,
    early_rate);

  EXPECT_TRUE(should_stop);
  EXPECT_GT(early_rate, 0.0);
}

TEST(EvolutionStoppingCriteriaTest, HighRecentRateComparedToEarlyRateDoesNotStop) {
  const std::vector<ImprovementLogEntry> improvements = {
    {1000, 1, 100.0, "Insert"},
    {2000, 2, 90.0, "Insert"},
    {3000, 3, 80.0, "Insert"}
  };
  const std::vector<IterationLogEntry> iterations = {
    {1, 1000, 100.0},
    {2, 2000, 90.0},
    {3, 3000, 80.0}
  };
  double early_rate = -1.0;

  const bool should_stop = evolutionary::stopping::sliding_window_improvement_rate_stop(
    3,
    improvements,
    iterations,
    2,
    2,
    0.2,
    100,
    early_rate);

  EXPECT_FALSE(should_stop);
  EXPECT_GT(early_rate, 0.0);
}

TEST(EvolutionStoppingCriteriaTest, RepeatedCallsReevaluateRecentWindowAndCanStop) {
  std::vector<ImprovementLogEntry> improvements = {
    {1000, 1, 100.0, "Insert"},
    {2000, 2, 90.0, "Insert"},
    {3000, 3, 85.0, "Insert"}
  };
  std::vector<IterationLogEntry> iterations = {
    {1, 1000, 100.0},
    {2, 2000, 90.0},
    {3, 3000, 85.0}
  };
  double early_rate = -1.0;

  const bool should_stop_first = evolutionary::stopping::sliding_window_improvement_rate_stop(
    3,
    improvements,
    iterations,
    2,
    2,
    0.2,
    1000,
    early_rate);

  EXPECT_FALSE(should_stop_first);
  ASSERT_GT(early_rate, 0.0);

  improvements.push_back({10000, 10, 84.0, "Insert"});
  iterations.push_back({10, 10000, 84.0});

  const bool should_stop_second = evolutionary::stopping::sliding_window_improvement_rate_stop(
    10,
    improvements,
    iterations,
    2,
    2,
    0.2,
    1000,
    early_rate);

  EXPECT_TRUE(should_stop_second);
}

TEST(EvolutionStoppingCriteriaTest, PlateauBoundaryEqualityDoesNotStopButStopsWhenAboveMaxIter) {
  const std::vector<ImprovementLogEntry> improvements = {
    {1000, 10, 100.0, "Insert"}
  };
  std::vector<IterationLogEntry> iterations = {
    {10, 1000, 100.0},
    {20, 2000, 101.0}
  };
  double early_rate = -1.0;

  const bool should_stop_first = evolutionary::stopping::sliding_window_improvement_rate_stop(
    10,
    improvements,
    iterations,
    2,
    2,
    0.2,
    10,
    early_rate);

  EXPECT_FALSE(should_stop_first);

  iterations.push_back({21,3000, 105.0});
  const bool should_stop_second = evolutionary::stopping::sliding_window_improvement_rate_stop(
    10,
    improvements,
    iterations,
    2,
    2,
    0.2,
    10,
    early_rate);
  EXPECT_TRUE(should_stop_second);
}

TEST(EvolutionStoppingCriteriaTest, NonPositiveEarlyRateDoesNotStop) {
  const std::vector<ImprovementLogEntry> improvements = {
    {1000, 1, 100.0, "Insert"},
    {2000, 2, 100.0, "Insert"},
    {3000, 3, 99.0, "Insert"}
  };
  const std::vector<IterationLogEntry> iterations = {
    {1, 1000, 100.0},
    {2, 2000, 100.0},
    {3, 3000, 99.0}
  };
  double early_rate = -1.0;

  const bool should_stop = evolutionary::stopping::sliding_window_improvement_rate_stop(
    3,
    improvements,
    iterations,
    2,
    2,
    0.2,
    100,
    early_rate);

  EXPECT_FALSE(should_stop);
  EXPECT_LE(early_rate, 0.0);
}

TEST(EvolutionStoppingCriteriaTest, MixedInitialAndInsertUsesOnlyEvoImprovements) {
  const std::vector<ImprovementLogEntry> improvements = {
    {1000, 1, 120.0, "Initial"},
    {2000, 2, 115.0, "Initial"},
    {3000, 3, 110.0, "Insert"}
  };
  const std::vector<IterationLogEntry> iterations = {
    {1, 1000, 120.0},
    {2, 2000, 115.0},
    {3, 3000, 110.0},
    {4, 4000, 110.0}
  };
  double early_rate = -1.0;

  const bool should_stop = evolutionary::stopping::sliding_window_improvement_rate_stop(
    3,
    improvements,
    iterations,
    2,
    2,
    0.2,
    100,
    early_rate);

  EXPECT_FALSE(should_stop);
  EXPECT_LT(early_rate, 0.0);
}

TEST(EvolutionStoppingCriteriaTest, RecentWindowLargerThanEvoImprovementsIsClamped) {
  const std::vector<ImprovementLogEntry> improvements = {
    {1000, 1, 100.0, "Insert"},
    {2000, 2, 90.0, "Insert"},
    {10000, 10, 89.0, "Insert"}
  };
  const std::vector<IterationLogEntry> iterations = {
    {1, 1000, 100.0},
    {2, 2000, 90.0},
    {10, 10000, 89.0}
  };
  double early_rate = -1.0;

  const bool should_stop = evolutionary::stopping::sliding_window_improvement_rate_stop(
    10,
    improvements,
    iterations,
    2,
    10,
    0.2,
    100,
    early_rate);

  EXPECT_TRUE(should_stop);
}

TEST(EvolutionStoppingCriteriaTest, RecentWindowOfOneStopsOnAnySlowdown) {
  const std::vector<ImprovementLogEntry> improvements = {
    {1000, 1, 100.0, "Insert"},
    {2000, 2, 90.0, "Insert"},
    {3000, 3, 80.0, "Insert"}
  };
  const std::vector<IterationLogEntry> iterations = {
    {1, 1000, 100.0},
    {2, 2000, 90.0},
    {3, 3000, 80.0}
  };
  double early_rate = -1.0;

  const bool should_stop = evolutionary::stopping::sliding_window_improvement_rate_stop(
    3,
    improvements,
    iterations,
    2,
    1,
    0.2,
    100,
    early_rate);

  EXPECT_TRUE(should_stop);
}

}  // namespace mt_kahypar
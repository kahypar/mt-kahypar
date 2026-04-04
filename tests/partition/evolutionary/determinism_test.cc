#include <utility>

#include "gmock/gmock.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/partition/evo_partitioner.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/partition/refinement/gains/km1/km1_gain_cache.h"

using ::testing::Test;

namespace mt_kahypar {

namespace {

using TypeTraits = StaticHypergraphTypeTraits;
using Hypergraph = typename TypeTraits::Hypergraph;
using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
using GainCache =
	typename GraphAndGainTypes<TypeTraits, Km1GainTypes>::GainCache;

}  // namespace

class EvoDeterminismTest : public Test {
 public:
  EvoDeterminismTest() : hypergraph(), partitioned_hypergraph(), context() {
	context.partition.graph_filename = "../tests/instances/powersim.mtx.hgr";
	context.partition.mode = Mode::direct;
	context.partition.preset_type = PresetType::deterministic;
	context.partition.instance_type = InstanceType::hypergraph;
	context.partition.partition_type = PartitionedHypergraph::TYPE;
	context.partition.epsilon = 0.3;
	context.partition.enable_logging = false;
	context.partition.k = 8;

	context.partition.objective = Objective::km1;
	context.partition.seed = 42;
	context.partition.deterministic = true;
	context.partition_evolutionary = true;
	context.partition.time_limit = 0;
  	context.evolutionary.max_iterations = 4;

	context.evolutionary.population_size = 4;
	context.evolutionary.replace_strategy = EvoReplaceStrategy::worst;
	context.evolutionary.batch_size = 1;
	context.partition.enable_benchmark_mode = true;

	context.shared_memory.num_threads = 4;
	context.evolutionary.num_threads_per_worker = 1;

	context.evolutionary.dynamic_population_size = false;
	context.evolutionary.meta_evo_mode = false;
	context.evolutionary.enable_modified_combine = false;
	context.evolutionary.modified_combine_mixed = false;

	context.evolutionary.iteration = 0;
	context.evolutionary.time_elapsed = std::chrono::milliseconds(0);

	hypergraph = io::readInputFile<Hypergraph>(context.partition.graph_filename,
											   FileFormat::hMetis, true);
	partitioned_hypergraph = PartitionedHypergraph(
		context.partition.k, hypergraph, parallel_tag_t());
  }

  Hypergraph hypergraph;
  PartitionedHypergraph partitioned_hypergraph;
  Context context;
  static constexpr size_t num_repetitions = 3;

  struct RunResult {
	std::vector<PartitionID> assignment;
	HyperedgeWeight quality;
	double imbalance;
	Context context;
  };

  RunResult runEvoOnce(const int seed, const size_t num_threads,
					   const size_t threads_per_worker) {
	Context run_context(context);
	run_context.partition.seed = seed;
	run_context.shared_memory.num_threads = num_threads;
	run_context.evolutionary.num_threads_per_worker = threads_per_worker;
	run_context.evolutionary.iteration = 0;
	run_context.evolutionary.time_elapsed = std::chrono::milliseconds(0);

	Hypergraph hg_copy = hypergraph.copy(parallel_tag_t{});
	PartitionedHypergraph phg =
		EvoPartitioner<TypeTraits>::partition(hg_copy, run_context, nullptr);

	RunResult result;
	result.assignment.resize(hg_copy.initialNumNodes());
	for (HypernodeID u : phg.nodes()) {
	  result.assignment[u] = phg.partID(u);
	}
	result.quality = metrics::quality(phg, run_context);
	result.imbalance = metrics::imbalance(phg, run_context);
	result.context = std::move(run_context);

	return result;
  }
};

TEST_F(EvoDeterminismTest, ExactRepeatSameSeedProducesIdenticalPartition) {
  const int seed = 42;
  const size_t num_threads = 8;
  const size_t threads_per_worker = 8;

  RunResult first = runEvoOnce(seed, num_threads, threads_per_worker);
  for (size_t i = 1; i < num_repetitions; ++i) {
	RunResult other = runEvoOnce(seed, num_threads, threads_per_worker);
	ASSERT_EQ(first.assignment, other.assignment);
	ASSERT_EQ(first.quality, other.quality);
	ASSERT_EQ(first.context.evolutionary.iteration,
			  other.context.evolutionary.iteration);
	ASSERT_DOUBLE_EQ(first.imbalance, other.imbalance);
  }
}

TEST_F(EvoDeterminismTest, SameSeedIsStableEvenAfterDifferentSeedRun) {
  const size_t num_threads = 8;
  const size_t threads_per_worker = 8;

  RunResult baseline = runEvoOnce(42, num_threads, threads_per_worker);
  RunResult ignored = runEvoOnce(777, num_threads, threads_per_worker);
  (void)ignored;
  RunResult again = runEvoOnce(42, num_threads, threads_per_worker);

  ASSERT_EQ(baseline.assignment, again.assignment);
  ASSERT_EQ(baseline.quality, again.quality);
  ASSERT_DOUBLE_EQ(baseline.imbalance, again.imbalance);
}

TEST_F(EvoDeterminismTest, CombineOnlyDeterministicRunIsStableForSameSeed) {
  context.evolutionary.mutation_chance = 0.0f;
  context.evolutionary.enable_modified_combine = false;

  const int seed = 42;
  const size_t num_threads = 8;
  const size_t threads_per_worker = 8;

  RunResult first = runEvoOnce(seed, num_threads, threads_per_worker);
  RunResult second = runEvoOnce(seed, num_threads, threads_per_worker);

  EXPECT_EQ(first.assignment, second.assignment);
  EXPECT_EQ(first.quality, second.quality);
  EXPECT_DOUBLE_EQ(first.imbalance, second.imbalance);
}

TEST_F(EvoDeterminismTest, MutationOnlyDeterministicRunIsStableForSameSeed) {
	context.evolutionary.mutation_chance = 1.0f;
	context.evolutionary.enable_modified_combine = false;

	const int seed = 42;
	const size_t num_threads = 8;
	const size_t threads_per_worker = 8;

	RunResult first = runEvoOnce(seed, num_threads, threads_per_worker);
	RunResult second = runEvoOnce(seed, num_threads, threads_per_worker);

	EXPECT_EQ(first.assignment, second.assignment);
	EXPECT_EQ(first.quality, second.quality);
	EXPECT_DOUBLE_EQ(first.imbalance, second.imbalance);
}

	TEST_F(EvoDeterminismTest, DifferentThreadDeterministicRunIsStableForSameSeed) {

	const int seed = 42;
	RunResult first = runEvoOnce(seed, 8, 2);
	RunResult second = runEvoOnce(seed, 8, 8);

	EXPECT_EQ(first.assignment, second.assignment);
	EXPECT_EQ(first.quality, second.quality);
	EXPECT_DOUBLE_EQ(first.imbalance, second.imbalance);
}

}  // namespace mt_kahypar


#include <chrono>
#include <random>

#include "gmock/gmock.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/partition/evo_partitioner.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/evolutionary/strategy_picker.h"

using ::testing::Test;

namespace mt_kahypar {

namespace {

using TypeTraits = StaticHypergraphTypeTraits;
using Hypergraph = TypeTraits::Hypergraph;
using PartitionedHypergraph = TypeTraits::PartitionedHypergraph;

}  // namespace

class EvoOperatorTest : public Test {
 public:
  EvoOperatorTest() : hypergraph(), context() {
	context.partition.graph_filename = "../tests/instances/powersim.mtx.hgr";
	context.partition.mode = Mode::direct;
	context.partition.preset_type = PresetType::deterministic;
	context.partition.instance_type = InstanceType::hypergraph;
	context.partition.partition_type = PartitionedHypergraph::TYPE;
	context.partition.objective = Objective::km1;
  context.coarsening.algorithm = CoarseningAlgorithm::multilevel_coarsener;
	context.partition.epsilon = 0.3;
	context.partition.k = 8;
	context.partition.seed = 42;
	context.partition.deterministic = false;
	context.partition.enable_logging = false;
	context.partition.enable_benchmark_mode = true;
	context.partition.time_limit = 0;
	context.partition_evolutionary = true;
  context.partition.num_vcycles = 1;

	context.shared_memory.num_threads = 4;

	context.evolutionary.population_size = 4;
	context.evolutionary.replace_strategy = EvoReplaceStrategy::worst;
	context.evolutionary.batch_size = 1;
	context.evolutionary.max_iterations = 1;
	context.evolutionary.num_threads_per_worker = 4;
	context.evolutionary.dynamic_population_size = false;
	context.evolutionary.meta_evo_mode = false;
	context.evolutionary.iteration = 0;
	context.evolutionary.time_elapsed = std::chrono::milliseconds(0);

	hypergraph = io::readInputFile<Hypergraph>(context.partition.graph_filename,
											   FileFormat::hMetis,
											   true);
  }

  Hypergraph hypergraph;
  Context context;

	void prepareEvolutionContext(Hypergraph& hg, Context& run_context) {
		EvoPartitioner<TypeTraits>::configurePreprocessing(hg, run_context);
		EvoPartitioner<TypeTraits>::setupContext(hg, run_context, nullptr);
	}

	void runSingleEvoIterationAndValidate(Context& run_context, Hypergraph hg) {
		prepareEvolutionContext(hg, run_context);

		Population population;
		EvoPartitioner<TypeTraits>::generateInitialPopulation(hg, run_context, nullptr, population);

		ASSERT_EQ(population.size(), run_context.evolutionary.population_size);
		const int iteration_before_evolution = run_context.evolutionary.iteration;
		const HyperedgeWeight best_before = population.bestFitness();

		EvoPartitioner<TypeTraits>::performEvolution(hg, run_context, nullptr, population);

		EXPECT_EQ(run_context.evolutionary.iteration, iteration_before_evolution + 1);
		EXPECT_EQ(population.size(), run_context.evolutionary.population_size);
		EXPECT_LE(population.bestFitness(), best_before);

		Hypergraph hg_copy = hg.copy(parallel_tag_t{});
		PartitionedHypergraph phg(run_context.partition.k, hg_copy, parallel_tag_t{});
		const std::vector<PartitionID> best_partition = population.bestPartitionCopySafe();
		ASSERT_EQ(best_partition.size(), hg_copy.initialNumNodes());

		for (const HypernodeID& hn : hg_copy.nodes()) {
			ASSERT_LT(best_partition[hn], run_context.partition.k);
			phg.setOnlyNodePart(hn, best_partition[hn]);
		}
		phg.initializePartition();
		EXPECT_EQ(population.bestFitness(), metrics::quality(phg, run_context));

		EXPECT_LE(metrics::imbalance(phg, run_context), run_context.partition.epsilon + 1e-9);
	}

	EvoDecision decideNextMove(const Context& decision_context, std::mt19937* rng) {
		return pick::decideNextMove(decision_context, rng);
	}

	Context modifyContext(Context old_context, const ContextModifierParameters &params) {
		return EvoPartitioner<TypeTraits>::modifyContext(old_context, params);;
	}
};

TEST_F(EvoOperatorTest, SingleModifiedCombineRandomPartitionsIterationKeepsPopulationAndPartitionValid) {
  Context run_context(context);
  run_context.evolutionary.mutation_chance = 0.0f;
  run_context.evolutionary.modified_combine_chance = 1.0f;
  run_context.evolutionary.enable_modified_combine = true;
  run_context.evolutionary.modified_combine_use_random_partitions = true;
  run_context.evolutionary.modified_combine_use_degree_sorted_partitions = false;
  run_context.evolutionary.modified_combine_mixed = false;

	Hypergraph hg = hypergraph.copy(parallel_tag_t{});
	runSingleEvoIterationAndValidate(run_context, std::move(hg));

}
TEST_F(EvoOperatorTest, SingleModifiedCombineDegreeSortedPartitionsIterationKeepsPopulationAndPartitionValid) {
	Context run_context(context);
	run_context.evolutionary.mutation_chance = 0.0f;
	run_context.evolutionary.modified_combine_chance = 1.0f;
	run_context.evolutionary.enable_modified_combine = true;
	run_context.evolutionary.modified_combine_use_random_partitions = false;
	run_context.evolutionary.modified_combine_use_degree_sorted_partitions = true;
	run_context.evolutionary.modified_combine_mixed = false;

	Hypergraph hg = hypergraph.copy(parallel_tag_t{});
	runSingleEvoIterationAndValidate(run_context, std::move(hg));

}
TEST_F(EvoOperatorTest, SingleModifiedCombineIterationKeepsPopulationAndPartitionValid) {
	Context run_context(context);
	run_context.evolutionary.mutation_chance = 0.0f;
	run_context.evolutionary.modified_combine_chance = 1.0f;
	run_context.evolutionary.enable_modified_combine = true;
	run_context.evolutionary.modified_combine_use_random_partitions = false;
	run_context.evolutionary.modified_combine_use_degree_sorted_partitions = false;
	run_context.evolutionary.modified_combine_mixed = false;

	Hypergraph hg = hypergraph.copy(parallel_tag_t{});
	runSingleEvoIterationAndValidate(run_context, std::move(hg));

}
TEST_F(EvoOperatorTest, SingleModifiedCombineMixedIterationKeepsPopulationAndPartitionValid) {
	Context run_context(context);
	run_context.evolutionary.mutation_chance = 0.0f;
	run_context.evolutionary.modified_combine_chance = 1.0f;
	run_context.evolutionary.enable_modified_combine = true;
	run_context.evolutionary.modified_combine_use_random_partitions = false;
	run_context.evolutionary.modified_combine_use_degree_sorted_partitions = false;
	run_context.evolutionary.modified_combine_mixed = true;

	Hypergraph hg = hypergraph.copy(parallel_tag_t{});
	runSingleEvoIterationAndValidate(run_context, std::move(hg));

}

TEST_F(EvoOperatorTest, SingleCombineIterationKeepsPopulationAndPartitionValid) {
	Context run_context(context);
	run_context.evolutionary.mutation_chance = 0.0f;
	run_context.evolutionary.enable_modified_combine = false;

	Hypergraph hg = hypergraph.copy(parallel_tag_t{});
	runSingleEvoIterationAndValidate(run_context, std::move(hg));
}

TEST_F(EvoOperatorTest, SingleMutationIterationKeepsPopulationAndPartitionValid) {
	Context run_context(context);
	run_context.evolutionary.mutation_chance = 1.0f;
	run_context.evolutionary.enable_modified_combine = false;

	Hypergraph hg = hypergraph.copy(parallel_tag_t{});
	runSingleEvoIterationAndValidate(run_context, std::move(hg));
}

TEST_F(EvoOperatorTest, DecideNextMoveSelectsExpectedEvoBranch) {
	Context run_context(context);
	run_context.partition.deterministic = true;

	std::mt19937 rng(0);

	run_context.evolutionary.mutation_chance = 1.0f;
	run_context.evolutionary.enable_modified_combine = false;
	EXPECT_EQ(decideNextMove(run_context, &rng),
	          EvoDecision::mutation);

	run_context.evolutionary.mutation_chance = 0.0f;
	run_context.evolutionary.enable_modified_combine = false;
	EXPECT_EQ(decideNextMove(run_context, &rng),
	          EvoDecision::combine);

	run_context.evolutionary.enable_modified_combine = true;
	run_context.evolutionary.modified_combine_chance = 1.0f;
	EXPECT_EQ(decideNextMove(run_context, &rng),
	          EvoDecision::modified_combine);
}

TEST_F(EvoOperatorTest, ModifyContextAppliesModifiedCombineParameters) {
	ContextModifierParameters params;
	params.k = 16;
	params.epsilon = 0.42;
	params.recursive_bipartitioning = true;

	const Context modified_context = modifyContext(context, params);

	EXPECT_EQ(modified_context.partition.k, params.k);
	EXPECT_DOUBLE_EQ(modified_context.partition.epsilon, params.epsilon);
	EXPECT_EQ(modified_context.partition.mode, Mode::recursive_bipartitioning);
}

}  // namespace mt_kahypar


#include <vector>

#include "gurobi_c++.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/io/sql_plottools_serializer.h"
#include "mt-kahypar/partition/partitioner.h"
#include "mt-kahypar/partition/refinement/ilp/ilp_scheduler.h"

using namespace mt_kahypar;

int main(int argc, char* argv[]) {
  Context context;
  processCommandLineInputForILPSolver(context, argc, argv);
  context.partition.is_ilp_solver = true;

  utils::Randomize::instance().setSeed(context.partition.seed);
  if ( context.shared_memory.use_localized_random_shuffle ) {
    utils::Randomize::instance().enableLocalizedParallelShuffle(
      context.shared_memory.shuffle_block_size);
  }

  size_t num_available_cpus = HardwareTopology::instance().num_cpus();
  if ( num_available_cpus < context.shared_memory.num_threads ) {
    WARNING("There are currently only" << num_available_cpus << "cpus available."
      << "Setting number of threads from" << context.shared_memory.num_threads
      << "to" << num_available_cpus);
    context.shared_memory.num_threads = num_available_cpus;
  }

  // Initialize TBB task arenas on numa nodes
  TBBNumaArena::instance(context.shared_memory.num_threads);

  // We set the membind policy to interleaved allocations in order to
  // distribute allocations evenly across NUMA nodes
  hwloc_cpuset_t cpuset = TBBNumaArena::instance().used_cpuset();
  parallel::HardwareTopology<>::instance().activate_interleaved_membind_policy(cpuset);
  hwloc_bitmap_free(cpuset);

  // Print Context Information
  if (context.partition.verbose_output) {
    io::printBanner();
    std::cout << "*******************************************************************************\n"
              << "*                            Partitioning Context                             *\n"
              << "*******************************************************************************\n"
              << context.partition
              << "-------------------------------------------------------------------------------\n"
              << context.refinement.ilp
              << "-------------------------------------------------------------------------------\n"
              << context.shared_memory
              << "-------------------------------------------------------------------------------\n"
              << std::endl;
  }

  // Read Hypergraph
  Hypergraph hypergraph = io::readHypergraphFile(
      context.partition.graph_filename,
      TBBNumaArena::GLOBAL_TASK_GROUP,
      context.preprocessing.stable_construction_of_incident_edges);
  context.setupPartWeights(hypergraph.totalWeight());
  if ( context.partition.verbose_output ) {
    io::printHypergraphInfo(hypergraph, "Input Hypergraph", false);
    io::printStripe();
  }

  // Read Partition
  std::vector<PartitionID> partition;
  io::readPartitionFile(context.partition.graph_partition_filename, partition);
  ASSERT(hypergraph.initialNumNodes() == ID(partition.size()));
  PartitionedHypergraph phg(
    context.partition.k, TBBNumaArena::GLOBAL_TASK_GROUP, hypergraph);
  phg.doParallelForAllNodes([&](const HypernodeID& hn) {
    phg.setOnlyNodePart(hn, partition[hn]);
  });
  phg.initializePartition(TBBNumaArena::GLOBAL_TASK_GROUP);
  if ( context.partition.verbose_output ) {
    io::printPartitioningResults(phg, context, "Input Partition");
    io::printStripe();
  }

  try {
    io::printILPBanner(context);
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    utils::Timer::instance().start_timer("ilp_solver", "ILP Solver");
    ILPScheduler ilp(phg, context);
    ilp.refine();
    utils::Timer::instance().stop_timer("ilp_solver");
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();

    // Print Stats
    auto elapsed_seconds = end - start;
    io::printPartitioningResults(phg, context, elapsed_seconds);

    if ( context.partition.sp_process_output ) {
      std::cout << io::serializer::serialize(phg, context, elapsed_seconds) << std::endl;
    }
  } catch(GRBException e) {
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
  } catch(...) {
    std::cout << "Exception during optimization" << std::endl;
  }

  mt_kahypar::TBBNumaArena::instance().terminate();
  return 0;
}

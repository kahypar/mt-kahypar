#include <vector>

#include "gurobi_c++.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/partitioner.h"
#include "mt-kahypar/partition/refinement/ilp/ilp_solver.h"

using namespace mt_kahypar;

Hypergraph generateRandomHypergraph(const HypernodeID num_hypernodes,
                                    const HyperedgeID num_hyperedges,
                                    const HypernodeID max_edge_size) {
  parallel::scalable_vector<parallel::scalable_vector<HypernodeID>> hyperedges;
  utils::Randomize& rand = utils::Randomize::instance();
  for ( size_t i = 0; i < num_hyperedges; ++i ) {
    parallel::scalable_vector<HypernodeID> net;
    const size_t edge_size = rand.getRandomInt(2, max_edge_size, sched_getcpu());
    for ( size_t i = 0; i < edge_size; ++i ) {
      const HypernodeID pin = rand.getRandomInt(0, num_hypernodes - 1, sched_getcpu());
      if ( std::find(net.begin(), net.end(), pin) == net.end() ) {
        net.push_back(pin);
      }
    }
    hyperedges.emplace_back(std::move(net));
  }
  return HypergraphFactory::construct(
    TBBNumaArena::GLOBAL_TASK_GROUP, num_hypernodes, num_hyperedges, hyperedges);
}

int main() {
  const HypernodeID num_nodes = 10000;
  const HyperedgeID num_edges = 10000;
  const HypernodeID max_edge_size = 20;

  const HypernodeID num_ilp_nodes = 1000;
  vec<HypernodeID> nodes;
  for ( HypernodeID hn = 0; hn < num_ilp_nodes; ++hn ) {
    nodes.push_back(hn);
  }

  Hypergraph hg = generateRandomHypergraph(num_nodes, num_edges, max_edge_size);
  const double epsilon = 0.03;
  const PartitionID k = 8;

  Context context;
  parseIniToContext(context, "/home/tobias/mt-kahypar/config/fast_preset.ini");
  context.partition.k = k;
  context.partition.epsilon = epsilon;
  context.partition.objective = kahypar::Objective::km1;
  context.partition.verbose_output = true;
  context.partition.enable_progress_bar = true;
  context.partition.show_detailed_timings = true;

  // Compute Initial Solution
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  PartitionedHypergraph phg = partition(hg, context);
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed_seconds(end - start);
  io::printPartitioningResults(phg, context, elapsed_seconds);

  try {
    GRBEnv env;
    ILPSolver solver(phg, context, env);
    utils::Timer::instance().start_timer("ilp_solver", "ILP Solver");
    solver.solve(nodes);
    utils::Timer::instance().stop_timer("ilp_solver");
  } catch(GRBException e) {
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
  } catch(...) {
    std::cout << "Exception during optimization" << std::endl;
  }
  end = std::chrono::high_resolution_clock::now();


  // Print Stats
  elapsed_seconds = end - start;
  io::printPartitioningResults(phg, context, elapsed_seconds);

  return 0;
}

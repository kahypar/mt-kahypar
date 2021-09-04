#include "mt-kahypar/partition/coarsening/multilevel_coarsener_base.h"
#include "mt-kahypar/partition/coarsening/nlevel_coarsener_base.h"

#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

#include "mt-kahypar/partition/refinement/rebalancing/rebalancer.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/io/partitioning_output.h"

namespace mt_kahypar {

  void MultilevelCoarsenerBase::finalize() {
    utils::Timer::instance().start_timer("finalize_multilevel_hierarchy", "Finalize Multilevel Hierarchy");
    // Free memory of temporary contraction buffer and
    // release coarsening memory in memory pool
    currentHypergraph().freeTmpContractionBuffer();
    if (_top_level) {
      parallel::MemoryPool::instance().release_mem_group("Coarsening");
    }
    // Construct top level partitioned hypergraph (memory is taken from memory pool)
    *_uncoarseningData.partitioned_hg = PartitionedHypergraph(
            _context.partition.k, _hg, parallel_tag_t());

    // Construct partitioned hypergraphs parallel
    tbb::task_group group;
    // Construct partitioned hypergraph for each coarsened hypergraph in the hierarchy
    for (size_t i = 0; i < _uncoarseningData.hierarchy.size(); ++i) {
      group.run([&, i] {
        (_uncoarseningData.hierarchy)[i].contractedPartitionedHypergraph() = PartitionedHypergraph(
                _context.partition.k, (_uncoarseningData.hierarchy)[i].contractedHypergraph(), parallel_tag_t());
      });
    }
    group.wait();

    // Set the representative partitioned hypergraph for each hypergraph
    // in the hierarchy
    if (_uncoarseningData.hierarchy.size() > 0) {
      (_uncoarseningData.hierarchy)[0].setRepresentativeHypergraph(&*_uncoarseningData.partitioned_hg);
      for (size_t i = 1; i < _uncoarseningData.hierarchy.size(); ++i) {
        (_uncoarseningData.hierarchy)[i].setRepresentativeHypergraph(&(_uncoarseningData.hierarchy)[i - 1].contractedPartitionedHypergraph());
      }
    }
    _uncoarseningData.is_finalized = true;
    utils::Timer::instance().stop_timer("finalize_multilevel_hierarchy");
  }

  void NLevelCoarsenerBase::finalize() {
    // Create compactified hypergraph containing only enabled vertices and hyperedges
    // with consecutive IDs => Less complexity in initial partitioning.
    utils::Timer::instance().start_timer("compactify_hypergraph", "Compactify Hypergraph");
    auto compactification = HypergraphFactory::compactify(_hg);
    *_uncoarseningData.compactified_hg = std::move(compactification.first);
    _uncoarseningData.compactified_hn_mapping = std::move(compactification.second);
    *_uncoarseningData.compactified_phg = PartitionedHypergraph(_context.partition.k, *_uncoarseningData.compactified_hg, parallel_tag_t());
    utils::Timer::instance().stop_timer("compactify_hypergraph");

    _uncoarseningData.is_finalized = true;
  }

  void MultilevelCoarsenerBase::performMultilevelContraction(
          parallel::scalable_vector<HypernodeID>&& communities,
          const HighResClockTimepoint& round_start) {
    ASSERT(!_uncoarseningData.is_finalized);
    Hypergraph& current_hg = currentHypergraph();
    ASSERT(current_hg.initialNumNodes() == communities.size());
    Hypergraph contracted_hg = current_hg.contract(communities);
    const HighResClockTimepoint round_end = std::chrono::high_resolution_clock::now();
    const double elapsed_time = std::chrono::duration<double>(round_end - round_start).count();
    _uncoarseningData.hierarchy.emplace_back(std::move(contracted_hg), std::move(communities), elapsed_time);
  }

}

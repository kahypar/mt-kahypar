#include "mt-kahypar/partition/multicandidate.h"

#include <memory>

#include "include/libmtkahypartypes.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/partition/coarsening/multilevel/multilevel_uncoarsener.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/hypergraph_statistics.h"
#include "mt-kahypar/utils/utilities.h"
#include "tbb/parallel_for.h"


#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/preprocessing/sparsification/degree_zero_hn_remover.h"
#include "mt-kahypar/partition/preprocessing/sparsification/large_he_remover.h"
#include "mt-kahypar/partition/initial_partitioning/pool_initial_partitioner.h"
#include "mt-kahypar/partition/recursive_bipartitioning.h"
#include "mt-kahypar/partition/deep_multilevel.h"
#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/partition/coarsening/nlevel/nlevel_uncoarsener.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/exception.h"

namespace mt_kahypar {

namespace {

void disableTimerAndStats(const Context& context) {
  if (context.type == ContextType::main) {
    utils::Utilities& utils = utils::Utilities::instance();
    parallel::MemoryPool::instance().deactivate_unused_memory_allocations();
    utils.getTimer(context.utility_id).disable();
    utils.getStats(context.utility_id).disable();
  }
}

void enableTimerAndStats(const Context& context) {
  if (context.type == ContextType::main) {
    utils::Utilities& utils = utils::Utilities::instance();
    parallel::MemoryPool::instance().activate_unused_memory_allocations();
    utils.getTimer(context.utility_id).enable();
    utils.getStats(context.utility_id).enable();
  }
}

  template <typename TypeTraits>
  typename TypeTraits::PartitionedHypergraph multicandidate_partitioning(
      typename TypeTraits::Hypergraph& hypergraph,
      const Context& context, 
      const TargetGraph* target_graph) {
    using Hypergraph = typename TypeTraits::Hypergraph;
    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
    // ################## COARSENING ##################
    mt_kahypar::io::printCoarseningBanner(context);
  
    const bool nlevel = context.isNLevelPartitioning();
    UncoarseningData<TypeTraits> uncoarseningData(nlevel, hypergraph, context);
    utils::Stats& stats = utils::Utilities::instance().getStats(context.utility_id);
    utils::Timer& timer = utils::Utilities::instance().getTimer(context.utility_id); 
    timer.start_timer("coarsening", "Coarsening");
    {
      std::unique_ptr<ICoarsener> coarsener = CoarsenerFactory::getInstance().createObject(
        context.coarsening.algorithm, utils::hypergraph_cast(hypergraph),
        context, uncoarsening::to_pointer(uncoarseningData));
      coarsener->coarsen();
      if (context.partition.verbose_output) {
        mt_kahypar_hypergraph_t coarsestHypergraph = coarsener->coarsestHypergraph();
        mt_kahypar::io::printHypergraphInfo(
          utils::cast<Hypergraph>(coarsestHypergraph), context,
          "Coarsened Hypergraph", context.partition.show_memory_consumption);
      }
    }
    timer.stop_timer("coarsening");
  
    // ################## INITIAL PARTITIONING ##################

    PartitionedHypergraph& partitioned_hg = uncoarseningData.coarsestPartitionedHypergraph();
    struct Partition {
      vec<PartitionID> partIDs;
      HyperedgeWeight quality;
      int initialLevel;
    };
    vec<Partition> partition_pool;
    
    auto uncoarsener = std::make_unique<MultilevelUncoarsener<TypeTraits>>(
        hypergraph, context, uncoarseningData, target_graph);

    auto replacePartition = [&partitioned_hg,
                             &uncoarsener](const Partition& partition) {
      partitioned_hg.resetData();
      partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {
        const PartitionID part_id = partition.partIDs[hn];
        partitioned_hg.setOnlyNodePart(hn, part_id);
      });
      partitioned_hg.initializePartition();
      uncoarsener->updateMetrics();
    };

    auto refine = [&](Partition& partition) {
      replacePartition(partition);
      GainCachePtr::resetGainCache(uncoarsener->getGainCache());
      timer.start_timer("refinement", "Refinement");
      uncoarsener->refine();
      partition.quality = metrics::quality(partitioned_hg, Objective::km1);
      timer.stop_timer("refinement");
      partition.partIDs.resize(partitioned_hg.initialNumNodes());
      partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {
        partition.partIDs[hn] = partitioned_hg.partID(hn);
      });
    };

    auto ip = [&](PartitionedHypergraph& phg, int current_level) {
      //partitioned_hg.resetData();
      //uncoarsener->updateMetrics();
      phg.resetData();

      Context ip_context(context);
      ip_context.type = ContextType::initial_partitioning;
      ip_context.refinement = context.initial_partitioning.refinement;
      if (context.initial_partitioning.mode == Mode::direct) {
        // The pool initial partitioner consist of several flat bipartitioning
        // techniques. This case runs as a base case (k = 2) within recursive
        // bipartitioning or the deep multilevel scheme.
        ip_context.partition.verbose_output = false;
        Pool<TypeTraits>::bipartition(phg, ip_context);
      } else if (context.initial_partitioning.mode ==
                 Mode::recursive_bipartitioning) {
        RecursiveBipartitioning<TypeTraits>::partition(
          phg, ip_context, target_graph);
      } else if (context.initial_partitioning.mode == Mode::deep_multilevel) {
        ASSERT(ip_context.partition.objective != Objective::steiner_tree);
        ip_context.partition.verbose_output = false;
        DeepMultilevel<TypeTraits>::partition(phg, ip_context);
      } else {
        throw InvalidParameterException(
            "Undefined initial partitioning algorithm");
      }
      static std::mutex mutex;
      {
        std::lock_guard<std::mutex> lock(mutex);
        Partition partition;
        partition.partIDs.resize(phg.initialNumNodes());
        phg.doParallelForAllNodes([&](const HypernodeID hn) {
          partition.partIDs[hn] = phg.partID(hn);
        });
        refine(partition);
        partition_pool.emplace_back(partition);
        partition_pool.back().initialLevel = current_level;
      }	
    };

    // ################## UNCOARSENING ##################

    // Functions

    auto isBetterThan = [](const Partition& l, const Partition& r) {
      return l.quality < r.quality;
    };

    auto projAndRefine = [&](Partition& partition) {
      replacePartition(partition);
      timer.start_timer("refinement", "Refinement");
      uncoarsener->projectToNextLevelAndRefine();
      partition.quality = metrics::quality(partitioned_hg, Objective::km1);
      timer.stop_timer("refinement");
      partition.partIDs.resize(partitioned_hg.initialNumNodes());
      partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {
        partition.partIDs[hn] = partitioned_hg.partID(hn);
      });
    };

    io::printInitialPartitioningBanner(context);

    uncoarsener->initialize();
    uncoarsener->stepNextLevel();
    int level = 3; //TODO: make this a parameter
    std::cout << "Current Level: " << uncoarsener->currentLevel() << std::endl;
    level = std::min(level, uncoarsener->currentLevel() + 2);
    while(!uncoarsener->isTopLevel()) {
      if (level > 0) {
        //#### Initial Partitioning ####
        timer.start_timer("initial_partitioning", "Initial Partitioning");
        DegreeZeroHypernodeRemover<TypeTraits> degree_zero_hn_remover(context);
        if (context.initial_partitioning.remove_degree_zero_hns_before_ip) {
          degree_zero_hn_remover.removeDegreeZeroHypernodes(partitioned_hg.hypergraph());
        }
        disableTimerAndStats(context);
        tbb::parallel_for(0, level * 12, [&](int i) {
          auto hg = partitioned_hg.hypergraph().copy();
          PartitionedHypergraph phg(context.partition.k, hg);
          ip(phg, level);
          std::cout << "Partition " << i << " done" << std::endl;
        });
        degree_zero_hn_remover.restoreDegreeZeroHypernodes(partitioned_hg);
        enableTimerAndStats(context);
        timer.stop_timer("initial_partitioning");
        // ####
        std::cout << "level " << level << std::endl;
        std::cout << "Number of partitions: " << partition_pool.size() << std::endl;

        // Project the first partition in the partition pool to the next level and refine
        projAndRefine(partition_pool[0]);

        // Project the rest of the partitions in the partition pool to the next level
        for(auto it = std::next(partition_pool.begin()); it != partition_pool.end(); it++) {
          vec<PartitionID> tmp_partition;
          tmp_partition.resize(partitioned_hg.initialNumNodes());
          std::copy(it->partIDs.begin(), it->partIDs.end(), tmp_partition.begin());

          it->partIDs.resize(partitioned_hg.initialNumNodes());
          partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {
            const HypernodeID coarse_hn = (uncoarseningData.hierarchy)[uncoarsener->currentLevel() + 1].mapToContractedHypergraph(hn);
            const PartitionID part_id = tmp_partition[coarse_hn];
            it->partIDs[hn] = part_id;
          });
        }
        // Refine all partitions
        for (auto it = std::next(partition_pool.begin()); it != partition_pool.end(); it++) {
          refine(*it);
        }

        // Do local tournament
        if(level > 1) {
          if(partition_pool.size() <= 2) {
            --level;
            continue;
          }
          // Sort the array of partitioned hypergraphs by quality
          std::sort(partition_pool.begin(), partition_pool.end(), isBetterThan);
          // first step: mark the winners of each match
          // mark the first 1 or 2 elements as winners depending on whether 
          // the number of partitions is even or odd
          std::vector<bool> winners(partition_pool.size(), true);
          int offset = partition_pool.size() % 2 ? 1 : 2;
          utils::Randomize::instance().parallelShuffleVector(partition_pool, offset, partition_pool.size());
          for(size_t i = offset; i < partition_pool.size(); i += 2) {
            if(isBetterThan(partition_pool[i], partition_pool[i + 1])) {
              winners[i + 1] = false;
            } else {
              winners[i] = false;
            }
          }
          //erase the losers
          partition_pool.erase(std::remove_if(partition_pool.begin(), partition_pool.end(),
                                              [&](const Partition& p) {
                                                return !winners[&p - &partition_pool[0]];
                                              }),
                              partition_pool.end());
        } else {
          int round = 0;
          // Do final tournament and choose winner and set the partitioned hypergraph to the winner
          while(partition_pool.size() > 1 && !uncoarsener->isTopLevel()) {
            std::cout << "Round " << round++ << std::endl;
            std::cout << "Number of partitions: " << partition_pool.size() << std::endl;
            // Print qualities of the partitions
            for(auto it = partition_pool.begin(); it != partition_pool.end(); it++) {
              std::cout << "Quality: " << it->quality << std::endl;
            }
            int offset = partition_pool.size() % 2 ? 1 : 0;
            utils::Randomize::instance().parallelShuffleVector(
                partition_pool, offset, partition_pool.size());

            std::vector<bool> winners(partition_pool.size(), true);
            for (size_t i = offset; i < partition_pool.size(); i += 2) {
              if (isBetterThan(partition_pool[i], partition_pool[i + 1])) {
                winners[i + 1] = false;
              } else {
                winners[i] = false;;
              }
            }
            // erase the losers
            partition_pool.erase(
                std::remove_if(partition_pool.begin(), partition_pool.end(),
                               [&](const Partition& p) {
                                 return !winners[&p - &partition_pool[0]];
                               }),
                partition_pool.end());
            // Project the partitions in the partition pool to the next level and refine
            projAndRefine(partition_pool[0]);
            for(auto it = std::next(partition_pool.begin()); it != partition_pool.end(); it++) {
              vec<PartitionID> tmp_partition;
              tmp_partition.resize(partitioned_hg.initialNumNodes());
              std::copy(it->partIDs.begin(), it->partIDs.end(), tmp_partition.begin());

              it->partIDs.resize(partitioned_hg.initialNumNodes());
              partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {
                const HypernodeID coarse_hn = (uncoarseningData.hierarchy)[uncoarsener->currentLevel() + 1].mapToContractedHypergraph(hn);
                const PartitionID part_id = tmp_partition[coarse_hn];
                it->partIDs[hn] = part_id;
              });
            }
            for(auto it = partition_pool.begin(); it != partition_pool.end(); it++) {
              refine(*it);
            }
          }
          std::sort(partition_pool.begin(), partition_pool.end(), isBetterThan);
          replacePartition(partition_pool[0]);
          std::cout << "Final partition quality: " << partition_pool[0].quality << std::endl;
        }
        --level;
      } else {
        timer.start_timer("refinement", "Refinement");
        uncoarsener->projectToNextLevelAndRefine();
        timer.stop_timer("refinement");
      }
    }
    if(partition_pool.size() == 0 || level == 1) {
      std::cout << "Partition pool is empty or level == 1" << std::endl;
      timer.start_timer("initial_partitioning", "Initial Partitioning");
      DegreeZeroHypernodeRemover<TypeTraits> degree_zero_hn_remover(context);
      if (context.initial_partitioning.remove_degree_zero_hns_before_ip) {
        degree_zero_hn_remover.removeDegreeZeroHypernodes(
            partitioned_hg.hypergraph());
      }
      disableTimerAndStats(context);
      tbb::parallel_for(0, 12, [&](int i) {
        auto hg = partitioned_hg.hypergraph().copy();
        PartitionedHypergraph phg(context.partition.k, hg);
        ip(phg, level);
        std::cout << "Partition " << i << " done" << std::endl;
      });
      degree_zero_hn_remover.restoreDegreeZeroHypernodes(partitioned_hg);
      enableTimerAndStats(context);
      timer.stop_timer("initial_partitioning");
    }
    std::sort(partition_pool.begin(), partition_pool.end(), isBetterThan);
    replacePartition(partition_pool[0]);
    uncoarsener->rebalancing();
    io::printPartitioningResults(partitioned_hg, context, "Local Search Results:");
    return uncoarsener->movePartitionedHypergraph();
  }
} // namespace

template<typename TypeTraits>
typename Multicandidate<TypeTraits>::PartitionedHypergraph
Multicandidate<TypeTraits>::partition(Hypergraph& hypergraph,
                                      const Context& context,
                                      const TargetGraph* target_graph) {
  PartitionedHypergraph partitioned_hg = 
    multicandidate_partitioning<TypeTraits>(hypergraph, context, target_graph);
  return partitioned_hg;
}

template<typename TypeTraits>
void Multicandidate<TypeTraits>::partition(
    PartitionedHypergraph& partitioned_hg, const Context& context,
    const TargetGraph* target_graph) {
  PartitionedHypergraph tmp_phg =
      partition(partitioned_hg.hypergraph(), context, target_graph);
}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(Multicandidate)

} // namespace mt_kahypar
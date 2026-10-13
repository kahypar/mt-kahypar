/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include "mt-kahypar/partition/recursive_bipartitioning.h"

#include <tbb/task_group.h>

#include <algorithm>
#include <vector>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/partition/bipartitioning_common.h"
#include "mt-kahypar/partition/multilevel.h"
#include "mt-kahypar/datastructures/fixed_vertex_support.h"
#include "mt-kahypar/partition/refinement/gains/bipartitioning_policy.h"
#ifdef KAHYPAR_ENABLE_STEINER_TREE_METRIC
#include "mt-kahypar/partition/mapping/initial_mapping.h"
#endif
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/utilities.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/weight/hypernode_weight_common.h"

#include "mt-kahypar/partition/metrics.h"

namespace mt_kahypar {


using OriginalHypergraphInfo = bipartitioning::OriginalHypergraphInfo;

namespace rb {

  static constexpr bool debug = false;

  // Sets the appropriate parameters for the recursive bipartitioning call
  Context setupRecursiveBipartitioningContext(const Context& context,
                                              const PartitionID k0, const PartitionID k1,
                                              const double degree_of_parallelism) {
    ASSERT((k1 - k0) >= 2);
    Context rb_context(context);
    rb_context.partition.k = k1 - k0;
    if (context.partition.mode == Mode::direct) {
      rb_context.type = ContextType::initial_partitioning;
    }

    rb_context.partition.perfect_balance_part_weights.changeSizeSequential(rb_context.partition.k);
    rb_context.partition.max_part_weights.changeSizeSequential(rb_context.partition.k);
    for ( PartitionID part_id = k0; part_id < k1; ++part_id ) {
      rb_context.partition.perfect_balance_part_weights[part_id - k0] =
              context.partition.perfect_balance_part_weights[part_id];
      rb_context.partition.max_part_weights[part_id - k0] =
              context.partition.max_part_weights[part_id];
    }

    rb_context.shared_memory.degree_of_parallelism *= degree_of_parallelism;

    return rb_context;
  }

  template<typename Hypergraph>
  void setupFixedVerticesForBipartitioning(Hypergraph& hg,
                                           const PartitionID k) {
    if ( hg.hasFixedVertices() ) {
      const PartitionID m = k / 2 + (k % 2);
      ds::FixedVertexSupport fixed_vertices(hg.initialNumNodes(), hg.dimension(), 2);
      hg.doParallelForAllNodes([&](const HypernodeID& hn) {
        if ( hg.isFixed(hn) ) {
          if ( hg.fixedVertexBlock(hn) < m ) {
            fixed_vertices.fixToBlock(hg, hn, 0);
          } else {
            fixed_vertices.fixToBlock(hg, hn, 1);
          }
        }
      });
      hg.addFixedVertexSupport(std::move(fixed_vertices));
    }
  }

  template<typename Hypergraph>
  void setupFixedVerticesForRecursion(const Hypergraph& input_hg,
                                      Hypergraph& extracted_hg,
                                      const vec<HypernodeID>& input2extracted,
                                      const PartitionID k0,
                                      const PartitionID k1) {
    if ( input_hg.hasFixedVertices() ) {
      ds::FixedVertexSupport fixed_vertices(
        extracted_hg.initialNumNodes(), extracted_hg.dimension(), k1 - k0);
      input_hg.doParallelForAllNodes([&](const HypernodeID& hn) {
        if ( input_hg.isFixed(hn) ) {
          const PartitionID block = input_hg.fixedVertexBlock(hn);
          if ( block >= k0 && block < k1 ) {
            fixed_vertices.fixToBlock(extracted_hg, input2extracted[hn], block - k0);
          }
        }
      });
      extracted_hg.addFixedVertexSupport(std::move(fixed_vertices));
    }
  }

  bool usesAdaptiveWeightOfNonCutEdges(const Context& context) {
    return BipartitioningPolicy::nonCutEdgeMultiplier(context.partition.gain_policy) != 1;
  }

  template<typename Hypergraph>
  void adaptWeightsOfNonCutEdges(Hypergraph& hg,
                                 const vec<uint8_t>& already_cut,
                                 const GainPolicy gain_policy,
                                 const bool revert) {
    const HyperedgeWeight multiplier = BipartitioningPolicy::nonCutEdgeMultiplier(gain_policy);
    if ( multiplier != 1 ) {
      ASSERT(static_cast<size_t>(hg.initialNumEdges()) == already_cut.size());
      hg.doParallelForAllEdges([&](const HyperedgeID& he) {
        if ( !already_cut[he] ) {
          hg.setEdgeWeight(he, static_cast<HyperedgeWeight>(( revert ? 1.0 / multiplier :
            static_cast<double>(multiplier) ) * hg.edgeWeight(he)));
        }
      });
    }
  }

  // Takes a hypergraph partitioned into two blocks as input and then recursively
  // partitions one block into (k1 - b0) blocks
  template<typename TypeTraits>
  void recursively_bipartition_block(typename TypeTraits::PartitionedHypergraph& phg,
                                     const Context& context,
                                     const PartitionID block, const PartitionID k0, const PartitionID k1,
                                     const OriginalHypergraphInfo& info,
                                     const vec<uint8_t>& already_cut,
                                     const double degree_of_parallism);

  // Uses multilevel recursive bipartitioning to partition the given hypergraph into (k1 - k0) blocks
  template<typename TypeTraits>
  void recursive_bipartitioning(typename TypeTraits::PartitionedHypergraph& phg,
                                const Context& context,
                                const PartitionID k0, const PartitionID k1,
                                const OriginalHypergraphInfo& info,
                                vec<uint8_t>& already_cut) {
    using Hypergraph = typename TypeTraits::Hypergraph;
    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
    if ( phg.initialNumNodes() > 0 ) {
      // Multilevel Bipartitioning
      const PartitionID k = (k1 - k0);
      Hypergraph& hg = phg.hypergraph();
      ds::FixedVertexSupport fixed_vertices = hg.copyOfFixedVertexSupport();
      Context b_context(context);
      bipartitioning::setupBipartitioningContext(
        b_context, context, info, 0, context.partition.k, hg.totalWeight(), Hypergraph::is_graph);
      setupFixedVerticesForBipartitioning(hg, k);
      adaptWeightsOfNonCutEdges(hg, already_cut, context.partition.gain_policy, false);

      DBG << "Multilevel Bipartitioning - Range = (" << k0 << "," << k1 << "), Epsilon =" << b_context.partition.epsilon;
      PartitionedHypergraph bipartitioned_hg = Multilevel<TypeTraits>::partition(hg, b_context);
      DBG << "Bipartitioning Result -"
          << "Objective =" << metrics::quality(bipartitioned_hg, b_context)
          << "Imbalance =" << metrics::imbalance(bipartitioned_hg, b_context)
          << "(Target Imbalance =" << b_context.partition.epsilon << ")";
      adaptWeightsOfNonCutEdges(hg, already_cut, context.partition.gain_policy, true);
      hg.addFixedVertexSupport(std::move(fixed_vertices));

      // Apply bipartition to the input hypergraph
      const PartitionID block_0 = 0;
      const PartitionID block_1 = k / 2 + (k % 2);
      phg.doParallelForAllNodes([&](const HypernodeID& hn) {
        PartitionID part_id = bipartitioned_hg.partID(hn);
        ASSERT(part_id != kInvalidPartition && part_id < phg.k());
        ASSERT(phg.partID(hn) == kInvalidPartition);
        if ( part_id == 0 ) {
          phg.setOnlyNodePart(hn, block_0);
        } else {
          phg.setOnlyNodePart(hn, block_1);
        }
      });
      phg.initializePartition();

      if ( usesAdaptiveWeightOfNonCutEdges(context) ) {
        // Update cut hyperedges
        phg.doParallelForAllEdges([&](const HyperedgeID& he) {
          already_cut[he] |= phg.connectivity(he) > 1;
        });
      }

      ASSERT(metrics::quality(bipartitioned_hg, context) ==
            metrics::quality(phg, context));

      ASSERT(context.partition.k >= 2);
      PartitionID rb_k0 = context.partition.k / 2 + context.partition.k % 2;
      PartitionID rb_k1 = context.partition.k / 2;
      if ( rb_k0 >= 2 && rb_k1 >= 2 ) {
        // Both blocks of the bipartition must to be further partitioned into at least two blocks.
        DBG << "Current k = " << context.partition.k << "\n"
            << "Block" << block_0 << "is further partitioned into k =" << rb_k0 << "blocks\n"
            << "Block" << block_1 << "is further partitioned into k =" << rb_k1 << "blocks\n";
        tbb::task_group tg;
        tg.run([&] { recursively_bipartition_block<TypeTraits>(phg, context, block_0, 0, rb_k0, info, already_cut, 0.5); });
        tg.run([&] { recursively_bipartition_block<TypeTraits>(phg, context, block_1, rb_k0, rb_k0 + rb_k1, info, already_cut, 0.5); });
        tg.wait();
      } else if ( rb_k0 >= 2 ) {
        ASSERT(rb_k1 < 2);
        // Only the first block needs to be further partitioned into at least two blocks.
        DBG << "Current k = " << context.partition.k << "\n"
            << "Block" << block_0 << "is further partitioned into k =" << rb_k0 << "blocks\n";
        recursively_bipartition_block<TypeTraits>(phg, context, block_0, 0, rb_k0, info, already_cut, 1.0);
      }
    }
  }
}

template<typename TypeTraits>
void rb::recursively_bipartition_block(typename TypeTraits::PartitionedHypergraph& phg,
                                        const Context& context,
                                        const PartitionID block, const PartitionID k0, const PartitionID k1,
                                        const OriginalHypergraphInfo& info,
                                        const vec<uint8_t>& already_cut,
                                        const double degree_of_parallism) {
  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  Context rb_context = setupRecursiveBipartitioningContext(context, k0, k1, degree_of_parallism);
  // Extracts the block of the hypergraph which we recursively want to partition
  const bool cut_net_splitting =
    BipartitioningPolicy::useCutNetSplitting(context.partition.gain_policy);
  auto extracted_block = phg.extract(block, !already_cut.empty() ? &already_cut : nullptr,
    cut_net_splitting, context.preprocessing.stable_construction_of_incident_edges);
  Hypergraph& rb_hg = extracted_block.hg;
  auto& mapping = extracted_block.hn_mapping;
  setupFixedVerticesForRecursion(phg.hypergraph(), rb_hg, mapping, k0, k1);

  if ( rb_hg.initialNumNodes() > 0 ) {
    // Recursively partition the given block into (k1 - k0) blocks
    PartitionedHypergraph rb_phg(rb_context.partition.k, rb_hg, parallel_tag_t());
    recursive_bipartitioning<TypeTraits>(rb_phg, rb_context,
      k0, k1, info, extracted_block.already_cut);

    ASSERT(phg.initialNumNodes() == mapping.size());
    // Apply k-way partition to the input hypergraph
    phg.doParallelForAllNodes([&](const HypernodeID& hn) {
      if ( phg.partID(hn) == block ) {
        ASSERT(hn < mapping.size());
        PartitionID to = block + rb_phg.partID(mapping[hn]);
        ASSERT(to != kInvalidPartition && to < phg.k());
        if ( block != to ) {
          phg.changeNodePartNoSync(hn, block, to, true);
        }
      }
    });
    phg.resetEdgeSynchronization();
    DBG << "Recursive Bipartitioning Result -"
        << "k =" << (k1 - k0)
        << "Objective =" << metrics::quality(phg, context)
        << "Imbalance =" << metrics::imbalance(phg, rb_context)
        << "(Target Imbalance =" << rb_context.partition.epsilon << ")";

  }
}

template<typename TypeTraits>
typename RecursiveBipartitioning<TypeTraits>::PartitionedHypergraph
RecursiveBipartitioning<TypeTraits>::partition(Hypergraph& hypergraph,
                                               const Context& context,
                                               const TargetGraph* target_graph) {
  PartitionedHypergraph partitioned_hypergraph(context.partition.k, hypergraph, parallel_tag_t());
  partition(partitioned_hypergraph, context, target_graph);
  return partitioned_hypergraph;
}

template<typename TypeTraits>
void RecursiveBipartitioning<TypeTraits>::partition(PartitionedHypergraph& hypergraph,
                                                    const Context& context,
                                                    const TargetGraph* target_graph) {
  unused(target_graph);
  utils::Utilities& utils = utils::Utilities::instance();
  if (context.partition.mode == Mode::recursive_bipartitioning) {
    utils.getTimer(context.utility_id).start_timer("rb", "Recursive Bipartitioning");
  }

  if (context.type == ContextType::main) {
    parallel::MemoryPool::instance().deactivate_unused_memory_allocations();
    utils.getTimer(context.utility_id).disable();
    utils.getStats(context.utility_id).disable();
  }

  Context rb_context(context);
  if ( rb_context.partition.objective == Objective::steiner_tree ) {
    // In RB mode, we optimize the km1 metric for the steiner tree metric and
    // apply the permutation computed in the target graph to the partition.
    rb_context.partition.objective = PartitionedHypergraph::is_graph ?
      Objective::cut : Objective::km1;
    rb_context.partition.gain_policy = PartitionedHypergraph::is_graph ?
      GainPolicy::cut_for_graphs : GainPolicy::km1;
  }
  if ( context.type == ContextType::initial_partitioning ) {
    rb_context.partition.enable_logging = false;
  }

  vec<uint8_t> already_cut(rb::usesAdaptiveWeightOfNonCutEdges(context) ?
    hypergraph.initialNumEdges() : 0, 0);
  rb::recursive_bipartitioning<TypeTraits>(hypergraph, rb_context, 0, rb_context.partition.k,
    OriginalHypergraphInfo { weight::sum(hypergraph.totalWeight()), rb_context.partition.k,
      rb_context.partition.epsilon }, already_cut);

  if (context.type == ContextType::main) {
    parallel::MemoryPool::instance().activate_unused_memory_allocations();
    utils.getTimer(context.utility_id).enable();
    utils.getStats(context.utility_id).enable();
  }

  #ifdef KAHYPAR_ENABLE_STEINER_TREE_METRIC
  if ( context.partition.objective == Objective::steiner_tree ) {
    ASSERT(target_graph);
    utils::Timer& timer = utils.getTimer(context.utility_id);
    const bool was_enabled = timer.isEnabled();
    timer.enable();
    timer.start_timer("one_to_one_mapping", "One-To-One Mapping");
    // Map partition onto target graph
    InitialMapping<TypeTraits>::mapToTargetGraph(
      hypergraph, *target_graph, context);
    timer.stop_timer("one_to_one_mapping");
    if ( !was_enabled ) {
      timer.disable();
    }
  }
  #endif

  if (context.partition.mode == Mode::recursive_bipartitioning) {
    utils.getTimer(context.utility_id).stop_timer("rb");
  }
}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(RecursiveBipartitioning)

} // namepace mt_kahypar
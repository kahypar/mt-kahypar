/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@kit.edu>
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#pragma once

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/coarsening/i_uncoarsener.h"
#include "mt-kahypar/partition/coarsening/uncoarsener_base.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/coarsening/coarsening_commons.h"
namespace mt_kahypar {

  class NLevelUncoarsener : public IUncoarsener,
                            private UncoarsenerBase {

  private:
    static constexpr bool enable_heavy_assert = false;

    using ParallelHyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<ParallelHyperedge>>;

  public:
    NLevelUncoarsener(Hypergraph& hypergraph,
                      const Context& context,
                      UncoarseningData& uncoarseningData) :
      UncoarsenerBase(hypergraph, context, uncoarseningData) {}

    NLevelUncoarsener(const NLevelUncoarsener&) = delete;
    NLevelUncoarsener(NLevelUncoarsener&&) = delete;
    NLevelUncoarsener & operator= (const NLevelUncoarsener &) = delete;
    NLevelUncoarsener & operator= (NLevelUncoarsener &&) = delete;

  private:
    PartitionedHypergraph&& doUncoarsen(std::unique_ptr<IRefiner>& label_propagation,
                                        std::unique_ptr<IRefiner>& fm);

    void localizedRefine(PartitionedHypergraph& partitioned_hypergraph,
                         const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                         std::unique_ptr<IRefiner>& label_propagation,
                         std::unique_ptr<IRefiner>& fm,
                         Metrics& current_metrics,
                         const bool force_measure_timings);

    void globalRefine(PartitionedHypergraph& partitioned_hypergraph,
                      std::unique_ptr<IRefiner>& fm,
                      std::unique_ptr<IRefiner>& flows,
                      Metrics& current_metrics,
                      const double time_limit);

  PartitionedHypergraph&& uncoarsenImpl(
      std::unique_ptr<IRefiner>& label_propagation,
      std::unique_ptr<IRefiner>& fm) override {
    initHierarchy();
    return doUncoarsen(label_propagation, fm);
  }
  void initHierarchy() {
    // Create n-level batch uncontraction hierarchy
    _timer.start_timer("create_batch_uncontraction_hierarchy", "Create n-Level Hierarchy");
    _hierarchy = _hg.createBatchUncontractionHierarchy(_context.refinement.max_batch_size);
    ASSERT(_uncoarseningData.removed_hyperedges_batches.size() == _hierarchy.size() - 1);
    _timer.stop_timer("create_batch_uncontraction_hierarchy");
  }

  // ! Represents the n-level hierarchy
  // ! A batch is vector of uncontractions/mementos that can be uncontracted in parallel
  // ! without conflicts. All batches of a specific version of the hypergraph are assembled
  // ! in a batch vector. Each time we perform single-pin and parallel net detection we create
  // ! a new version (simply increment a counter) of the hypergraph. Once a batch vector is
  // ! completly processed single-pin and parallel nets have to be restored.
  VersionedBatchVector _hierarchy;
  };
}

/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
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

#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/partition/refinement/spectral/algebraic_wrappers/matrix.h"
#include "mt-kahypar/partition/refinement/spectral/algebraic_wrappers/vector.h"


namespace mt_kahypar {
template <typename GraphAndGainTypes>
class SpectralRefiner final : public IRefiner {
 private:
  using Hypergraph = typename GraphAndGainTypes::Hypergraph;
  using PartitionedHypergraph = typename GraphAndGainTypes::PartitionedHypergraph;
  using GainCache = typename GraphAndGainTypes::GainCache;

  static constexpr bool debug = true;
  static constexpr bool enable_heavy_assert = true;

 public:
  explicit SpectralRefiner(const HypernodeID num_hypernodes,
                           const HyperedgeID num_hyperedges,
                           const Context& context,
                           GainCache& gain_cache) :
    _context(context),
    _gain_cache(gain_cache) {
      unused(num_hypernodes);
      unused(num_hyperedges);
    }

  explicit SpectralRefiner(const HypernodeID num_hypernodes,
                                   const HyperedgeID num_hyperedges,
                                   const Context& context,
                                   gain_cache_t gain_cache) :
    SpectralRefiner(num_hypernodes, num_hyperedges, context,
      GainCachePtr::cast<GainCache>(gain_cache)) { }

  SpectralRefiner(const SpectralRefiner&) = delete;
  SpectralRefiner(SpectralRefiner&&) = delete;

  SpectralRefiner & operator= (const SpectralRefiner &) = delete;
  SpectralRefiner & operator= (SpectralRefiner &&) = delete;

 private:
  bool refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                  const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                  Metrics& best_metrics,
                  double) final;

  void initializeImpl(mt_kahypar_partitioned_hypergraph_t&) final;

  void kSpecPartAlgorithm(PartitionedHypergraph& hypergraph);

  void dehyperizeToLaplacian(Hypergraph& hypergraph, spectral::Matrix& target);

  void buildWeightBalanceGraphLaplacian(Hypergraph& hypergraph, spectral::Matrix& target);

  void generate2WayVertexEmbedding(spectral::Matrix& baseBalance, spectral::Matrix& graphLaplacian, PartitionedHypergraph& hintSolution, vec<spectral::Vector>& target);

  void generateHintGraphLaplacian(PartitionedHypergraph& hintSolution, spectral::Matrix& target);


  void resizeDataStructuresForCurrentK() {
    /* TODO to be implemented far in the future */
  }

  const Context& _context;
  GainCache& _gain_cache;
};

}  // namespace kahypar

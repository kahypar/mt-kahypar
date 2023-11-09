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

#include "mt-kahypar/partition/refinement/spectral/spectral_refiner.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"

namespace mt_kahypar {

  template <typename GraphAndGainTypes>
  bool SpectralRefiner<GraphAndGainTypes>::refineImpl(mt_kahypar_partitioned_hypergraph_t& phg,
                                                      const vec<HypernodeID>& refinement_nodes,
                                                      Metrics& best_metrics,
                                                      const double)  {
    ASSERT(refinement_nodes.empty());  // these are always empty for your case
    unused(refinement_nodes);

    PartitionedHypergraph& hypergraph = utils::cast<PartitionedHypergraph>(phg);
    Gain old_quality = best_metrics.quality;
    // resizeDataStructuresForCurrentK();

    // disable the gain cache for spectral partitioning
    _gain_cache.reset();


    // implementation goes here
    DBG << "Spectral Refiner called";
    kSpecPartAlgorithm(hypergraph);
    if constexpr (Hypergraph::is_graph) {

    }

    // recalculate metrics


    HEAVY_REFINEMENT_ASSERT(hypergraph.checkTrackedPartitionInformation(_gain_cache));
    HEAVY_REFINEMENT_ASSERT(best_metrics.quality ==
      metrics::quality(hypergraph, _context,
        !_context.refinement.label_propagation.execute_sequential),
      V(best_metrics.quality) << V(metrics::quality(hypergraph, _context,
          !_context.refinement.label_propagation.execute_sequential)));

    // Update metrics statistics
    Gain delta = old_quality - best_metrics.quality;
    ASSERT(delta >= 0, "Refiner worsened solution quality");
    utils::Utilities::instance().getStats(_context.utility_id).update_stat("spectral_improvement", delta);
    return delta > 0;
  }


  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::kSpecPartAlgorithm(PartitionedHypergraph& inputSolution) {
    Hypergraph& inputHypergraph /* = TODO */;

    // dehyperisation
    SimpleGraph inputGraph; /* TODO allocation?? */
    dehyperize(inputHypergraph, inputGraph);
    /*auto inputLaplacian /* = TODO(inputGraph) */; /* TODO even possible as general Laplacian? see (8) in paper */

    // weight-balance graph construction
    SimpleGraph weightBalanceGraph; /* TODO allocation?? */
    buildWeightBalanceGraph(inputHypergraph, weightBalanceGraph);
    /*auto baseBalance /* = TODO(weightBalanceGraph) */; /* TODO same as l74 with (11) in paper */

    // actual refinement
    /* vec<PartitionedHypergraph>*///auto candidateSolutions;
    /* candidateSolutions->reserve("beta"); */
    /* candidateSolutions->pushBack(inputSolution); */
    for (int i = 0; i < 1 /* TODO parameter "beta" */; i++) {
      if (/* "k" == 2 */true) {
        /* TODO returntype/result, header file entry, implementation *///solveGEVP(inputLaplacian, baseBalance, candidateSolutions->back());
      } else {
        /* TODO */
      }
    }
  }

  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::dehyperize(Hypergraph& hypergraph, SimpleGraph& target) {
    // dehyperisation via clique expansion graph

    /* TODO maybe directly/only compute laplacian */
  }
  
  
  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::buildWeightBalanceGraph(Hypergraph& hypergraph, SimpleGraph& target) {
    /* TODO maybe directly/only compute laplacian */
  }


  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::initializeImpl(mt_kahypar_partitioned_hypergraph_t& phg) {
    unused(phg);
  }



  namespace {
  #define SPECTRAL_REFINER(X) SpectralRefiner<X>
  }

  // explicitly instantiate so the compiler can generate them when compiling this cpp file
  INSTANTIATE_CLASS_WITH_VALID_TRAITS(SPECTRAL_REFINER)
}

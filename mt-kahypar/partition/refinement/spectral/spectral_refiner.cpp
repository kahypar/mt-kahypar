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
// #include "mt-kahypar/partition/refinement/spectral/solvers/gevp.h"
#include "mt-kahypar/partition/refinement/spectral/solvers/slepc_gevp.cpp" /* TODO only makeshift */

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"

namespace mt_kahypar {

  using namespace spectral;

  template <typename GraphAndGainTypes>
  bool SpectralRefiner<GraphAndGainTypes>::refineImpl(mt_kahypar_partitioned_hypergraph_t& phg,
                                                      const vec<HypernodeID>& refinement_nodes,
                                                      Metrics& best_metrics,
                                                      const double)  {
    ASSERT(refinement_nodes.empty());  // these are always empty for your case
    unused(refinement_nodes);

    PartitionedHypergraph& hypergraph = utils::cast<PartitionedHypergraph>(phg);
    Gain old_quality = best_metrics.quality;
    resizeDataStructuresForCurrentK();

    // disable the gain cache for spectral partitioning
    _gain_cache.reset();


    // implementation goes here
    DBG << "Spectral Refiner called";
    partition(hypergraph);
    DBG << "Spectral Refiner finished kSpecPart";
    if constexpr (Hypergraph::is_graph) {

    }

    // recalculate metrics
    /* TODO */


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
  void SpectralRefiner<GraphAndGainTypes>::partition(PartitionedHypergraph& inputSolution) {
    Hypergraph& inputHypergraph  = inputSolution.hypergraph();
    PartitionID k = inputSolution.k(); /* TODO extract from argv */
    
    numNodes = inputHypergraph.initialNumNodes();

    // dehyperisation
    Operator inputGraphLaplacian(numNodes);
    dehyperizeToLaplacian(inputHypergraph, inputGraphLaplacian);

    // weight-balance graph construction
    Operator weightBalanceLaplacian(numNodes);
    buildWeightBalanceGraphLaplacian(inputHypergraph, weightBalanceLaplacian);

    // actual refinement

    vec<PartitionedHypergraph*> candidateSolutions; 
    candidateSolutions.reserve(1/* TODO argv "beta" */ + 1);
    candidateSolutions.push_back(&inputSolution);

    for (int i = 0; i < 1 /* TODO argv "beta" */; i++) {
      vec<spectral::Vector> embedding; /* TODO type alias */
      if (k == 2) {
        generate2WayVertexEmbedding(weightBalanceLaplacian, inputGraphLaplacian, *candidateSolutions.back(), embedding);
      } else {
        /* TODO */
      }
    }
  }

  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::dehyperizeToLaplacian(Hypergraph& hypergraph, spectral::Operator& target) {
    // dehyperisation via clique expansion graph

    target.ctx = (void *) &hypergraph;

    target.effect = [](Operator *self, Vector& operand, Vector& target_vector) {
      size_t n = operand.dimension();
      Hypergraph *hg = (Hypergraph *) self->ctx;

      spectral::Skalar clique_edge_weight_sum = 0.0;
      spectral::Vector edges_superposition(n);
      for (const HyperedgeID& he : hg->edges()) { /* TODO use parallel?? */
        HypernodeID edge_size = hg->edgeSize(he);

        clique_edge_weight_sum += 1.0 / (edge_size - 1.0);

        spectral::Skalar operand_norm_he = 0.0;
        for (const HypernodeID& pin : hg->pins(he)) {
          operand_norm_he += operand.get(pin /* TODO calculate index */);
        }
        operand_norm_he /= edge_size;

        // accumulate edges superposition vector
        for (const HypernodeID& pin : hg->pins(he)) {
          size_t index = pin; /* TODO calculate index */
          edges_superposition.set(index, edges_superposition.get(index) + operand_norm_he);
        }
      }
      
      // calculate result
      for (size_t i = 0; i < target_vector.dimension(); i++) {
        target_vector.set(i, clique_edge_weight_sum * operand.get(i) - edges_superposition.get(i));
      }
    };

    target.calc_diagonal = [] (Operator *self, Vector& target_vector) {
      Hypergraph *hg = (Hypergraph *) self->ctx;
      for (const HypernodeID& node : hg->nodes()) {
        for (const HyperedgeID& edge : hg->incidentEdges(node)) {
          size_t index = node; /* TODO calculate index */
          target_vector.set(node, target_vector.get(node) + hg->edgeWeight(edge));
        }
      }
    };
  }
  
  
  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::buildWeightBalanceGraphLaplacian(Hypergraph& hypergraph, spectral::Operator& target) {
    target.ctx = (void *) &hypergraph;

    target.effect = [](Operator *self, Vector& operand, Vector& target_vector) {
      Hypergraph *hg = (Hypergraph *) self->ctx;
      auto dimension = operand.dimension();

      spectral::Skalar sum_x = 0.0;
      for (size_t i = 0; i < dimension; i++) {
        sum_x += operand.get(i);
      }
      
      for (const HypernodeID& node : hg->nodes()) {
        size_t index = node; /* TODO calculate index */
        target_vector.set(index, hg->nodeWeight(node) * (operand.get(index) - sum_x / dimension));
      }
    };

    target.calc_diagonal = [] (Operator *self, Vector& target_vector) {
      Hypergraph *hg = (Hypergraph *) self->ctx;

      spectral::Skalar total_node_weight = 0.0; // TODO use Hypergraph::computeAndSetTotalNodeWeight(parallel_tag_t) ???
      for (const HypernodeID node : hg->nodes()) {
        total_node_weight += hg->nodeWeight(node);
      }

      for (const HypernodeID node : hg->nodes()) {
        size_t index = node; /* TODO calculate index */
        spectral::Skalar weight = hg->nodeWeight(node);
        target_vector.set(index, weight * (total_node_weight - weight));
      }
    };
  }


  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::generate2WayVertexEmbedding(spectral::Operator& baseBalance, spectral::Operator& graphLaplacian, PartitionedHypergraph& hintSolution, vec<spectral::Vector>& target) {
    // hint graph
    Operator hintGraphLaplacian(numNodes);
    generateHintGraphLaplacian(hintSolution, hintGraphLaplacian);
    
    spectral::SLEPcGEVPSolver solver; /* TODO get gevp variant otherwise */
    solver.setProblem(graphLaplacian, baseBalance /*+ hintGraphLaplacian*/);

    /* TODO */
    // to see something:
    spectral::Skalar a;
    spectral::Vector v(numNodes);
    // v.set(0, 5);
    // graphLaplacian.getDiagonal(v);
    //DBG << v.get(0);
    solver.nextEigenpair(a, v);
    v.get_all();
  }


  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::generateHintGraphLaplacian(PartitionedHypergraph& hintSolution, spectral::Operator& target) {
    target.ctx = (void *) &(hintSolution.hypergraph());
    /* TODO */
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

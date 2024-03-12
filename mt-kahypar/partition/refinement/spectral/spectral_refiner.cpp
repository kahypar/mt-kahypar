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

#include <iostream>

namespace mt_kahypar {

  using namespace spectral;

  template <typename GraphAndGainTypes>
  bool SpectralRefiner<GraphAndGainTypes>::refineImpl(mt_kahypar_partitioned_hypergraph_t& phg,
                                                      const vec<HypernodeID>& refinement_nodes,
                                                      Metrics& best_metrics,
                                                      const double)  {
    ASSERT(refinement_nodes.empty());  // these are always empty for your case
    unused(refinement_nodes);

    PartitionedHypergraph& partionedHypergraph = utils::cast<PartitionedHypergraph>(phg);
    Gain old_quality = best_metrics.quality;
    resizeDataStructuresForCurrentK();

    // disable the gain cache for spectral partitioning
    _gain_cache.reset();


    // implementation goes here
    DBG << "Spectral Refiner called";
    utils::Timer& timer = utils::Utilities::instance().getTimer(_context.utility_id);
    timer.start_timer("partition_sp", "Partition");
    bool found_new_partition = partition(partionedHypergraph, best_metrics);
    timer.stop_timer("partition_sp");
    if (found_new_partition) {
      DBG << "found new partitioning solution";
    }
    DBG << "Spectral Refiner finished partitioning";

    // recalculate metrics
    
    /* TODO */


    HEAVY_REFINEMENT_ASSERT(partionedHypergraph.checkTrackedPartitionInformation(_gain_cache));
    HEAVY_REFINEMENT_ASSERT(best_metrics.quality ==
      metrics::quality(partionedHypergraph, _context,
        !_context.refinement.label_propagation.execute_sequential),
      V(best_metrics.quality) << V(metrics::quality(partionedHypergraph, _context,
          !_context.refinement.label_propagation.execute_sequential)));

    // Update metrics statistics
    Gain delta = old_quality - best_metrics.quality;
    ASSERT(delta >= 0, "Refiner worsened solution quality");
    utils::Utilities::instance().getStats(_context.utility_id).update_stat("spectral_improvement", delta);
    return delta > 0;
  }


  template <typename GraphAndGainTypes>
  bool SpectralRefiner<GraphAndGainTypes>::partition(PartitionedHypergraph& phg, Metrics &best_metrics) {
    if constexpr (Hypergraph::is_graph) {
      /* TODO */
    }

    Hypergraph& inputHypergraph  = phg.hypergraph();
    const PartitionID k = phg.k(); /* TODO extract from argv */
    
    numNodes = inputHypergraph.initialNumNodes();

    vec<PartitionID> inputPartition;
    inputPartition.reserve(numNodes);
    for (const HypernodeID node : phg.nodes()) {
      inputPartition.push_back(phg.partID(node));
    }

    // dehyperisation
    Operator inputGraphLaplacian(numNodes);
    dehyperizeToLaplacian(inputHypergraph, inputGraphLaplacian);

    /* vec<spectral::Vector> inputGraphLaplacianMatrix;
    spectral::Vector diagonal(numNodes);
    inputGraphLaplacian.getDiagonal(diagonal);
    inputGraphLaplacian.getMatrix(inputGraphLaplacianMatrix);
    inputGraphLaplacianMatrix.push_back(diagonal);
    for(spectral::Vector row : inputGraphLaplacianMatrix) {
      std::ostringstream row_str;
      for (size_t i = 0; i < numNodes; i++) {
        char buf[100];
        sprintf(buf, " %+.2f", row[i]);
        row_str << buf;
      }
      DBG << row_str.str();
    } */

    // weight-balance graph construction
    Operator weightBalanceLaplacian(numNodes);
    buildWeightBalanceGraphLaplacian(inputHypergraph, weightBalanceLaplacian);

    // actual refinement

    vec<vec<PartitionID>> candidateSolutions; /* TODO alias */
    candidateSolutions.reserve(1/* TODO argv "beta" */ + 1);
    Gain best_cutsize;
    size_t best_index;

    for (int i = 0; i < 1 /* TODO argv "beta" */; i++) {
      vec<spectral::Vector> embedding; /* TODO type alias */
      if (k == 2) {
        generate2WayVertexEmbedding(weightBalanceLaplacian, inputGraphLaplacian, phg, embedding);
      } else {
        /* TODO */
      }

      vec<PartitionID> newSolution;
      generateSolution(phg, embedding, newSolution);
      candidateSolutions.push_back(newSolution);

      // calulate metrics
      Gain new_cutsize = metrics::quality(phg, _context, !_context.refinement.label_propagation.execute_sequential);
      if (candidateSolutions.size() == 1 || new_cutsize < best_cutsize) {
        best_cutsize = new_cutsize;
        best_index = candidateSolutions.size() - 1;
      }
      
    }

    bool found_valid_solution = best_cutsize <= best_metrics.quality;
    best_metrics.quality = found_valid_solution ? best_cutsize : best_metrics.quality;

    DBG << V(best_cutsize);

    if (found_valid_solution && best_index != candidateSolutions.size() - 1) {
      setPartition(phg, candidateSolutions[best_index]);
    } else {
      setPartition(phg, inputPartition);
    }

    return found_valid_solution;
  }

  template <typename GraphAndGainTypes>
  template <typename Collection>
  void SpectralRefiner<GraphAndGainTypes>::setPartition(PartitionedHypergraph &phg, Collection &partition) {
    phg.resetPartition(); /* TODO change changed nodes only? */
    for (size_t i = 0; i < numNodes; i++) {
      HypernodeID node = i; /* TODO calculate index */
      phg.setNodePart(node, partition[i]);
    }
  }

  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::dehyperizeToLaplacian(Hypergraph& hypergraph, spectral::Operator& target) {
    // dehyperisation via clique expansion graph

    target.ctx = (void *) &hypergraph;

    target.effects[0] = [](Operator *self, Vector& operand, Vector& target_vector) {
      size_t n = operand.dimension();
      Hypergraph *hg = (Hypergraph *) self->ctx;

      spectral::Vector factor(n);
      spectral::Vector subtrahend(n);
      for (const HyperedgeID& he : hg->edges()) {
        HypernodeID edge_size = hg->edgeSize(he);

        // isolated vertex
        if (edge_size == 1) {
          continue;
        }
        
        spectral::Skalar operand_dot_e = 0.0;
        for (const HypernodeID& pin : hg->pins(he)) {
          operand_dot_e += operand[pin];
        }

        spectral::Skalar clique_edge_weight = ((spectral::Skalar) hg->edgeWeight(he)) / (-1.0 + (spectral::Skalar) edge_size);

        for (const HypernodeID& pin : hg->pins(he)) {
          factor.set(pin, factor[pin] + clique_edge_weight * edge_size);
          subtrahend.set(pin, subtrahend[pin] + clique_edge_weight * operand_dot_e);
        }
      }

      // isolated vertices
      spectral::Skalar sum_operand = INT_MAX;
      spectral::Skalar subtrahend_iso_default = 0.0;
      spectral::Skalar zero = 0.0;
      vec<spectral::Skalar*> subtrahend_iso;
      subtrahend_iso.resize(n, &subtrahend_iso_default);
      size_t isolated_vertices = 0;
      for (const HypernodeID &v : hg->nodes()) {
        if (hg->nodeDegree(v) > 0) {
          continue;
        }

        isolated_vertices++;

        if (sum_operand == INT_MAX) {
            sum_operand = 0.0;
            for (size_t i = 0; i < n; i++) {
              sum_operand += operand[i];
            }
          }
          
          subtrahend_iso_default += operand[v];
          subtrahend_iso[v] = &zero;
          subtrahend.set(v, subtrahend[v] + sum_operand - ((spectral::Skalar) n) * operand[v]);
      }

      // calculate result
      target_vector.reset();
      for (size_t i = 0; i < target_vector.dimension(); i++) {
        target_vector.set(i, target_vector[i] + operand[i] * factor[i] - subtrahend[i] - *subtrahend_iso[i]);
      }
    };

    target.calc_diagonal_ops[0] = [] (Operator *self, Vector& target_vector) {
      Hypergraph *hg = (Hypergraph *) self->ctx;
      for (const HypernodeID& node : hg->nodes()) {
        if (hg->nodeDegree(node) == 0) {
          target_vector.set(node, self->dimension() - 1);
          continue;
        }
        for (const HyperedgeID& edge : hg->incidentEdges(node)) {
          target_vector.set(node, target_vector[node] + hg->edgeWeight(edge));
        }
      }
    };
  }
  
  
  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::buildWeightBalanceGraphLaplacian(Hypergraph& hypergraph, spectral::Operator& target) {
    target.ctx = (void *) &hypergraph;

    target.effects.push_back([](Operator *self, Vector& operand, Vector& target_vector) {
      Hypergraph *hg = (Hypergraph *) self->ctx;
      auto dimension = operand.dimension();

      spectral::Skalar sum_x = 0.0;
      for (size_t i = 0; i < dimension; i++) {
        sum_x += operand[i];
      }
      
      for (const HypernodeID& node : hg->nodes()) {
        size_t index = node; /* TODO calculate index */
        target_vector.set(index, target_vector[index] + hg->nodeWeight(node) * (operand[index] - sum_x / dimension));
      }
    });

    target.calc_diagonal_ops.push_back([] (Operator *self, Vector& target_vector) {
      Hypergraph *hg = (Hypergraph *) self->ctx;

      spectral::Skalar total_node_weight = 0.0; // TODO use Hypergraph::computeAndSetTotalNodeWeight(parallel_tag_t) ???
      for (const HypernodeID node : hg->nodes()) {
        total_node_weight += hg->nodeWeight(node);
      }

      for (const HypernodeID node : hg->nodes()) {
        size_t index = node; /* TODO calculate index */
        spectral::Skalar weight = hg->nodeWeight(node);
        target_vector.set(index, target_vector[index] + weight * (total_node_weight - weight));
      }
    });
  }


  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::generate2WayVertexEmbedding(spectral::Operator& baseBalance, spectral::Operator& graphLaplacian, PartitionedHypergraph& hintSolution, vec<spectral::Vector>& target) {
    // hint graph
    Operator hintGraphLaplacian(numNodes);
    generateHintGraphLaplacian(hintSolution, hintGraphLaplacian);
    
    spectral::SLEPcGEVPSolver solver; /* TODO get gevp variant otherwise */
    spectral::Vector one(numNodes, 1.0);
    Skalar zero = 0.0;
    spectral::Operator dummy(numNodes);
    solver.setProblem(graphLaplacian, dummy, one, zero);//baseBalance /*+ hintGraphLaplacian*/);

    spectral::Skalar a;
    spectral::Vector v(numNodes);

    size_t num_evecs = 2; /* only fiedler */
    target.push_back(one);
    while (target.size() < num_evecs) {
      solver.nextEigenpair(a, v); /* TODO assert return value is positive */
      target.push_back(v);
    }  
  }

  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::generateSolution(PartitionedHypergraph &phg, vec<spectral::Vector> &embedding, vec<PartitionID> &target) { /* TODO aliase */
    // definitions
    spectral::Skalar max_part_weight = round(0.5 * (_context.partition.epsilon + 1) * phg.totalWeight());
    auto objective_delta = [&](const SynchronizedEdgeUpdate& sync_update) {
      _gain.computeDeltaForHyperedge(sync_update);
    };    
    spectral::Skalar split_index = 0;
    spectral::Vector &fiedler = embedding[1];
    vec<HypernodeID> indices(phg.nodes().begin(), phg.nodes().end());
    auto move_node = [&](bool use_gain, bool dest) {
      target[indices[split_index]] = dest ? 1 : 0;
      if (use_gain) {
        phg.changeNodePart(indices[split_index], dest ? 0 : 1, dest ? 1 : 0, objective_delta);
      } else {
        phg.changeNodePart(indices[split_index], dest ? 0 : 1, dest ? 1 : 0);
      }
    };

    // sort by fiedler entry
    std::sort(indices.begin(), indices.end(), [&](const HypernodeID& u, const HypernodeID& v) { return fiedler[u] < fiedler[v]; });

    // initialize with all nodes belonging to partition 1
    target.resize(numNodes, 1);
    setPartition(phg, target);

    // move into range
    for (; phg.partWeight(1) > max_part_weight; move_node(false, false)) {
      split_index++;
    }
    
    // find optimum in range
    _gain.reset();
    size_t best_index = split_index;
    Gain best_delta = _gain.localDelta();
    while (phg.partWeight(0) + phg.nodeWeight(indices[split_index + 1]) <= max_part_weight) {
      split_index++;
      move_node(true, false);

      if (_gain.localDelta() <= best_delta) {
        best_index = split_index;
        best_delta =_gain.localDelta();
      }  
    }
    
    // undo bad moves
    for (; best_index < split_index; split_index--) {
      move_node(false, true);
    }
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

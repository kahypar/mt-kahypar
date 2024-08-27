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

/* TODO thats not how cpp includes work */
// #include "mt-kahypar/partition/refinement/spectral/solvers/gevp.h"
// #include "mt-kahypar/partition/refinement/spectral/solvers/slepc_gevp.cpp"
#include "mt-kahypar/partition/refinement/spectral/solvers/julia_gevp.cpp"


#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/preprocessing/sparsification/degree_zero_hn_remover.h"

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
    _gain.reset();


    // implementation goes here
    LOG << "Spectral Refiner called";
    utils::Timer& timer = utils::Utilities::instance().getTimer(_context.utility_id);
    timer.start_timer("partition_sp", "Partition");
    bool found_new_partition = partition(partionedHypergraph, best_metrics);
    timer.stop_timer("partition_sp");

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
    
    DBG << "Spectral Refiner finished";
    return delta > 0;
  }


  template <typename GraphAndGainTypes>
  bool SpectralRefiner<GraphAndGainTypes>::partition(PartitionedHypergraph& phg, Metrics &best_metrics) {
    // remove single pins
    auto dzhr = DegreeZeroHypernodeRemover<GraphAndGainTypes>(_context);
    size_t numRemovedSinglePins = dzhr.removeDegreeZeroHypernodes(phg.hypergraph(), true);
    LOG << "removed single pins: " << numRemovedSinglePins;

    DBG << "setting up...";

    if constexpr (Hypergraph::is_graph) {
      /* TODO */
    }

    const PartitionID k = phg.k(); /* TODO extract from argv? */    
    numNodes = 0;
    for (auto np = phg.nodes().begin(); *np < ~np; *(++np) & ++numNodes) {}

    HypernodeWeight total_weight = phg.totalWeight() - phg.hypergraph().weightOfRemovedDegreeZeroVertices();
    max_part_weight = floor(0.5 * (_context.partition.epsilon + 1) * total_weight);
    HypernodeWeight max_optimal_weight = ceil(total_weight / float(k));
    max_part_weight = max_part_weight < max_optimal_weight ? max_optimal_weight : max_part_weight;

    DBG << "node count calculated";

    vec<PartitionID> inputPartition;
    inputPartition.reserve(numNodes);
    for (const HypernodeID node : phg.nodes()) {
      inputPartition.push_back(phg.partID(node));
    }

    DBG << "hint read";

    // dehyperisation
    Operator inputGraphLaplacian(numNodes);
    dehyperizeToLaplacian(phg.hypergraph(), inputGraphLaplacian);
    // DBG << "---- Laplace";
    // inputGraphLaplacian.printMatrix([&](std::string s){DBG<<s;});

    DBG << "laplacian set up";

    // weight-balance graph construction
    Operator weightBalanceLaplacian(numNodes);
    buildWeightBalanceGraphLaplacian(phg.hypergraph(), weightBalanceLaplacian);

    DBG << "problem operators set";

    // actual refinement

    vec<vec<PartitionID>> candidateSolutions; /* TODO alias */
    candidateSolutions.reserve(1/* TODO argv "beta" */ + 1);
    candidateSolutions.push_back(inputPartition);
    vec<Gain> cut_sizes;
    cut_sizes.push_back(best_metrics.quality);
    Gain best_cutsize = -1;
    size_t best_index = -1;

    DBG << "solution vec prepared";

    for (int i = 0; i < params.numCandidates; i++) {
      DBG << "generating embedding...";
      vec<spectral::Vector> embedding; /* TODO type alias */
      if (k == 2) {
        generate2WayVertexEmbedding(phg.hypergraph(), weightBalanceLaplacian, inputGraphLaplacian, candidateSolutions.back(), cut_sizes.back(), embedding);
      } else {
        /* TODO */
      }

      /* the following is heavily adapted to julia not just solving the gevp but also generating a solution */
      // read solution from embedding, recalculate if not valid
      vec<PartitionID> newSolution;
      phg.resetPartition();
      Skalar eps = 1.0e-9;
      bool is_partition = true;
      for (auto nptr = phg.nodes().begin(), i = 0UL; i < numNodes; *(++nptr) & ++i) {
        if (!is_partition) {
          LOG << "recomputing solution...";
          newSolution.clear();
          generateSolution(phg, embedding, newSolution);
          break;
        }
        Skalar v = embedding.back()[i];
        if (v < eps && v > -eps) {
          newSolution.push_back(0);
        } else if (v < 1.0 + eps && v > 1.0 - eps) {
          newSolution.push_back(1);
        } else {
          is_partition = false;
          continue;
        }
        phg.setNodePart(*nptr, newSolution[i]);
        if (phg.partWeight(newSolution[i]) > max_part_weight) {
          is_partition = false;
        }
      }
      candidateSolutions.push_back(newSolution);


      // calulate metrics
      cut_sizes.push_back(metrics::quality(phg, _context, !_context.refinement.label_propagation.execute_sequential));
      if (best_index == -1 || cut_sizes.back() < best_cutsize) {
        best_cutsize = cut_sizes.back();
        best_index = candidateSolutions.size() - 1;
      }
    }

    LOG << "finished partitioning";

    bool found_new_partition = best_cutsize == best_metrics.quality && inputPartition != candidateSolutions[best_index];
    bool found_valid_solution = best_cutsize <= best_metrics.quality;

    LOG << "spectral results: "
      << (found_valid_solution ? "found valid solution, " : "")
      << (found_new_partition ? "found alternative partition, " : "")
      << V(best_cutsize) << ", "
      << V(best_metrics.quality);

    best_metrics.quality = found_valid_solution ? best_cutsize : best_metrics.quality;

    if (found_valid_solution && best_index != candidateSolutions.size() - 1) {
      setPartition(phg, candidateSolutions[best_index]);
    } else {
      setPartition(phg, inputPartition);
    }

    dzhr.restoreDegreeZeroHypernodes(phg);

    return found_valid_solution;
  }

  template <typename GraphAndGainTypes>
  template <typename Collection>
  void SpectralRefiner<GraphAndGainTypes>::setPartition(PartitionedHypergraph &phg, Collection &partition) {
    phg.resetPartition(); /* TODO change changed nodes only? */
    for (auto nptr = phg.nodes().begin(), i = 0UL; i < numNodes; *(++nptr) & ++i) {
      phg.setNodePart(*nptr, partition[i]);
    }
  }

  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::dehyperizeToLaplacian(Hypergraph& hypergraph, spectral::Operator& target) {
    // dehyperisation via clique expansion graph

    target.ctx[0] = (void *) &hypergraph;

    // matrix vector multiplication
    target.effects[0] = [](void *ctx, Vector& operand, Vector& target_vector) {
      size_t n = operand.dimension();
      Hypergraph *hg = (Hypergraph *) ctx;

      spectral::Skalar weight_factor = 1.0;

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

        spectral::Skalar clique_edge_weight = weight_factor * ((spectral::Skalar) hg->edgeWeight(he)) / (-1.0 + (spectral::Skalar) edge_size);

        for (const HypernodeID& pin : hg->pins(he)) {
          factor.set(pin, factor[pin] + clique_edge_weight * edge_size);
          subtrahend.set(pin, subtrahend[pin] + clique_edge_weight * operand_dot_e);
        }
      }

      // calculate result
      target_vector.reset();
      for (size_t i = 0; i < target_vector.dimension(); i++) {
        target_vector.set(i, target_vector[i] + operand[i] * factor[i] - subtrahend[i]);
      }
    };
    
    // matrix diagonal
    target.calc_diagonal_ops[0] = [] (void *ctx, Vector& target_vector) {
      Hypergraph *hg = (Hypergraph *) ctx;
      for (auto nptr = hg->nodes().begin(), i = 0UL; i < target_vector.dimension(); *(++nptr) & ++i) {
        for (const HyperedgeID& edge : hg->incidentEdges(*nptr)) {
          target_vector.set(i, target_vector[i] + hg->edgeWeight(edge));
        }
      }
    };

    // export hypergraph
    target.ctx_exporter[0] = [] (void *ctx, vec<size_t> &result) {
      Hypergraph *hg = (Hypergraph *) ctx;

      // format: n, m, node weights, edge weights, pin list indices, pin lists

      size_t n = 0;
      for (auto np = hg->nodes().begin(); *np < ~np; *(++np) & ++n) {}
      size_t m = hg->initialNumEdges();

      if (hg->is_graph) {
        // graph edges are stored directed
        m /= 2;
      }

      result.push_back(n);
      result.push_back(m);

      for (const HypernodeID &n : hg->nodes()) {
        result.push_back(hg->nodeWeight(n));
      }

      vec<size_t> pin_indices;
      vec<size_t> pin_lists;
      
      pin_indices.push_back(0);

      spectral::Vector node_index_map(hg->initialNumNodes());
      for (auto nptr = hg->nodes().begin(), i = 0UL; i < n; *(++nptr) & ++i) {
        node_index_map.set(*nptr, i);
      }

      for (const HyperedgeID &e : hg->edges()) {
        if (hg->is_graph) {
          auto pin_iter = hg->pins(e).begin();
          if (*pin_iter > *(++pin_iter)) {
            // backwards edge
            continue;
          }
        }

        result.push_back(hg->edgeWeight(e));
        pin_indices.push_back(pin_indices.back() + hg->edgeSize(e));
        for (const HypernodeID &pin : hg->pins(e)) {
          pin_lists.push_back(node_index_map[pin]);
        }
      }
      
      result.insert(result.end(), pin_indices.begin(), pin_indices.end());
      result.insert(result.end(), pin_lists.begin(), pin_lists.end());
    };
  }
  
  
  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::buildWeightBalanceGraphLaplacian(Hypergraph& hypergraph, spectral::Operator& target) {
    target.ctx[0] = (void *) &hypergraph;

    // matrix vector multiplication
    target.effects[0] = [](void *ctx, Vector& operand, Vector& target_vector) {
      Hypergraph *hg = (Hypergraph *) ctx;
      size_t dimension = operand.dimension();

      spectral::Skalar w_dot_x = 0.0;
      for (size_t i = 0; i < dimension; i++) {
        w_dot_x += ((spectral::Skalar) hg->nodeWeight(i)) * operand[i];
      }
      
      HypernodeWeight current_total_weight = hg->totalWeight() - hg->weightOfRemovedDegreeZeroVertices();
      for (auto nptr = hg->nodes().begin(), i = 0UL; i < dimension; *(++nptr) & ++i) {
        target_vector.set(i, target_vector[i] + hg->nodeWeight(*nptr) * (((spectral::Skalar) current_total_weight) * operand[i] - w_dot_x));
      }
    };

    // matrix diagonal
    target.calc_diagonal_ops[0] = [] (void *ctx, Vector& target_vector) {
      Hypergraph *hg = (Hypergraph *) ctx;

      HypernodeWeight current_total_weight = hg->totalWeight() - hg->weightOfRemovedDegreeZeroVertices();
      for (auto nptr = hg->nodes().begin(), i = 0UL; i < target_vector.dimension(); *(++nptr) & ++i) {
        spectral::Skalar weight = hg->nodeWeight(*nptr);
        target_vector.set(i, target_vector[i] + weight * (current_total_weight - weight));
      }
    };
  }


  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::generate2WayVertexEmbedding(Hypergraph &hypergraph, spectral::Operator& weightBalance, spectral::Operator& graphLaplacian, vec<PartitionID>& hintSolution, Gain &hint_quality, vec<spectral::Vector>& target) {
    /* this method is heavily adapted to julia not just solving the gevp but also generating a solution */
    
    // hint graph
    Operator balanceOperator(numNodes);
    generateHintGraphLaplacian(hintSolution, balanceOperator);
    // DBG << "---- hint";
    // balanceOperator.printMatrix([&](std::string s){DBG<<s;});
    // DBG << "---- weight";
    // weightBalance.printMatrix([&](std::string s){DBG<<s;});
    balanceOperator += weightBalance;
    // DBG << "---- balance";
    // balanceOperator.printMatrix([&](std::string s){DBG<<s;});
    
    /* TODO use GEVPSolver interface! */
    spectral::JuliaGEVPSolver solver;
    // spectral::SLEPcGEVPSolver solver;
    vec<spectral::Vector> known_evecs;
    vec<spectral::Skalar> known_evals;
    // trivial 1 0 epair
    known_evecs.push_back(spectral::Vector(numNodes, 1.0));
    known_evals.push_back(0.0);
    // target vector
    /* spectral::Vector hint(numNodes);
    for (size_t i = 0; i < numNodes; i++) {
      hint.set(i, hintSolution[i] == (PartitionID) 0 ? 1 : -1);
    }
    known_evecs.push_back(hint);
    known_evals.push_back(hint_quality); */
    // spectral::Operator dummy(numNodes); TODO flag
    solver.setProblem(graphLaplacian, balanceOperator, known_evecs, known_evals, known_evals.size()/*  - 1 */);

    DBG << "eigenproblem set";

    spectral::Skalar a;
    spectral::Vector v(numNodes);

    size_t num_evecs = 1; /* only fiedler */
    while (target.size() < num_evecs) {
      int info = solver.nextEigenpair(a, v, false); /* TODO assert return value is positive */
      target.push_back(v);
      // if (info == 0) {
      //   DBG << "---- Laplacian";
      //   graphLaplacian.printMatrix([&](std::string s){DBG<<s;});
      //   DBG << "---- balance";
      //   balanceOperator.printMatrix([&](std::string s){DBG<<s;});
      // }
      
    }
  }

  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::generateSolution(PartitionedHypergraph &phg, vec<spectral::Vector> &embedding, vec<PartitionID> &target) { /* TODO aliase */
    /* fielder distilling */
    
    // definitions
    spectral::Skalar split_index = 0;
    spectral::Vector &fiedler = embedding.back();
    vec<size_t> node_indices_fiedler;
    vec<HypernodeID> index_node_map(phg.nodes().begin(), phg.nodes().end());

    const auto objective_delta = [&](const SynchronizedEdgeUpdate& sync_update) {
      _gain.computeDeltaForHyperedge(sync_update);
    };    
    const auto move_next_node = [&](bool use_gain, bool dest) {
      !dest && split_index++;
      size_t node_index = node_indices_fiedler[split_index];
      HypernodeID node = index_node_map[node_index];
      target[node_index] = dest ? 1 : 0;
      if (use_gain) {
        phg.changeNodePart(node, dest ? 0 : 1, dest ? 1 : 0, objective_delta);
      } else {
        phg.changeNodePart(node, dest ? 0 : 1, dest ? 1 : 0);
      }
      dest && split_index --;
    };

    // sort by fiedler entry
    for (size_t i = 0; i < numNodes; node_indices_fiedler.push_back(i++)) {}
    std::sort(node_indices_fiedler.begin(), node_indices_fiedler.end(),
      [&](const size_t& u, const size_t& v) { return fiedler[u] < fiedler[v]; });

    // initialize with all nodes belonging to partition 1 except for first
    target.clear();
    target.resize(numNodes, 1);
    target[node_indices_fiedler[0]] = 0;
    setPartition(phg, target);

    // move into range
    while (phg.partWeight(1) > max_part_weight) {
      move_next_node(false, false);
    }

    if (numNodes < 20 && debug) {
      for (size_t i=0; i < fiedler.dimension(); i++) {
        DBG << target[i];
      }
    }
    
    // find optimum in range
    _gain.reset();
    size_t best_index = split_index;
    Gain best_delta = _gain.localDelta();
    while (phg.partWeight(0) + phg.nodeWeight(index_node_map[node_indices_fiedler[split_index + 1]]) <= max_part_weight) {
      move_next_node(true, false);

      if (_gain.localDelta() <= best_delta) {
        best_index = split_index;
        best_delta =_gain.localDelta();
      }  
    }
    
    // undo bad moves
    while (best_index < split_index) {
      move_next_node(false, true);
    }

    if (numNodes < 20 && debug) {
      for (const auto p : target) {
        DBG << p;
      }
    }
  }


  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::generateHintGraphLaplacian(vec<PartitionID> &hintSolution, spectral::Operator& target) {
    target.ctx[0] = (void *) &(hintSolution);

    // matrix vector multiplication
    target.effects[0] = [](void *ctx, Vector& operand, Vector& target_vector) {
      vec<PartitionID> *hint = (vec<PartitionID> *) ctx;
      
      // partition counts
      size_t n_0 = 0;
      size_t n_1;

      // partial sums of x
      spectral::Skalar x_0 = 0.0;
      spectral::Skalar x_1 = 0.0;

      for (size_t i = 0; i < hint->size(); i++) {
        if (hint->at(i) == (PartitionID) 0) {
          n_0++;
          x_0 += operand[i];
        } else {
          x_1 += operand[i];
        }
      }

      n_1 = operand.dimension() - n_0;

      for (size_t i = 0; i < target_vector.dimension(); i++) {
        bool p0 = hint->at(i) == (PartitionID) 0;
        target_vector.set(i, target_vector[i] + operand[i] * (p0 ? n_1 : n_0) - (p0 ? x_1 : x_0));
      }
    };

    // matrix diagonal
    target.calc_diagonal_ops[0] = [] (void *ctx, Vector& target_vector) {
      vec<PartitionID> *hint = (vec<PartitionID> *) ctx;
      
      // partition counts
      size_t n_0 = 0;
      size_t n_1;

      for (size_t i = 0; i < hint->size(); i++) {
        if (hint->at(i) == (PartitionID) 0) {
          n_0++;
        }
      }

      n_1 = hint->size() - n_0;

      for (size_t i = 0; i < target_vector.dimension(); i++) {
        bool p0 = hint->at(i) == (PartitionID) 0;
        target_vector.set(i, target_vector[i] + (p0 ? n_1 : n_0));
      }
    };

    // export hint partition
    target.ctx_exporter[0] = [] (void *ctx, vec<size_t> &target_vector) {
      vec<PartitionID> *hint = (vec<PartitionID> *) ctx;

      target_vector.insert(target_vector.begin(), hint->begin(), hint->end());
    };
  }

  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::readConfigFile() {
    /* TODO due to procastination, the config file param is currently interpreted as beta */
    if (_context.refinement.spectral.config_path.length() > 0) {
      params.numCandidates = (size_t) std::stoi(_context.refinement.spectral.config_path, nullptr);
    }
  }


  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::initializeImpl(mt_kahypar_partitioned_hypergraph_t& phg) {
    unused(phg);
  }

  template <typename GraphAndGainTypes>
  void SpectralRefiner<GraphAndGainTypes>::initializeSolver() {
    JuliaGEVPSolver::initialize();
  }


  template <typename GraphAndGainTypes>
  SpectralRefiner<GraphAndGainTypes>::~SpectralRefiner() {
    JuliaGEVPSolver::exit();
  }

  namespace {
  #define SPECTRAL_REFINER(X) SpectralRefiner<X>
  }

  // explicitly instantiate so the compiler can generate them when compiling this cpp file
  INSTANTIATE_CLASS_WITH_VALID_TRAITS(SPECTRAL_REFINER)
}

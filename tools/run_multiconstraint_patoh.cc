/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2025 Nikolai Maas <nikolai.maas@kit.edu>
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

// Note: the file was created with the help of GPT o3


#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#include <boost/program_options.hpp>

#include "patoh/patoh.h"

#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_graph.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/weight/hypernode_weight_common.h"

using namespace mt_kahypar;
namespace po = boost::program_options;

template <class TypeTraits>
std::pair<HyperedgeWeight, HypernodeWeightArray> computeObjectiveAndWeight(typename TypeTraits::Hypergraph& hg,
                                                                           const std::vector<int>& parts,
                                                                           int k,
                                                                           Objective objective) {
  using PHG = typename TypeTraits::PartitionedHypergraph;
  ASSERT(parts.size() == hg.initialNumNodes());

  PHG phg(k, hg, parallel_tag_t{});
  hg.doParallelForAllNodes([&](const HypernodeID hn) {
    phg.setOnlyNodePart(hn, parts[hn]);
  });
  phg.initializePartition();

  HyperedgeWeight obj = metrics::quality(phg, objective);
  std::cout << "Cut: " << metrics::quality(phg, Objective::cut) << '\n';
  std::cout << "Connectivity: " << metrics::quality(phg, Objective::km1) << '\n';
  HypernodeWeightArray weights = phg.partWeights().copy();
  return {obj, std::move(weights)};
}

/*********************************************************************
 *  transformToPaToHInput
 *
 *  Fills        _c, _n, _nconst
 *  allocates    cwghts  (size = _c * _nconst)
 *               nwghts  (size = _n)
 *               xpins   (size = _n + 1)
 *               pins    (size = ∑ |net|)
 *
 *  Memory must be released by the caller with  free().
 ********************************************************************/
template <class Hypergraph>
void transformToPaToHInput(const Hypergraph& hg,
                           int*        _c,        int*        _n,
                           int*        _nconst,   int**       cwghts,
                           int**       nwghts,    int**       xpins,
                           int**       pins)
{
  /* ------------------------------------------------------------
   * 1.  Basic sizes
   * ----------------------------------------------------------*/
  const int c        = static_cast<int>(hg.initialNumNodes());
  const int n        = static_cast<int>(Hypergraph::is_graph ? hg.initialNumEdges() / 2 : hg.initialNumEdges());
  const int n_const  = static_cast<int>(hg.dimension());   // = #constraints

  *_c      = c;
  *_n      = n;
  *_nconst = n_const ? n_const : 1;                        // PaToH needs >=1

  /* ------------------------------------------------------------
   * 2.  Allocate PaToH arrays
   * ----------------------------------------------------------*/
  *cwghts  = static_cast<int*>(malloc(sizeof(int) * c * *_nconst));
  *nwghts  = static_cast<int*>(malloc(sizeof(int) * n));
  *xpins   = static_cast<int*>(malloc(sizeof(int) * (n + 1)));

  // total number of pins (hypergraph incidence entries)
  std::size_t tot_pins = hg.initialNumPins();

  *pins    = static_cast<int*>(malloc(sizeof(int) * tot_pins));

  /* ------------------------------------------------------------
   * 3.  Fill cwghts   (node weights)
   * ----------------------------------------------------------*/
  for (HypernodeID v : hg.nodes()) {
    const auto& w = hg.nodeWeight(v);          // weight is n_const-dim vector
    for (Dimension d = 0; d < static_cast<Dimension>(*_nconst); ++d) {
      (*cwghts)[v * *_nconst + d] = static_cast<int>(w.at(d));
    }
  }

  /* ------------------------------------------------------------
   * 4.  Fill nwghts, xpins, pins   (edge weights & incidence)
   * ----------------------------------------------------------*/
  std::size_t pin_idx = 0;
  (*xpins)[0] = 0;

  auto to_uniq_id = [&](const HyperedgeID he) {
    if constexpr (Hypergraph::is_graph) return hg.uniqueEdgeID(he);
    else return he;
  };
  for (HyperedgeID he : hg.edges()) {
    if constexpr (Hypergraph::is_graph) {
      if (hg.edgeSource(he) > hg.edgeTarget(he)) continue;
    }

    const int e = static_cast<int>(he);                      // contiguous id
    (*nwghts)[to_uniq_id(e)] = static_cast<int>(hg.edgeWeight(he));

    for (HypernodeID v : hg.pins(he)) {
      (*pins)[pin_idx++] = static_cast<int>(v);              // 0-based
    }
    (*xpins)[to_uniq_id(e) + 1] = static_cast<int>(pin_idx);   // next start pos
  }
  ALWAYS_ASSERT(pin_idx == tot_pins);
}

/*--------------------------------------------------------------------------*/
/* Small helper to print weight statistics – adapted from PaToH example     */
/*--------------------------------------------------------------------------*/
void print_info(int k, const std::vector<int>& partweights, int nconst) {
  std::vector<int> avg(nconst, 0.0);
  std::vector<double> maxi(nconst, 0.0);

  for (int p = 0; p < k; ++p)
    for (int j = 0; j < nconst; ++j)
      avg[j] += partweights[p * nconst + j];

  for (int j = 0; j < nconst; ++j) {
    avg[j] = std::ceil(avg[j] /= static_cast<double>(k));
  }

  std::cout << "\n------------------------------------------------------------\n";
  std::cout << " AvgWeights: ";
  for (int j = 0; j < nconst; ++j) {
    std::cout << std::setw(10) << avg[j] << ' ';
  }
  std::cout << '\n';

  std::cout << " PartWeights:\n";
  for (int p = 0; p < k; ++p) {
    std::cout << "  P" << std::setw(3) << p << ": ";
    for (int j = 0; j < nconst; ++j) {
      double imbalance =
          (static_cast<double>(partweights[p * nconst + j]) - avg[j]) / avg[j];
      maxi[j] = std::max(maxi[j], imbalance);
      std::cout << std::setw(10) << partweights[p * nconst + j] << ' ';
    }
    std::cout << '\n';
  }
  std::cout << " Max. imbalance (percent): ";
  for (int j = 0; j < nconst; ++j)
    std::cout << std::fixed << std::setprecision(3) << 100.0 * maxi[j] << ' ';
  std::cout << "\n------------------------------------------------------------\n";
}

using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

int main(int argc, char* argv[]) try {
  /*------------------------------------------------------------------*/
  /* 1. Parse command-line arguments (already declared in stub)        */
  /*------------------------------------------------------------------*/
  std::string hgr_filename, preset, objective;
  int k, seed;
  double epsilon;
  bool verbose = false;
  bool is_metis = false;

  po::options_description opts("Options");
  opts.add_options()
    ("hypergraph,h", po::value<std::string>(&hgr_filename)->required(),
     "Hypergraph filename (.hgr/.hg)")
    ("preset,p",     po::value<std::string>(&preset)->required(),
     "PaToH preset: default | quality | speed")
    ("objective,o",  po::value<std::string>(&objective)->required(),
     "Objective: cut | km1")
    ("epsilon,e",    po::value<double>(&epsilon)->required(), "Allowed imbalance")
    ("n_parts,k",    po::value<int>(&k)->required(), "Number of parts")
    ("seed,s",       po::value<int>(&seed)->required(), "Random seed")
    ("verbose,v",    po::bool_switch(&verbose), "Verbose output")
    ("is_metis",     po::value<bool>(&is_metis), "Metis input file");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, opts), vm);
  po::notify(vm);

  /*------------------------------------------------------------------*/
  /* 2. Read hypergraph                                               */
  /*------------------------------------------------------------------*/
  int _c = 0, _n = 0, _nconst = 0;
  int *cwghts = nullptr, *nwghts = nullptr, *xpins = nullptr, *pins = nullptr;

  if (is_metis) {
    auto graph = io::readInputFile<ds::StaticGraph>(hgr_filename, FileFormat::Metis, true);
    transformToPaToHInput(graph, &_c, &_n, &_nconst, &cwghts, &nwghts, &xpins, &pins);
  } else {
    auto hypergraph = io::readInputFile<ds::StaticHypergraph>(hgr_filename, FileFormat::hMetis, true);
    transformToPaToHInput(hypergraph, &_c, &_n, &_nconst, &cwghts, &nwghts, &xpins, &pins);
  }

  /*------------------------------------------------------------------*/
  /* 3. Initialise PaToH parameters                                   */
  /*------------------------------------------------------------------*/
  PaToH_Parameters args;
  ALWAYS_ASSERT(objective == "km1" || objective == "cut");
  int cuttype = (objective == "km1") ? PATOH_CONPART : PATOH_CUTPART;

  ALWAYS_ASSERT(preset == "default" || preset == "quality" || preset == "speed");
  PaToH_Initialize_Parameters(&args, cuttype,
    (preset == "speed")   ? PATOH_SUGPARAM_SPEED :
    (preset == "quality") ? PATOH_SUGPARAM_QUALITY :
                            PATOH_SUGPARAM_DEFAULT);

  args._k    = k;
  args.seed  = seed;
  args.final_imbal = epsilon;
  args.MemMul_CellNet = 100;
  args.MemMul_Pins = 100;
  args.MemMul_General = 100;
  args.balance = PATOH_BALANCE_STRICT;
  if (!verbose) args.outputdetail = PATOH_OD_LOW;

  /* You may fine-tune further parameters here, e.g. balance, coarsening … */

  /*------------------------------------------------------------------*/
  /* 4. Allocate PaToH internal memory                                */
  /*------------------------------------------------------------------*/
  if (PaToH_Alloc(&args, _c, _n, _nconst,
                  cwghts, nwghts, xpins, pins) != 0) {
    std::cerr << "ERROR: PaToH_Alloc failed\n";
    return EXIT_FAILURE;
  }

  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();

  /*------------------------------------------------------------------*/
  /* 5. Partition                                                     */
  /*------------------------------------------------------------------*/
  std::vector<int> partvec(_c, 0);
  std::vector<int> partweights(k * _nconst, 0);
  int cut = 0;

  /* PaToH_Part is the unified entry; if _nconst>1 it automatically    */
  /* calls the multi-constraint variant inside the library [1].        */
  if (PaToH_Part(&args, _c, _n, _nconst, /*useFixCells*/ 0,
                 cwghts, nwghts, xpins, pins,
                 /*targetweights*/ nullptr,
                 partvec.data(), partweights.data(), &cut) != 0) {
    std::cerr << "ERROR: PaToH_Part failed\n";
    PaToH_Free();
    return EXIT_FAILURE;
  }

  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds(end - start);

  /*------------------------------------------------------------------*/
  /* 6. Verify result                                                 */
  /*------------------------------------------------------------------*/
  int recomputed_objective = 0;
  HypernodeWeightArray recomputed_weights;

  if (is_metis) {
    auto graph = io::readInputFile<ds::StaticGraph>(hgr_filename, FileFormat::Metis, true);
    std::tie(recomputed_objective, recomputed_weights) = computeObjectiveAndWeight<StaticGraphTypeTraits>(
      graph, partvec, k, (objective == "km1") ? Objective::km1 : Objective::cut);
  } else {
    auto hypergraph = io::readInputFile<ds::StaticHypergraph>(hgr_filename, FileFormat::hMetis, true);
    std::tie(recomputed_objective, recomputed_weights) = computeObjectiveAndWeight<StaticHypergraphTypeTraits>(
      hypergraph, partvec, k, (objective == "km1") ? Objective::km1 : Objective::cut);
  }
  ALWAYS_ASSERT(cut == recomputed_objective, V(cut) << V(recomputed_objective));
  for (int part = 0; part < k; ++part) {
    for (Dimension d = 0; d < static_cast<Dimension>(_nconst); ++d) {
      ALWAYS_ASSERT(partweights[part * _nconst + d] == recomputed_weights[part].at(d));
    }
  }

  /*------------------------------------------------------------------*/
  /* 7. Output result                                                 */
  /*------------------------------------------------------------------*/
  std::cout << "Time: " << elapsed_seconds.count() << '\n';
  print_info(k, partweights, _nconst);

  /*------------------------------------------------------------------*/
  /* 8. Clean up                                                      */
  /*------------------------------------------------------------------*/
  free(cwghts);  free(nwghts);  free(xpins);  free(pins);
  PaToH_Free();

  return EXIT_SUCCESS;
}
/*--------------------------------------------------------------------*/
catch (const std::exception& ex) {
  std::cerr << "FATAL: " << ex.what() << '\n';
  return EXIT_FAILURE;
}

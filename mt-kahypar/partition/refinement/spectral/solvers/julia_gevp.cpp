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

#include "mt-kahypar/partition/refinement/spectral/solvers/julia_gevp.h"
#include "mt-kahypar/partition/refinement/spectral/algebraic_wrappers/vector.cpp" /* TODO imports should work otherwise... */
#include "mt-kahypar/partition/refinement/spectral/algebraic_wrappers/operator.cpp" /* TODO imports should work otherwise... */

#include "mt-kahypar/partition/refinement/spectral/datatypes.h"

#include <filesystem>

#include <julia.h>
// JULIA_DEFINE_FAST_TLS
// #include <jluna.hpp> doesnt work with cpp17 :(


namespace mt_kahypar {
namespace spectral {

// using namespace jluna;    :/

bool JuliaGEVPSolver::julia_initialized = false;

void JuliaGEVPSolver::setProblem(Operator& a, Operator& b) {
  if (!julia_initialized) {
    

    julia_initialized = true;
  }
  

  op_a = &a;
  op_b = &b;

  solved = false;
  epairs_found = 0;
}

void JuliaGEVPSolver::setProblem(Operator& a, Operator& b, vec<Vector>& known_evecs, vec<Skalar> &known_evals, size_t deflation_epairs) {
  setProblem(a, b);

  evecs.insert(evecs.end(), known_evecs.begin(), known_evecs.end());
  evals.insert(evals.end(), known_evals.begin(), known_evals.end());

  num_deflation_epairs = deflation_epairs;
}

int JuliaGEVPSolver::nextEigenpair(Skalar& eval, Vector& evec, bool try_from_above) {
  if (!solved) {
    /* TODO try from above */
    solve();
  }
  eval = evals.back();
  evec = evecs.back();
  return 0;
}

void JuliaGEVPSolver::solve() {
  bool debug = false;

  DBG << "exporting contexts...";

  vec<uint64_t> hgr; /* TODO type */
  op_a->exportContext(0, hgr);
  size_t n = hgr[0];
  size_t m = hgr[1];

  DBG << "laplacian exported";

  vec<uint64_t> hint;
  op_b->exportContext(0, hint);

  DBG << "hint exported";

  vec<double> deflation_evecs;
  for (size_t i = 0; i < num_deflation_epairs; i++) {
    deflation_evecs.insert(deflation_evecs.end(), evecs[i].get_all(), evecs[i].get_all() + n);
  }

  DBG << "deflation space prepared";

  // initialize array objects and their types
  
  jl_value_t *node_array_type = jl_apply_array_type((jl_value_t *) jl_uint64_type, 1); /* TODO check sizeof HypernodeID */
  jl_value_t *double_array_type = jl_apply_array_type((jl_value_t *) jl_float64_type, 1); 
  
  jl_array_t *hgr_jl = jl_ptr_to_array_1d(node_array_type, hgr.data(), hgr.size(), 0);
  jl_array_t *hint_jl = jl_ptr_to_array_1d(node_array_type, hint.data(), hint.size(), 0);
  jl_array_t *constraints_jl = jl_ptr_to_array_1d(double_array_type, deflation_evecs.data(), num_deflation_epairs * n, 0);

  jl_function_t *solve = (jl_function_t *) jl_eval_string("main_auto");

  DBG << "launching Julia code...";

  jl_array_t *evecs_jl = (jl_array_t *) jl_call3(solve, (jl_value_t *) hgr_jl, (jl_value_t *) hint_jl, (jl_value_t *) constraints_jl);

  DBG << "returning result...";

  JL_GC_PUSH1(evecs_jl);

  Vector result(n);
  result.set_all((Skalar *) jl_array_data(evecs_jl));
  evecs.push_back(result);

  JL_GC_POP();
}

JuliaGEVPSolver::~JuliaGEVPSolver() {
  if (julia_initialized) {
    // jl_atexit_hook(0);
  }
}

void JuliaGEVPSolver::initialize() {
  jl_init();
  std::filesystem::path file_path = __FILE__;
  file_path = file_path.parent_path() / "julia";
  jl_eval_string(("cd(\"" + file_path.string() + "\");include(\"main.jl\")").c_str());
}

void JuliaGEVPSolver::exit() {
  // jl_atexit_hook(0);
}

}
}

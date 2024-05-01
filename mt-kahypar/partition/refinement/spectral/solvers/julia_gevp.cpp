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

#include <julia.h>
// #include <jluna.hpp>


namespace mt_kahypar {
namespace spectral {

// using namespace jluna;

bool JuliaGEVPSolver::julia_initialized = false;

void JuliaGEVPSolver::setProblem(Operator& a, Operator& b) {
  if (!julia_initialized) {
    jl_init();
    /* TODO path */
    jl_eval_string("cd(\"/home/julian/Dokumente/Studium/BA/mt-kahypar/mt-kahypar/partition/refinement/spectral/solvers/julia/\");include(\"lobpcg.jl\")");

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
  bool debug = true;

  vec<uint64_t> hgr; /* TODO type */
  op_a->exportContext(0, hgr);
  size_t n = hgr[0];
  size_t m = hgr[1];

  vec<uint64_t> hint;
  op_b->exportContext(0, hint);

  vec<double> deflation_evecs;
  for (size_t i = 0; i < num_deflation_epairs; i++) {
    deflation_evecs.insert(deflation_evecs.end(), evecs[i].get_all(), evecs[i].get_all() + n);
  }


  jl_value_t *node_array_type = jl_apply_array_type((jl_value_t *) jl_uint64_type, 1); /* TODO check sizeof HypernodeID */
  jl_value_t *double_array_type = jl_apply_array_type((jl_value_t *) jl_float64_type, 1); 

  /* jl_array_t *node_weights = jl_ptr_to_array_1d(node_array_type, hgr.data() + 2, n, 0);
  jl_array_t *edge_weights = jl_ptr_to_array_1d(node_array_type, hgr.data() + 2 + n, m, 0);
  jl_array_t *pin_indices = jl_ptr_to_array_1d(node_array_type, hgr.data() + 2 + n + m, m + 1, 0);
  jl_array_t *pin_lists = jl_ptr_to_array_1d(node_array_type, hgr.data() + 2 + n + m + (m + 1), hgr.size() - (2 + n + m + (m + 1)), 0); */
  
  jl_array_t *hgr_jl = jl_ptr_to_array_1d(node_array_type, hgr.data(), hgr.size(), 0);
  jl_array_t *hint_jl = jl_ptr_to_array_1d(node_array_type, hint.data(), hint.size(), 0);
  jl_array_t *constraints_jl = jl_ptr_to_array_1d(double_array_type, deflation_evecs.data(), num_deflation_epairs * n, 0);

  // jl_value_t *res = jl_eval_string("typeof(include('')) == Module ? sqrt(2.) : 0.");
  // DBG << (jl_typeis(res, jl_float64_type) ? jl_unbox_float64(res) : -1.);
  // jl_function_t *func = (jl_function_t *) jl_get_function(jl_main_module, "pwd");
  // jl_value_t *res = (jl_value_t *) jl_call0(func);
  // DBG << ((const char *) jl_string_ptr(res));
  
  // jl_function_t *print = (jl_function_t *) jl_eval_string("print");
  // jl_call2(print, (jl_value_t *) hgr_jl, (jl_value_t *) hint_jl);

  // jl_module_t *GEVP = (jl_module_t *) jl_eval_string("GEVP");
  // jl_function_t *func = (jl_function_t *) jl_get_function(GEVP, "test_julia_from_c");//"solve_lobpcg");
  // jl_function_t *test = (jl_function_t *) jl_eval_string("test_julia_from_c");
  // jl_call0(test);
  // jl_eval_string("using GraphSignals");

  jl_function_t *solve = (jl_function_t *) jl_eval_string("solve_lobpcg");
  DBG << "launching Julia code...";
  jl_array_t *evecs_jl = (jl_array_t *) jl_call3(solve, (jl_value_t *) hgr_jl, (jl_value_t *) hint_jl, (jl_value_t *) constraints_jl);

  Vector result(n);
  result.set_all((Skalar *) jl_array_data(evecs_jl));
  evecs.push_back(result);
}

JuliaGEVPSolver::~JuliaGEVPSolver() {
  if (julia_initialized) {
    // jl_atexit_hook(0);
  }
}

}
}

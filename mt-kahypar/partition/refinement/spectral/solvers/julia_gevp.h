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

#include "mt-kahypar/partition/refinement/spectral/solvers/gevp.h"

namespace mt_kahypar {
namespace spectral {

class JuliaGEVPSolver : public GEVPSolver {
 public:
  void setProblem(Operator& a, Operator& b) final;

  void setProblem(Operator& a, Operator& b, vec<Vector>& evecs, vec<Skalar> &evals, size_t num_deflation_epairs) final;

  int nextEigenpair(Skalar& eval, Vector& evec, bool try_from_above) final;

  ~JuliaGEVPSolver() final;
  
 private:
 
  static bool julia_initialized;

  // attributes
  Operator *op_a = nullptr;
  Operator *op_b = nullptr;

  bool solved = false;
  vec<Skalar> evals;
  vec<Vector> evecs;
  size_t epairs_found;
  size_t num_deflation_epairs;

  // methods
  void solve();
};

}
}
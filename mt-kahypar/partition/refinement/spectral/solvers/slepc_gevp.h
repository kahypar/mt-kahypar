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


#include <slepceps.h>

#include "mt-kahypar/partition/refinement/spectral/solvers/gevp.h"

#define CallPetsc(...) PetscCallAbort(PETSC_COMM_WORLD, __VA_ARGS__)

namespace mt_kahypar {
namespace spectral {

class SLEPcGEVPSolver : public spectral::GEVPSolver {
 public:
  void initialize(spectral::Matrix& a, spectral::Matrix& b) final;

  bool nextEigenpair(spectral::Skalar& eval, spectral::Vector& evec) final;

  SLEPcGEVPSolver();
  ~SLEPcGEVPSolver() final;

 private:
  spectral::Matrix *_a;
  spectral::Matrix *_b;

  bool _solved = false;

  EPS _context;

  void solve();
};

}
}

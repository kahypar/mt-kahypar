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

class SLEPcGEVPSolver : public GEVPSolver {
 public:
  void setProblem(Operator& a, Operator& b) final;

  bool nextEigenpair(Skalar& eval, Vector& evec) final;

  ~SLEPcGEVPSolver() final;

 private:
  // attributes

  Operator *op_a;
  Operator *op_b;

  Mat mat_A;
  Mat mat_B;

  EPS eps;
  bool solved = false;

  // constants

  static constexpr MPI_Comm& GLOBAL_COMMUNICATOR = PETSC_COMM_WORLD;

  // matrix storage size calculation strategys
  static size_t (*GET_N_NAIVE) (size_t);
  /* ... */

  // problem type strategys
  static EPSProblemType (*GET_PBT_NAIVE) (Operator&, Operator&);
  /* ... */

  // eps type strategys
  static EPSType (*GET_TYPE_NAIVE) (Operator&, Operator&);

  // methods

  void solve();

  static void vector2Vec(Vector& vector, Vec& vec);
  static void vec2vector(Vec& vec, Vector& vector);

  static void getMatContext(Mat& M, Operator *op);

  
  // lambdas
  static size_t (*getN) (size_t);
  static EPSProblemType (*getProblemType) (Operator&, Operator&);
  static EPSType (*getEpsType) (Operator&, Operator&);

  static int (*matMult) (Mat& M, Vec& x, Vec& y);
  static int (*matMultTranspose) (Mat& M, Vec& x, Vec& y);
};

}
}

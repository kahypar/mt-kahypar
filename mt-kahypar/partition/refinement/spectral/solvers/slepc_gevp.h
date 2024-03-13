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
#define Range(length) PetscInt range[length]; for (PetscInt i = 0; i < length; i++) { range[i] = i; }

namespace mt_kahypar {
namespace spectral {

class SLEPcGEVPSolver : public GEVPSolver {
 public:
  void setProblem(Operator& a, Operator& b) final;

  void setProblem(Operator& a, Operator& b, vec<Vector>& trivial_evecs, vec<Skalar> &trivial_evals) final;

  int nextEigenpair(Skalar& eval, Vector& evec) final;

  ~SLEPcGEVPSolver() final;

 private:
  // attributes

  bool slepc_running = false;
  bool solved = false;
  bool tried_from_above = false;
  vec<PetscScalar> evals;
  vec<Vec> evecs;

  // custom matrices
  Operator *op_a = nullptr;
  Operator *op_b = nullptr;

  // petsc matrices
  Mat mat_A = nullptr;
  Mat mat_B = nullptr;

  // petsc eigen problem solver
  EPS eps;


  // methods

  void start_slepc();
  void solve();
  void reset_matrices();
  void end_slepc();

  // conversion
  /* TODO maybe not needed? */
  static void vector2Vec(Vector& vector, Vec vec);
  static void vec2vector(Vec vec, Vector& vector);

  static void getMatContext(Mat M, Operator **op);

  // algebraic operations
  static PetscErrorCode matMult(Mat M, Vec x, Vec y);
  static PetscErrorCode matMultTranspose(Mat M, Vec x, Vec y);
  static PetscErrorCode getDiagonal(Mat M, Vec diag);

  // lambdas
  size_t (*getN) (size_t) = [](size_t n) { throw "unassigned"; return n; };
  EPSProblemType (*getProblemType) (Operator&, Operator&) = [](Operator& a, Operator& b) { throw "unassigned"; return EPS_HEP; };
  EPSType (*getEpsType) (Operator&, Operator&) = [](Operator& a, Operator& b) { throw "unassigned"; return EPSLOBPCG; };


  // constants

  static constexpr MPI_Comm& GLOBAL_COMMUNICATOR = PETSC_COMM_WORLD;

  // matrix storage size calculation strategys
  static constexpr size_t (*GET_N_NAIVE) (size_t) = [](size_t n) { return n/*n*/; };
  /* ... */

  // problem type strategys
  static constexpr EPSProblemType (*GET_PBT_NAIVE) (Operator&, Operator&) = [](Operator& a, Operator& b) {
    return a.isSymmetric() && b.isSymmetric() ? EPS_GHEP : EPS_GNHEP;
  };
  /* ... */

  // eps type strategys
  static constexpr EPSType (*GET_TYPE_NAIVE) (Operator&, Operator&) = [](Operator& a, Operator& b) {
    return EPSLOBPCG;
  };
};

}
}

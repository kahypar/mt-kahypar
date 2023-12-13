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

#include "mt-kahypar/partition/refinement/spectral/solvers/slepc_gevp.h"
#include "mt-kahypar/partition/refinement/spectral/algebraic_wrappers/vector.cpp" /* TODO imports should work otherwise... */
#include "mt-kahypar/partition/refinement/spectral/algebraic_wrappers/operator.cpp" /* TODO imports should work otherwise... */

namespace mt_kahypar {
namespace spectral {

SLEPcGEVPSolver::SLEPcGEVPSolver() {
  PetscFunctionBeginUser;

  /* TODO parse args */

  CallPetsc(SlepcInitialize(0, nullptr, (char*)0, nullptr)); /* TODO cl args */
  CallPetsc(EPSCreate(GLOBAL_COMMUNICATOR, &eps));

  SLEPcGEVPSolver::getN = GET_N_NAIVE;
  SLEPcGEVPSolver::getProblemType = GET_PBT_NAIVE;
}

void SLEPcGEVPSolver::initialize(Operator& a, Operator& b) {
  PetscFunctionBeginUser;

  // validation

  /*
  assert a->dimension() == b->dimension()
  assert a->isSymmetric() && b->isSymmetric()
  */

  // definitions
  
  op_a = &a;
  op_b = &b;
  solved = false;

  size_t n = op_a->dimension();

  // create matrices

  /* CallPetsc(MatCreateShell(GLOBAL_COMMUNICATOR, SLEPcGEVPSolver::getN(n), SLEPcGEVPSolver::getN(n), SLEPcGEVPSolver::getN(n), SLEPcGEVPSolver::getN(n), &n, &mat_A));
  CallPetsc(MatShellSetContext(mat_A, (void *) op_a));
  CallPetsc(MatShellSetOperation(mat_A, MATOP_MULT, (void(*)(void)) &matMult));
  CallPetsc(MatShellSetOperation(mat_A, MATOP_MULT_TRANSPOSE,(void(*)(void)) &matMultTranspose));
  // CallPetsc(MatShellSetOperation(mat_A,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_Laplacian2D));

  CallPetsc(MatCreateShell(GLOBAL_COMMUNICATOR, SLEPcGEVPSolver::getN(n), SLEPcGEVPSolver::getN(n), SLEPcGEVPSolver::getN(n), SLEPcGEVPSolver::getN(n), &n, &mat_B));
  CallPetsc(MatShellSetContext(mat_B, (void *) op_b));
  CallPetsc(MatShellSetOperation(mat_B, MATOP_MULT, (void(*)(void)) &matMult));
  CallPetsc(MatShellSetOperation(mat_B, MATOP_MULT_TRANSPOSE,(void(*)(void)) &matMultTranspose));
  // CallPetsc(MatShellSetOperation(mat_B,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_Laplacian2D));*/

  // set options

  CallPetsc(EPSSetType(eps, SLEPcGEVPSolver::getEpsType(*op_a, *op_b)));
  CallPetsc(EPSSetProblemType(eps, SLEPcGEVPSolver::getProblemType(*op_a, *op_b)));

  /*CallPetsc(EPSSetOperators(eps, mat_A, mat_B));*/

  CallPetsc(EPSSetFromOptions(eps)); /* TODO really needed? */ 
}

bool SLEPcGEVPSolver::nextEigenpair(Skalar& eval, Vector& evec) {
  PetscFunctionBeginUser;

  if (!solved) {
    solve();
  }

  // ...
  return true;
}

SLEPcGEVPSolver::~SLEPcGEVPSolver() {
  PetscFunctionBeginUser;

  CallPetsc(EPSDestroy(&eps));
  CallPetsc(SlepcFinalize());
}

void SLEPcGEVPSolver::solve() {
  PetscFunctionBeginUser;

  // solve

  /* CallPetsc(EPSSolve(eps)); */
  solved = true;
}

void SLEPcGEVPSolver::vector2Vec(Vector& vector, Vec& vec) {
  /* TODO */
}

void SLEPcGEVPSolver::vec2vector(Vec& vec, Vector& vector) {
  /* TODO */
}

int (*SLEPcGEVPSolver::matMult) (Mat&, Vec&, Vec&) = [](Mat& M, Vec& x, Vec& y) {
  PetscFunctionBeginUser;

  Operator *op_m;
  getMatContext(M, op_m);

  Vector sp_x(op_m->dimension());
  Vector sp_y(op_m->dimension());
  SLEPcGEVPSolver::vec2vector(x, sp_x);
  op_m->apply(sp_x, sp_y);
  SLEPcGEVPSolver::vector2Vec(sp_y, y);
  return 0; /* TODO error code */
};

int (*SLEPcGEVPSolver::matMultTranspose) (Mat&, Vec&, Vec&) = [](Mat& M, Vec& x, Vec& y) {
  PetscFunctionBeginUser;

  Operator *op_m;
  getMatContext(M, op_m);

  if (op_m->isSymmetric()) {
    return SLEPcGEVPSolver::matMult(M, x, y);
  } else {
    return -1; /* TODO */
  }
};

 void SLEPcGEVPSolver::getMatContext(Mat& M, Operator *op) {
  PetscFunctionBeginUser;

  void *ctx;
  CallPetsc(MatShellGetContext(M, &ctx));
  *op = *((Operator *)ctx);
}


size_t (*SLEPcGEVPSolver::GET_N_NAIVE) (size_t) = [](size_t n) { return n*n; };

EPSProblemType (*SLEPcGEVPSolver::GET_PBT_NAIVE) (Operator&, Operator&) = [](Operator& a, Operator& b) {
  return a.isSymmetric() && b.isSymmetric() ? EPS_GHEP : EPS_GNHEP;
};

EPSType (*SLEPcGEVPSolver::GET_TYPE_NAIVE) (Operator&, Operator&) = [](Operator& a, Operator& b) {
  return EPSLOBPCG;
};

size_t (*SLEPcGEVPSolver::getN) (size_t) = GET_N_NAIVE;
EPSProblemType (*SLEPcGEVPSolver::getProblemType) (Operator&, Operator&) = GET_PBT_NAIVE;
EPSType (*SLEPcGEVPSolver::getEpsType) (Operator&, Operator&) = GET_TYPE_NAIVE;

}
}

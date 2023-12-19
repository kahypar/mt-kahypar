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

void SLEPcGEVPSolver::start_slepc() {
  PetscFunctionBeginUser;

  /* TODO parse cl args */

  PetscBool slepc_running_but_with_pets_datatype;
  SlepcInitialized(&slepc_running_but_with_pets_datatype);

  if (slepc_running_but_with_pets_datatype == PETSC_FALSE) {
    CallPetsc(SlepcInitialize(0, nullptr, (char*)0, nullptr)); /* TODO pass cl args? */
  } else {
    /* TODO
    if(EPSIsCreated(&eps)) {
    CallPetsc(EPSDestroy(&eps));
    }*/
  }

  CallPetsc(EPSCreate(GLOBAL_COMMUNICATOR, &eps));

  getN = GET_N_NAIVE;
  getProblemType = GET_PBT_NAIVE;
  getEpsType = GET_TYPE_NAIVE;

  slepc_running = true;

  PetscFunctionReturnVoid();
}

void SLEPcGEVPSolver::setProblem(Operator& a, Operator& b) {
  PetscFunctionBeginUser;

  // initialize

  if (!slepc_running) {
    start_slepc();
  }

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

  reset_matrices();

  // A
  CallPetsc(MatCreateShell(GLOBAL_COMMUNICATOR, getN(n), getN(n), getN(n), getN(n), (void *) op_a, &mat_A));
  CallPetsc(MatShellSetOperation(mat_A, MATOP_MULT, (void(*)(void)) &matMult));
  CallPetsc(MatShellSetOperation(mat_A, MATOP_MULT_TRANSPOSE,(void(*)(void)) &matMultTranspose));
  CallPetsc(MatShellSetOperation(mat_A,MATOP_GET_DIAGONAL,(void(*)(void)) &getDiagonal));

  // B
  CallPetsc(MatCreateShell(GLOBAL_COMMUNICATOR, getN(n), getN(n), getN(n), getN(n), (void *) op_b, &mat_B));
  CallPetsc(MatShellSetOperation(mat_B, MATOP_MULT, (void(*)(void)) &matMult));
  CallPetsc(MatShellSetOperation(mat_B, MATOP_MULT_TRANSPOSE,(void(*)(void)) &matMultTranspose));
  CallPetsc(MatShellSetOperation(mat_B,MATOP_GET_DIAGONAL,(void(*)(void)) &getDiagonal));

  // set options

  CallPetsc(EPSSetType(eps, getEpsType(*op_a, *op_b)));
  CallPetsc(EPSSetProblemType(eps, getProblemType(*op_a, *op_b)));

  CallPetsc(EPSSetOperators(eps, mat_A, mat_B));

  CallPetsc(EPSSetFromOptions(eps)); /* TODO really needed? */

  PetscFunctionReturnVoid();
}

bool SLEPcGEVPSolver::nextEigenpair(Skalar& eval, Vector& evec) {
  PetscFunctionBeginUser;

  if (!solved) {
    solve();
  }

  // ...
  PetscFunctionReturn(true);
}

void SLEPcGEVPSolver::reset_matrices() {
  PetscFunctionBeginUser;

  if (mat_A != nullptr) {
    MatDestroy(&mat_A);
  }
  if (mat_B != nullptr) {
    MatDestroy(&mat_B);
  }
  
  PetscFunctionReturnVoid();
}

void SLEPcGEVPSolver::end_slepc() {
  PetscFunctionBeginUser;

  reset_matrices();
  
  CallPetsc(EPSDestroy(&eps));

  //CallPetsc(SlepcFinalize());

  slepc_running = false;

  PetscFunctionReturnVoid();
}

SLEPcGEVPSolver::~SLEPcGEVPSolver() {
  if (slepc_running) {
    end_slepc();
  }
}

void SLEPcGEVPSolver::solve() {
  PetscFunctionBeginUser;

  // solve

  // CallPetsc(EPSSolve(eps));
  solved = true;

  PetscFunctionReturnVoid();
}

void SLEPcGEVPSolver::vector2Vec(Vector& vector, Vec vec) {
  /* TODO */
}

void SLEPcGEVPSolver::vec2vector(Vec vec, Vector& vector) {
  /* TODO */
}

PetscErrorCode SLEPcGEVPSolver::matMult(Mat M, Vec x, Vec y) {
  PetscFunctionBeginUser;

  Operator *op_m;
  getMatContext(M, &op_m);

  Vector sp_x(op_m->dimension());
  Vector sp_y(op_m->dimension());
  vec2vector(x, sp_x);
  op_m->apply(sp_x, sp_y);
  vector2Vec(sp_y, y);

  PetscFunctionReturn(0); /* TODO error code */
}

PetscErrorCode SLEPcGEVPSolver::matMultTranspose(Mat M, Vec x, Vec y) {
  PetscFunctionBeginUser;

  Operator *op_m;
  getMatContext(M, &op_m);

  if (op_m->isSymmetric()) {
    PetscFunctionReturn(matMult(M, x, y));
  } else {
    PetscFunctionReturn(-1); /* TODO */
  }
}

PetscErrorCode SLEPcGEVPSolver::getDiagonal(Mat M, Vec diag) {
  PetscFunctionBeginUser;
  
  Operator *op_m;
  getMatContext(M, &op_m);

  Vector d(op_m->dimension());
  op_m->getDiagonal(d);
  vector2Vec(d, diag);

  PetscFunctionReturn(0); /* TODO error code */
}

 void SLEPcGEVPSolver::getMatContext(Mat M, Operator **op) {
  PetscFunctionBeginUser;

  void *ctx;
  CallPetsc(MatShellGetContext(M, &ctx));
  *op = (Operator *) ctx;

  PetscFunctionReturnVoid();
}

}
}

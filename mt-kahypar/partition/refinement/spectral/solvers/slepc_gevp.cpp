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

#include "mt-kahypar/partition/refinement/spectral/datatypes.h"

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
  epairs_found = 0;

  size_t n = op_a->dimension();

  // create matrices

  reset_matrices();

  /* TODO eliminate duplicate code below 8S */
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
  CallPetsc(EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL));

  /* CallPetsc(EPSSetTolerances(eps, PETSC_DEFAULT, 500)); */

  CallPetsc(EPSSetOperators(eps, mat_A, mat_B));

  CallPetsc(EPSSetFromOptions(eps)); /* TODO really needed? */

  PetscFunctionReturnVoid();
}

void SLEPcGEVPSolver::setProblem(Operator& a, Operator& b, vec<Vector>& known_evecs, vec<Skalar> &known_evals, size_t deflation_epairs) {
  PetscFunctionBeginUser;

  setProblem(a, b);

  for (Vector& tevec : known_evecs) {
    Vec v;
    CallPetsc(VecCreateSeq(GLOBAL_COMMUNICATOR, tevec.dimension(), &v));
    vector2Vec(tevec, v);
    evecs.push_back(v);
  }
  evals.insert(evals.end(), known_evals.begin(), known_evals.end());

  num_deflation_epairs = deflation_epairs;

  PetscFunctionReturnVoid();
}

int SLEPcGEVPSolver::nextEigenpair(Skalar& eval, Vector& evec, bool try_from_above) {
  PetscFunctionBeginUser;

  // call solver
  solved = false;
  solve();
  if (!solved) {
    // try to search largest evec instead
    if (tried_from_above) {
      eval = evals.back();
      vec2vector(evecs.back(), evec);
      return 0;
    } else {
      tried_from_above = true;
      CallPetsc(EPSSetOperators(eps, mat_B, mat_A));
      CallPetsc(EPSSetWhichEigenpairs(eps, EPS_LARGEST_REAL));
      return nextEigenpair(eval, evec, try_from_above);
    }
  }

  PetscScalar kr, ki;
  Vec xr, xi;

  // debug output
  /* TODO opt off for performance reasons? */
  CallPetsc(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO_DETAIL));
  CallPetsc(EPSErrorView(eps, EPS_ERROR_RELATIVE, PETSC_VIEWER_STDOUT_WORLD));
  CallPetsc(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));

  // get eigenpair
  CallPetsc(VecCreateSeq(GLOBAL_COMMUNICATOR, evec.dimension(), &xr));
  CallPetsc(VecCreateSeq(GLOBAL_COMMUNICATOR, evec.dimension(), &xi));
  CallPetsc(EPSGetEigenpair(eps, 0, &kr, &ki, xr, xi));

  eval = tried_from_above ? 1.0 / kr : kr;
  vec2vector(xr, evec);
  
  evals.push_back(kr);
  evecs.push_back(xr);
  epairs_found++;

  if (tried_from_above) {
    tried_from_above = false;
    PetscFunctionReturn(-1);
  } else {
    PetscFunctionReturn(1);
  }
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

  /* CallPetsc(SlepcFinalize()); TODO execute this somehow */

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

  // deflation
  CallPetsc(EPSSetDeflationSpace(eps, num_deflation_epairs, evecs.data()));

  // initial vectors
  if (num_deflation_epairs + epairs_found < evecs.size()) {
    CallPetsc(EPSSetInitialSpace(eps, evecs.size() - num_deflation_epairs - epairs_found, evecs.data() + num_deflation_epairs));
  }

  // preconditioning
  /* TODO ??? */

  // solve
  PetscErrorCode err = EPSSolve(eps);
  if (err > 0) {
    solved = false;
  } else {
    PetscInt nconv;
    CallPetsc(EPSGetConverged(eps, &nconv));
    solved = nconv > 0;
  }

  // debug output
  /* TODO opt off for performance reasons? */
  CallPetsc(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO_DETAIL));
  CallPetsc(EPSConvergedReasonView(eps, PETSC_VIEWER_STDOUT_WORLD));
  CallPetsc(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));

  PetscFunctionReturnVoid();
}

void SLEPcGEVPSolver::vector2Vec(Vector& vector, Vec vec) {
  PetscFunctionBeginUser;

  Range(vector.dimension());
  
  CallPetsc(VecSetValues(vec, vector.dimension(), range, vector.get_all(), INSERT_VALUES));
  CallPetsc(VecAssemblyBegin(vec));
  CallPetsc(VecAssemblyEnd(vec));
  
  PetscFunctionReturnVoid();
}

void SLEPcGEVPSolver::vec2vector(Vec slepc_vec, Vector& vector) {
  PetscFunctionBeginUser;

  Range(vector.dimension());

  CallPetsc(VecGetValues(slepc_vec, vector.dimension(), range, vector.get_all()));

  /* vector.setGetter([&](size_t i){
    PetscScalar value;
    CallPetsc(VecGetValues(slepc_vec, (PetscInt) 1, (const PetscInt *) &i, &value));
    return value;
  }); */
  
  PetscFunctionReturnVoid();
}

PetscErrorCode SLEPcGEVPSolver::matMult(Mat M, Vec x, Vec y) {
  PetscFunctionBeginUser;

  Operator *op_m;
  getMatContext(M, &op_m);

  // vector objects in spectral namespace
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
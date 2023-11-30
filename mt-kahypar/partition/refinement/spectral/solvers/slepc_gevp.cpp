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

namespace mt_kahypar {
namespace spectral {

SLEPcGEVPSolver::SLEPcGEVPSolver() {
  PetscFunctionBeginUser;

  CallPetsc(SlepcInitialize(0, nullptr, (char*)0, nullptr));
  CallPetsc(EPSCreate(PETSC_COMM_WORLD, &_context));
}

void SLEPcGEVPSolver::initialize(spectral::Matrix& a, spectral::Matrix& b) {
  _a = &a;
  _b = &b;
  _solved = false;
}

bool SLEPcGEVPSolver::nextEigenpair(spectral::Skalar& eval, spectral::Vector& evec) {
  PetscFunctionBeginUser;

  if (!_solved) {
    solve();
  }

  // ...
  return true;
}

SLEPcGEVPSolver::~SLEPcGEVPSolver() {
  PetscFunctionBeginUser;

  CallPetsc(EPSDestroy(&_context));
  CallPetsc(SlepcFinalize());
}

void SLEPcGEVPSolver::solve() {
  PetscFunctionBeginUser;

  CallPetsc(PetscPrintf(PETSC_COMM_WORLD, "\n\nhello world from SLEPc!\n\n"));

  // /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //     Create the operator matrix that defines the eigensystem, Ax=kx
  //   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  // CallPetsc(MatCreateShell(PETSC_COMM_WORLD,N,N,N,N,&n,&A));
  // CallPetsc(MatShellSetOperation(A,MATOP_MULT,(void(*)(void))MatMult_Laplacian2D));
  // CallPetsc(MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_Laplacian2D));
  // CallPetsc(MatShellSetOperation(A,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_Laplacian2D));

  // /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //               Create the eigensolver and set various options
  //   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  // /*
  //   Create eigensolver context
  // */
  // CallPetsc(EPSCreate(PETSC_COMM_WORLD,&eps));

  // /*
  //   Set operators. In this case, it is a standard eigenvalue problem
  // */
  // CallPetsc(EPSSetOperators(eps,A,NULL));
  // CallPetsc(EPSSetProblemType(eps,EPS_HEP));

  // /*
  //   Set solver parameters at runtime
  // */
  // CallPetsc(EPSSetFromOptions(eps));

  // /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //                     Solve the eigensystem
  //   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  // CallPetsc(EPSSolve(eps));

  _solved = true;
}

}
}

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

#pragma once

#include <cstddef>

#include "mt-kahypar/partition/refinement/spectral/algebraic_wrappers/vector.h"
#include "mt-kahypar/partition/refinement/spectral/algebraic_wrappers/operator.h"
#include "mt-kahypar/partition/refinement/spectral/datatypes.h"

namespace mt_kahypar {
namespace spectral {

class GEVPSolver {
 public:
  virtual void setProblem(Operator& a, Operator& b) = 0;

  virtual void setProblem(Operator& a, Operator& b, vec<Vector>& trivial_evecs, vec<Skalar> &trivial_evals) = 0;

  // ! positive (negative) return value: ascending (descending), zero: no new pair found, returning last pair
  virtual int nextEigenpair(Skalar& eval, Vector& evec) = 0; /* TODO return rather all pairs, an iterator or YAGNI? */

  virtual ~GEVPSolver() = default;
};

}
}

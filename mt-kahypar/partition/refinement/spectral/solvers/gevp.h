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
#include "mt-kahypar/partition/refinement/spectral/algebraic_wrappers/matrix.h"
#include "mt-kahypar/partition/refinement/spectral/datatypes.h"

namespace mt_kahypar {
namespace spectral {

class GEVPSolver {
 public:
  virtual void initialize(spectral::Matrix& a, spectral::Matrix& b) = 0;

  virtual bool nextEigenpair(spectral::Skalar& eval, spectral::Vector& evec) = 0; /* TODO return rather all pairs, an iterator or YAGNI? */

  virtual ~GEVPSolver() = default;
};

}
}

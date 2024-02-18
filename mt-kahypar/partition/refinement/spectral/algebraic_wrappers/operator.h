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

#include "mt-kahypar/partition/refinement/spectral/datatypes.h"
#include "mt-kahypar/partition/refinement/spectral/algebraic_wrappers/vector.h"


namespace mt_kahypar {
namespace spectral {

class Operator {
 public:
  Operator(size_t dimension);

  size_t dimension();

  void apply(Vector& operand, Vector& target);

  void getDiagonal(Vector& target);

  void getMatrix(vec<Vector> &target);

  bool isSymmetric();

  /* only makeshift */
  void (*effect) (Operator*, Vector&, Vector&) = [](Operator *self, Vector& operand, Vector& target) { throw "unassigned"; };
  void (*calc_diagonal) (Operator*, Vector&) = [](Operator *self, Vector& target) { throw "unassigned"; };
  void *ctx;

 private:
  size_t dim;
};

Operator operator+ (Operator& a, Operator& b);

}
}

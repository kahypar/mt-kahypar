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

#include "mt-kahypar/partition/refinement/spectral/algebraic_wrappers/operator.h"

namespace mt_kahypar {
namespace spectral {

/* make getters more elegant? */

Operator::Operator(size_t dimension) {
  dim = dimension;
  // initialize with id
  effects.push_back([](Operator *self, Vector& operand, Vector& target_vector) {
    for (size_t i = 0; i < operand.dimension(); i++) {
      target_vector.set(i, operand[i]);
    }
  });
  calc_diagonal_ops.push_back([] (Operator *self, Vector& target_vector) {
    for (size_t i = 0; i < target_vector.dimension(); i++) {
      target_vector.set(i, 1.0);
    }
  });
}

size_t Operator::dimension() {
  return dim;
}

void Operator::apply(Vector& operand, Vector& target) {
  for (auto effect: effects) {
    effect(this, operand, target);
  }
}

void Operator::getDiagonal(Vector& target) {
  for (auto op: calc_diagonal_ops) {
    op(this, target);
  }
}

void Operator::getMatrix(vec<Vector> &target) {
  applyOnMatrix([&](Vector &col) { target.push_back(col); });
}

template <typename P>
void Operator::printMatrix(P printer) {
  applyOnMatrix([&] (Vector &row) {
    std::ostringstream row_str;
    for (size_t i = 0; i < dimension(); i++) {
      char buf[100];
      sprintf(buf, " %+.2f", row[i]);
      row_str << buf;
    }
    printer(row_str.str());
  });
}

template <typename F>
void Operator::applyOnMatrix(F fun) {
  Vector filter(dim);
  for (size_t i = 0; i < dim; i++) {
    filter.set(i, 1);
    Vector column(dim);
    apply(filter, column);
    fun(column);
    filter.set(i, 0);
  }
}

bool Operator::isSymmetric() {
  return true; /* TODO */
}

Operator& Operator::operator+=(Operator& op) { /* TODO not that clean algebraically nor softwaretechnically */
  effects.insert(effects.end(), op.effects.begin(), op.effects.end());
  return *this;
}

}
}

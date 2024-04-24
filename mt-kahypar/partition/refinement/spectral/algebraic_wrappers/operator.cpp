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
  effects.push_back([](void *ctx, Vector& operand, Vector& target_vector) {
    unused(ctx);
    for (size_t i = 0; i < operand.dimension(); i++) {
      target_vector.set(i, operand[i]);
    }
  });
  calc_diagonal_ops.push_back([] (void *ctx, Vector& target_vector) {
    unused(ctx);
    for (size_t i = 0; i < target_vector.dimension(); i++) {
      target_vector.set(i, 1.0);
    }
  });
  ctx_exporter.push_back([] (void *ctx, vec<size_t> &target_vector) {
    unused(ctx);
    unused(target_vector);
  });
  ctx.push_back(nullptr);
}

size_t Operator::dimension() {
  return dim;
}

void Operator::apply(Vector& operand, Vector& target) {
  for (size_t i = 0; i < effects.size(); i++) {
    effects[i](ctx[i], operand, target);
  }
}

void Operator::getDiagonal(Vector& target) {
  for (size_t i = 0; i < calc_diagonal_ops.size(); i++) {
    calc_diagonal_ops[i](ctx[i], target);
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

void Operator::exportContext(size_t index, vec<size_t> &target) {
  ctx_exporter[index](ctx[index], target);
}

Operator& Operator::operator+=(Operator& op) { /* TODO not that clean algebraically nor softwaretechnically */
  effects.insert(effects.end(), op.effects.begin(), op.effects.end());
  calc_diagonal_ops.insert(calc_diagonal_ops.end(), op.calc_diagonal_ops.begin(), op.calc_diagonal_ops.end());
  ctx.insert(ctx.end(), op.ctx.begin(), op.ctx.end());
  return *this;
}

}
}

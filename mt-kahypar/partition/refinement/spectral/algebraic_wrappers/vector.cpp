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

#include "mt-kahypar/partition/refinement/spectral/algebraic_wrappers/vector.h"

namespace mt_kahypar {
namespace spectral {

/* make getters more elegant? */

Vector::Vector(size_t dimension, Skalar default_value) {
  _default_value = default_value;
  dim = dimension;
  data.reserve(dim);
}

size_t Vector::dimension() {
  return dim;
}

Skalar Vector::get(size_t index) {
  if (use_getter) {
    return custom_getter(index);
  }

  if (index >= data.size()) {
    data.resize(index + 1, _default_value);
  }
  return data[index];
}

template <typename F>
void Vector::set_generalized(size_t index, F &&set_op) {
  disableGetter();
  if (index >= data.size()) { /* TODO copied code 8S */
    data.resize(index + 1, _default_value);
  }
  set_op();
}

void Vector::set(size_t index, const Skalar &&value) {
  set_generalized(index, [&](){data[index] = value;});
}

void Vector::set(size_t index, const Skalar &value) {
  set_generalized(index, [&](){data[index] = value;});
}


const Skalar *Vector::get_all() {
  if (use_getter) {
    for (size_t i = 0; i < dimension(); i++) {
      set(i, get(i));
    }
  }
  get(dimension() - 1);
  return data.data();
}

void Vector::set_all(const Skalar *data_array) {
  for (size_t i = 0; i < dimension(); i++) {
    set(i, data_array[i]);
  }
  
}

void Vector::setGetter(Skalar (*getter) (size_t)) {
  custom_getter = getter;
  use_getter = true;
}

void Vector::disableGetter() {
  use_getter = false;
}

Vector& Vector::operator+=(Vector& v) {
  for (size_t i = 0; i < dimension(); i++) {
    set(i, get(i) + v[i]);
  }
  return *this;
}

Skalar Vector::operator[](size_t index) {
  return get(index);
}

}
}

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


namespace mt_kahypar {
namespace spectral {

class Vector {
 public:
  Vector(size_t dimension) : Vector(dimension, 0.0) {};
  Vector(size_t dimension, Skalar default_value);
  size_t dimension();

  Skalar get(size_t index);
  void set(size_t index, const Skalar &&value);
  void set(size_t index, const Skalar &value);

  template <typename F>
  void setGetter(F getter);
  
  void disableGetter();

  Skalar operator[](size_t index);

  Skalar *get_all();
  void set_all(const Skalar *data_array);

  Skalar get_default();
  void set_default(Skalar val);

  void reset();

  Vector& operator+=(Vector& v);
 private:
  size_t dim;
  Skalar _default_value;
  bool use_getter = false;
  std::function<Skalar(size_t)> custom_getter;
  vec<Skalar> data;

  template <typename F>
  void set_generalized(size_t index, F &&set_op);
};
}
}

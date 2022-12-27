/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <array>
#include <string>
#include <utility>
#include <vector>


#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {

class IHypergraphSparsifier {

 public:
  IHypergraphSparsifier(const IHypergraphSparsifier&) = delete;
  IHypergraphSparsifier(IHypergraphSparsifier&&) = delete;
  virtual ~IHypergraphSparsifier() = default;
  IHypergraphSparsifier & operator= (const IHypergraphSparsifier &) = delete;
  IHypergraphSparsifier & operator= (IHypergraphSparsifier &&) = delete;

  bool isSparsified() const {
    return _is_sparsified;
  }

  Hypergraph& sparsifiedHypergraph() {
    ASSERT(_is_sparsified);
    return sparsifiedHypergraphImpl();
  }

  PartitionedHypergraph& sparsifiedPartitionedHypergraph() {
    ASSERT(_is_sparsified);
    return sparsifiedPartitionedHypergraphImpl();
  }

  void sparsify(const Hypergraph& hypergraph) {
    ASSERT(!_is_sparsified);
    sparsifyImpl(hypergraph);
    _is_sparsified = true;
  }

  void undoSparsification(PartitionedHypergraph& hypergraph) {
    ASSERT(_is_sparsified);
    undoSparsificationImpl(hypergraph);
  }

 protected:
  IHypergraphSparsifier() :
    _is_sparsified(false) { }

 private:
  virtual Hypergraph& sparsifiedHypergraphImpl() = 0;
  virtual PartitionedHypergraph& sparsifiedPartitionedHypergraphImpl() = 0;
  virtual void sparsifyImpl(const Hypergraph& hypergraph) = 0;
  virtual void undoSparsificationImpl(PartitionedHypergraph& hypergraph) = 0;

  bool _is_sparsified;
};

}  // namespace mt_kahypar

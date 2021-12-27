/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "kahypar/partition/metrics.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {

class IHypergraphSparsifier {

 public:
  IHypergraphSparsifier(const IHypergraphSparsifier&) = delete;
  IHypergraphSparsifier(IHypergraphSparsifier&&) = delete;
  IHypergraphSparsifier & operator= (const IHypergraphSparsifier &) = delete;
  IHypergraphSparsifier & operator= (IHypergraphSparsifier &&) = delete;

  virtual ~IHypergraphSparsifier() = default;

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

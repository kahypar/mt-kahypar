/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
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
template<typename TypeTraits>
class IHypergraphSparsifierT {
  using HyperGraph = typename TypeTraits::HyperGraph;
  using PartitionedHyperGraph = typename TypeTraits::PartitionedHyperGraph;

 public:
  IHypergraphSparsifierT(const IHypergraphSparsifierT&) = delete;
  IHypergraphSparsifierT(IHypergraphSparsifierT&&) = delete;
  IHypergraphSparsifierT & operator= (const IHypergraphSparsifierT &) = delete;
  IHypergraphSparsifierT & operator= (IHypergraphSparsifierT &&) = delete;

  virtual ~IHypergraphSparsifierT() = default;

  bool isSparsified() const {
    return _is_sparsified;
  }

  HyperGraph& sparsifiedHypergraph() {
    ASSERT(_is_sparsified);
    return sparsifiedHypergraphImpl();
  }

  PartitionedHyperGraph& sparsifiedPartitionedHypergraph() {
    ASSERT(_is_sparsified);
    return sparsifiedPartitionedHypergraphImpl();
  }

  void sparsify(const HyperGraph& hypergraph) {
    ASSERT(!_is_sparsified);
    sparsifyImpl(hypergraph);
    _is_sparsified = true;
  }

  void undoSparsification(PartitionedHyperGraph& hypergraph) {
    ASSERT(_is_sparsified);
    undoSparsificationImpl(hypergraph);
  }

 protected:
  IHypergraphSparsifierT() :
    _is_sparsified(false) { }

 private:
  virtual HyperGraph& sparsifiedHypergraphImpl() = 0;
  virtual PartitionedHyperGraph& sparsifiedPartitionedHypergraphImpl() = 0;
  virtual void sparsifyImpl(const HyperGraph& hypergraph) = 0;
  virtual void undoSparsificationImpl(PartitionedHyperGraph& hypergraph) = 0;

  bool _is_sparsified;
};

using IHypergraphSparsifier = IHypergraphSparsifierT<GlobalTypeTraits>;

}  // namespace mt_kahypar

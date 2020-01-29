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

#include <string>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"

namespace mt_kahypar {
template<typename TypeTraits>
class ICoarsenerT {

  using HyperGraph = typename TypeTraits::HyperGraph;
  using Refiner = IRefinerT<TypeTraits>;

 public:
  ICoarsenerT(const ICoarsenerT&) = delete;
  ICoarsenerT(ICoarsenerT&&) = delete;
  ICoarsenerT & operator= (const ICoarsenerT &) = delete;
  ICoarsenerT & operator= (ICoarsenerT &&) = delete;

  void coarsen() {
    coarsenImpl();
  }

  bool uncoarsen(std::unique_ptr<Refiner>& label_propagation) {
    return uncoarsenImpl(label_propagation);
  }

  HyperGraph& coarsestHypergraph() {
    return coarsestHypergraphImpl();
  }

  virtual ~ICoarsenerT() = default;

 protected:
  ICoarsenerT() = default;

 private:
  virtual void coarsenImpl() = 0;
  virtual bool uncoarsenImpl(std::unique_ptr<Refiner>& label_propagation) = 0;
  virtual HyperGraph& coarsestHypergraphImpl() = 0;
};

using ICoarsener = ICoarsenerT<GlobalTypeTraits>;
}  // namespace kahypar

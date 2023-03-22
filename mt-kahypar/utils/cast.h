/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <string>

#include "include/libmtkahypartypes.h"

namespace mt_kahypar::utils {

namespace {

std::string typeToString(const mt_kahypar_hypergraph_type_t type) {
  switch ( type ) {
    case STATIC_GRAPH: return "STATIC_GRAPH";
    case DYNAMIC_GRAPH: return "DYNAMIC_GRAPH";
    case STATIC_HYPERGRAPH: return "STATIC_HYPERGRAPH";
    case DYNAMIC_HYPERGRAPH: return "DYNAMIC_HYPERGRAPH";
    case NULLPTR: return "NULLPTR";
  }
  return "UNDEFINED";
}

} // namespace

template<typename Hypergraph>
Hypergraph& cast(mt_kahypar_hypergraph_t hypergraph) {
  if ( Hypergraph::TYPE != hypergraph.type ) {
    ERR("Cannot cast" << typeToString(hypergraph.type) << "to" << typeToString(Hypergraph::TYPE));
  }
  return *reinterpret_cast<Hypergraph*>(hypergraph.hypergraph);
}

template<typename Hypergraph>
const Hypergraph& cast_const(const mt_kahypar_hypergraph_t hypergraph) {
  if ( Hypergraph::TYPE != hypergraph.type ) {
    ERR("Cannot cast" << typeToString(hypergraph.type) << "to" << typeToString(Hypergraph::TYPE));
  }
  return *reinterpret_cast<const Hypergraph*>(hypergraph.hypergraph);
}

}  // namespace mt_kahypar
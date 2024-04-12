/*******************************************************************************
 * MIT License
 *
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include <vector>

#include "include/libmtkahypartypes.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/flows/deterministic/deterministic_quotient_graph.h"

namespace mt_kahypar {

template<typename TypeTraits>
class IDeterministicBlockSchedule {

public:
    IDeterministicBlockSchedule(const IDeterministicBlockSchedule&) = delete;
    IDeterministicBlockSchedule(IDeterministicBlockSchedule&&) = delete;
    IDeterministicBlockSchedule& operator= (const IDeterministicBlockSchedule&) = delete;
    IDeterministicBlockSchedule& operator= (IDeterministicBlockSchedule&&) = delete;

    virtual ~IDeterministicBlockSchedule() = default;

    void initialize(mt_kahypar_partitioned_hypergraph_t& hypergraph, const DeterministicQuotientGraph<TypeTraits>& qg) {
        initializeImpl(hypergraph, qg);
    }

    vec<BlockPair> getNextMatching(const DeterministicQuotientGraph<TypeTraits>& qg) {
        return getNextMatchingImpl(qg);
    }

protected:
    IDeterministicBlockSchedule() = default;

private:
    virtual void initializeImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph, const DeterministicQuotientGraph<TypeTraits>& qg) = 0;

    virtual vec<BlockPair> getNextMatchingImpl(const DeterministicQuotientGraph<TypeTraits>& qg) = 0;
};

}  // namespace mt_kahypar

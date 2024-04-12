/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "tbb/concurrent_queue.h"
#include "tbb/concurrent_vector.h"
#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/flows/refiner_adapter.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/randomize.h"
#include <tbb/parallel_sort.h>

namespace mt_kahypar {

class DeterministicQuotientGraphEdge {
public:
    DeterministicQuotientGraphEdge() : cut_hyperedges(), num_cut_hyperedges(0UL), cut_hyperedge_weight(0), total_improvement(0) {}

    void addHyperedge(const HyperedgeID he, const HyperedgeWeight weight) {
        cut_hyperedges.push_back(he);
        cut_hyperedge_weight += weight;
        num_cut_hyperedges++;
    }

    void reset() {
        cut_hyperedges.clear();
        cut_hyperedge_weight = 0;
        num_cut_hyperedges = 0;
    }
    vec<HyperedgeID>& getCutEdges() {
        return cut_hyperedges;
    }

    vec<HyperedgeID> cut_hyperedges;
    size_t num_cut_hyperedges;
    HyperedgeWeight cut_hyperedge_weight;
    HyperedgeWeight total_improvement;
};

template<typename TypeTraits>
class DeterministicQuotientGraph {

    static constexpr bool debug = false;
    static constexpr bool enable_heavy_assert = false;

    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

public:

    explicit DeterministicQuotientGraph(const Context& context) :
        _edges(context.partition.k, vec<DeterministicQuotientGraphEdge>(context.partition.k)),
        _k(context.partition.k),
        _dirty(false) {}

    DeterministicQuotientGraph(const DeterministicQuotientGraph&) = delete;
    DeterministicQuotientGraph(DeterministicQuotientGraph&&) = delete;

    DeterministicQuotientGraph& operator= (const DeterministicQuotientGraph&) = delete;
    DeterministicQuotientGraph& operator= (DeterministicQuotientGraph&&) = delete;

    void initialize(const PartitionedHypergraph& phg) {
        if (_dirty) {
            reset();
        }
        // TOOD: Parallelize? We might have to sort edges though
        for (const HyperedgeID he : phg.edges()) {
            for (const PartitionID i : phg.connectivitySet(he)) {
                for (const PartitionID j : phg.connectivitySet(he)) {
                    if (i < j) {
                        _edges[i][j].addHyperedge(he, phg.edgeWeight(he));
                    }
                }
            }
        }
        _dirty = true;
    }

    void reset() {
        for (PartitionID i = 0; i < _k; ++i) {
            for (PartitionID j = i + 1; j < _k; ++j) {
                _edges[i][j].reset();
            }
        }
    }

    DeterministicQuotientGraphEdge& getEdge(const PartitionedHypergraph& phg, const PartitionID a, const PartitionID b) {
        const PartitionID i = std::min(a, b);
        const PartitionID j = std::max(a, b);
        DeterministicQuotientGraphEdge& quotientEdge = _edges[i][j];
        vec<HyperedgeID>& edges = quotientEdge.cut_hyperedges;
        for (size_t i = 0; i < edges.size(); ++i) {
            const HyperedgeID he = edges[i];
            if (phg.pinCountInPart(he, i) == 0 || phg.pinCountInPart(he, j) == 0) {
                edges[i] = edges.back();
                edges.pop_back();
                quotientEdge.num_cut_hyperedges--;
                quotientEdge.cut_hyperedge_weight -= phg.edgeWeight(he);
                --i;
            }
        }
        // For determinism
        tbb::parallel_sort(edges.begin(), edges.end());
        return quotientEdge;
    }

private:
    vec<vec<DeterministicQuotientGraphEdge>> _edges;
    PartitionID _k;
    bool _dirty;
};

}  // namespace kahypar

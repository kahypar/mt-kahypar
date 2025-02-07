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
#include <tbb/parallel_sort.h>

namespace mt_kahypar {

class DeterministicQuotientGraphEdge {
public:
    DeterministicQuotientGraphEdge() : cut_hyperedges(), num_cut_hyperedges(0UL), cut_hyperedge_weight(0), total_improvement(0), previously_scheduled(false) {}

    void addHyperedge(const HyperedgeID he, const HyperedgeWeight weight) {
        cut_hyperedges.push_back(he);
        cut_hyperedge_weight += weight;
        num_cut_hyperedges++;
    }

    void reset() {
        cut_hyperedges.clear();
        cut_hyperedge_weight.store(0, std::memory_order_relaxed);
        num_cut_hyperedges.store(0, std::memory_order_relaxed);
        previously_scheduled.store(false, std::memory_order_relaxed);
    }

    tbb::concurrent_vector<HyperedgeID> cut_hyperedges; // reset
    CAtomic<size_t> num_cut_hyperedges; // reset
    CAtomic<HyperedgeWeight> cut_hyperedge_weight; // reset
    CAtomic<HyperedgeWeight> total_improvement; // NOT reset but ok
    CAtomic<bool> previously_scheduled;
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

    DeterministicQuotientGraphEdge& getEdgeFiltered(const PartitionedHypergraph& phg, const PartitionID i, const PartitionID j) {
        assert(i < j);
        DeterministicQuotientGraphEdge& quotientEdge = _edges[i][j];
        vec<HyperedgeID>& edges = quotientEdge.cut_hyperedges;
        for (size_t idx = 0; idx < edges.size(); ++idx) {
            const HyperedgeID he = edges[idx];
            if constexpr (debug) {
                size_t count0 = 0;
                size_t count1 = 0;
                for (const HypernodeID pin : phg.pins(he)) {
                    if (phg.partID(pin) == i) count0++;
                    if (phg.partID(pin) == j) count1++;
                }
                assert(count0 == phg.pinCountInPart(he, i));
                assert(count1 == phg.pinCountInPart(he, j));
            }
            if (phg.pinCountInPart(he, i) == 0 || phg.pinCountInPart(he, j) == 0) {
                edges[idx] = edges.back();
                edges.pop_back();
                quotientEdge.num_cut_hyperedges--;
                quotientEdge.cut_hyperedge_weight -= phg.edgeWeight(he);
                --idx;
            }
        }
        // For determinism
        tbb::parallel_sort(edges.begin(), edges.end());
        return quotientEdge;
    }
    template<typename F>
    void doForAllCutHyperedgesOfPair(const PartitionedHypergraph& phg, const PartitionID i, const PartitionID j, const F& f) {
        auto& edges = _edges[i][j].cut_hyperedges;
        const size_t num_cut_hes = _edges[i][j].num_cut_hyperedges;
        tbb::parallel_sort(edges.begin(), edges.begin() + num_cut_hes);
        HyperedgeID prev_he = kInvalidHyperedge;
        for (size_t k = 0; k < num_cut_hes; ++k) {
            const HyperedgeID he = edges[k];
            if (phg.pinCountInPart(he, i) > 0 && phg.pinCountInPart(he, j) > 0 && he != prev_he) {
                f(he);
            }
            prev_he = he;
        }
    }

    HyperedgeWeight getCutWeight(const PartitionID i, const PartitionID j)const {
        assert(i < j);
        return _edges[i][j].cut_hyperedge_weight;
    }

    HyperedgeWeight getImprovement(const PartitionID i, const PartitionID j)const {
        assert(i < j);
        return _edges[i][j].total_improvement;
    }

    bool wasAlreadyScheduled(const PartitionID i, const PartitionID j)const {
        assert(i < j);
        return _edges[i][j].previously_scheduled;
    }

    void addNewCutHyperedge(const HyperedgeID he, const PartitionID block, const PartitionedHypergraph& phg) {
        ASSERT(phg.pinCountInPart(he, block) > 0);
        // Add hyperedge he as a cut hyperedge to each block pair that contains 'block'
        for (const PartitionID& other_block : phg.connectivitySet(he)) {
            if (other_block != block) {
                _edges[std::min(block, other_block)][std::max(block, other_block)].addHyperedge(he, phg.edgeWeight(he));
            }
        }
    }

    void reportImprovement(const PartitionID i, const PartitionID j, const HyperedgeWeight improvement) {
        _edges[i][j].total_improvement += improvement;
        _edges[i][j].previously_scheduled.store(true);
    }

private:
    vec<vec<DeterministicQuotientGraphEdge>> _edges;
    PartitionID _k;
    bool _dirty;
};

}  // namespace kahypar

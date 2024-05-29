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

#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/stl/scalable_queue.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/partition/refinement/flows/deterministic/deterministic_quotient_graph.h"


namespace mt_kahypar {

template<typename TypeTraits>
class DeterministicProblemConstruction {

    static constexpr bool debug = false;

    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

    /**
     * Contains data required to grow two region around
     * the cut of two blocks of the partition.
     */
    struct BFSData {
        explicit BFSData(const HypernodeID num_nodes,
            const HyperedgeID num_edges,
            const PartitionID k) :
            current_distance(0),
            queue(),
            next_queue(),
            visited_hn(num_nodes, false),
            visited_he(num_edges, false),
            contained_hes(num_edges, false),
            locked_blocks(k, false),
            queue_weight_block_0(0),
            queue_weight_block_1(0),
            lock_queue(false) {}

        void clearQueue();

        void reset();

        HypernodeID pop_hypernode();

        void add_pins_of_hyperedge_to_queue(const HyperedgeID& he,
            const PartitionedHypergraph& phg,
            const size_t max_bfs_distance,
            const HypernodeWeight max_weight_block_0,
            const HypernodeWeight max_weight_block_1);

        bool is_empty() const {
            return queue.empty();
        }

        bool is_next_empty() const {
            return next_queue.empty();
        }

        void swap_with_next_queue() {
            if (!is_next_empty()) {
                std::swap(queue, next_queue);
                ++current_distance;
            }
        }

        PartitionID block0; // reset
        PartitionID block1; // reset
        size_t current_distance; // reset
        parallel::scalable_queue<HypernodeID> queue; // reset
        parallel::scalable_queue<HypernodeID> next_queue; // reset
        vec<bool> visited_hn; // reset
        vec<bool> visited_he; // reset
        vec<bool> contained_hes; // reset
        vec<bool> locked_blocks; // reset
        HypernodeWeight queue_weight_block_0; // reset
        HypernodeWeight queue_weight_block_1; // reset
        bool lock_queue; // reset
    };

public:
    explicit DeterministicProblemConstruction(const HypernodeID num_hypernodes,
        const HyperedgeID num_hyperedges,
        const Context& context) :
        _context(context),
        _scaling(1.0 + _context.refinement.flows.alpha *
            std::min(0.05, _context.partition.epsilon)),
        _num_hypernodes(num_hypernodes),
        _num_hyperedges(num_hyperedges),
        _bfs(_num_hypernodes, _num_hyperedges, _context.partition.k) {}

    DeterministicProblemConstruction(const DeterministicProblemConstruction&) = delete;
    DeterministicProblemConstruction(DeterministicProblemConstruction&&) = delete;

    DeterministicProblemConstruction& operator= (const DeterministicProblemConstruction&) = delete;
    DeterministicProblemConstruction& operator= (DeterministicProblemConstruction&&) = delete;

    Subhypergraph construct(const PartitionedHypergraph& phg,
        DeterministicQuotientGraph<TypeTraits>& quotient_graph,
        const PartitionID block0,
        const PartitionID block1);

private:

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool isMaximumProblemSizeReached(
        const Subhypergraph& sub_hg,
        const HypernodeWeight max_weight_block_0,
        const HypernodeWeight max_weight_block_1,
        vec<bool>& locked_blocks) const;

    const Context& _context;
    double _scaling; // constatn
    HypernodeID _num_hypernodes; // nur inital einmal verwendet
    HyperedgeID _num_hyperedges; // nur inistial einmal verwerndet

    // ! Contains data required for BFS construction algorithm
    BFSData _bfs; // reset
};

}  // namespace kahypar
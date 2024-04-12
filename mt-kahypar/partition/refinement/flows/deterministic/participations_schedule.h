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
#include "mt-kahypar/partition/refinement/flows/deterministic/i_deterministic_block_schedule.h"

namespace mt_kahypar {

template<typename TypeTraits>
class ParticipationsSchedule final : public IDeterministicBlockSchedule<TypeTraits> {


    ParticipationsSchedule(const ParticipationsSchedule&) = delete;
    ParticipationsSchedule(ParticipationsSchedule&&) = delete;
    ParticipationsSchedule& operator= (const ParticipationsSchedule&) = delete;
    ParticipationsSchedule& operator= (ParticipationsSchedule&&) = delete;
public:
    explicit ParticipationsSchedule(const Context& context) :
        _active_blocks(context.partition.k, true),
        _active_blocks_next_round(context.partition.k, false),
        _participations(context.partition.k, 0),
        _processed(context.partition.k, vec<uint8_t>(context.partition.k, false)),
        _partitions_sorted_by_participations(),
        _scheduled(context.partition.k, false),
        _round(0),
        _k(context.partition.k) {}

    void resetForNewRound(const DeterministicQuotientGraph<TypeTraits>& qg) {
        _scheduled.assign(_scheduled.size(), false);
        _active_blocks.assign(_active_blocks.size(), true);
        _active_blocks_next_round.assign(_active_blocks_next_round.size(), false);
        _participations.assign(_participations.size(), 0);
        for (auto& v : _processed) {
            v.assign(v.size(), false);
        }
        for (PartitionID i = 0; i < _k - 1; ++i) {
            for (PartitionID j = i + 1; j < _k; ++j) {
                if (isEligible(i, j, qg)) {
                    _participations[i]++;
                    _participations[j]++;
                }
            }
        }
        _partitions_sorted_by_participations.reserve(_k);
        for (PartitionID i = 0; i < _k; ++i) {
            if (_participations[i] > 0) {
                _partitions_sorted_by_participations.push_back(i);
            }
        }
        std::sort(_partitions_sorted_by_participations.begin(), _partitions_sorted_by_participations.end(), [&](const PartitionID i, const PartitionID j) {
            return _participations[i] > _participations[j] || (_participations[i] == _participations[j] && i < j);
        });
        _round++;
    }


    void reportResults(const PartitionID block0, const PartitionID block1, const MoveSequence& sequence) {
        if (sequence.moves.size() > 0) {
            // There was improvement
            _active_blocks_next_round[block0] = true;
            _active_blocks_next_round[block1] = true;
        }
    }

private:
    void initializeImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph, const DeterministicQuotientGraph<TypeTraits>& qg) {
        unused(hypergraph);
        resetForNewRound(qg);
        _round = 0;
    }

    vec<BlockPair> getNextMatchingImpl(const DeterministicQuotientGraph<TypeTraits>& qg) {
        vec<BlockPair> tasks;
        tasks.reserve(_k / 2);
        for (size_t i = 0; i < _partitions_sorted_by_participations.size() - 1; ++i) {
            const PartitionID block0 = _partitions_sorted_by_participations[i];
            if (_scheduled[block0]) continue;
            for (size_t j = i + 1; j < _partitions_sorted_by_participations.size(); ++j) {
                const PartitionID block1 = _partitions_sorted_by_participations[j];
                if (_scheduled[block1] || _processed[block0][block1]) continue;

                if (isEligible(block0, block1, qg)) {
                    addBlockPair(block0, block1);
                    tasks.push_back({ block0, block1 });
                    --i;
                    break;
                }
            }
            if (tasks.size() == size_t(_k) / 2) { break; }
        }
        return tasks;
    }

    void addBlockPair(const PartitionID i, const PartitionID j) {
        _processed[i][j] = true;
        _scheduled[i] = true;
        _scheduled[j] = true;
        _participations[i]--;
        fixOrdering(i);
        _participations[j]--;
        fixOrdering(j);
    }

    void fixOrdering(const PartitionID id) {
        for (size_t i = 0; i < _partitions_sorted_by_participations.size(); ++i) {
            if (_partitions_sorted_by_participations[i] == id) {
                size_t swapPosition = i + 1;
                while (swapPosition < _partitions_sorted_by_participations.size() && _participations[_partitions_sorted_by_participations[swapPosition]] > _participations[_partitions_sorted_by_participations[i]]) { swapPosition++; }
                swapPosition--;
                assert(swapPosition < _partitions_sorted_by_participations.size());
                std::swap(_partitions_sorted_by_participations[i], _partitions_sorted_by_participations[swapPosition]);
                break;
            }
        }
        if (_participations[_partitions_sorted_by_participations[_partitions_sorted_by_participations.size() - 1]] == 0) { _partitions_sorted_by_participations.pop_back(); }
    }

    bool isEligible(const PartitionID i, const PartitionID j, const DeterministicQuotientGraph<TypeTraits>& qg) {
        assert(i < j);
        const HyperedgeWeight weight = qg.getCutWeight(i, j);
        const HyperedgeWeight improvement = qg.getImprovement(i, j);
        return weight > 0                               // cut between blocks
            && (_active_blocks[i] || _active_blocks[j]) // at least one block active
            && !_scheduled[i] && !_scheduled[j]         // none of the blocks is already scheduled
            && !_processed[i][j]                        // the pairing has not been processed yet
            && (_round < 2 || improvement > 0);
    }

    vec<uint8_t> _active_blocks;
    vec<uint8_t> _active_blocks_next_round;
    vec<size_t> _participations;
    vec<vec<uint8_t>> _processed;
    vec<PartitionID> _partitions_sorted_by_participations;
    vec<bool> _scheduled;
    size_t _round;
    PartitionID _k;
};

}  // namespace mt_kahypar

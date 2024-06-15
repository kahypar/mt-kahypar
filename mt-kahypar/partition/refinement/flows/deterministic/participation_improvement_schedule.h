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
#include "tbb/parallel_for.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/flows/deterministic/i_deterministic_block_schedule.h"

namespace mt_kahypar {

template<typename TypeTraits>
class ParticipationImprovementSchedule final : public IDeterministicBlockSchedule<TypeTraits> {
    static constexpr bool debug = false;

    ParticipationImprovementSchedule(const ParticipationImprovementSchedule&) = delete;
    ParticipationImprovementSchedule(ParticipationImprovementSchedule&&) = delete;
    ParticipationImprovementSchedule& operator= (const ParticipationImprovementSchedule&) = delete;
    ParticipationImprovementSchedule& operator= (ParticipationImprovementSchedule&&) = delete;
public:
    explicit ParticipationImprovementSchedule(const Context& context) :
        _active_blocks(context.partition.k, true),
        _active_blocks_next_round(context.partition.k, false),
        _participations(context.partition.k, 0),
        _processed(context.partition.k, vec<uint8_t>(context.partition.k, false)),
        _partitions_sorted_by_participations(),
        _scheduled(context.partition.k, false),
        _round(0),
        _k(context.partition.k),
        _rng(context.partition.seed),
        _active_block_pairs(_k) {}

    void resetForNewRound(const DeterministicQuotientGraph<TypeTraits>& qg) {
        _scheduled.assign(_scheduled.size(), false);
        std::swap(_active_blocks, _active_blocks_next_round);
        _active_blocks_next_round.assign(_active_blocks_next_round.size(), false);
        _participations.assign(_participations.size(), 0);
        for (auto& v : _processed) {
            v.assign(v.size(), false);
        }
        for (auto& active_pairs : _active_block_pairs) {
            active_pairs.clear();
            active_pairs.reserve(_k - 1);
        }
        for (PartitionID i = 0; i < _k - 1; ++i) {
            for (PartitionID j = i + 1; j < _k; ++j) {
                if (isEligible(i, j, qg)) {
                    _participations[i]++;
                    _participations[j]++;
                    _active_block_pairs[i].push_back(j);
                    _active_block_pairs[j].push_back(i);
                }
            }
        }
        _partitions_sorted_by_participations.clear();
        _partitions_sorted_by_participations.reserve(_k);
        for (PartitionID i = 0; i < _k; ++i) {
            if (_participations[i] > 0) {
                _partitions_sorted_by_participations.push_back(i);
            }
        }
        std::sort(_partitions_sorted_by_participations.begin(), _partitions_sorted_by_participations.end(), [&](const PartitionID i, const PartitionID j) {
            return _participations[i] > _participations[j] || (_participations[i] == _participations[j] && i < j);
        });
        tbb::parallel_for(PartitionID(0), _k, [&](const size_t i) {
            auto& active_pairs = _active_block_pairs[i];
            std::sort(active_pairs.begin(), active_pairs.end(), [&](const PartitionID& lhs, const PartitionID& rhs) {
                const HyperedgeWeight lImprove = qg.getImprovement(i, lhs);
                const HyperedgeWeight rImprove = qg.getImprovement(i, rhs);
                const HyperedgeWeight lCut = qg.getCutWeight(i, lhs);
                const HyperedgeWeight rCut = qg.getCutWeight(i, rhs);

                return std::tie(lImprove, lCut, lhs) > std::tie(rImprove, rCut, rhs);
            });
        });

        _round++;
    }

    bool hasActiveBlocks() {
        return _partitions_sorted_by_participations.size() > 0;
    }


    void reportResults(const PartitionID block0, const PartitionID block1, const MoveSequence& sequence) {
        if (sequence.expected_improvement > 0) {
            _active_blocks_next_round[block0] = true;
            _active_blocks_next_round[block1] = true;
        }
    }

    void setTopLevelFlag() {
        _top_level = true;
    }

private:
    void initializeImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph, const DeterministicQuotientGraph<TypeTraits>& qg) {
        unused(hypergraph);
        _active_blocks_next_round.assign(_active_blocks_next_round.size(), true);
        resetForNewRound(qg);
        _round = 0;
    }

    size_t getNextMatchingImpl(tbb::concurrent_queue<ScheduledPair>& tasks, const DeterministicQuotientGraph<TypeTraits>& qg) {
        _scheduled.assign(_scheduled.size(), false);
        assert(_scheduled.size() == size_t(_k));
        assert(_partitions_sorted_by_participations.size() <= size_t(_k));
        size_t taskSize = 0;
        if (_partitions_sorted_by_participations.size() > 0) {
            for (size_t i = 0; i < _partitions_sorted_by_participations.size(); ++i) {
                const PartitionID block0 = _partitions_sorted_by_participations[i];
                DBG << V(block0) << ", " << V(_participations[block0]);
                assert(block0 < _k);
                if (_scheduled[block0]) continue;
                for (PartitionID block1 : _active_block_pairs[block0]) {
                    assert(block1 < _k);
                    if (_scheduled[block1] || _processed[block0][block1] || block0 == block1) continue;
                    const PartitionID smaller = std::min(block0, block1);
                    const PartitionID larger = std::max(block0, block1);
                    if (isEligible(smaller, larger, qg)) {
                        addBlockPair(smaller, larger);
                        tasks.push({ {smaller, larger}, _rng() });
                        taskSize++;
                        i--;
                        break;
                    }
                }
                if (taskSize == size_t(_k) / 2) { break; }
            }
        }
        if constexpr (debug) {
            for (auto it = tasks.unsafe_begin(); it != tasks.unsafe_end(); ++it) {
                auto t = *it;
                DBG << "Scheduling: (" << t.bp.i << ", " << t.bp.j << ", " << t.seed << ") in round " << _round;
            }
        }
        return taskSize;
    }

    void addBlockPair(const PartitionID i, const PartitionID j) {
        //DBG << V(i) << ", " << V(j);
        assert(i < _k);
        assert(j < _k);
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
                while (swapPosition < _partitions_sorted_by_participations.size()
                    && (
                        _participations[_partitions_sorted_by_participations[swapPosition]] > _participations[_partitions_sorted_by_participations[i]]
                        || (_participations[_partitions_sorted_by_participations[swapPosition]] == _participations[_partitions_sorted_by_participations[i]] && _partitions_sorted_by_participations[swapPosition] < _partitions_sorted_by_participations[i])
                        )) {
                    swapPosition++;
                }
                swapPosition--;
                assert(swapPosition < _partitions_sorted_by_participations.size());
                std::swap(_partitions_sorted_by_participations[i], _partitions_sorted_by_participations[swapPosition]);
                break;
            }
        }
        if (_participations[_partitions_sorted_by_participations[_partitions_sorted_by_participations.size() - 1]] == 0) { _partitions_sorted_by_participations.pop_back(); }
    }

    bool isEligible(const PartitionID i, const PartitionID j, const DeterministicQuotientGraph<TypeTraits>& qg) {
        //DBG << "isEligible: " << V(i) << ", " << V(j);
        assert(i < j);
        const HyperedgeWeight weight = qg.getCutWeight(i, j);
        const bool skip_small_cuts = !_top_level;
        const bool contains_enough_cut_hes =
            (skip_small_cuts && weight > 10) ||
            (!skip_small_cuts && weight > 0);
        return (_active_blocks[i] || _active_blocks[j]) // at least one block active
            && !_scheduled[i] && !_scheduled[j]         // none of the blocks is already scheduled
            && !_processed[i][j]                        // the pairing has not been processed yet
            && (!qg.wasAlreadyScheduled(i, j) || qg.getImprovement(i, j) > 0)
            && contains_enough_cut_hes;
    }

    vec<uint8_t> _active_blocks;
    vec<uint8_t> _active_blocks_next_round;
    vec<size_t> _participations;
    vec<vec<uint8_t>> _processed;
    vec<PartitionID> _partitions_sorted_by_participations;
    vec<bool> _scheduled;
    size_t _round;
    PartitionID _k;
    std::mt19937 _rng;
    bool _top_level = false;
    vec<vec<PartitionID>> _active_block_pairs;
};

}  // namespace mt_kahypar

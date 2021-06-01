//
// Created by mlaupichler on 26.05.21.
//

#pragma once

namespace mt_kahypar::ds {

    template<typename IndexType>
    class ThreadSafeFlagArray {

    private:
        static constexpr bool debug = false;
        static constexpr bool enable_heavy_assert = true;

    public:
        explicit ThreadSafeFlagArray(IndexType size) : _data(size, CAtomic<bool>(false)) {}
        ThreadSafeFlagArray(const ThreadSafeFlagArray&) = delete;
        ThreadSafeFlagArray& operator= (const ThreadSafeFlagArray&) = delete;

        ThreadSafeFlagArray(ThreadSafeFlagArray&&) = delete;
        ThreadSafeFlagArray& operator= (ThreadSafeFlagArray&&) = delete;

        ~ThreadSafeFlagArray() {
            HEAVY_REFINEMENT_ASSERT(checkAllFalse());
        }

        IndexType size() {
            return _data.size();
        }

        bool compare_and_set_to_true(IndexType idx) {
            ASSERT(idx < size());
            bool expected = _data[idx].load(std::memory_order_relaxed);
            bool desired = true;
            return expected != desired && _data[idx].compare_exchange_strong(expected, desired);
        }

        bool compare_and_set_to_false(IndexType idx) {
            ASSERT(idx < size());
            bool expected = _data[idx].load(std::memory_order_relaxed);
            bool desired = false;
            return expected != desired && _data[idx].compare_exchange_strong(expected, desired);
        }

        bool loadAt(IndexType idx) {
            ASSERT(idx < size());
            return _data[idx].load(std::memory_order_relaxed);
        }

        bool operator[](IndexType idx) {
            return loadAt(idx);
        }

        // ! Only for testing
        bool checkAllFalse() {
            return std::all_of(_data.begin(),_data.end(), [&](const CAtomic<bool>& flag){return !flag.load(std::memory_order_relaxed);});
        }

    private:
        Array<CAtomic<bool>> _data;
    };

} // namespace mt_kahypar::ds
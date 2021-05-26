//
// Created by mlaupichler on 26.05.21.
//

#pragma once

namespace mt_kahypar::ds {

    template<typename IndexType>
    class ThreadSafeFlagArray {

    public:
        ThreadSafeFlagArray() : _data(), _initialized(false) {}
        ThreadSafeFlagArray(const ThreadSafeFlagArray&) = delete;
        ThreadSafeFlagArray& operator= (const ThreadSafeFlagArray&) = delete;

        ThreadSafeFlagArray(ThreadSafeFlagArray&&) = delete;
        ThreadSafeFlagArray& operator= (ThreadSafeFlagArray&&) = delete;

        void initialize(IndexType size) {
            ASSERT(!isInitialized());
            _data = parallel::scalable_vector<CAtomic<bool>>(size,CAtomic(false));
            _initialized = true;
        }

        void reset() {
            ASSERT(isInitialized());
            _initialized = false;
            _data.clear();
        }

        bool isInitialized() {
            return _initialized;
        }

        IndexType size() {
            ASSERT(isInitialized());
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
            ASSERT(isInitialized());
            return std::all_of(_data.begin(),_data.end(), [&](const CAtomic<bool>& flag){return !flag.load(std::memory_order_relaxed);});
        }

    private:
        parallel::scalable_vector<CAtomic<bool>> _data;
        bool _initialized;
    };

} // namespace mt_kahypar::ds
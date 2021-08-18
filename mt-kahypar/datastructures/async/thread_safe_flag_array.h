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

        using UnderlyingAtomic = CAtomic<bool>;

    public:
        explicit ThreadSafeFlagArray(IndexType size) : _size(size), _data(std::make_unique<UnderlyingAtomic[]>(size)) {
            init();
        }
        ThreadSafeFlagArray(const ThreadSafeFlagArray&) = delete;
        ThreadSafeFlagArray& operator= (const ThreadSafeFlagArray&) = delete;

        ThreadSafeFlagArray(ThreadSafeFlagArray&&) = delete;
        ThreadSafeFlagArray& operator= (ThreadSafeFlagArray&&) = delete;

        ~ThreadSafeFlagArray() {
            HEAVY_REFINEMENT_ASSERT(checkAllFalse());
        }

        void init(const bool initialiser = false) {
            for ( size_t i = 0; i < _size; ++i ) {
                _data[i].store(initialiser, std::memory_order_relaxed);
            }
        }

        IndexType size() {
            return _size;
        }

        bool compare_and_set_to_true(const IndexType& idx) {
            ASSERT(idx < size());
            auto expected = _data[idx].load(std::memory_order_relaxed);
            auto desired = true;
            return expected != desired && _data[idx].compare_exchange_strong(expected, desired,
                                                      /* memory order for write on success: */ std::memory_order_relaxed,
                                                      /* memory order for load on fail: */ std::memory_order_relaxed);
        }

        bool compare_and_set_to_false(const IndexType& idx) {
            ASSERT(idx < size());
            auto expected = _data[idx].load(std::memory_order_relaxed);
            auto desired = false;
            return expected != desired && _data[idx].compare_exchange_strong(expected, desired,
                    /* memory order for write on success: */ std::memory_order_relaxed,
                    /* memory order for load on fail: */ std::memory_order_relaxed);
        }

        bool isSet(IndexType idx) {
            ASSERT(idx < size());
            return _data[idx].load(std::memory_order_relaxed);
        }

//        bool operator[](IndexType idx) {
//            return isSet(idx);
//        }

        // ! Only for testing
        bool checkAllFalse() {
            for (IndexType i = 0; i < _size; ++i) {
                if (_data[i].load(std::memory_order_relaxed)) return false;
            }
            return true;
        }

    private:
        const IndexType _size;
        std::unique_ptr<UnderlyingAtomic[]> _data;
    };

} // namespace mt_kahypar::ds
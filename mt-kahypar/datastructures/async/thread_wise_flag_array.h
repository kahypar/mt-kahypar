//
// Created by mlaupichler on 18.08.21.
//

#pragma once


namespace mt_kahypar::ds {

    template<typename UnsafeBlock, typename IndexType = uint32_t>
    struct BlockThreadWiseFlagArray {

        static_assert(std::is_integral<UnsafeBlock>::value && std::is_unsigned<UnsafeBlock>::value,
                      "Block for BlockThreadWiseFlagArray has to be an unsigned integral type!");
        static_assert(std::is_integral<IndexType>::value && std::is_unsigned<IndexType>::value,
                      "IndexType for BlockThreadWiseFlagArray has to be an unsigned integral type!");

    private:

        using Block = CAtomic<UnsafeBlock>;

    public:

        BlockThreadWiseFlagArray(const size_t num_elements, const size_t num_threads) :
            _num_elements(num_elements),
            _num_threads(num_threads),
            _blocks(_num_elements, CAtomic<UnsafeBlock>(0)) {
          if (_num_threads > BLOCK_SIZE) {
            ERROR("Blocks of size " << BLOCK_SIZE << " cannot accommodate one bit each for " << _num_threads
                                    << " threads!");
          }
        }

        size_t size_in_bytes() const {
          return _blocks.size() * sizeof(Block);
        }

        bool is_set(const IndexType idx, const size_t tid) const {
          ASSERT(idx < _num_elements);
          ASSERT(tid < _num_threads);
          UnsafeBlock mask = UnsafeBlock(1) << tid;
          const Block &block = _blocks[idx];
          UnsafeBlock val = block.load(std::memory_order_relaxed);
          return val & mask;
        }

        bool any_set(const IndexType idx) const {
          ASSERT(idx < _num_elements);
          return _blocks[idx].load(std::memory_order_relaxed) != 0;
        }

        bool any_set_except_thread(const IndexType idx, const size_t except_tid) const {
          ASSERT(idx < _num_elements);
          ASSERT(except_tid < _num_threads);
          UnsafeBlock mask = ~(UnsafeBlock(1) << except_tid);
          const Block &block = _blocks[idx];
          const UnsafeBlock block_except_thread = block.load(std::memory_order_relaxed) & mask;
          return block_except_thread != 0;
        }

        std::pair<bool, bool> any_set_with_and_without_thread(const IndexType idx, const size_t except_tid) const {
          ASSERT(idx < _num_elements);
          ASSERT(except_tid < _num_threads);
          UnsafeBlock mask = ~(UnsafeBlock(1) << except_tid);
          const Block &block = _blocks[idx];
          const UnsafeBlock unsafe_block = block.load(std::memory_order_relaxed);
          const UnsafeBlock block_except_thread = unsafe_block & mask;
          return std::make_pair(unsafe_block != 0, block_except_thread != 0);
        }

        // ! Guarantees that when this call returns, the entry is set to true. Returns true if this call was the one to
        // ! set the bit and false if it was already set to true or a different call to this function set it to true concurrently.
        bool set_true(const IndexType idx, const size_t tid, bool &block_was_zero) {
          ASSERT(idx < _num_elements);
          ASSERT(tid < _num_threads);
          block_was_zero = false;
          UnsafeBlock mask = UnsafeBlock(1) << tid;
          Block &block = _blocks[idx];
          UnsafeBlock expected = block.load(std::memory_order_relaxed);
          UnsafeBlock desired = expected | mask;
          if (expected == desired) return false;
          while (!block.compare_exchange_strong(expected, desired, std::memory_order_relaxed)) {
            // If block changed in the meantime, adapt desired to changes in other bits
            desired = expected | mask;
            if (expected == desired) return false;
          }
          if (expected == 0) {
            block_was_zero = true;
          }
          return true;
        }

        // ! Guarantees that when this call returns, the entry is set to false. Returns true if this call was the one to
        // ! set the bit to false and true if it was already set to false or a different call to this function set it to false concurrently.
        bool set_false(const IndexType idx, const size_t tid, bool &set_block_to_zero) {
          ASSERT(idx < _num_elements);
          ASSERT(tid < _num_threads);
          set_block_to_zero = false;
          UnsafeBlock mask = ~(UnsafeBlock(1) << tid);
          Block &block = _blocks[idx];
          UnsafeBlock expected = block.load(std::memory_order_relaxed);
          UnsafeBlock desired = expected & mask;
          if (expected == desired) return false;
          while (!block.compare_exchange_strong(expected, desired, std::memory_order_relaxed)) {
            // If block changed in the meantime, adapt desired to changes in other bits
            desired = expected & mask;
            if (expected == desired) return false;
          }
          if (desired == 0) {
            set_block_to_zero = true;
          }
          return true;
        }

    private:

        static constexpr int BLOCK_SIZE = std::numeric_limits<UnsafeBlock>::digits;
        const size_t _num_elements;
        const size_t _num_threads;
        Array <Block> _blocks;

    };

    template<typename IndexType = uint32_t>
    struct MutexThreadWiseFlagArray {

        static_assert(std::is_integral<IndexType>::value && std::is_unsigned<IndexType>::value,
                      "IndexType for MutexThreadWiseFlagArray has to be an unsigned integral type!");

    private:

        using UnsafeBlock = uint32_t;
        using BlockIndex = size_t;

    public:

        MutexThreadWiseFlagArray(const size_t num_elements, const size_t num_threads) :
            _num_elements(num_elements),
            _num_threads(num_threads),
            _blocks_per_element(num_threads / BLOCK_SIZE + (num_threads % BLOCK_SIZE == 0 ? 0 : 1)),
            _blocks(_num_elements * _blocks_per_element, UnsafeBlock(0)),
            _spinlocks(_num_elements) {
        }

        size_t size_in_bytes() const {
          return _blocks.size() * sizeof(UnsafeBlock) + _spinlocks.size() * sizeof(SpinLock);
        }

        bool any_set(const IndexType idx) {
          ASSERT(idx < _num_elements);
          auto block_range = IteratorRange(_blocks.begin() + idx * _blocks_per_element, _blocks.begin() + (idx + 1) *_blocks_per_element);
          _spinlocks[idx].lock();
          for (const UnsafeBlock& block : block_range) {
            if (block != 0) {
              _spinlocks[idx].unlock();
              return true;
            }
          }
          _spinlocks[idx].unlock();
          return false;
        }

        bool is_set(const IndexType idx, const size_t tid) {
          ASSERT(idx < _num_elements);
          ASSERT(tid < _num_threads);
          BlockIndex block_idx = idx * _blocks_per_element + tid / BLOCK_SIZE;
          size_t offset = tid % BLOCK_SIZE;
          UnsafeBlock mask = UnsafeBlock(1) << offset;
          _spinlocks[idx].lock();
          bool set = _blocks[block_idx] & mask;
          _spinlocks[idx].unlock();
          return set;
        }

        bool any_set_except_thread(const IndexType idx, const size_t except_tid) {
          ASSERT(idx < _num_elements);
          ASSERT(except_tid < _num_threads);
          const BlockIndex first_block = idx * _blocks_per_element;
          const BlockIndex first_block_of_next = (idx + 1) * _blocks_per_element;
          const BlockIndex except_block_idx = idx * _blocks_per_element + except_tid / BLOCK_SIZE;
          size_t except_offset = except_tid % BLOCK_SIZE;
          UnsafeBlock mask = ~(UnsafeBlock(1) << except_offset);

          _spinlocks[idx].lock();
          for (BlockIndex bidx = first_block; bidx < first_block_of_next; ++bidx) {
            UnsafeBlock block = _blocks[bidx];
            if (bidx == except_block_idx) {
              block &= mask;
            }
            if (block != 0) {
              _spinlocks[idx].unlock();
              return true;
            }
          }
          _spinlocks[idx].unlock();
          return false;
        }

        std::pair<bool, bool> any_set_with_and_without_thread(const IndexType idx, const size_t except_tid) {
          ASSERT(idx < _num_elements);
          ASSERT(except_tid < _num_threads);
          const BlockIndex first_block = idx * _blocks_per_element;
          const BlockIndex first_block_of_next = (idx + 1) * _blocks_per_element;
          const BlockIndex except_block_idx = idx * _blocks_per_element + except_tid / BLOCK_SIZE;
          size_t except_offset = except_tid % BLOCK_SIZE;
          UnsafeBlock mask = ~(UnsafeBlock(1) << except_offset);

          auto result = std::make_pair(false, false);

          _spinlocks[idx].lock();
          for (BlockIndex bidx = first_block; bidx < first_block_of_next; ++bidx) {
            UnsafeBlock block = _blocks[bidx];
            if (bidx == except_block_idx) {
              if (block != 0) {
                result.first = true;
              }
              const UnsafeBlock block_without_thread = block & mask;
              if (block_without_thread != 0) {
                result.second = true;
              }
            } else {
              if (block != 0) {
                result.first = true;
                result.second = true;
              }
            }

            if (result.second) {
              _spinlocks[idx].unlock();
              return result;
            }
          }
          _spinlocks[idx].unlock();
          return result;
        }

        // ! Guarantees that when this call returns, the entry is set to true. Returns true if this call was the one to
        // ! set the bit and false if it was already set to true or a different call to this function set it to true concurrently.
        bool set_true(const IndexType idx, const size_t tid, bool &block_was_zero) {
          ASSERT(idx < _num_elements);
          ASSERT(tid < _num_threads);
          BlockIndex block_idx = idx * _blocks_per_element + tid / BLOCK_SIZE;
          size_t offset = tid % BLOCK_SIZE;
          UnsafeBlock mask = UnsafeBlock(1) << offset;
          UnsafeBlock &block = _blocks[block_idx];
          _spinlocks[idx].lock();
          const UnsafeBlock before = block;
          block |= mask;
          bool changed = before != block;
          _spinlocks[idx].unlock();
          block_was_zero = (before == 0);
          return changed;
        }

        // ! Guarantees that when this call returns, the entry is set to false. Returns true if this call was the one to
        // ! set the bit to false and true if it was already set to false or a different call to this function set it to false concurrently.
        bool set_false(const IndexType idx, const size_t tid, bool &set_block_to_zero) {
          ASSERT(idx < _num_elements);
          ASSERT(tid < _num_threads);
          BlockIndex block_idx = idx * _blocks_per_element + tid / BLOCK_SIZE;
          size_t offset = tid % BLOCK_SIZE;
          UnsafeBlock mask = ~(UnsafeBlock(1) << offset);
          UnsafeBlock &block = _blocks[block_idx];
          _spinlocks[idx].lock();
          const UnsafeBlock before = block;
          block &= mask;
          bool changed = before != block;
          set_block_to_zero = (changed && block == 0);
          _spinlocks[idx].unlock();
          return changed;
        }

    private:

        static constexpr int BLOCK_SIZE = std::numeric_limits<UnsafeBlock>::digits;

        const size_t _num_elements;
        const size_t _num_threads;
        const size_t _blocks_per_element;
        Array <UnsafeBlock> _blocks;
        Array <SpinLock> _spinlocks;

    };

    // Ugliest facade ever but due to the fact that the number of threads is not statically known and we can use no interfaces due to slow vtable lookup, I don't see a way around this.
    template<typename IndexType = uint32_t>
    class ThreadWiseFlagArray {

        static_assert(std::is_integral<IndexType>::value && std::is_unsigned<IndexType>::value,
                      "IndexType for ThreadWiseFlagArray has to be an unsigned integral type!");

    private:
        enum Type {
            eight_bit_block,
            sixteen_bit_block,
            thirtytwo_bit_block,
            sixtyfour_bit_block,
            mutex
        };

    public:

        ThreadWiseFlagArray(const IndexType num_elements, const size_t num_threads) : _type(get_type(num_threads)), _num_elements(num_elements) {
          make_array(num_elements, num_threads);
        }

        IndexType numElements() const {
          return _num_elements;
        }

        size_t size_in_bytes() const {
          switch (_type) {
            case eight_bit_block :
              return _eight_bit_block_array->size_in_bytes();
            case sixteen_bit_block :
              return _sixteen_bit_block_array->size_in_bytes();
            case thirtytwo_bit_block :
              return _thirtytwo_bit_block_array->size_in_bytes();
            case sixtyfour_bit_block :
              return _sixtyfour_bit_block_array->size_in_bytes();
            case mutex :
              return _mutex_array->size_in_bytes();
          }
          ERROR("Type not set in ThreadWiseFlagArray!");
        }

        bool is_set(const IndexType idx, const size_t tid) {
          switch (_type) {
            case eight_bit_block :
              return _eight_bit_block_array->is_set(idx, tid);
            case sixteen_bit_block :
              return _sixteen_bit_block_array->is_set(idx, tid);
            case thirtytwo_bit_block :
              return _thirtytwo_bit_block_array->is_set(idx, tid);
            case sixtyfour_bit_block :
              return _sixtyfour_bit_block_array->is_set(idx, tid);
            case mutex :
              return _mutex_array->is_set(idx, tid);
          }
          ERROR("Type not set in ThreadWiseFlagArray!");
        }


        bool any_set(const IndexType idx) {
          switch (_type) {
            case eight_bit_block :
              return _eight_bit_block_array->any_set(idx);
            case sixteen_bit_block :
              return _sixteen_bit_block_array->any_set(idx);
            case thirtytwo_bit_block :
              return _thirtytwo_bit_block_array->any_set(idx);
            case sixtyfour_bit_block :
              return _sixtyfour_bit_block_array->any_set(idx);
            case mutex :
              return _mutex_array->any_set(idx);
          }
          ERROR("Type not set in ThreadWiseFlagArray!");
        }

        bool any_set_except_thread(const IndexType idx, const size_t except_tid) {
          switch (_type) {
            case eight_bit_block :
              return _eight_bit_block_array->any_set_except_thread(idx, except_tid);
            case sixteen_bit_block :
              return _sixteen_bit_block_array->any_set_except_thread(idx, except_tid);
            case thirtytwo_bit_block :
              return _thirtytwo_bit_block_array->any_set_except_thread(idx, except_tid);
            case sixtyfour_bit_block :
              return _sixtyfour_bit_block_array->any_set_except_thread(idx, except_tid);
            case mutex :
              return _mutex_array->any_set_except_thread(idx, except_tid);
          }
          ERROR("Type not set in ThreadWiseFlagArray!");
        }

        std::pair<bool, bool> any_set_with_and_without_thread(const IndexType idx, const size_t except_tid) {
          switch (_type) {
            case eight_bit_block :
              return _eight_bit_block_array->any_set_with_and_without_thread(idx, except_tid);
            case sixteen_bit_block :
              return _sixteen_bit_block_array->any_set_with_and_without_thread(idx, except_tid);
            case thirtytwo_bit_block :
              return _thirtytwo_bit_block_array->any_set_with_and_without_thread(idx, except_tid);
            case sixtyfour_bit_block :
              return _sixtyfour_bit_block_array->any_set_with_and_without_thread(idx, except_tid);
            case mutex :
              return _mutex_array->any_set_with_and_without_thread(idx, except_tid);
          }
          ERROR("Type not set in ThreadWiseFlagArray!");
        }

        // ! Guarantees that when this call returns, the entry is set to true. Returns true if this call was the one to
        // ! set the bit and false if it was already set to true or a different call to this function set it to true concurrently.
        bool set_true(const IndexType idx, const size_t tid, bool &block_was_zero) {
          switch (_type) {
            case eight_bit_block :
              return _eight_bit_block_array->set_true(idx, tid, block_was_zero);
            case sixteen_bit_block :
              return _sixteen_bit_block_array->set_true(idx, tid, block_was_zero);
            case thirtytwo_bit_block :
              return _thirtytwo_bit_block_array->set_true(idx, tid, block_was_zero);
            case sixtyfour_bit_block :
              return _sixtyfour_bit_block_array->set_true(idx, tid, block_was_zero);
            case mutex :
              return _mutex_array->set_true(idx, tid, block_was_zero);
          }
          ERROR("Type not set in ThreadWiseFlagArray!");
        }

        // ! Guarantees that when this call returns, the entry is set to false. Returns true if this call was the one to
        // ! set the bit to false and true if it was already set to false or a different call to this function set it to false concurrently.
        bool set_false(const IndexType idx, const size_t tid, bool &set_block_to_zero) {
          switch (_type) {
            case eight_bit_block :
              return _eight_bit_block_array->set_false(idx, tid, set_block_to_zero);
            case sixteen_bit_block :
              return _sixteen_bit_block_array->set_false(idx, tid, set_block_to_zero);
            case thirtytwo_bit_block :
              return _thirtytwo_bit_block_array->set_false(idx, tid, set_block_to_zero);
            case sixtyfour_bit_block :
              return _sixtyfour_bit_block_array->set_false(idx, tid, set_block_to_zero);
            case mutex :
              return _mutex_array->set_false(idx, tid, set_block_to_zero);
          }
          ERROR("Type not set in ThreadWiseFlagArray!");
        }

    private:

        void make_array(const size_t num_elements, const size_t num_threads) {
          switch (_type) {
            case eight_bit_block :
              _eight_bit_block_array = std::make_unique<BlockThreadWiseFlagArray<uint8_t, IndexType>>(num_elements,
                                                                                                      num_threads);
            case sixteen_bit_block :
              _sixteen_bit_block_array = std::make_unique<BlockThreadWiseFlagArray<uint16_t, IndexType>>(num_elements,
                                                                                                         num_threads);
            case thirtytwo_bit_block :
              _thirtytwo_bit_block_array = std::make_unique<BlockThreadWiseFlagArray<uint32_t, IndexType>>(num_elements,
                                                                                                           num_threads);
            case sixtyfour_bit_block :
              _sixtyfour_bit_block_array = std::make_unique<BlockThreadWiseFlagArray<uint64_t, IndexType>>(num_elements,
                                                                                                           num_threads);
            case mutex :
              _mutex_array = std::make_unique<MutexThreadWiseFlagArray<IndexType>>(num_elements, num_threads);
          }
        }

        static Type get_type(const size_t num_threads) {
          if (num_threads <= 8) return eight_bit_block;
          if (num_threads <= 16) return sixteen_bit_block;
          if (num_threads <= 32) return thirtytwo_bit_block;
          if (num_threads <= 64) return sixtyfour_bit_block;
          return mutex;
        }

        const Type _type;
        const IndexType _num_elements;

        std::unique_ptr <BlockThreadWiseFlagArray<uint8_t, IndexType>> _eight_bit_block_array;
        std::unique_ptr <BlockThreadWiseFlagArray<uint16_t, IndexType>> _sixteen_bit_block_array;
        std::unique_ptr <BlockThreadWiseFlagArray<uint32_t, IndexType>> _thirtytwo_bit_block_array;
        std::unique_ptr <BlockThreadWiseFlagArray<uint64_t, IndexType>> _sixtyfour_bit_block_array;
        std::unique_ptr <MutexThreadWiseFlagArray<IndexType>> _mutex_array;
    };

} // end namespace
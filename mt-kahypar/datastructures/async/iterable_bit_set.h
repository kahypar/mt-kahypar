//
// Created by mlaupichler on 07.06.21.
//

#pragma once

#include "mt-kahypar/datastructures/array.h"

namespace mt_kahypar::ds {

/// Non-thread-safe bit-set that is memory efficient and provides an iterator over indices that are set to true.
/// Changes to the bits will have iterators act according to the changes, e.g. if bits i and i+1 are set and an iterator
/// has reached i but then bit i+1 is set to false, i+1 will not be considered by the iterator on increment. (Changes
/// to bits i and smaller will not affect an iterator at bit i.)
template<typename IndexType>
class IterableBitSet {

private:
    using Block = uint64_t;
    using BlockIterator = Array<Block>::const_iterator;
    static constexpr uint32_t BITS_PER_BLOCK = std::numeric_limits<Block>::digits;

    class BitReference {
    public:
        BitReference(Block* const block_ptr, const IndexType bit_offset) :
                _block_ptr(block_ptr),
                _bit_offset(bit_offset) {}

        BitReference(const BitReference&) = delete;
        BitReference(BitReference&&) = delete;
        BitReference& operator=(const BitReference&) = delete;
        BitReference operator=(BitReference&&) = delete;

        // Implicit conversion to bool operator
        operator bool() const {return (*_block_ptr) & (Block(1) << _bit_offset);} // NOLINT(google-explicit-constructor)
        BitReference& operator=(const bool bit) {
            setBit(bit);
            return *this;
        }

    private:

        void setBit(bool val) {
            if (val) {
                (*_block_ptr) |= (Block(1) << _bit_offset);
            } else {
                (*_block_ptr) &= ~(Block(1) << _bit_offset);
            }
        }

        Block* const _block_ptr;
        const IndexType _bit_offset;
    };

public:

    explicit IterableBitSet(IndexType n) :
            _n(n),
            _blocks("Refinement", "bit_set", n / BITS_PER_BLOCK + (n % BITS_PER_BLOCK == 0 ? 0 : 1), true, false) {}

    IterableBitSet() : _n(0), _blocks("Refinement", "bit_set", 0) {}

    void resize(IndexType n) {
      _n = n;
      _blocks.resize("Refinement", "bit_set", n / BITS_PER_BLOCK + (n % BITS_PER_BLOCK == 0 ? 0 : 1), true, false);
    }

    IndexType size() {
        return _n;
    }

    size_t size_in_bytes() const {
      return static_cast<size_t>(sizeof(Block)) * static_cast<size_t>(_blocks.size());
    }

    bool isSet(const IndexType i) const {
        ASSERT(i < _n);
        const IndexType div = i / BITS_PER_BLOCK;
        const IndexType rem = i % BITS_PER_BLOCK;
        ASSERT(div < _blocks.size());
        return _blocks[div] & (Block(1) << rem);
    }

    BitReference operator[](const IndexType i) {
        ASSERT(i < _n);
        const IndexType div = i / BITS_PER_BLOCK;
        const IndexType rem = i % BITS_PER_BLOCK;
        ASSERT(div < _blocks.size());
        return BitReference(&_blocks[div], rem);
    }

    void set_true(const IndexType i) {
        ASSERT(i < _n);
        const IndexType div = i / BITS_PER_BLOCK;
        const IndexType rem = i % BITS_PER_BLOCK;
        ASSERT(div < _blocks.size());
        _blocks[div] |= (Block(1) << rem);
    }

    void set_false(const IndexType i) {
        ASSERT(i < _n);
        const IndexType div = i / BITS_PER_BLOCK;
        const IndexType rem = i % BITS_PER_BLOCK;
        ASSERT(div < _blocks.size());
        _blocks[div] &= ~(Block(1) << rem);
    }

    class Iterator : public std::iterator<std::forward_iterator_tag, IndexType, std::ptrdiff_t, const IndexType*, IndexType> {
    public:
        Iterator(const BlockIterator& first, IndexType i, const IndexType n) : currentElement(i), _n(n), firstBlockIt(first) {
            findNextBit();
        }

        IndexType operator*() const {
            return currentElement;
        }

        Iterator& operator++() {
            findNextBit();
            return *this;
        }

        Iterator operator++(int ) {
            const Iterator res = *this;
            findNextBit();
            return res;
        }

        bool operator==(const Iterator& o) const {
            return currentElement == o.currentElement && firstBlockIt == o.firstBlockIt;
        }

        bool operator!=(const Iterator& o) const {
            return !operator==(o);
        }

    private:

        IndexType currentElement;
        const IndexType _n;
        const BlockIterator firstBlockIt;

        void findNextBit() {
            ++currentElement;
            Block b = *currentBlock();
            while (b >> (currentElement % BITS_PER_BLOCK) == 0 && currentElement < _n) {
                currentElement += (BITS_PER_BLOCK - (currentElement % BITS_PER_BLOCK));   // skip rest of block
                b = *currentBlock();
            }
            if (currentElement < _n) {
                currentElement += utils::lowest_set_bit_64(b >> (currentElement % BITS_PER_BLOCK));
            } else {
                currentElement = _n;
            }
        }

        BlockIterator currentBlock() const {
            return firstBlockIt + (currentElement / BITS_PER_BLOCK);
        }
    };

    Iterator begin() const {
        return Iterator(_blocks.cbegin(), -1, _n);
    }

    Iterator end() const {
        return Iterator(_blocks.cbegin(), _n-1, _n);
    }

private:

    IndexType _n;
    Array<Block> _blocks;

};


using PartitionBitSet = IterableBitSet<PartitionID>;

} // end namespace mt_kahypar::ds
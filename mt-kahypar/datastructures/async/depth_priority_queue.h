//
// Created by mlaupichler on 05.07.21.
//

#ifndef KAHYPAR_DEPTH_PRIORITY_QUEUE_H
#define KAHYPAR_DEPTH_PRIORITY_QUEUE_H

#include <tbb/concurrent_vector.h>

namespace mt_kahypar::ds {

    class DepthPriorityQueue {

        using depth_type = uint32_t;

    public:

        DepthPriorityQueue(const depth_type num_depths, std::vector<ContractionGroupID>&& max_elements_per_depth) :
            _num_depths(num_depths),
            _total_elements_per_depth(std::move(max_elements_per_depth)),
            _min_non_empty_depth(_num_depths),
            _cur_elements_per_depth(_num_depths, CAtomic<int64_t>(0)),
            _offsets(_num_depths, CAtomic<ContractionGroupID>(0)),
            _vecs(_num_depths, tbb::concurrent_vector<ContractionGroupID>()) {
          ASSERT(_num_depths == _total_elements_per_depth.size());
          // Reserve twice as much memory as elements per depth to accommodate re-insertions. (If that memory will be
          // exceeded further memory will be allocated on the fly)
          for (ContractionGroupID i = 0; i < _num_depths; ++i) {
            _vecs[i].reserve(_total_elements_per_depth[i] * 2);
          }
        }

        // ! Calls to push are expected to have a depth equal to or larger than the current smallest non-completed depth
        void push(const ContractionGroupID id, const depth_type depth) {
          ASSERT(depth < _num_depths);
          _vecs[depth].push_back(id);
          _cur_elements_per_depth[depth].fetch_add(1, std::memory_order_relaxed);

          // Update min non-empty if necessary
          depth_type current_min = _min_non_empty_depth.load(std::memory_order_relaxed);
          while (depth < current_min && !_min_non_empty_depth.compare_exchange_strong(current_min, depth)) {
            /* Attempt to update min depth if the depth of this push is smaller than the current min. Use CAS until
             * update works or a different thread has set min depth to an equal or smaller value */
          }
        }

        bool try_pop(ContractionGroupID& dest) {

          depth_type min_ptr = _min_non_empty_depth.load(std::memory_order_relaxed);

          for (depth_type depth = min_ptr; depth < _num_depths; ++depth) {

            if (_cur_elements_per_depth[depth].load(std::memory_order_relaxed) <= 0) continue;

            int64_t num_before_sub = _cur_elements_per_depth[depth].fetch_sub(1, std::memory_order_relaxed);
            if (num_before_sub > 0) {
              // Success: an element in this depth is available
              ContractionGroupID idx = _offsets[depth].fetch_add(1, std::memory_order_relaxed);
              ASSERT(idx < _vecs[depth].size());
              dest = _vecs[depth][idx];

              // Check if the popped element was the last one for this depth and if it is the min depth. If so increment
              // the min depth to a non-empty depth
              if (depth == min_ptr) {
                depth_type expected = min_ptr;
                while (expected < _num_depths
                        && _offsets[expected] == _vecs[expected].size()
                        && _min_non_empty_depth.compare_exchange_strong(expected, expected + 1, std::memory_order_relaxed)) {
                  ++expected;
                }
              }

              return true;
            } else {
              // Failure: Another thread took the last remaining element in this depth. Repair number of elements by
              // adding 1 back and continue with next depth
              _cur_elements_per_depth[depth].fetch_add(1, std::memory_order_relaxed);
              continue;
            }
          }

          return false;
        }

        ContractionGroupID unsafe_size() const {
          ContractionGroupID size = 0;
          for(depth_type i = 0; i < _num_depths; ++i) {
            size += _vecs[i].size() - _offsets[i].load(std::memory_order_relaxed);
          }
          return size;
        }

        bool unsafe_empty() const {
          return (unsafe_size() == 0);
        }

        // ! Clears the bucket contents while not changing number of depths or maximum number of elements per depth.
        // ! Maintains memory reserved for buckets.
        void reset() {
          for(depth_type i = 0; i < _num_depths; ++i) {
            _vecs[i].clear();
            _cur_elements_per_depth[i].store(0, std::memory_order_relaxed);
            _offsets[i].store(0, std::memory_order_relaxed);
          }
        }

    private:

        const depth_type _num_depths;
        const std::vector<ContractionGroupID> _total_elements_per_depth;

        CAtomic<depth_type> _min_non_empty_depth;
        Array<CAtomic<int64_t>> _cur_elements_per_depth;
        Array<CAtomic<ContractionGroupID>> _offsets;

        // This has to be a std::vector instead of an Array from array.h as constructing the array fails for some
        // tbb internal reason
        std::vector<tbb::concurrent_vector<ContractionGroupID>> _vecs;

    };

} // end namespace


#endif //KAHYPAR_DEPTH_PRIORITY_QUEUE_H

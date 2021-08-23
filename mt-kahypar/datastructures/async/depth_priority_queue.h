//
// Created by mlaupichler on 05.07.21.
//

#ifndef KAHYPAR_DEPTH_PRIORITY_QUEUE_H
#define KAHYPAR_DEPTH_PRIORITY_QUEUE_H

#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_queue.h>

namespace mt_kahypar::ds {

//    class DepthPriorityQueue {
//
//        using depth_type = uint32_t;
//
//    public:
//
//        DepthPriorityQueue(const depth_type num_depths, std::vector<ContractionGroupID>&& max_elements_per_depth) :
//            _num_depths(num_depths),
//            _total_elements_per_depth(std::move(max_elements_per_depth)),
//            _min_non_completed_depth(_num_depths),
//            _cur_elements_per_depth(_num_depths, CAtomic<int64_t>(0)),
//            _offsets(_num_depths, CAtomic<ContractionGroupID>(0)),
//            _vecs(_num_depths, tbb::concurrent_vector<ContractionGroupID>()) {
//          ASSERT(_num_depths == _total_elements_per_depth.size());
//          // Reserve twice as much memory as elements per depth to accommodate re-insertions. (If that memory will be
//          // exceeded further memory will be allocated on the fly)
//          for (ContractionGroupID i = 0; i < _num_depths; ++i) {
//            _vecs[i].reserve(_total_elements_per_depth[i] * 2);
//          }
//        }
//
//        // ! Calls to push are expected to have a depth equal to or larger than the current smallest non-completed depth
//        void push(const ContractionGroupID id, const depth_type depth) {
//          ASSERT(depth < _num_depths);
//          _vecs[depth].push_back(id);
//          _cur_elements_per_depth[depth].fetch_add(1, std::memory_order_relaxed);
//
//          // Update min non-empty if necessary
//          depth_type current_min = _min_non_completed_depth.load(std::memory_order_relaxed);
//          while (depth < current_min && !_min_non_completed_depth.compare_exchange_strong(current_min, depth)) {
//            /* Attempt to update min depth if the depth of this push is smaller than the current min. Use CAS until
//             * update works or a different thread has set min depth to an equal or smaller value */
//          }
//        }
//
//        bool try_pop(ContractionGroupID& dest) {
//
//          depth_type min_ptr = _min_non_completed_depth.load(std::memory_order_relaxed);
//
//          for (depth_type depth = min_ptr; depth < _num_depths; ++depth) {
//
//            if (_cur_elements_per_depth[depth].load(std::memory_order_relaxed) <= 0) continue;
//
//            int64_t num_before_sub = _cur_elements_per_depth[depth].fetch_sub(1, std::memory_order_relaxed);
//            if (num_before_sub > 0) {
//              // Success: an element in this depth is available
//              ContractionGroupID idx = _offsets[depth].fetch_add(1, std::memory_order_relaxed);
//              ASSERT(idx < _vecs[depth].size());
//              dest = _vecs[depth][idx];
//
//              // Check if the popped element was the last one for this depth and if it is the min depth. If so increment
//              // the min depth to a non-empty depth
//              if (depth == min_ptr) {
//                depth_type expected = min_ptr;
//                while (expected < _num_depths
//                        && _offsets[expected] == _vecs[expected].size()
//                        && _min_non_completed_depth.compare_exchange_strong(expected, expected + 1, std::memory_order_relaxed)) {
//                  ++expected;
//                }
//              }
//
//              return true;
//            } else {
//              // Failure: Another thread took the last remaining element in this depth. Repair number of elements by
//              // adding 1 back and continue with next depth
//              _cur_elements_per_depth[depth].fetch_add(1, std::memory_order_relaxed);
//              continue;
//            }
//          }
//
//          return false;
//        }
//
//        ContractionGroupID unsafe_size() const {
//          ContractionGroupID size = 0;
//          for(depth_type i = 0; i < _num_depths; ++i) {
//            size += _vecs[i].size() - _offsets[i].load(std::memory_order_relaxed);
//          }
//          return size;
//        }
//
//        bool unsafe_empty() const {
//          return (unsafe_size() == 0);
//        }
//
//        // ! Clears the bucket contents while not changing number of depths or maximum number of elements per depth.
//        // ! Maintains memory reserved for buckets.
//        void reset() {
//          for(depth_type i = 0; i < _num_depths; ++i) {
//            _vecs[i].clear();
//            _cur_elements_per_depth[i].store(0, std::memory_order_relaxed);
//            _offsets[i].store(0, std::memory_order_relaxed);
//          }
//        }
//
//    private:
//
//        const depth_type _num_depths;
//        const std::vector<ContractionGroupID> _total_elements_per_depth;
//
//        CAtomic<depth_type> _min_non_completed_depth;
//        Array<CAtomic<int64_t>> _cur_elements_per_depth;
//        Array<CAtomic<ContractionGroupID>> _offsets;
//
//        // This has to be a std::vector instead of an Array from array.h as constructing the array fails for some
//        // tbb internal reason
//        std::vector<tbb::concurrent_vector<ContractionGroupID>> _vecs;
//
//    };

    class DepthPriorityQueue {

        using depth_type = uint32_t;

    public:

        DepthPriorityQueue(const depth_type num_depths, std::vector<ContractionGroupID>&& max_elements_per_depth) :
            _num_depths(num_depths),
            _total_elements_per_depth(std::move(max_elements_per_depth)),
            _min_non_completed_depth(0),
            _num_finished_per_depth(_num_depths, CAtomic<ContractionGroupID>(0)),
            _completed(_num_depths, CAtomic<uint8_t>(uint8_t(false))),
            _queues(_num_depths, tbb::concurrent_bounded_queue<ContractionGroupID>()) {
          ASSERT(_num_depths == _total_elements_per_depth.size());
          // Reserve as much memory as elements per depth
          depth_type min_non_completed = 0;
          for (ContractionGroupID i = 0; i < _num_depths; ++i) {
            _queues[i].set_capacity(_total_elements_per_depth[i]);
            if (_total_elements_per_depth[i] == 0) {
              _completed[i].store(uint8_t(true), std::memory_order_relaxed);
              if (min_non_completed == i) {
                ++min_non_completed;
              }
            }
          }
          _min_non_completed_depth.store(min_non_completed, std::memory_order_relaxed);
        }


        void push(const ContractionGroupID id, const depth_type depth) {
          if (depth >= _num_depths || _completed[depth].load(std::memory_order_relaxed)) {
            ERROR("depth is unexpectedly too large or already completed!");
          }
          ASSERT(depth < _num_depths);
          ASSERT(!_completed[depth].load(std::memory_order_relaxed));
          if (id == invalidGroupID) {
            ERROR("Trying to push invalidGroupID to DepthPriorityQueue!");
          }
          bool pushed = _queues[depth].try_push(id);
          unused(pushed);
          ASSERT(pushed);
          if (!pushed) {
            ERROR("Pushes into depth PQ need to succeed!");
          }
        }

        bool try_pop(ContractionGroupID& dest) {

          if (_completed.size() != _num_depths || _queues.size() != _num_depths) {
            ERROR("Vectors in depth PQ have wrong size." << V(_completed.size()) << V(_queues.size()) << V(_num_depths));
          }

          depth_type current_min = _min_non_completed_depth.load(std::memory_order_relaxed);
          ASSERT(current_min <= _num_depths);
          if (current_min == _num_depths) {
            return false;
          }

          while (current_min < _num_depths) {

            if (_completed.size() != _num_depths || _queues.size() != _num_depths) {
              ERROR("Vectors in depth PQ have wrong size.");
            }

            bool depth_completed = _completed[current_min].load(std::memory_order_relaxed);
            if (!depth_completed) {
              if (_queues[current_min].try_pop(dest)) {
                // Success, dest holds a value that was popped
                ASSERT(dest != invalidGroupID);
                if (dest == invalidGroupID) {
                  ERROR("invalid group picked!");
                }
                return true;
              }
            }

            ++current_min;
          }

          // Failure, no value found
          return false;
        }

        void increment_finished(const depth_type depth) {
          bool completed_depth = false;
          increment_finished(depth, completed_depth);
        }

        void increment_finished(const depth_type depth, bool& completed_depth) {
          if (depth >= _num_depths) {
            ERROR("depth is unexpectedly too large!");
          }
          ASSERT(depth < _num_depths);
          ASSERT(!_completed[depth].load(std::memory_order_relaxed));
          ContractionGroupID num_finished = _num_finished_per_depth[depth].add_fetch(1, std::memory_order_relaxed);
          ASSERT(num_finished <= _total_elements_per_depth[depth]);
          if (num_finished == _total_elements_per_depth[depth]) {
            completed_depth = true;
            auto exp = uint8_t(false);
            bool changed_completed_flag = _completed[depth].compare_exchange_strong(exp, uint8_t(true), std::memory_order_relaxed);
            ASSERT(!exp);
            if (!changed_completed_flag) {
              ERROR("Multiple threads setting completed flag of same depth in DepthPriorityQueue!");
            }

            // Update min ptr if min was just completed
            if (depth == _min_non_completed_depth.load(std::memory_order_relaxed)) {
              depth_type next_uncompleted = depth + 1;
              while (next_uncompleted < _num_depths && _completed[next_uncompleted].load(std::memory_order_relaxed)) {
                ++next_uncompleted;
              }
              depth_type current_min = depth;
              bool changed_min = _min_non_completed_depth.compare_exchange_strong(current_min, next_uncompleted, std::memory_order_relaxed);
              unused(changed_min);
              ASSERT(changed_min && (current_min == depth));
              if (!changed_min) {
                ERROR("Multiple threads updating min pointer in DepthPriorityQueue simultaneously!");
              }
            }
          } else {
            completed_depth = false;
          }
        }

        int64_t unsafe_size() const {
          int64_t size = 0;
          for(depth_type i = 0; i < _num_depths; ++i) {
            size += _queues[i].size();
          }
          return size;
        }

        bool unsafe_empty() const {
          return (unsafe_size() == 0);
        }

        // ! Clears the bucket contents while not changing number of depths or maximum number of elements per depth.
        void reset() {
          depth_type min_non_completed = 0;
          for(depth_type i = 0; i < _num_depths; ++i) {
            _queues[i].clear();
            if (_total_elements_per_depth[i] != 0) {
              _completed[i].store(uint8_t(false), std::memory_order_relaxed);
            } else if (min_non_completed == i) {
              ++min_non_completed;
            }
          }
          _min_non_completed_depth.store(min_non_completed, std::memory_order_relaxed);
        }

        bool allDepthsCompleted() const {
          return _min_non_completed_depth.load(std::memory_order_relaxed) == _num_depths;
        }

        void memoryConsumption(utils::MemoryTreeNode* parent) const {
          ASSERT(parent);

          parent->addChild("Total Elements per Depth", _total_elements_per_depth.size() * sizeof(ContractionGroupID));
          parent->addChild("Number Finished per Depth", _num_finished_per_depth.size() * sizeof(CAtomic<ContractionGroupID>));
          parent->addChild("Depth Completed Flags", _completed.size() * sizeof(CAtomic<uint8_t>));
        }

    private:

        const depth_type _num_depths;
        const std::vector<ContractionGroupID> _total_elements_per_depth;

        CAtomic<depth_type> _min_non_completed_depth;
        std::vector<CAtomic<ContractionGroupID>> _num_finished_per_depth;
        std::vector<CAtomic<uint8_t>> _completed;

        // This has to be a std::vector instead of an Array from array.h as constructing the array fails for some
        // tbb internal reason
        std::vector<tbb::concurrent_bounded_queue<ContractionGroupID>> _queues;

    };

} // end namespace


#endif //KAHYPAR_DEPTH_PRIORITY_QUEUE_H

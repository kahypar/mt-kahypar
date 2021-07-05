//
// Created by mlaupichler on 05.07.21.
//

#ifndef KAHYPAR_DEPTH_PRIORITY_QUEUE_H
#define KAHYPAR_DEPTH_PRIORITY_QUEUE_H

namespace mt_kahypar::ds {

    class DepthPriorityQueue {

        using depth_type = uint32_t;

    public:

        DepthPriorityQueue(const depth_type num_depths, Array<ContractionGroupID>&& num_elements_per_depth) :
          _num_depths(num_depths),
          _num_elements_per_depth(std::move(num_elements_per_depth)),
          _offsets(_num_depths, CAtomic<ContractionGroupID>(0)),
          _vecs(_num_depths) {
          for (ContractionGroupID i = 0; i < _num_depths; ++i) {
            _vecs[i].reserve(_num_elements_per_depth[i]);
          }
        }

        void push(const ContractionGroupID id, const depth_type depth) {
          ASSERT(depth < _num_depths);
          _vecs[depth].push_back(id);
        }

        bool tryPop(const depth_type depth, ContractionGroupID& destination) {
          ASSERT(depth < _num_depths);

        }

    private:

        const depth_type _num_depths;
        const Array<ContractionGroupID> _num_elements_per_depth;
        Array<CAtomic<ContractionGroupID>> _offsets;
        Array<tbb::concurrent_vector<ContractionGroupID>> _vecs;

    };

} // end namespace


#endif //KAHYPAR_DEPTH_PRIORITY_QUEUE_H

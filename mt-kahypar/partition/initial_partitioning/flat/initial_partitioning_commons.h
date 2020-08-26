#pragma once

#include "kahypar/datastructure/kway_priority_queue.h"
#include "tbb/enumerable_thread_specific.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "kahypar/datastructure/fast_reset_flag_array.h"

namespace mt_kahypar {

using KWayPriorityQueue = kahypar::ds::KWayPriorityQueue<HypernodeID, Gain, std::numeric_limits<Gain>, true>;
using ThreadLocalKWayPriorityQueue = tbb::enumerable_thread_specific<KWayPriorityQueue>;

using ThreadLocalFastResetFlagArray = tbb::enumerable_thread_specific<kahypar::ds::FastResetFlagArray<> >;

}
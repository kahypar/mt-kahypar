//
// Created by mlaupichler on 19.04.21.
//

#include "asynch_contraction_pool.h"
#include "mock_group_hierarchy.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace mt_kahypar::ds {

template<typename GroupHierarchy>
uint32_t mt_kahypar::ds::SequentialContractionGroupPool<GroupHierarchy>::getNumActive() const {
    return _active_ids.size();
}

template<typename GroupHierarchy>
const mt_kahypar::ds::ContractionGroup &
mt_kahypar::ds::SequentialContractionGroupPool<GroupHierarchy>::group(mt_kahypar::ds::ContractionGroupID id) const {
    return _hierarchy->group(id);
}

template<typename GroupHierarchy>
mt_kahypar::ds::ContractionGroupID mt_kahypar::ds::SequentialContractionGroupPool<GroupHierarchy>::pickAnyActiveID() {
    ASSERT(!_active_ids.empty());
    auto picked = _active_ids.front();
    ASSERT(isActive(picked));
    _active[picked] = false;
    _active_ids.pop();
    return picked;
}

template<typename GroupHierarchy>
void mt_kahypar::ds::SequentialContractionGroupPool<GroupHierarchy>::activateSuccessors(mt_kahypar::ds::ContractionGroupID id) {
    ASSERT(!isActive(id));
    ASSERT(!_successors_activated[id]);

    auto succs = _hierarchy->successors(id);
    for (auto s : succs) {
        activate(s);
    }

    _successors_activated[id] = true;
}

template<typename GroupHierarchy>
bool mt_kahypar::ds::SequentialContractionGroupPool<GroupHierarchy>::hasActive() const {
    return !(_active_ids.empty());
}

template<typename GroupHierarchy>
mt_kahypar::ds::BlockedGroupIDIterator mt_kahypar::ds::SequentialContractionGroupPool<GroupHierarchy>::all() const {
    return _hierarchy->all();
}

template<typename GroupHierarchy>
uint32_t mt_kahypar::ds::SequentialContractionGroupPool<GroupHierarchy>::getNumTotal() const {
    return _hierarchy->getNumGroups();
}

template<typename GroupHierarchy>
size_t mt_kahypar::ds::SequentialContractionGroupPool<GroupHierarchy>::getVersion() const {
    return _hierarchy->getVersion();
}

template<typename GroupHierarchy>
void mt_kahypar::ds::SequentialContractionGroupPool<GroupHierarchy>::reactivate(mt_kahypar::ds::ContractionGroupID id) {

    activate(id);

//    // Shuffle on reactivate to break LIFO picking in order to provide other possible picked ID
//    utils::Randomize::instance().shuffleVector(
//            _active_ids, 0UL, _active_ids.size(), sched_getcpu());
}

template<typename GroupHierarchy>
bool mt_kahypar::ds::SequentialContractionGroupPool<GroupHierarchy>::isActive(ContractionGroupID id) const {
    ASSERT(id < getNumTotal());
    return _active[id];
}

template<typename GroupHierarchy>
void mt_kahypar::ds::SequentialContractionGroupPool<GroupHierarchy>::activate(mt_kahypar::ds::ContractionGroupID id) {
    ASSERT(!isActive(id));
    ASSERT(!_successors_activated[id]);
    _active[id] = true;
    _active_ids.push(id);
}

// explicitly instantiate so the compiler can generate them when compiling this cpp file
template class SequentialContractionGroupPool<UncontractionGroupTree>;
template class SequentialContractionGroupPool<MockGroupHierarchy>;

}



//
// Created by mlaupichler on 19.04.21.
//

#include "asynch_contraction_pool.h"

uint32_t mt_kahypar::ds::SequentialContractionGroupPool::getNumActive() const {
    return _active_ids.size();
}

const mt_kahypar::ds::ContractionGroup &
mt_kahypar::ds::SequentialContractionGroupPool::group(mt_kahypar::ds::ContractionGroupID id) const {
    return _hierarchy->group(id);
}

mt_kahypar::ds::ContractionGroupID mt_kahypar::ds::SequentialContractionGroupPool::pickAnyActiveID() {
    ASSERT(!_active_ids.empty());
    // This implementation is LIFO
    auto picked = _active_ids.back();
    ASSERT(isActive(picked));
    _active[picked] = false;
    _active_ids.pop_back();
    return picked;
}

void mt_kahypar::ds::SequentialContractionGroupPool::activateSuccessors(mt_kahypar::ds::ContractionGroupID id) {
    ASSERT(!isActive(id));
    ASSERT(!_successors_activated[id]);

    auto succs = _hierarchy->successors(id);
    for (auto s : succs) {
        activate(s);
    }

    _successors_activated[id] = true;
}

bool mt_kahypar::ds::SequentialContractionGroupPool::hasActive() const {
    return !(_active_ids.empty());
}

mt_kahypar::ds::BlockedGroupIDIterator mt_kahypar::ds::SequentialContractionGroupPool::all() const {
    return _hierarchy->all();
}

uint32_t mt_kahypar::ds::SequentialContractionGroupPool::getNumTotal() const {
    return _hierarchy->getNumGroups();
}

size_t mt_kahypar::ds::SequentialContractionGroupPool::getVersion() const {
    return _hierarchy->getVersion();
}

void mt_kahypar::ds::SequentialContractionGroupPool::reactivate(mt_kahypar::ds::ContractionGroupID id) {

    activate(id);

    // Shuffle on reactivate to break LIFO picking in order to provide other possible picked ID
    utils::Randomize::instance().shuffleVector(
            _active_ids, 0UL, _active_ids.size(), sched_getcpu());
}

bool mt_kahypar::ds::SequentialContractionGroupPool::isActive(ContractionGroupID id) const {
    ASSERT(id < getNumTotal());
    return _active[id];
}

void mt_kahypar::ds::SequentialContractionGroupPool::activate(mt_kahypar::ds::ContractionGroupID id) {
    ASSERT(!isActive(id));
    ASSERT(!_successors_activated[id]);
    _active[id] = true;
    _active_ids.push_back(id);
}


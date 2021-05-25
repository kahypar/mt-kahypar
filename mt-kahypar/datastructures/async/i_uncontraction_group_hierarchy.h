//
// Created by mlaupichler on 28.04.21.
//

#ifndef KAHYPAR_I_UNCONTRACTION_GROUP_HIERARCHY_H
#define KAHYPAR_I_UNCONTRACTION_GROUP_HIERARCHY_H

#include "gmock/gmock.h"

#include "async_common.h"

namespace mt_kahypar::ds {

    /// Pure virtual interface to define functions for a hierarchy of UncontractionGroups in a version.
    /// In this hierarchy an UncontractionGroup is supposed to be uncontracted only once its predecessor has been uncontracted.
    /// This expresses the minimal order between uncontraction groups based on the order of contractions.
    class IUncontractionGroupHierarchy {

    public:

        static constexpr size_t kInvalidVersion = std::numeric_limits<size_t>::max();

        virtual size_t getVersion() const = 0;
        virtual uint32_t getNumGroups() const = 0;
        virtual const ContractionGroup& group(ContractionGroupID id) const = 0;
        virtual ContractionGroupID predecessor(ContractionGroupID id) const = 0;
        virtual ContractionGroupIDIteratorRange successors(ContractionGroupID id) const = 0;
        virtual ContractionGroupIDIteratorRange roots() const = 0;

        virtual BlockedGroupIDIterator all() const = 0;

        virtual ~IUncontractionGroupHierarchy() = default;
    };





}

#endif //KAHYPAR_I_UNCONTRACTION_GROUP_HIERARCHY_H

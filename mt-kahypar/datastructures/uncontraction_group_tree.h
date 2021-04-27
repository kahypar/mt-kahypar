//
// Created by mlaupichler on 24.04.21.
//

#ifndef KAHYPAR_UNCONTRACTION_GROUP_TREE_H
#define KAHYPAR_UNCONTRACTION_GROUP_TREE_H

#include <utility>
#include <boost/range/counting_range.hpp>

#include "contraction_tree.h"

namespace mt_kahypar::ds {

    using ContractionGroupID = HypernodeID;
    using ContractionGroupIDIteratorRange = IteratorRange<parallel::scalable_vector<ContractionGroupID>::const_iterator>;
    using BlockedGroupIDIterator = tbb::blocked_range<ContractionGroupID>;
    using Contraction = Memento;
    typedef std::vector<Contraction>::iterator ContractionIterator;

    /// Represents a group of (un-)contractions that have been contracted simultaneously with the same representative.
    /// This means the contractions in a group have to be uncontracted simultaneously as well in the sense that no local
    /// refinement is performed in between uncontractions of the same group as those states may be inconsistent.
    /// Groups are expected to be small (constant against the number of uncontractions in a graph version).
    class ContractionGroup {

    private:

        std::vector<Contraction> _contractions;
        HypernodeID _representative;

        /// Check that group is not empty and that all contractions have the same representative
        bool sanityCheck() {
            if (empty()) return false;
            auto repr = _contractions.at(0).u;
            return std::all_of(_contractions.begin(),_contractions.end(),[repr](Contraction c){return c.u == repr;});
        };

    public:

        bool empty() const {return _contractions.empty();};
        uint64_t size() const {return _contractions.size();};

        Contraction at(int i) const {return _contractions.at(i);};

//        ContractionGroup() = default;
        ContractionGroup(std::initializer_list<Contraction> init) : _contractions(init) {
            ASSERT(sanityCheck());
            _representative = _contractions.front().u;
        };

        explicit ContractionGroup(std::vector<Contraction> init) : _contractions(std::move(init)) {
            ASSERT(sanityCheck());
            _representative = _contractions.front().u;
        }

        ContractionGroup(ContractionGroup& other) = default;
        ContractionGroup(const ContractionGroup& other) = default;

        HypernodeID getRepresentative() const;

        /// This function is linear in the number of contractions.
        bool contains(Contraction contraction) const;

        /// This function is quadratic in the group size.
        bool operator==(const ContractionGroup &rhs) const;

        /// This function is quadratic in the group size.
        bool operator!=(const ContractionGroup &rhs) const;

        ContractionIterator begin() {return _contractions.begin();};
        ContractionIterator end() {return _contractions.end();};
        auto cbegin() const {return _contractions.begin();};
        auto cend() const {return _contractions.end();};
        auto begin() const {return _contractions.begin();};
        auto end() const {return _contractions.end();};

        void debugPrint() const;
    };

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

        virtual ~IUncontractionGroupHierarchy() = 0;
    };



/// Data structure that contains a tree specifying the order in which contraction groups of a particular version have to
/// be uncontracted due to contraction order between parent and child in the contraction tree
/// ("vertical order") and due to contraction order between siblings in the contraction tree ("horizontal order"). A group
/// can have at most as many vertical children as it has members (for each group member the group of its children that was
/// contracted last). A group can have at most one horizontal child (the group of siblings in the contraction tree that
/// was contracted right before it).
/// A contraction group is a set of contractions with the same representative that were contracted simultaneously
/// (i.e. their contraction intervals have transitive overlap) which means they have to be uncontracted simultaneously
/// as well. Starting with the root groups for a version, a group can only be uncontracted once its parent group in the
/// GroupUncontractionForest has been uncontracted.
    class UncontractionGroupTree : public IUncontractionGroupHierarchy {

        using ContractionInterval = typename ContractionTree::Interval;

    private:
        struct GroupNode {
        public:
            GroupNode(ContractionGroup& group, size_t version, ContractionGroupID parent) : _parentGroup(parent),
                                                                                    _version(version), _group(group) {}

            GroupNode(GroupNode& other) = default;
            GroupNode(const GroupNode& other) = default;

            GroupNode& operator=(const GroupNode& other) = default;

            ContractionGroupID getParentGroup() const {
                return _parentGroup;
            }

            size_t getVersion() const {
                return _version;
            }

            const ContractionGroup &getGroup() const {
                return _group;
            }

        private:
            // GroupNodeID of the nodes parent group
            ContractionGroupID _parentGroup;
            // version of the group (has to be the same for all GroupNodes in a tree)
            size_t _version;
            // the contained group
            ContractionGroup _group;
        };

    public:

        UncontractionGroupTree(ContractionTree &contractionTree, size_t version);

        ~UncontractionGroupTree() override {
            freeInternalData();
        }

        const ContractionGroup &group(ContractionGroupID id) const override {
            ASSERT(id < _num_group_nodes);
            return _tree[id].getGroup();
        }

        ContractionGroupIDIteratorRange successors(ContractionGroupID id) const override{
            ASSERT(id < _num_group_nodes);
            return ContractionGroupIDIteratorRange(
                    _incidence_array.cbegin() + _out_degrees[id],
                    _incidence_array.cbegin() + _out_degrees[id + 1]);
        }

        ContractionGroupID predecessor(ContractionGroupID id) const override {
            ASSERT(id < _num_group_nodes);
            return _tree[id].getParentGroup();
        }

//        size_t version(GroupNodeID id) const {
//            ASSERT(id < _num_group_nodes);
//            return _tree[id].getVersion();
//        }

        size_t getVersion() const override {
            return _version;
        }

        uint32_t getNumGroups() const override {
            return _num_group_nodes;
        }

        ContractionGroupIDIteratorRange roots() const override {
            return mt_kahypar::ds::ContractionGroupIDIteratorRange(_roots.cbegin(), _roots.cend());
        }

        BlockedGroupIDIterator all() const override {
            ASSERT(_num_group_nodes == _tree.size());
            return BlockedGroupIDIterator (0,_num_group_nodes);
        }

    private:

        void freeInternalData();

        void insertRootBranchesForVersion();

        void insertHorizontalBranch(IteratorRange<ContractionTree::ChildIterator> childrenIt, ContractionGroupID parentGroup,
                                    HypernodeID parentHypernode);

        void insertRootHorizontalBranch(IteratorRange<ContractionTree::ChildIterator> childrenIt,
                                        HypernodeID parentHypernode);

        /**
         * Inserts the given ContractionGroup into the UncontractionGroupTree as a child of the group at GroupNodeID parent
         * @param group group to insert
         * @param parent the parent GroupNode
         * @param hasHorizontalChild indicates whether the group to insert has horizontal children,
         *  i.e. siblings in the ContractionTree that were contracted earlier than this group but still within the same version
         * @return The GroupNodeID of the newly inserted group
         */
        ContractionGroupID insertGroup(ContractionGroup group, ContractionGroupID parent, bool hasHorizontalChild);

        bool doIntervalsIntersect(const ContractionInterval &i1, const ContractionInterval &i2);

        // Contraction tree that this is based on
        ContractionTree &_contraction_tree;
        uint32_t _num_group_nodes;

        parallel::scalable_vector <GroupNode> _tree;
        parallel::scalable_vector <ContractionGroupID> _roots;

        parallel::scalable_vector <uint32_t> _current_child_offsets;
        parallel::scalable_vector <uint32_t> _out_degrees;
        parallel::scalable_vector <ContractionGroupID> _incidence_array;

        // The hypergraph version that this UncontractionGroupTree is for
        size_t _version;

    };

} // namespace mt_kahypar::ds
#endif //KAHYPAR_UNCONTRACTION_GROUP_TREE_H

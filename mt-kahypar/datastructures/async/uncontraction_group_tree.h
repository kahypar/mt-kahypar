//
// Created by mlaupichler on 24.04.21.
//

#ifndef KAHYPAR_UNCONTRACTION_GROUP_TREE_H
#define KAHYPAR_UNCONTRACTION_GROUP_TREE_H

#include <mt-kahypar/datastructures/sparse_map.h>
#include "mt-kahypar/datastructures/contraction_tree.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/async/async_common.h"

namespace mt_kahypar::ds {



/// Data structure that contains a tree specifying the order in which contraction groups of a particular version have to
/// be uncontracted due to contraction order between parent and child in the contraction tree
/// ("vertical order") and due to contraction order between siblings in the contraction tree ("horizontal order"). A group
/// can have at most as many vertical children as it has members (for each group member the group of its children that was
/// contracted last). A group can have at most one horizontal child (the group of siblings in the contraction tree that
/// was contracted right before it).
/// A contraction group is a set of contractions with the same representative that were contracted simultaneously
/// (i.e. their contraction intervals have transitive overlap) which means they have to be uncontracted simultaneously
/// as well. Starting with the root groups for a version, a group can only be uncontracted once its parent group in the
/// GroupUncontractionTree has been uncontracted.
    class UncontractionGroupTree {

        using ContractionInterval = typename ContractionTree::Interval;

    private:
        struct GroupNode {
        public:
            GroupNode(std::vector<Contraction>& contractions, size_t version, ContractionGroupID parent, size_t depth) : _parentGroup(parent),
                                                                                    _version(version), _group(contractions), _depth(depth) {}

            GroupNode(GroupNode& other) :
                _parentGroup(other._parentGroup),
                _version(other._version),
                _group(std::vector<Contraction>(other._group.begin(),other._group.end())),
                _depth(other._depth){}

            GroupNode(const GroupNode& other) :
                    _parentGroup(other._parentGroup),
                    _version(other._version),
                    _group(std::vector<Contraction>(other._group.begin(),other._group.end())),
                    _depth(other._depth) {}

            GroupNode(GroupNode&& other) = delete;
            GroupNode& operator=(const GroupNode& other) = delete;

            ContractionGroupID getParentGroup() const {
                return _parentGroup;
            }

            size_t getVersion() const {
                return _version;
            }

            const ContractionGroup &getGroup() const {
                return _group;
            }

            const HypernodeID& getDepth() const {
                return _depth;
            }

        private:
            // GroupNodeID of the nodes parent group
            ContractionGroupID _parentGroup;
            // version of the group (has to be the same for all GroupNodes in a tree)
            size_t _version;
            // the contained group
            ContractionGroup _group;
            // depth in the UncontractionGroupTree (where roots have depth 0)
            HypernodeID _depth;
        };

    public:

        UncontractionGroupTree(ContractionTree &contractionTree, size_t version);

        ~UncontractionGroupTree() {
            freeInternalData();
        }

        const ContractionGroup &group(ContractionGroupID id) const {
            ASSERT(id < _num_group_nodes);
            return _tree[id].getGroup();
        }

        const HypernodeID& depth(ContractionGroupID id) const {
            ASSERT(id < _num_group_nodes);
            return _tree[id].getDepth();
        }
//
//        const size_t& node_depth(HypernodeID node_id) const {
//            ASSERT(node_id < _contraction_tree.num_hypernodes());
//            ContractionGroupID group_id = (*_node_to_group_map)[node_id];
//            ASSERT(group_id < _num_group_nodes);
//            return _tree[group_id].getDepth();
//        }

        ContractionGroupIDIteratorRange successors(ContractionGroupID id) const {
            ASSERT(id < _num_group_nodes);
            return ContractionGroupIDIteratorRange(
                    _incidence_array.cbegin() + _out_degrees[id],
                    _incidence_array.cbegin() + _out_degrees[id + 1]);
        }

        ContractionGroupID numSuccessors(ContractionGroupID id) const {
            return _out_degrees[id+1] - _out_degrees[id];
        }

        ContractionGroupID predecessor(ContractionGroupID id) const {
            ASSERT(id < _num_group_nodes);
            return _tree[id].getParentGroup();
        }

        size_t getVersion() const {
            return _version;
        }

        uint32_t getNumGroups() const {
            return _num_group_nodes;
        }

        ContractionGroupIDIteratorRange roots() const {
            return mt_kahypar::ds::ContractionGroupIDIteratorRange(_roots.cbegin(), _roots.cend());
        }

        const BlockedGroupIDIterator all() const {
            return BlockedGroupIDIterator(0, _num_group_nodes);
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
         * @param contractions group to insert
         * @param parent the parent GroupNode
         * @param hasHorizontalChild indicates whether the group to insert has horizontal children,
         *  i.e. siblings in the ContractionTree that were contracted earlier than this group but still within the same version
         * @return The GroupNodeID of the newly inserted group
         */
        ContractionGroupID insertGroup(std::vector<Contraction> &contractions, ContractionGroupID parent, bool hasHorizontalChild);

        bool doIntervalsIntersect(const ContractionInterval &i1, const ContractionInterval &i2);

        // Contraction tree that this is based on
        ContractionTree &_contraction_tree;
        ContractionGroupID _num_group_nodes;
        HypernodeID _num_contained_contracted_nodes;

        parallel::scalable_vector <GroupNode> _tree;
        parallel::scalable_vector <ContractionGroupID> _roots;

        parallel::scalable_vector <uint32_t> _current_child_offsets;
        parallel::scalable_vector <uint32_t> _out_degrees;
        parallel::scalable_vector <ContractionGroupID> _incidence_array;

        // The hypergraph version that this UncontractionGroupTree is for
        size_t _version;

        // ! Mapping from HypernodeID to the ContractionGroupID of the group that the contraction (Memento) of the
        // ! respective Hypernode is in. Used to find the depth of the contraction of a  node in an UncontractionGroupTree.
        // ! (ContractionGroupIDs are only unique within a version and the respective hierarchy but combined with
        // ! information about which version a node has been contracted in, the right hierarchy and depth can be found.)
//        std::unique_ptr<parallel::scalable_vector<ContractionGroupID>> _node_to_group_map;


    };

} // namespace mt_kahypar::ds
#endif //KAHYPAR_UNCONTRACTION_GROUP_TREE_H

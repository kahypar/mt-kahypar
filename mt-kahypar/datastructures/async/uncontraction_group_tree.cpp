//
// Created by mlaupichler on 24.04.21.
//

#include "uncontraction_group_tree.h"

namespace mt_kahypar::ds {

    // Static initialization
//    std::unique_ptr<parallel::scalable_vector<ContractionGroupID>> UncontractionGroupTree::_node_to_group_map;

    void UncontractionGroupTree::freeInternalData() {
        if (_num_group_nodes > 0 ) {
            parallel::parallel_free(_tree, _roots, _out_degrees,_current_child_offsets, _incidence_array);
        }
        _num_group_nodes = 0;
    }

    UncontractionGroupTree::UncontractionGroupTree(const ContractionTree &contractionTree, const size_t version)
            : _contraction_tree(contractionTree),
              _num_group_nodes(0),
              _num_contained_contracted_nodes(0),
              _version(version),
              _last_uncontraction_group_in_version(contractionTree.num_hypernodes(), invalidGroupID) {
        ASSERT(_contraction_tree.isFinalized());

        _out_degrees.push_back(0);

        insertRootBranchesForVersion();

        if(_num_group_nodes == 0) return;

        ContractionGroupID currentParentID = 0;
        while (currentParentID != _num_group_nodes) {
            const GroupNode& currentParentNode = _tree[currentParentID];
            // explore one step vertically and then the horizontal branch
            for (const auto& member : currentParentNode.getGroup()) {
                insertHorizontalBranch(_contraction_tree.childs(member.v),currentParentID,member.v);
            }
            ++currentParentID;
        }

//        size_t max_depth = 0;
//        for (const auto& node : _tree) {
//            if (node.getDepth() > max_depth) max_depth = node.getDepth();
//        }

//        std::cout << "Version " << _version << " has a max depth of " << max_depth << std::endl;

//        _node_to_group_map = std::make_unique<parallel::scalable_vector<ContractionGroupID>>(_contraction_tree.num_hypernodes());
//        for (ContractionGroupID id = 0; id < _num_group_nodes; ++id) {
//            for (const auto& memento : group(id)) {
//                (*_node_to_group_map)[memento.v] = id;
//            }
//        }
    }

    void UncontractionGroupTree::insertRootBranchesForVersion() {

        auto versionRoots = _contraction_tree.roots_of_version(_version);

        for (auto root: versionRoots) {
            auto it = _contraction_tree.childs(root);
            insertRootHorizontalBranch(it,root);
        }
    }

    bool UncontractionGroupTree::doIntervalsIntersect(const ContractionInterval& i1, const ContractionInterval& i2) {
        if (i1.start == kInvalidHypernode || i2.start == kInvalidHypernode) {
            return false;
        }
        return (i1.start <= i2.end && i1.end >= i2.end) ||
               (i2.start <= i1.end && i2.end >= i1.end);
    }

    ContractionGroupID UncontractionGroupTree::insertGroup(std::vector<Contraction> &&contractions, const ContractionGroupID parent, const bool hasHorizontalChild) {

        ContractionGroupID newID = _num_group_nodes;
        ++_num_group_nodes;
        _num_contained_contracted_nodes += contractions.size();

        // Add new node into the tree at newID (i.e. its ID is newID)
        size_t depth = (newID == parent) ? 0 : _tree[parent].getDepth() + 1;
        _tree.emplace_back(std::move(contractions), _version, parent, depth);
        const auto& newGroupNode = _tree.back();
        _current_child_offsets.emplace_back(0);
        ASSERT(_tree.size() == newID + 1);
        ASSERT(_current_child_offsets.size() == newID + 1);

        // Add it into the incidence array of the parent at the right offset and increase the child offset by 1.
        // Roots do not have parents though so ignore them.
        if (newID != parent) {
            ASSERT(_out_degrees.at(parent)+_current_child_offsets.at(parent) < _out_degrees.at(parent+1));
            _incidence_array[_out_degrees.at(parent)+_current_child_offsets.at(parent)] = newID;
            ++_current_child_offsets[parent];
        }

        auto numVertical = std::count_if(newGroupNode.getGroup().begin(), newGroupNode.getGroup().end(), [&](Memento member) {
            // Search for any children in the contraction tree that have this version. If any exist for member.v then
            // member.v contributes to the number of vertical children this group has with one vertical child group
            // (the group of the child in the contraction tree that was contracted last)
            auto childrenInContractionTree = _contraction_tree.childs(member.v);
            return std::any_of(childrenInContractionTree.begin(),childrenInContractionTree.end(),[&](auto child) {
                return _contraction_tree.version(child) == _version;
            });
        });
        auto numChildren = numVertical + (hasHorizontalChild? 1 : 0);
        _out_degrees.push_back(_out_degrees.back() + numChildren);
        ASSERT(_out_degrees.size() == newID + 2);
        _incidence_array.resize(_incidence_array.size() + numChildren); // REVIEW do one resize at the end?

        // If new is its own parent then it is a root
        if(newID == parent) {
            _roots.push_back(newID);
        }

        // If this group does not have a horizontal child, it is the last group for this representative in the version
        if (!hasHorizontalChild) {
            const auto& representative = newGroupNode.getGroup().getRepresentative();
            ASSERT(representative < _last_uncontraction_group_in_version.size());
            _last_uncontraction_group_in_version[representative] = newID;
        }

        return newID;
    }

    void UncontractionGroupTree::insertHorizontalBranch(IteratorRange<ContractionTree::ChildIterator> childrenIt,
                                                        const ContractionGroupID parentGroup,
                                                        const HypernodeID parentHypernode) {

        auto current = childrenIt.begin();
        auto end = childrenIt.end();
        while ( current != end && _contraction_tree.version(*current) != _version ) {
            ++current;
        }
        // Return if no children of this version in the given range
        if (current == end) return;

        auto moreChildrenInThisVersion = [&](ContractionTree::ChildIterator& it) {
            return it != end && _contraction_tree.version(*it) == _version;
        };

        std::vector<Contraction> inCurrentGroup;
        ContractionInterval current_ival = _contraction_tree.interval(*current);
        inCurrentGroup.push_back(Contraction {parentHypernode, *current});
        auto currentParentGroup = parentGroup;
        ++current;
        while (moreChildrenInThisVersion(current)) {
            auto sibling = *current;
            ContractionInterval sibling_ival = _contraction_tree.interval(sibling);

            ASSERT(_contraction_tree.parent(sibling) == parentHypernode);
            if (doIntervalsIntersect(current_ival,sibling_ival)) {
                // Add sibling to group and continue with next sibling (in the contraction tree)
                inCurrentGroup.push_back(Contraction {parentHypernode, sibling});
                current_ival.start = std::min(current_ival.start, sibling_ival.start);
                current_ival.end = std::max(current_ival.end, sibling_ival.end);
            } else {
                // Group is finished as interval does not intersect anymore (and intervals are ordered)
                // Add it to the UncontractionGroupTree and continue with its next sibling in the contraction tree
                // (its horizontal child in the UncontractionGroupTree)
                auto idOfInserted = insertGroup(std::move(inCurrentGroup),currentParentGroup,true);

                // And reset current group/parent
                current_ival = sibling_ival;
                // inCurrentGroup is in valid but unspecified state after std::move, so clear() which has no
                // preconditions is safe to use and afterwards it becomes a regular empty vector
                inCurrentGroup.clear();
                inCurrentGroup.push_back(Contraction {parentHypernode, sibling});
                currentParentGroup = idOfInserted;
            }
            ++current;
        }

        // End of children in this version (that intersect with last contraction group) has been reached so finish up by inserting remaining group
        ASSERT(!moreChildrenInThisVersion(current));
        if (!inCurrentGroup.empty()) {
            insertGroup(std::move(inCurrentGroup),currentParentGroup,false);
        }

    }

    void UncontractionGroupTree::insertRootHorizontalBranch(IteratorRange<ContractionTree::ChildIterator> childrenIt, const HypernodeID parentHypernode) {
        // Same as any horizontal branch, except the parentGroupID is the next free available ID, so the first group in the horizontal
        // branch becomes its own parent (i.e. a root)
        auto num_roots = _roots.size();
        unused(num_roots);
        insertHorizontalBranch(childrenIt,_num_group_nodes,parentHypernode);
        ASSERT(_roots.size() == num_roots+1);
    }

} // namespace mt_kahypar::ds
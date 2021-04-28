//
// Created by mlaupichler on 24.04.21.
//

#include "uncontraction_group_tree.h"

namespace mt_kahypar::ds {

    void UncontractionGroupTree::freeInternalData() {
        if (_num_group_nodes > 0 ) {
            parallel::parallel_free(_tree, _roots, _out_degrees,_current_child_offsets, _incidence_array);
        }
        _num_group_nodes = 0;
    }

    UncontractionGroupTree::UncontractionGroupTree(ContractionTree &contractionTree, size_t version)
            : _contraction_tree(contractionTree),
              _num_group_nodes(0),
              _version(version) {
        ASSERT(_contraction_tree.isFinalized());

        _out_degrees.push_back(0);

        insertRootBranchesForVersion();

        if(_num_group_nodes == 0) return;

        ContractionGroupID currentParentID = 0;
        while (currentParentID != _num_group_nodes) {
            GroupNode currentParentNode = _tree[currentParentID];
            // explore one step vertically and then the horizontal branch
            for (auto member : currentParentNode.getGroup()) {
                insertHorizontalBranch(_contraction_tree.childs(member.v),currentParentID,member.v);
            }
            ++currentParentID;
        }

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

    ContractionGroupID UncontractionGroupTree::insertGroup(ContractionGroup group, ContractionGroupID parent, bool hasHorizontalChild) {

        ContractionGroupID newID = _num_group_nodes;
        ++_num_group_nodes;

        // Add new node into the tree at newID (i.e. its ID is newID)
        _tree.emplace_back(group,_version,parent);
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

        auto numVertical = std::count_if(group.begin(),group.end(),[&](Memento member) {
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
        _incidence_array.resize(_incidence_array.size() + numChildren);

        // If new is its own parent then it is a root
        if(newID == parent) {
            _roots.push_back(newID);
        }

        return newID;
    }

    void UncontractionGroupTree::insertHorizontalBranch(IteratorRange<ContractionTree::ChildIterator> childrenIt,
                                                        ContractionGroupID parentGroup,
                                                        HypernodeID parentHypernode) {

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
                auto group = ContractionGroup(inCurrentGroup);
                auto idOfInserted = insertGroup(group,currentParentGroup,true);

                // And reset current group/parent
                current_ival = sibling_ival;
                inCurrentGroup.clear();
                inCurrentGroup.push_back(Contraction {parentHypernode, sibling});
                currentParentGroup = idOfInserted;
            }
            ++current;
        }

        // End of children in this version (that intersect with last contraction group) has been reached so finish up by inserting remaining group
        ASSERT(!moreChildrenInThisVersion(current));
        if (!inCurrentGroup.empty()) {
            auto group = ContractionGroup(inCurrentGroup);
            insertGroup(group,currentParentGroup,false);
        }

    }

    void UncontractionGroupTree::insertRootHorizontalBranch(IteratorRange<ContractionTree::ChildIterator> childrenIt,  HypernodeID parentHypernode) {
        // Same as any horizontal branch, except the parentGroupID is the next free available ID, so the first group in the horizontal
        // branch becomes its own parent (i.e. a root)
        auto num_roots = _roots.size();
        insertHorizontalBranch(childrenIt,_num_group_nodes,parentHypernode);
        ASSERT(_roots.size() == num_roots+1);
    }

} // namespace mt_kahypar::ds
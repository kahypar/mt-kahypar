/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <map>
#include <memory>
#include <string>
#include <functional>

namespace mt_kahypar::utils {

enum class OutputType : uint8_t {
  BYTES = 0,
  KILOBYTE = 1,
  MEGABYTE = 2,
  PERCENTAGE = 3
};

class MemoryTreeNode {

 using map_type = std::map<std::string, std::unique_ptr<MemoryTreeNode>>;

 public:
  static constexpr int MAX_LINE_LENGTH = 45;
  static char LINE_PREFIX[];
  static constexpr size_t LINE_PREFIX_LENGTH = 3;

  MemoryTreeNode(const std::string& name, const OutputType& output_type = OutputType::MEGABYTE) :
    _name(name),
    _size_in_bytes(0),
    _output_type(output_type),
    _children() { }

  MemoryTreeNode* addChild(const std::string& name, const size_t size_in_bytes = 0) {
    auto child_iter = _children.find(name);
    if ( child_iter == _children.end() ) {
      MemoryTreeNode* child = new MemoryTreeNode(name, _output_type);
      child->_size_in_bytes = size_in_bytes;
      _children[name] = std::unique_ptr<MemoryTreeNode>(child);
      return child;
    } else {
      return (*child_iter).second.get();
    }
  }

  void updateSize(const size_t delta) {
    _size_in_bytes += delta;
  }

  void finalize() {
    for ( auto& child : _children ) {
      child.second->finalize();
    }

    // Aggregate size of childs
    for ( auto& child : _children ) {
      _size_in_bytes += child.second->_size_in_bytes;
    }
  }

  friend std::ostream & operator<< (std::ostream& str, const MemoryTreeNode& memory_tree_node);

 private:
  std::string _name;
  size_t _size_in_bytes;
  OutputType _output_type;
  map_type _children;
};

static std::string serialize_in_bytes(const size_t size_in_bytes) {
  std::stringstream ss;
  ss << size_in_bytes << " bytes";
  return ss.str();
}

static std::string serialize_in_kilobytes(const size_t size_in_bytes) {
  std::stringstream ss;
  const double size_in_kb = static_cast<double>(size_in_bytes) / 1000.0;
  ss << std::fixed << std::setprecision(3) << size_in_kb << " KB";
  return ss.str();
}

static std::string serialize_in_megabytes(const size_t size_in_bytes) {
  std::stringstream ss;
  const double size_in_mb = static_cast<double>(size_in_bytes) / 1000000.0;
  ss << std::fixed << std::setprecision(3) << size_in_mb << " MB";
  return ss.str();
}

static std::string serialize_in_percentage(const size_t parent_size_in_bytes,
                                           const size_t size_in_bytes) {
  if ( parent_size_in_bytes > 0 ) {
    std::stringstream ss;
    const double percentage = ( static_cast<double>(size_in_bytes) /
        static_cast<double>(parent_size_in_bytes) ) * 100.0;
    ss << std::fixed << std::setprecision(2) << percentage << "%";
    return ss.str();
  } else {
    return serialize_in_megabytes(size_in_bytes);
  }
}

static std::string serialize_metric(const OutputType& type,
                                    const size_t parent_size_in_bytes,
                                    const size_t size_in_bytes) {
  switch(type) {
    case OutputType::BYTES:
      return serialize_in_bytes(size_in_bytes);
    case OutputType::KILOBYTE:
      return serialize_in_kilobytes(size_in_bytes);
    case OutputType::MEGABYTE:
      return serialize_in_megabytes(size_in_bytes);
    case OutputType::PERCENTAGE:
      return serialize_in_percentage(parent_size_in_bytes, size_in_bytes);
  }
  return "";
}

inline char MemoryTreeNode::LINE_PREFIX[] = " + ";

inline std::ostream & operator<< (std::ostream& str, const MemoryTreeNode& memory_tree_node) {

  auto print = [&](std::ostream& str,
                   const size_t parent_size_in_bytes,
                   const MemoryTreeNode& memory_tree_node,
                   int level) {
                 std::string prefix = "";
                 prefix += level == 0 ? std::string(MemoryTreeNode::LINE_PREFIX, MemoryTreeNode::LINE_PREFIX_LENGTH) :
                           std::string(MemoryTreeNode::LINE_PREFIX_LENGTH, ' ');
                 prefix += level > 0 ? std::string(MemoryTreeNode::LINE_PREFIX_LENGTH * (level - 1), ' ') : "";
                 prefix += level > 0 ? std::string(MemoryTreeNode::LINE_PREFIX, MemoryTreeNode::LINE_PREFIX_LENGTH) : "";
                 size_t length = prefix.size() + memory_tree_node._name.size();
                 str << prefix
                     << memory_tree_node._name;
                 if (length < MemoryTreeNode::MAX_LINE_LENGTH) {
                   str << std::string(MemoryTreeNode::MAX_LINE_LENGTH - length, ' ');
                 }
                 str << " = "
                     << serialize_metric(memory_tree_node._output_type,
                          parent_size_in_bytes, memory_tree_node._size_in_bytes) << "\n";
               };


  std::function<void(std::ostream&, const size_t, const MemoryTreeNode&, int)> dfs =
    [&](std::ostream& str,
        const size_t parent_size_in_bytes,
        const MemoryTreeNode& parent,
        int level) {
      if ( parent._size_in_bytes > 0 ) {
        print(str, parent_size_in_bytes, parent, level);
        for (const auto& child : parent._children) {
          dfs(str, parent._size_in_bytes, *child.second.get(), level + 1);
        }
      }
    };

  dfs(str, 0UL, memory_tree_node, 0);

  return str;
}

}  // namespace mt_kahypar
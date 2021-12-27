/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "mt-kahypar/utils/memory_tree.h"



#include <functional>
#include <iomanip>

namespace mt_kahypar::utils {

  MemoryTreeNode::MemoryTreeNode(const std::string& name, const OutputType& output_type) :
          _name(name),
          _size_in_bytes(0),
          _output_type(output_type),
          _children() { }

  MemoryTreeNode* MemoryTreeNode::addChild(const std::string& name, const size_t size_in_bytes) {
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

  void MemoryTreeNode::finalize() {
    for ( auto& child : _children ) {
      child.second->finalize();
    }

    // Aggregate size of childs
    for ( auto& child : _children ) {
      _size_in_bytes += child.second->_size_in_bytes;
    }
  }


  std::string serialize_in_bytes(const size_t size_in_bytes) {
    std::stringstream ss;
    ss << size_in_bytes << " bytes";
    return ss.str();
  }

  std::string serialize_in_kilobytes(const size_t size_in_bytes) {
    std::stringstream ss;
    const double size_in_kb = static_cast<double>(size_in_bytes) / 1000.0;
    ss << std::fixed << std::setprecision(3) << size_in_kb << " KB";
    return ss.str();
  }

  std::string serialize_in_megabytes(const size_t size_in_bytes) {
    std::stringstream ss;
    const double size_in_mb = static_cast<double>(size_in_bytes) / 1000000.0;
    ss << std::fixed << std::setprecision(3) << size_in_mb << " MB";
    return ss.str();
  }

  std::string serialize_in_percentage(const size_t parent_size_in_bytes,
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

  std::string serialize_metric(const OutputType& type,
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


  void MemoryTreeNode::print(std::ostream& str, const size_t parent_size_in_bytes, int level) const {

    constexpr int MAX_LINE_LENGTH = 45;
    constexpr size_t LINE_PREFIX_LENGTH = 3;

    std::string prefix;
    prefix += level == 0 ? " + " :
              std::string(LINE_PREFIX_LENGTH, ' ');
    prefix += level > 0 ? std::string(LINE_PREFIX_LENGTH * (level - 1), ' ') : "";
    prefix += level > 0 ? " + " : "";
    size_t length = prefix.size() + _name.size();
    str << prefix
        << _name;
    if (length < MAX_LINE_LENGTH) {
      str << std::string(MAX_LINE_LENGTH - length, ' ');
    }
    str << " = "
        << serialize_metric(_output_type,
                            parent_size_in_bytes, _size_in_bytes) << "\n";

  }

  void MemoryTreeNode::dfs(std::ostream& str, const size_t parent_size_in_bytes, int level) const {
    if ( _size_in_bytes > 0 ) {
      print(str, parent_size_in_bytes, level);
      for (const auto& child : _children) {
        child.second->dfs(str, parent_size_in_bytes + _size_in_bytes, level + 1);
      }
    }
  }

  std::ostream & operator<< (std::ostream& str, const MemoryTreeNode& root) {
    root.dfs(str, 0UL, 0);
    return str;
  }


}
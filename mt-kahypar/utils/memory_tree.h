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

#pragma once

#include <map>
#include <string>
#include <memory>

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
  MemoryTreeNode(const std::string& name, const OutputType& output_type = OutputType::MEGABYTE);

  MemoryTreeNode* addChild(const std::string& name, const size_t size_in_bytes = 0);

  void updateSize(const size_t delta) {
    _size_in_bytes += delta;
  }

  void finalize();

 private:

  void dfs(std::ostream& str, const size_t parent_size_in_bytes, int level) const ;
  void print(std::ostream& str, const size_t parent_size_in_bytes, int level) const ;

  friend std::ostream& operator<<(std::ostream& str, const MemoryTreeNode& root);

  std::string _name;
  size_t _size_in_bytes;
  OutputType _output_type;
  map_type _children;
};


std::ostream & operator<< (std::ostream& str, const MemoryTreeNode& root);

}  // namespace mt_kahypar
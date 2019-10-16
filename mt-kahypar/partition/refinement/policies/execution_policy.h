/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2018 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2018 Tobias Heuer <tobias.heuer@live.com>
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

#include <vector>

#include "kahypar/meta/policy_registry.h"
#include "kahypar/meta/typelist.h"

#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {

template <class Derived = Mandatory>
class ExecutionPolicy : public kahypar::meta::PolicyBase {
 public:
  ExecutionPolicy() :
    _execution_levels(),
    _alpha(2.0) { }

  ExecutionPolicy(const double alpha) :
    _execution_levels(),
    _alpha(alpha) { }

  template < typename HyperGraph >
  void initialize(const HyperGraph& hg, const HypernodeID current_num_nodes) {
    static_cast<Derived*>(this)->initializeImpl(hg, current_num_nodes);
  }

  bool execute(const size_t level) {
    if (_execution_levels.size() == 0) {
      return false;
    }

    if (level >= _execution_levels.back()) {
      _execution_levels.pop_back();
      return true;
    } else {
      return false;
    }
  }

 protected:
  std::vector<size_t> _execution_levels;
  double _alpha;
};


class MultilevelExecutionPolicy : public ExecutionPolicy<MultilevelExecutionPolicy>{
 public:
  MultilevelExecutionPolicy() :
    ExecutionPolicy() { }

  MultilevelExecutionPolicy(const double alpha) :
    ExecutionPolicy(alpha) { }

  template < typename HyperGraph >
  void initialize(const HyperGraph& hg, const HypernodeID current_num_nodes) {
    ASSERT(_alpha >= 1.0);
    size_t last_level = 0;
    for (size_t i = 0; hg.initialNumNodes() / std::pow(_alpha, i) >= current_num_nodes; ++i) {
      size_t next_level = hg.initialNumNodes() / std::pow(_alpha, i) - current_num_nodes;
      if ( last_level != next_level ) {
        _execution_levels.push_back(next_level);
        last_level = next_level;
      }
    }
  }

 private:
  using ExecutionPolicy::_execution_levels;
  using ExecutionPolicy::_alpha;
};

class ExponentialExecutionPolicy : public ExecutionPolicy<ExponentialExecutionPolicy>{
 public:
  ExponentialExecutionPolicy() :
    ExecutionPolicy() { }

  ExponentialExecutionPolicy(const double alpha) :
    ExecutionPolicy(alpha) { }

  template < typename HyperGraph >
  void initialize(const HyperGraph& hg, const HypernodeID current_num_nodes) {
    ASSERT(_alpha >= 1.0);
    size_t last_level = 0;
    for (size_t i = 0; current_num_nodes + std::pow(_alpha, i) < hg.initialNumNodes(); ++i) {
      size_t next_level = std::pow(_alpha, i);
      if ( last_level != next_level ) {
        _execution_levels.push_back(next_level);
        last_level = next_level;
      }
    }
    _execution_levels.push_back(hg.initialNumNodes() - current_num_nodes);
    std::reverse(_execution_levels.begin(), _execution_levels.end());
  }

 private:
  using ExecutionPolicy::_execution_levels;
  using ExecutionPolicy::_alpha;
};

class ConstantExecutionPolicy : public ExecutionPolicy<ConstantExecutionPolicy>{
 public:
  ConstantExecutionPolicy() :
    ExecutionPolicy() { }

  ConstantExecutionPolicy(const double alpha) :
    ExecutionPolicy(alpha) { }

  template < typename HyperGraph >
  void initialize(const HyperGraph& hg, const HypernodeID current_num_nodes) {
    ASSERT(_alpha >= 1.0);
    size_t constant = _alpha;
    size_t last_level = 0;
    for (size_t i = 0; current_num_nodes + i * constant < hg.initialNumNodes(); ++i) {
      size_t next_level = std::max( i * constant, 1UL );
      if ( last_level != next_level ) {
        _execution_levels.push_back(next_level);
        last_level = next_level;
      }
    }
    _execution_levels.push_back(hg.initialNumNodes() - current_num_nodes);
    std::reverse(_execution_levels.begin(), _execution_levels.end());
  }

 private:
  using ExecutionPolicy::_execution_levels;
  using ExecutionPolicy::_alpha;
};

using ExecutionPolicyClasses = meta::Typelist<ExponentialExecutionPolicy, MultilevelExecutionPolicy, ConstantExecutionPolicy>;

}  // namespace mt_kahypar

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

#include <mutex>
#include <string>
#include <unordered_map>
#include <functional>

#include "kahypar/macros.h"

namespace mt_kahypar {
namespace utils {

class Timer {

 static constexpr bool debug = false;

 public:
  static constexpr int MAX_LINE_LENGTH = 40;
  static char TOP_LEVEL_PREFIX[];
  static constexpr size_t TOP_LEVEL_PREFIX_LENGTH = 3;
  static char SUB_LEVEL_PREFIX[];
  static constexpr size_t SUB_LEVEL_PREFIX_LENGTH = 3;

  enum class Type : uint8_t {
    IMPORT = 0,
    PREPROCESSING = 1,
    COARSENING = 2,
    INITIAL_PARTITIONING = 3,
    REFINEMENT = 4
  };

 private:

  using Key = std::pair<std::string, std::string>;

  struct PairHasher {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const
    {
      return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
  };

  class Timing {

    public:
      Timing(const std::string& name,
             const std::string& description,
             const std::string& parent,
             const Type& type,
             const int order) :
        _name(name),
        _description(description),
        _parent(parent),
        _type(type),
        _order(order),
        _timing(0.0) { }

    std::string name() const {
      return _name;
    }

    std::string description() const {
      return _description;
    }

    std::string parent() const {
      return _parent;
    }

    Type type() const {
      return _type;
    }

    int order() const {
      return _order;
    }

    bool is_root() const {
      return (_parent == "");
    }

    double timing() const {
      return _timing;
    }

    void add_timing(const double timing) {
      _timing += timing;
    }

    private:
      std::string _name;
      std::string _description;
      std::string _parent;
      Type _type;
      int _order;
      double _timing;
  };

 public:

  static Timer& instance(bool show_detailed_timings = false) {
    if ( _instance == nullptr ) {
      std::lock_guard<std::mutex> _lock(_mutex);
      if ( _instance == nullptr ) {
        _instance = new Timer(show_detailed_timings);
      }
    }
    _instance->_show_detailed_timings = show_detailed_timings;
    return *_instance;
  }

  void add_timing( const std::string& name, const std::string& description,
                   const std::string& parent, const Type& type, const int order,
                   const double timing ) {
    std::lock_guard<std::mutex> lock(_timing_mutex);
    Key key = std::make_pair(parent, name);
    if ( _timings.find(key) == _timings.end() ) {
      _timings.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(key),
        std::forward_as_tuple(name, description, parent, type, order));
        _timings.at(key).add_timing(timing);
    }
  }

  void update_timing( const std::string& name, const std::string& description,
                      const std::string& parent, const Type& type, const int order,
                      const double timing ) {
    std::lock_guard<std::mutex> lock(_timing_mutex);
    Key key = std::make_pair(parent, name);
    if ( _timings.find(key) == _timings.end() ) {
      _timings.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(key),
        std::forward_as_tuple(name, description, parent, type, order));
    }
    _timings.at(key).add_timing(timing);
  }

  void serialize(std::ostream& str) {
    for ( const auto& timing : _timings ) {
      const Timing& time_point = timing.second;
      if ( _show_detailed_timings || time_point.parent() == "" ) {
        str << " " << time_point.name() << "=" << time_point.timing();
      }
    }
  }

  friend std::ostream& operator<<(std::ostream& str, const Timer& timer);

 private:
  explicit Timer(const bool show_detailed_timings) :
    _timing_mutex(),
    _timings(),
    _show_detailed_timings(show_detailed_timings) { }

  static std::mutex _mutex;
  static Timer* _instance;

  std::mutex _timing_mutex;
  std::unordered_map<Key, Timing, PairHasher> _timings;
  bool _show_detailed_timings;
};

Timer* Timer::_instance { nullptr };
std::mutex Timer::_mutex;

char Timer::TOP_LEVEL_PREFIX[] = " + ";
char Timer::SUB_LEVEL_PREFIX[] = " + ";


std::ostream& operator<<(std::ostream& str, const Timer& timer) {
  std::vector<Timer::Timing> timings;
  for ( const auto& timing : timer._timings ) {
    timings.emplace_back(timing.second);
  }
  std::sort(timings.begin(), timings.end(),
    [&](const Timer::Timing& lhs, const Timer::Timing& rhs) {
      return lhs.type() < rhs.type() || (lhs.type() == rhs.type() && lhs.order() < rhs.order());
    });

  auto print = [&](std::ostream& str, const Timer::Timing& timing, int level) {
    std::string prefix = "";
    prefix += level == 0 ? std::string(Timer::TOP_LEVEL_PREFIX, Timer::TOP_LEVEL_PREFIX_LENGTH) :
                           std::string(Timer::TOP_LEVEL_PREFIX_LENGTH, ' ');
    prefix += level > 0 ? std::string(Timer::SUB_LEVEL_PREFIX_LENGTH * (level - 1), ' ' ) : "";
    prefix += level > 0 ? std::string(Timer::SUB_LEVEL_PREFIX, Timer::SUB_LEVEL_PREFIX_LENGTH) : "";
    size_t length = prefix.size() + timing.description().size();
    str << prefix
        << timing.description();
    if ( length < Timer::MAX_LINE_LENGTH ) {
      str << std::string(Timer::MAX_LINE_LENGTH - length, ' ');
    }
    str << " = " << timing.timing() << " s\n";
  };

  std::function<void(std::ostream&,const Timer::Timing&,int)> dfs =
    [&](std::ostream& str, const Timer::Timing& parent, int level) {
    for ( const Timer::Timing& timing : timings ) {
      if ( timing.parent() == parent.name() ) {
        print(str, timing, level);
        dfs(str, timing, level + 1);
      }
    }
  };

  for ( const Timer::Timing& timing : timings ) {
    if ( timing.is_root() ) {
      print(str, timing, 0);
      if ( timer._show_detailed_timings ) {
        dfs(str, timing, 1);
      }
    }
  }

  return str;
}

} // namespace utils
} // namespace mt_kahypar
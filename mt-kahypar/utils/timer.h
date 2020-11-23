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

#include <atomic>
#include <chrono>
#include <functional>
#include <mutex>
#include <string>
#include <unordered_map>

#include <tbb/enumerable_thread_specific.h>

#include "mt-kahypar/macros.h"

namespace mt_kahypar {
namespace utils {
class Timer {
  static constexpr bool debug = false;

  using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

 public:
  static constexpr int MAX_LINE_LENGTH = 45;
  static char TOP_LEVEL_PREFIX[];
  static constexpr size_t TOP_LEVEL_PREFIX_LENGTH = 3;
  static char SUB_LEVEL_PREFIX[];
  static constexpr size_t SUB_LEVEL_PREFIX_LENGTH = 3;

 private:
  struct Key {
    std::string parent;
    std::string key;
  };

  struct KeyHasher {
    std::size_t operator() (const Key& key) const {
      return std::hash<std::string>()(key.parent) ^ std::hash<std::string>()(key.key);
    }
  };

  struct KeyEqual {
    bool operator() (const Key& lhs, const Key& rhs) const {
      return lhs.parent == rhs.parent && lhs.key == rhs.key;
    }
  };

  class ActiveTiming {
   public:
    ActiveTiming() :
      _key(""),
      _description(""),
      _start() { }

    ActiveTiming(const std::string& key,
                 const std::string& description,
                 const HighResClockTimepoint& start) :
      _key(key),
      _description(description),
      _start(start) { }

    std::string key() const {
      return _key;
    }

    std::string description() const {
      return _description;
    }

    HighResClockTimepoint start() const {
      return _start;
    }

   private:
    std::string _key;
    std::string _description;
    HighResClockTimepoint _start;
  };

  class Timing {
   public:
    Timing(const std::string& key,
           const std::string& description,
           const std::string& parent,
           const int order) :
      _key(key),
      _description(description),
      _parent(parent),
      _order(order),
      _timing(0.0) { }

    std::string key() const {
      return _key;
    }

    std::string description() const {
      return _description;
    }

    std::string parent() const {
      return _parent;
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
    std::string _key;
    std::string _description;
    std::string _parent;
    int _order;
    double _timing;
  };

  using ActiveTimingStack = std::vector<ActiveTiming>;
  using LocalActiveTimingStack = tbb::enumerable_thread_specific<ActiveTimingStack>;

 public:
  Timer(const Timer&) = delete;
  Timer & operator= (const Timer &) = delete;

  Timer(Timer&&) = delete;
  Timer & operator= (Timer &&) = delete;

  static Timer& instance(bool show_detailed_timings = false) {
    static Timer instance;
    instance._show_detailed_timings = show_detailed_timings;
    return instance;
  }

  void setMaximumOutputDepth(const size_t max_output_depth) {
    _max_output_depth = max_output_depth;
  }

  bool isEnabled() const {
    return _is_enabled;
  }

  void enable() {
    std::lock_guard<std::mutex> lock(_timing_mutex);
    _is_enabled = true;
  }

  void disable() {
    std::lock_guard<std::mutex> lock(_timing_mutex);
    _is_enabled = false;
  }

  void clear() {
    std::lock_guard<std::mutex> lock(_timing_mutex);
    _timings.clear();
    _active_timings.clear();
    _index = 0;
  }

  void start_timer(const std::string& key,
                   const std::string& description,
                   bool is_parallel_context = false,
                   bool force = false) {
    if (_is_enabled || force) {
      std::lock_guard<std::mutex> lock(_timing_mutex);
      if (force || is_parallel_context) {
        _local_active_timings.local().emplace_back(key, description, std::chrono::high_resolution_clock::now());
      } else {
        _active_timings.emplace_back(key, description, std::chrono::high_resolution_clock::now());
      }
    }
  }

  void stop_timer(const std::string& key, bool force = false) {
    unused(key);
    if (_is_enabled || force) {
      std::lock_guard<std::mutex> lock(_timing_mutex);
      HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
      ASSERT(!force || !_local_active_timings.local().empty());
      ActiveTiming current_timing;
      // First check if there are some active timings on the local stack
      // (in that case we are in a parallel context) and if there are
      // no active timings we pop from global stack
      if (!_local_active_timings.local().empty()) {
        ASSERT(_local_active_timings.local().back().key() == key);
        current_timing = _local_active_timings.local().back();
        _local_active_timings.local().pop_back();
      } else {
        ASSERT(!_active_timings.empty());
        ASSERT(_active_timings.back().key() == key, V(_active_timings.back().key()) << V(key));
        current_timing = _active_timings.back();
        _active_timings.pop_back();
      }

      // Parent is either the last element on the local stack and
      // if there are no timings on the local stack the parent is
      // on the global stack. If there are no elements on the global
      // stack the timing represents a root.
      std::string parent = "";
      if (!_local_active_timings.local().empty()) {
        parent = _local_active_timings.local().back().key();
      } else if (!_active_timings.empty()) {
        parent = _active_timings.back().key();
      }

      Key timing_key { parent, current_timing.key() };
      if (_timings.find(timing_key) == _timings.end()) {
        _timings.emplace(
          std::piecewise_construct,
          std::forward_as_tuple(timing_key),
          std::forward_as_tuple(current_timing.key(),
                                current_timing.description(), parent, _index++));
      }
      double time = std::chrono::duration<double>(end - current_timing.start()).count();
      _timings.at(timing_key).add_timing(time);
    }
  }

  void serialize(std::ostream& str) {
    std::vector<Timing> timings;
    for (const auto& timing : _timings) {
      timings.emplace_back(timing.second);
    }
    std::sort(timings.begin(), timings.end(),
              [&](const Timing& lhs, const Timing& rhs) {
          return lhs.key() < rhs.key();
        });

    if (!timings.empty()) {
      auto print = [&](const std::string& key, const double time, const bool is_root) {
                     if (_show_detailed_timings || is_root) {
                       str << " " << key << "=" << time;
                     }
                   };

      // We aggregate timings in case there multiple timings with same key
      // but different parent
      std::string last_key = timings[0].key();
      double time = timings[0].timing();
      bool is_root = timings[0].is_root();
      for (size_t i = 1; i < timings.size(); ++i) {
        if (last_key == timings[i].key()) {
          time += timings[i].timing();
          is_root |= timings[i].is_root();
        } else {
          print(last_key, time, is_root);
          last_key = timings[i].key();
          time = timings[i].timing();
          is_root = timings[i].is_root();
        }
      }
      print(last_key, time, is_root);
    }
  }

  friend std::ostream & operator<< (std::ostream& str, const Timer& timer);

  double get(std::string key) const {
    for (const auto& x : _timings) {
      // unfortunately it has to be linear search because the parent (which we can't lookup at this stage) is part of the map key
      if (x.first.key == key) {
        return x.second.timing();
      }
    }
    return 0.0;
  }

 private:
  explicit Timer() :
    _timing_mutex(),
    _timings(),
    _active_timings(),
    _local_active_timings(),
    _index(0),
    _is_enabled(true),
    _show_detailed_timings(false),
    _max_output_depth(std::numeric_limits<size_t>::max()) { }

  std::mutex _timing_mutex;
  std::unordered_map<Key, Timing, KeyHasher, KeyEqual> _timings;
  // Global Active Timing Stack
  // Timings are pushed to the global stack
  // if we are in a sequential context
  ActiveTimingStack _active_timings;
  // Thread Local Timing Stack
  // Timings are pushed to the thread local stack
  // if we are in a parallel context
  LocalActiveTimingStack _local_active_timings;
  std::atomic<int> _index;
  bool _is_enabled;
  bool _show_detailed_timings;
  size_t _max_output_depth;
};

inline char Timer::TOP_LEVEL_PREFIX[] = " + ";
inline char Timer::SUB_LEVEL_PREFIX[] = " + ";

inline std::ostream & operator<< (std::ostream& str, const Timer& timer) {
  std::vector<Timer::Timing> timings;
  for (const auto& timing : timer._timings) {
    timings.emplace_back(timing.second);
  }
  std::sort(timings.begin(), timings.end(),
            [&](const Timer::Timing& lhs, const Timer::Timing& rhs) {
        return lhs.order() < rhs.order();
      });

  auto print = [&](std::ostream& str, const Timer::Timing& timing, int level) {
                 std::string prefix = "";
                 prefix += level == 0 ? std::string(Timer::TOP_LEVEL_PREFIX, Timer::TOP_LEVEL_PREFIX_LENGTH) :
                           std::string(Timer::TOP_LEVEL_PREFIX_LENGTH, ' ');
                 prefix += level > 0 ? std::string(Timer::SUB_LEVEL_PREFIX_LENGTH * (level - 1), ' ') : "";
                 prefix += level > 0 ? std::string(Timer::SUB_LEVEL_PREFIX, Timer::SUB_LEVEL_PREFIX_LENGTH) : "";
                 size_t length = prefix.size() + timing.description().size();
                 str << prefix
                     << timing.description();
                 if (length < Timer::MAX_LINE_LENGTH) {
                   str << std::string(Timer::MAX_LINE_LENGTH - length, ' ');
                 }
                 str << " = " << timing.timing() << " s\n";
               };

  std::function<void(std::ostream&, const Timer::Timing&, int)> dfs =
    [&](std::ostream& str, const Timer::Timing& parent, int level) {
      if ( level <= timer._max_output_depth ) {
        for (const Timer::Timing& timing : timings) {
          if (timing.parent() == parent.key()) {
            print(str, timing, level);
            dfs(str, timing, level + 1);
          }
        }
      }
    };

  for (const Timer::Timing& timing : timings) {
    if (timing.is_root()) {
      print(str, timing, 0);
      if (timer._show_detailed_timings) {
        dfs(str, timing, 1);
      }
    }
  }

  return str;
}
}  // namespace utils
}  // namespace mt_kahypar

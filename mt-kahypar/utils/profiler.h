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

#include <thread>
#include <chrono>
#include <vector>
#include <fstream>

#ifdef __linux__
#include <cstring>
#include <cstdint>
#include <sys/ioctl.h>
#include <asm/unistd.h>
#include <linux/perf_event.h>
#else
#error "profiler.h can only be used on linux systems"
#endif

#include "mt-kahypar/macros.h"

namespace mt_kahypar {
namespace utils {

class Profiler {

  using sys_time = unsigned long long;

  struct CPUStat {
    int cpu_id;
    sys_time user;
    sys_time user_low;
    sys_time system;
    sys_time idle;

    double computeUtilization(const CPUStat& last) const {
      if ( user < last.user || user_low < last.user_low ||
           system < last.system || idle < last.idle ) {
        // Overflow detected => Invalid Value
        return -1.0;
      } else {
        double busy = (user - last.user) + (user_low - last.user_low) + (system - last.system);
        double total = busy + (idle - last.idle);
        return busy / total;
      }
    }
  };

  struct CPUSnapshot {
    CPUSnapshot() :
      _stats(std::thread::hardware_concurrency()) { }

    std::vector<CPUStat> _stats;
  };

  class PerfHardwareEvent {

   public:
    explicit PerfHardwareEvent(const perf_hw_id config) :
      _hardware_event(),
      _event_fd(-1) {
      _hardware_event.type = PERF_TYPE_HARDWARE;
      _hardware_event.config = config;
      _hardware_event.disabled = 1;
      _hardware_event.exclude_kernel = 1;
      _hardware_event.exclude_hv = 1;
    }

    ~PerfHardwareEvent() {
      stop();
    }

    void start() {
      _event_fd = syscall(__NR_perf_event_open, &_hardware_event, 0, -1, -1, 0);
      if ( _event_fd != -1 ) {
        ioctl(_event_fd, PERF_EVENT_IOC_RESET, 0);
        ioctl(_event_fd, PERF_EVENT_IOC_ENABLE, 0);
      }
    }

    void stop() {
      if ( _event_fd != -1 ) {
        ioctl(_event_fd, PERF_EVENT_IOC_DISABLE, 0);
        close(_event_fd);
      }
    }

    long long readValue() {
      long long value = 0LL;
      if ( _event_fd != -1 ) {
        size_t ret = read(_event_fd, &value, sizeof(long long));
        unused(ret);
      }
      return value;
    }

   private:
    perf_event_attr _hardware_event;
    int _event_fd;
  };

  struct PerfStats {
    PerfStats() :
      cache_misses(0),
      cache_references(0) { }

    long long cache_misses;
    long long cache_references;
  };

  class PerfEventManager {

   public:
    PerfEventManager() :
      _cache_miss_event(PERF_COUNT_HW_CACHE_MISSES),
      _cache_reference_event(PERF_COUNT_HW_CACHE_REFERENCES),
      _last_stats(),
      _total_stats() {
      _cache_miss_event.start();
      _cache_reference_event.start();
    }

    void start() {
      _last_stats.cache_misses = _cache_miss_event.readValue();
      _last_stats.cache_references = _cache_reference_event.readValue();
    }

    void stop() {
      long long cache_misses = _cache_miss_event.readValue();
      long long cache_references = _cache_reference_event.readValue();
      _total_stats.cache_misses += (cache_misses - _last_stats.cache_misses);
      _total_stats.cache_references += (cache_references - _last_stats.cache_references);
    }

    const PerfStats& perf_stats() const {
      return _total_stats;
    }

   private:
    PerfHardwareEvent _cache_miss_event;
    PerfHardwareEvent _cache_reference_event;
    PerfStats _last_stats;
    PerfStats _total_stats;
  };

 public:
  class StatContainer {

   public:
    StatContainer(const std::string& description) :
      _is_active(true),
      _description(description),
      _cpu_snapshots(),
      _memory(),
      _perf_manager() {
      _perf_manager.start();
    }

    bool is_active() const {
      return _is_active;
    }

    void activate() {
      ASSERT(!_is_active);
      _is_active = true;
      _perf_manager.start();
    }

    void deactivate() {
      ASSERT(_is_active);
      _is_active = false;
      _perf_manager.stop();
    }

    const std::string& description() const {
      return _description;
    }

    void add_snapshot(const CPUSnapshot& snapshot, const int memory_in_kb) {
      _cpu_snapshots.push_back(snapshot);
      _memory.push_back(memory_in_kb);
    }

    double cacheMissRatio() const {
      const PerfStats& perf_stats = _perf_manager.perf_stats();
      return static_cast<double>(perf_stats.cache_misses) /
        static_cast<double>(perf_stats.cache_references);
    }

  friend std::ostream & operator<< (std::ostream& str, const StatContainer& stat);

   private:
    bool _is_active;
    std::string _description;
    std::vector<CPUSnapshot> _cpu_snapshots;
    std::vector<int> _memory;
    PerfEventManager _perf_manager;
  };

  Profiler(const Profiler&) = delete;
  Profiler & operator= (const Profiler &) = delete;

  Profiler(Profiler&&) = delete;
  Profiler & operator= (Profiler &&) = delete;

  static Profiler& instance(const int snapshot_interval = std::numeric_limits<int>::max()) {
    static Profiler instance(snapshot_interval);
    return instance;
  }

  void start() {
    ASSERT(!_profiler && !_terminate);
    _profiler = std::make_unique<std::thread>([&] {
      while ( !_terminate ) {
        {
          std::lock_guard<std::mutex> lock(_profiler_mutex);
          const CPUSnapshot snapshot = readCPUSnapshot();
          const int memory_in_kb = readPhysicalMemoryUsage();
          for ( StatContainer& container : _stats ) {
            if ( container.is_active() ) {
              container.add_snapshot(snapshot, memory_in_kb);
            }
          }
        }
        std::this_thread::sleep_for(_snapshot_interval);
      }
      _stats[0].deactivate();
    });
  }

  void stop() {
    _terminate = true;
    if ( _profiler ) {
      _profiler->join();
    }
  }

  void activate(const std::string& description) {
    if ( _profiler ) {
      std::lock_guard<std::mutex> lock(_profiler_mutex);
      bool found = false;
      for ( StatContainer& container : _stats ) {
        if ( description == container.description() ) {
          container.activate();
          found = true;
          break;
        }
      }

      if ( !found ) {
        _stats.emplace_back(description);
      }
    }
  }

  void deactivate(const std::string& description) {
    if ( _profiler ) {
      std::lock_guard<std::mutex> lock(_profiler_mutex);
      for ( StatContainer& container : _stats ) {
        if ( description == container.description() ) {
          container.deactivate();
          break;
        }
      }
    }
  }

  friend std::ostream & operator<< (std::ostream& str, const Profiler& profiler);

 private:
  explicit Profiler(const int snapshot_interval) :
    _terminate(false),
    _profiler_mutex(),
    _profiler(nullptr),
    _snapshot_interval(snapshot_interval),
    _stats() {
    _stats.emplace_back("Overall");
  }

  CPUSnapshot readCPUSnapshot() {
    CPUSnapshot snapshot;

    std::string stat_line;
    std::string cpu_attr;
    std::ifstream stat_stream("/proc/stat");
    std::getline(stat_stream, stat_line); // Skip overall cpu summary
    for ( uint32_t cpu_id = 0; cpu_id < std::thread::hardware_concurrency(); ++cpu_id ) {
      std::getline(stat_stream, stat_line);
      std::istringstream ss(stat_line);
      snapshot._stats[cpu_id].cpu_id = cpu_id;
      ss >> cpu_attr >> snapshot._stats[cpu_id].user >> snapshot._stats[cpu_id].user_low
         >> snapshot._stats[cpu_id].system >> snapshot._stats[cpu_id].idle;
    }
    stat_stream.close();
    return snapshot;
  }

  int readPhysicalMemoryUsage() {
    std::string stat_line;
    std::ifstream stat_stream("/proc/self/status");
    while ( getline(stat_stream, stat_line) ) {
      if ( strncmp( stat_line.c_str(), "VmRSS", 5 ) == 0 ) {
        return parsePhysicalMemoryUsage(stat_line);
      }
    }
    return 0;
  }

  int parsePhysicalMemoryUsage(const std::string& stat_line) {
    // This assumes that a digit will be found and the line ends in " Kb".
    const int n = stat_line.length();
    char* line = new char[n + 1];
    strcpy(line, stat_line.c_str());
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    free(line);
    return i;
  }

  bool _terminate;
  std::mutex _profiler_mutex;
  std::unique_ptr<std::thread> _profiler;
  std::chrono::milliseconds _snapshot_interval;

  std::vector<StatContainer> _stats;
};

std::ostream & operator<< (std::ostream& str, const Profiler::StatContainer& stats) {
  if ( stats._cpu_snapshots.size() > 1 ) {
    str << "\n" << BOLD << stats._description << ":" << END << std::endl;
    // Memory Consumption
    size_t min_mem = std::numeric_limits<size_t>::max();
    size_t avg_mem = 0;
    size_t max_mem = std::numeric_limits<size_t>::min();
    for ( const size_t mem : stats._memory ) {
      min_mem = std::min(mem, min_mem);
      avg_mem += mem;
      max_mem = std::max(mem, max_mem);
    }
    avg_mem /= static_cast<int>(stats._memory.size());
    str << "  Physical Memory Usage: [ "
        << "Min = " << ( min_mem / 1000 ) << " MB, "
        << "Avg = " << ( avg_mem / 1000 ) << " MB, "
        << "Max = " << ( max_mem / 1000 ) << " MB ]" << std::endl;

    // Cache Miss Ratio
    str << "  L3 Cache Miss Ratio: ";
    const size_t cache_miss_ratio = ( stats.cacheMissRatio() * 100.0 );
    if ( cache_miss_ratio < 25 ) {
      str << GREEN;
    } else if ( cache_miss_ratio < 50 ) {
      str << YELLOW;
    } else {
      str << RED;
    }
    str << cache_miss_ratio << "%" << END << std::endl;

    // CPU Utilization for each processor
    for ( uint32_t cpu_id = 0; cpu_id < std::thread::hardware_concurrency(); ++cpu_id ) {
      double avg_utilization = 0.0;
      size_t num_snapshots = 0;
      for ( size_t i = 1; i < stats._cpu_snapshots.size(); ++i ) {
        double snapshot_utilization =
          stats._cpu_snapshots[i]._stats[cpu_id].computeUtilization(
            stats._cpu_snapshots[i - 1]._stats[cpu_id]);
        if ( snapshot_utilization >= 0 ) {
          avg_utilization += snapshot_utilization;
          ++num_snapshots;
        }
      }
      avg_utilization /= static_cast<double>(num_snapshots);
      const size_t utilization_percent = avg_utilization * 100.0;
      str << "  CPU " << cpu_id << " [ ";
      if ( utilization_percent > 90 ) {
        std::cout << GREEN;
      } else if ( utilization_percent > 80 ) {
        std::cout << YELLOW;
      } else {
        std::cout << RED;
      }
      for ( size_t j = 0; j < utilization_percent / 2UL; ++j ) {
        std::cout << "#";
      }
      for ( size_t j = 0; j < 50UL - ( utilization_percent / 2UL ); ++j ) {
        std::cout << " ";
      }
      std::cout << END << " ] " << utilization_percent << "%" << std::endl;
    }
  }
  return str;
}

std::ostream & operator<< (std::ostream& str, const Profiler& profiler) {
  str << "\n********************************************************************************\n";
  str << "*                              Profiler Results                                *\n";
  str << "********************************************************************************\n";
  for ( const auto& stat : profiler._stats ) {
    str << stat << std::endl;
  }
  return str;
}
}  // namespace utils
}  // namespace mt_kahypar

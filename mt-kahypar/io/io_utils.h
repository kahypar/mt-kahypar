/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#pragma once

#include <fcntl.h>
#include <string>
#include <thread>
#include <sys/stat.h>

#if _WIN32
#include <windows.h>
#include <process.h>
#include <memoryapi.h>
#else
#include <sys/mman.h>
#include <unistd.h>
#endif

#include <tbb/parallel_for.h>

#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/exception.h"

namespace mt_kahypar {
namespace io {

#if _WIN32
struct FileHandle {
  HANDLE hFile;
  HANDLE hMem;
  char* mapped_file;
  size_t length;

  void closeHandle() {
    CloseHandle(hFile);
    CloseHandle(hMem);
  }
};
#else
struct FileHandle {
  int fd;
  char* mapped_file;
  size_t length;

  void closeHandle() {
    close(fd);
  }
};
#endif

size_t file_size(const std::string& filename) {
  struct stat stat_buf;
  const int res = stat( filename.c_str(), &stat_buf);
  if (res < 0) {
    throw InvalidInputException("Could not open: " + filename);
  }
  return static_cast<size_t>(stat_buf.st_size);
}

FileHandle mmap_file(const std::string& filename) {
  FileHandle handle;
  handle.length = file_size(filename);

  #ifdef _WIN32
    PSECURITY_DESCRIPTOR pSD;
    SECURITY_ATTRIBUTES  sa;

    /* create security descriptor (needed for Windows NT) */
    pSD = (PSECURITY_DESCRIPTOR) malloc( SECURITY_DESCRIPTOR_MIN_LENGTH );
    if( pSD == NULL ) {
      throw SystemException("Error while creating security descriptor!");
    }

    InitializeSecurityDescriptor(pSD, SECURITY_DESCRIPTOR_REVISION);
    SetSecurityDescriptorDacl(pSD, TRUE, (PACL) NULL, FALSE);

    sa.nLength = sizeof(sa);
    sa.lpSecurityDescriptor = pSD;
    sa.bInheritHandle = TRUE;

    // open file
    handle.hFile = CreateFile ( filename.c_str(), GENERIC_READ, FILE_SHARE_READ,
      &sa, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);

    if (handle.hFile == INVALID_HANDLE_VALUE) {
      free( pSD);
      throw InvalidInputException("Invalid file handle when opening: " + filename);
    }

    // Create file mapping
    handle.hMem = CreateFileMapping( handle.hFile, &sa, PAGE_READONLY, 0, handle.length, NULL);
    free(pSD);
    if (handle.hMem == NULL) {
      throw InvalidInputException("Invalid file mapping when opening: " + filename);
    }

    // map file to memory
    handle.mapped_file = (char*) MapViewOfFile(handle.hMem, FILE_MAP_READ, 0, 0, 0);
    if ( handle.mapped_file == NULL ) {
      throw SystemException("Failed to map file to main memory:" + filename);
    }
  #else
    handle.fd = open(filename.c_str(), O_RDONLY);
    if ( handle.fd < -1 ) {
      throw InvalidInputException("Could not open: " + filename);
    }
    handle.mapped_file = (char*) mmap(0, handle.length, PROT_READ, MAP_PRIVATE, handle.fd, 0);
    if ( handle.mapped_file == MAP_FAILED ) {
      close(handle.fd);
      throw SystemException("Error while mapping file to memory");
    }
  #endif

  return handle;
}

void munmap_file(FileHandle& handle) {
  #ifdef _WIN32
  UnmapViewOfFile(handle.mapped_file);
  #else
  munmap(handle.mapped_file, handle.length);
  #endif
  handle.closeHandle();
}


MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
bool is_line_ending(char* mapped_file, size_t pos) {
  return mapped_file[pos] == '\r' || mapped_file[pos] == '\n' || mapped_file[pos] == '\0';
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
void do_line_ending(char* mapped_file, size_t& pos) {
  ASSERT(is_line_ending(mapped_file, pos));
  if (mapped_file[pos] != '\0') {
    if (mapped_file[pos] == '\r') {     // windows line ending
      ++pos;
      ASSERT(mapped_file[pos] == '\n');
    }
    ++pos;
  }
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
void goto_next_line(char* mapped_file, size_t& pos, const size_t /*length*/) {
  for ( ; ; ++pos ) {
    if ( is_line_ending(mapped_file, pos) ) {
      do_line_ending(mapped_file, pos);
      break;
    }
  }
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
int64_t read_number(char* mapped_file, size_t& pos, const size_t length) {
  int64_t number = 0;
  while ( mapped_file[pos] == ' ' ) {
    ++pos;
  }
  for ( ; pos < length; ++pos ) {
    if ( mapped_file[pos] == ' ' || is_line_ending(mapped_file, pos) ) {
      while ( mapped_file[pos] == ' ' ) {
        ++pos;
      }
      break;
    }
    ASSERT(mapped_file[pos] >= '0' && mapped_file[pos] <= '9');
    number = number * 10 + (mapped_file[pos] - '0');
  }
  return number;
}

struct LineRange {
  const size_t start;
  const size_t end;
  const size_t line_start_index;
  const size_t num_lines;
};

size_t num_parallel_ranges() {
  return 2 * std::thread::hardware_concurrency();
}

vec<LineRange> split_lines(char* mapped_file,
                           size_t& pos,
                           const size_t length,
                           const size_t expected_num_lines) {
  const size_t initial_pos = pos;
  vec<size_t> line_positions;
  line_positions.reserve(expected_num_lines + 1);

  while ( line_positions.size() < expected_num_lines ) {
    // skip comments
    while ( mapped_file[pos] == '%' ) {
      goto_next_line(mapped_file, pos, length);
      ASSERT(pos < length); // TODO: pos >= length
    }

    line_positions.push_back(pos);
    goto_next_line(mapped_file, pos, length);
  }
  line_positions.push_back(pos);

  // split the lines into ranges that can be processed in parallel
  const size_t num_ranges = num_parallel_ranges();
  const size_t bytes_per_range = std::max((pos - initial_pos) / num_ranges, UL(1));

  vec<LineRange> result;
  result.reserve(num_ranges);

  size_t line_index = 0;
  while (line_index + 1 < line_positions.size()) {
    const size_t start_index = line_index;
    const size_t start_pos = line_positions[start_index];
    size_t delta_pos = 0;
    do {
      ++line_index;
      delta_pos = line_positions[line_index] - start_pos;
    } while (line_index + 1 < line_positions.size() && delta_pos < bytes_per_range);
    result.push_back(LineRange{ start_pos, start_pos + delta_pos, start_index, line_index - start_index });
  }
  ASSERT(result.back().end == line_positions.back());
  return result;
}

template<typename T>
void copy_to_global_list(const vec<vec<T>>& source, vec<T>& target) {
  vec<size_t> prefix_sum;
  prefix_sum.reserve(source.size() + 1);
  prefix_sum.push_back(0);
  size_t sum = 0;
  for (const auto& local_list: source) {
    sum += local_list.size();
    prefix_sum.push_back(sum);
  }

  tbb::parallel_for(UL(0), source.size(), [&](const size_t i) {
    const size_t start = prefix_sum[i];
    const size_t end = prefix_sum[i + 1];
    for (size_t j = 0; start + j < end; ++j) {
      ASSERT(start + j < target.size() && j < source[i].size());
      target[start + j] = source[i][j];
    }
  });
}

}  // namespace io
}  // namespace mt_kahypar

/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2025 Patrick Steil <patrick@steil.dev>
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

#include <cstdlib>

#ifdef _WIN32
#include <io.h>
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <unistd.h>
#endif

namespace mt_kahypar::color {

namespace internal {
inline bool stdoutSupportsColor() {
  // https://no-color.org
  const char* no_color = std::getenv("NO_COLOR");
  if ( no_color != nullptr && no_color[0] != '\0' ) {
    return false;
  }
#ifdef _WIN32
  if ( !_isatty(_fileno(stdout)) ) {
    return false;
  }
  // ANSI escape sequences have to be enabled explicitly on Windows
  HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
  DWORD mode = 0;
  return handle != INVALID_HANDLE_VALUE && GetConsoleMode(handle, &mode) &&
         SetConsoleMode(handle, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
#else
  return isatty(STDOUT_FILENO);
#endif
}

inline bool& enabledFlag() {
  static bool enabled = stdoutSupportsColor();
  return enabled;
}
}  // namespace internal

inline bool enabled() {
  return internal::enabledFlag();
}

inline void setEnabled(const bool enabled) {
  internal::enabledFlag() = enabled;
}

// Returns the escape sequence if colored output is enabled and an empty string otherwise
inline const char* code(const char* escape_sequence) {
  return enabled() ? escape_sequence : "";
}

}  // namespace mt_kahypar::color

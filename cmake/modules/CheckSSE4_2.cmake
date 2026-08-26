# Check if the CPU provides fast operations
# for popcount, leftmost and rightmost bit

set(BUILTIN_POPCNT 0)
# Check if we are on a Linux system
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
	# Use /proc/cpuinfo to get the information
	file(STRINGS "/proc/cpuinfo" _cpuinfo)
	if(_cpuinfo MATCHES "(sse4_2)|(sse4a)")
		set(BUILTIN_POPCNT 1)
	endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    execute_process(
        COMMAND powershell.exe
            -NoProfile
            -NonInteractive
            -ExecutionPolicy Bypass
            -File "${CMAKE_CURRENT_LIST_DIR}/../check-ise.ps1"
        OUTPUT_VARIABLE CPU_FEATURES
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    # The MSVC STL assumes that popcount is available if AVX is supported.
    # Since we use BUILTIN_POPCOUNT for general SSE4.2 availability detection,
    # only set it if both features are supported.
    if(("sse4_2" IN_LIST CPU_FEATURES) AND ("avx" IN_LIST CPU_FEATURES))
        set(BUILTIN_POPCNT 1)
    endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
#  handle MacOs
execute_process(COMMAND sysctl -n machdep.cpu.features
                OUTPUT_VARIABLE _cpuinfo OUTPUT_STRIP_TRAILING_WHITESPACE)
	if(_cpuinfo MATCHES "SSE4.2")
		set(BUILTIN_POPCNT 1)
	endif()
endif()

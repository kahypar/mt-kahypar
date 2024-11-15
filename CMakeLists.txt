cmake_minimum_required(VERSION 3.16)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

project(MtKaHyPar
        DESCRIPTION "Mt-KaHyPar: Multi-Threaded Karlsruhe Hypergraph Partitioning"
        LANGUAGES CXX C)
set(PROJECT_VENDOR "Tobias Heuer")
set(PROJECT_CONTACT "tobias.heuer@kit.edu")
set(PROJECT_URL "https://github.com/kahypar/mt-kahypar")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


#################################################################
## Sanity checks and policies                                  ##
#################################################################

if(NOT CMAKE_BUILD_TYPE)
  message(FATAL_ERROR "CMAKE_BUILD_TYPE must be set. Options are: Debug, Release, RelWithDebInfo")
endif()
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "Build directory must be different from source directory.")
endif()
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  message(FATAL_ERROR "MSVC not supported at the moment")
endif()

if (POLICY CMP0135)  # policy for download timestamps
  cmake_policy(SET CMP0135 NEW)
endif()


#################################################################
## Meta targets                                               ##
#################################################################

# meta targets for flags
add_library(MtKaHyPar-BuildFlags INTERFACE)
set_target_properties(MtKaHyPar-BuildFlags PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)

# meta targets for includes and source files
add_library(MtKaHyPar-Include INTERFACE)
target_include_directories(MtKaHyPar-Include INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(MtKaHyPar-Sources INTERFACE)

add_library(MtKaHyPar-ToolsSources INTERFACE)

# meta target for CLI/test build
add_library(MtKaHyPar-BuildSources OBJECT "")
target_link_libraries(MtKaHyPar-BuildSources PRIVATE MtKaHyPar-Sources)
target_link_libraries(MtKaHyPar-BuildSources PUBLIC MtKaHyPar-Include MtKaHyPar-BuildFlags)

# meta target for library build, which must be built with all features enabled
add_library(MtKaHyPar-LibraryBuildSources OBJECT "")
target_link_libraries(MtKaHyPar-LibraryBuildSources PRIVATE MtKaHyPar-Sources)
target_link_libraries(MtKaHyPar-LibraryBuildSources PUBLIC MtKaHyPar-Include MtKaHyPar-BuildFlags)
target_compile_definitions(MtKaHyPar-LibraryBuildSources PUBLIC MT_KAHYPAR_LIBRARY_MODE
                           KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
                           KAHYPAR_ENABLE_SOED_METRIC KAHYPAR_ENABLE_STEINER_TREE_METRIC)
# building shared libraries requires position independent code
if(BUILD_SHARED_LIBS)
  set_target_properties(MtKaHyPar-LibraryBuildSources PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

#################################################################
## Options                                                     ##
#################################################################

option(KAHYPAR_ENABLE_TESTING
  "Enables tests, which requires dowloading googletest." OFF)

option(KAHYPAR_DOWNLOAD_BOOST
  "Download boost automatically and compile required libraries." OFF)

option(KAHYPAR_DOWNLOAD_TBB
  "Download TBB automatically." OFF)

option(KAHYPAR_USE_GCOV
  "Compile and run tests with gcov for coverage analysis." OFF)

option(KAHYPAR_DISABLE_ASSERTIONS
  "Disable KaHyPar's internal assertions." OFF)

option(KAHYPAR_USE_STANDARD_ASSERTIONS
  "Use standard C++ asserts instead of custom assertions." OFF)

option(KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS
  "Enable heavy assertions in preprocessing phase." OFF)

option(KAHYPAR_ENABLE_HEAVY_COARSENING_ASSERTIONS
  "Enable heavy assertions in coarsening phase." OFF)

option(KAHYPAR_ENABLE_HEAVY_INITIAL_PARTITIONING_ASSERTIONS
  "Enable heavy assertions in initial partitioning phase." OFF)

option(KAHYPAR_ENABLE_HEAVY_REFINEMENT_ASSERTIONS
  "Enable heavy assertions in refinement phase." OFF)

option(KAHYPAR_USE_64_BIT_IDS
  "Enables 64-bit vertex and hyperedge IDs." OFF)

option(KAHYPAR_CI_BUILD
  "Indicate that this build is executed on GitHub Actions." OFF)

option(KAHYPAR_ADD_ADDRESS_SANITIZER
  "Adds address sanitizer to compile options." ON)

option(KAHYPAR_ENABLE_THREAD_PINNING
  "Enables thread pinning in Mt-KaHyPar." ON)

option(KAHYPAR_ENABLE_EXTENDED_INSTRUCTIONS
  "Depeding on detected compiler features, allows instructions that might not be perfectly portable: `-mcx16 -msse4.2 -mcrc32`" ON)

option(KAHYPAR_ENABLE_ARCH_COMPILE_OPTIMIZATIONS
  "Adds the compile flags `-mtune=native -march=native`" ON)

option(KAHYPAR_ENABLE_EXPERIMENTAL_FEATURES
  "Enables some experimental features. Can be turned off for faster compilation." OFF)

option(KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
  "Enables graph partitioning features. Can be turned off for faster compilation." OFF)

option(KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES
  "Enables highest quality preset features. Can be turned off for faster compilation." OFF)

option(KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
  "Enables large k partitioning features. Can be turned off for faster compilation." OFF)

option(KAHYPAR_ENABLE_SOED_METRIC
  "Enables the sum-of-external-degree metric. Can be turned off for faster compilation." OFF)

option(KAHYPAR_ENABLE_STEINER_TREE_METRIC
  "Enables the Steiner tree metric. Can be turned off for faster compilation." OFF)

option(KAHYPAR_PYTHON
  "Build the Python interface. Can be turned off in case Python is not available." OFF)

option(MT_KAHYPAR_DISABLE_BOOST
  "Whether to exclude components requiring Boost::program_options. Will result in no binary target and the C and Python interface not being able to load configuration files." OFF)


if(KAHYPAR_ENABLE_STEINER_TREE_METRIC AND NOT KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES)
  message(FATAL_ERROR "Steiner tree metric requires graph features. Add -DKAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES=On to your cmake command")
endif()


#################################################################
## Some options translate into compile flags                   ##
#################################################################

if(KAHYPAR_DISABLE_ASSERTIONS)
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE KAHYPAR_DISABLE_ASSERTIONS)
endif(KAHYPAR_DISABLE_ASSERTIONS)

if(KAHYPAR_USE_STANDARD_ASSERTIONS)
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE KAHYPAR_USE_STANDARD_ASSERTIONS)
endif(KAHYPAR_USE_STANDARD_ASSERTIONS)

if(KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS)
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS)
endif(KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS)

if(KAHYPAR_ENABLE_HEAVY_COARSENING_ASSERTIONS)
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE KAHYPAR_ENABLE_HEAVY_COARSENING_ASSERTIONS)
endif(KAHYPAR_ENABLE_HEAVY_COARSENING_ASSERTIONS)

if(KAHYPAR_ENABLE_HEAVY_INITIAL_PARTITIONING_ASSERTIONS)
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE KAHYPAR_ENABLE_HEAVY_INITIAL_PARTITIONING_ASSERTIONS)
endif(KAHYPAR_ENABLE_HEAVY_INITIAL_PARTITIONING_ASSERTIONS)

if(KAHYPAR_ENABLE_HEAVY_REFINEMENT_ASSERTIONS)
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE KAHYPAR_ENABLE_HEAVY_REFINEMENT_ASSERTIONS)
endif(KAHYPAR_ENABLE_HEAVY_REFINEMENT_ASSERTIONS)

if(KAHYPAR_USE_64_BIT_IDS)
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE KAHYPAR_USE_64_BIT_IDS)
endif(KAHYPAR_USE_64_BIT_IDS)

if(KAHYPAR_ENABLE_THREAD_PINNING)
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE KAHYPAR_ENABLE_THREAD_PINNING)
endif(KAHYPAR_ENABLE_THREAD_PINNING)

if(KAHYPAR_ENABLE_EXPERIMENTAL_FEATURES)
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE KAHYPAR_ENABLE_EXPERIMENTAL_FEATURES)
endif(KAHYPAR_ENABLE_EXPERIMENTAL_FEATURES)

if(KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES)
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES)
endif(KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES)

if(KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES)
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES)
endif(KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES)

if(KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES)
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES)
endif(KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES)

if(KAHYPAR_ENABLE_SOED_METRIC)
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE KAHYPAR_ENABLE_SOED_METRIC)
endif(KAHYPAR_ENABLE_SOED_METRIC)

if(KAHYPAR_ENABLE_STEINER_TREE_METRIC)
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE KAHYPAR_ENABLE_STEINER_TREE_METRIC)
endif(KAHYPAR_ENABLE_STEINER_TREE_METRIC)

if(MT_KAHYPAR_DISABLE_BOOST)
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE MT_KAHYPAR_DISABLE_BOOST)
endif()


#################################################################
## Setup of source files and dependencies                      ##
#################################################################

include(FetchContent)


if (KAHYPAR_ENABLE_TESTING)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.15.2
  )
  FetchContent_MakeAvailable(googletest)
  include(gmock)
  enable_testing()

  add_library(MtKaHyPar-Test INTERFACE)
  target_link_libraries(MtKaHyPar-Test INTERFACE gmock gtest gtest_main)
endif()

# Include Submodules
target_include_directories(MtKaHyPar-Include INTERFACE external_tools/kahypar-shared-resources)
target_include_directories(MtKaHyPar-Include INTERFACE external_tools/growt)
target_include_directories(MtKaHyPar-Include INTERFACE external_tools/WHFC)

if(NOT MT_KAHYPAR_DISABLE_BOOST)
  if(KAHYPAR_DOWNLOAD_BOOST)
    # Download Boost program options
    set(BOOST_INCLUDE_LIBRARIES program_options)
    set(BOOST_ENABLE_CMAKE ON)
    set(FETCHCONTENT_QUIET FALSE)
    FetchContent_Declare(
      Boost EXCLUDE_FROM_ALL SYSTEM
      URL https://github.com/boostorg/boost/releases/download/boost-1.86.0/boost-1.86.0-cmake.tar.gz
    )
    set(FETCHCONTENT_QUIET TRUE)
    FetchContent_MakeAvailable(Boost)
  else()
    # Find Boost program options system library
    find_package(Boost 1.69 REQUIRED COMPONENTS program_options)
    if(NOT Boost_FOUND)
      message(FATAL_ERROR "
        Boost not found. Install Boost on your system or
        add -DKAHYPAR_DOWNLOAD_BOOST=On to the cmake build command.")
    endif()
    message(STATUS "Boost Include: ${Boost_INCLUDE_DIRS}, Boost Library: ${Boost_LIBRARY_DIRS}")
  endif()
  target_link_libraries(MtKaHyPar-Include INTERFACE Boost::program_options)
endif()


if(KAHYPAR_DOWNLOAD_TBB)
  # Download TBB library
  FetchContent_Declare(
    TBB EXCLUDE_FROM_ALL SYSTEM
    GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
    GIT_TAG        v2022.0.0
    GIT_SHALLOW    FALSE  # TBB seems to assume that a git repo is present
  )
  FetchContent_MakeAvailable(TBB)
else()
  # Find system TBB library
  find_package(TBB 2021.5 COMPONENTS tbb tbbmalloc)
  if(NOT TBB_FOUND)
    message(FATAL_ERROR "
      TBB library not found or current TBB version is too old. Install TBB on your system
      or add -DKAHYPAR_DOWNLOAD_TBB=On to the cmake build command.")
  endif()
  get_target_property(KAHYPAR_TBB_INCLUDE_DIRS TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "TBB Version: ${TBB_VERSION_MAJOR}.${TBB_VERSION_MINOR}, TBB Include: ${KAHYPAR_TBB_INCLUDE_DIRS}")
endif()
target_link_libraries(MtKaHyPar-Include INTERFACE TBB::tbb TBB::tbbmalloc)


# Find HWLOC Library
FIND_PATH(HWLOCK_INCLUDE_DIR NAME hwloc.h
  HINTS $ENV{HOME}/local/include /opt/local/include /usr/local/include /usr/include)
FIND_LIBRARY(HWLOCK_LIBRARY NAME hwloc
  HINTS $ENV{HOME}/local/lib64 $ENV{HOME}/local/lib /usr/local/lib64 /usr/local/lib /opt/local/lib64 /opt/local/lib /usr/lib64 /usr/lib
)

IF (HWLOCK_INCLUDE_DIR AND HWLOCK_LIBRARY)
  message(STATUS "Found hwlock library: inc=${HWLOCK_INCLUDE_DIR}, lib=${HWLOCK_LIBRARY}")
  target_include_directories(MtKaHyPar-Include INTERFACE ${HWLOCK_INCLUDE_DIR})
  target_link_libraries(MtKaHyPar-Sources INTERFACE ${HWLOCK_LIBRARY})
ELSE ()
  message(FATAL_ERROR "
    HwLoc library not found. Install HwLoc on your system.")
ENDIF ()


#################################################################
## Compile options and specialized builds                      ##
#################################################################

if(KAHYPAR_USE_GCOV)
  # Add targets for code coverage analysis
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SetupCoverage.cmake)
endif(KAHYPAR_USE_GCOV)

# Get current commit hash
find_package(Git)
include(GetGitRevisionDescription)
get_git_head_revision(MT_KAHYPAR_VERSION_GIT_REFSPEC MT_KAHYPAR_VERSION_GIT_SHA1)
if(MT_KAHYPAR_VERSION_GIT_REFSPEC)
  message(STATUS "Detected git refspec ${MT_KAHYPAR_VERSION_GIT_REFSPEC} sha ${MT_KAHYPAR_VERSION_GIT_SHA1}")
  configure_file(${PROJECT_SOURCE_DIR}/mt-kahypar/application/git_revision.txt.in ${PROJECT_BINARY_DIR}/mt-kahypar/application/git_head.txt)
endif(MT_KAHYPAR_VERSION_GIT_REFSPEC)

if (CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
    set(KAHYPAR_X86 TRUE)
else()
    set(KAHYPAR_X86 FALSE)
endif()

# Add compile flags that enable warnings
target_compile_options(MtKaHyPar-BuildFlags INTERFACE
  -W -Wall -Wextra -Wunused -Wuninitialized -Wfatal-errors -Wcast-qual -Woverloaded-virtual
  -Wredundant-decls -Wno-unused-function -Winit-self -pedantic -DPARANOID -Wno-unused-function)

# lld is way faster than ld. If you have it, use it!
find_program(LLD_BIN lld)
if (LLD_BIN AND KAHYPAR_X86)
  message(STATUS "Found and will use LLVM linker " ${LLD_BIN})
  target_compile_options(MtKaHyPar-BuildFlags INTERFACE -fuse-ld=lld)
else()
  message(STATUS "Default linker")
endif()

target_compile_options(MtKaHyPar-BuildFlags INTERFACE $<$<CONFIG:RELEASE>:-O3>
                                                      $<$<CONFIG:RELWITHDEBINFO>:-g3 -UNDEBUG>  # keep assertions activated
                                                      $<$<CONFIG:DEBUG>:-g3 -fno-omit-frame-pointer>)
target_link_options(MtKaHyPar-BuildFlags INTERFACE $<$<CONFIG:DEBUG>:-fno-omit-frame-pointer>)

if(UNIX AND NOT WIN32)
  target_compile_options(MtKaHyPar-BuildFlags INTERFACE $<$<CONFIG:DEBUG>:-fsanitize=undefined>)
  target_link_options(MtKaHyPar-BuildFlags INTERFACE $<$<CONFIG:DEBUG>:-fsanitize=undefined>)
endif()

if (KAHYPAR_ENABLE_EXTENDED_INSTRUCTIONS)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-mcrc32 KAHYPAR_HAS_CRC32)
  include(CheckSSE4_2)

  if(KAHYPAR_X86)
    target_compile_options(MtKaHyPar-BuildFlags INTERFACE -mcx16)
  endif()
  if(BUILTIN_POPCNT)
    target_compile_options(MtKaHyPar-BuildFlags INTERFACE -msse4.2)
  endif()
  if(KAHYPAR_HAS_CRC32 AND KAHYPAR_X86)
    target_compile_options(MtKaHyPar-BuildFlags INTERFACE -mcrc32)
  endif()
endif()

if(KAHYPAR_ENABLE_ARCH_COMPILE_OPTIMIZATIONS)
  target_compile_options(MtKaHyPar-BuildFlags INTERFACE $<$<CONFIG:RELEASE>:-mtune=native -march=native>
                                                        $<$<CONFIG:RELWITHDEBINFO>:-mtune=native -march=native>)
endif()

if(KAHYPAR_ADD_ADDRESS_SANITIZER AND (NOT KAHYPAR_CI_BUILD))
  target_compile_options(MtKaHyPar-BuildFlags INTERFACE $<$<CONFIG:DEBUG>:-fsanitize=address>)
  target_link_options(MtKaHyPar-BuildFlags INTERFACE $<$<CONFIG:DEBUG>:-fsanitize=address>)
endif()


message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CMAKE_CXX_FLAGS_RELWITHDEBINFO: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")


#################################################################
## Include the source code and targets via subdirectories      ##
#################################################################

if(NOT MT_KAHYPAR_DISABLE_BOOST)
  add_subdirectory(mt-kahypar/application)
endif()

if(KAHYPAR_PYTHON)
  if (NOT BUILD_SHARED_LIBS)
    message(SEND_ERROR "Python interface must be built as shared library.
                        Add -DBUILD_SHARED_LIBS=On to your cmake command.")
  endif()
  add_subdirectory(python)
endif()

if (KAHYPAR_ENABLE_TESTING)
  add_subdirectory(tests)
endif()

add_subdirectory(tools)
add_subdirectory(lib)
add_subdirectory(mt-kahypar)

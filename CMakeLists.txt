cmake_minimum_required(VERSION 3.26)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

# define the project version
set(MT_KAHYPAR_VERSION "1.5.1")
string(SUBSTRING "${MT_KAHYPAR_VERSION}" 0 1 MT_KAHYPAR_SO_VERSION)
set(MTKAHYPAR_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

project(MtKaHyPar
        VERSION ${MT_KAHYPAR_VERSION}
        DESCRIPTION "Mt-KaHyPar: Multi-Threaded Karlsruhe Hypergraph Partitioning"
        LANGUAGES CXX C)
set(PROJECT_URL "https://github.com/kahypar/mt-kahypar")


#################################################################
## Sanity checks and policies                                  ##
#################################################################

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
  message(FATAL_ERROR "CMAKE_BUILD_TYPE must be set. Options are: Debug, Release, RelWithDebInfo")
endif()
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "Build directory must be different from source directory.")
endif()
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  message(FATAL_ERROR "MSVC not supported at the moment")
endif()

if (POLICY CMP0135)  # download timestamps
  cmake_policy(SET CMP0135 NEW)
endif()
if (POLICY CMP0144)  # uppercase find_package variables
  cmake_policy(SET CMP0144 NEW)
endif()
if (POLICY CMP0167)  # new FindBoost behavior
  cmake_policy(SET CMP0167 NEW)
endif()

if (NOT DEFINED BUILD_SHARED_LIBS AND NOT WIN32)
  # since TBB does not support static linking, we build with dynamic linking by default (except on Windows)
  set(BUILD_SHARED_LIBS ON)
endif()

#################################################################
## Meta targets                                               ##
#################################################################

# meta targets for flags
add_library(MtKaHyPar-BuildFlags INTERFACE)
target_compile_features(MtKaHyPar-BuildFlags INTERFACE cxx_std_17)

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

# general build options
option(KAHYPAR_PYTHON "Include the Python interface in the build." OFF)
option(KAHYPAR_ENABLE_TESTING "Enables tests, which requires dowloading googletest." OFF)
option(KAHYPAR_STATIC_LINK_DEPENDENCIES "In static build, also link dependencies (other than TBB) statically." OFF)
option(KAHYPAR_STATIC_LINK_TBB "In static build, also link TBB statically. Note that this is not officially supported!" OFF)
option(KAHYPAR_INSTALL_CLI "Provide a target to install an executable binary `mtkahypar`." OFF)
option(KAHYPAR_BUILD_DEBIAN_PACKAGE "Provide a target to build a debian package." OFF)

# dependencies
option(KAHYPAR_DOWNLOAD_BOOST "Download boost automatically and compile required libraries." OFF)
option(KAHYPAR_DOWNLOAD_TBB "Download TBB automatically." OFF)
option(KAHYPAR_DISABLE_HWLOC "Exclude components requiring hwloc. Might impact running time on NUMA platforms." OFF)
option(KAHYPAR_DISABLE_GROWT "Remove the growt dependency. Might significantly impact running time of the Steiner tree metric." OFF)
option(KAHYPAR_DISABLE_PARLAY "Remove the parlay dependency. Might impact running time of deterministic mode on some instances." OFF)

# specific compile features/build properties
option(KAHYPAR_USE_64_BIT_IDS "Enables 64-bit vertex and hyperedge IDs." OFF)
option(KAHYPAR_USE_ADDRESS_SANITIZER "Adds address sanitizer to compile options." OFF)
option(KAHYPAR_ENABLE_EXTENDED_INSTRUCTIONS "Allows instructions that might not be fully portable: `-mcx16 -msse4.2 -mcrc32`" OFF)
option(KAHYPAR_ENABLE_ARCH_COMPILE_OPTIMIZATIONS "Adds the compile flags `-mtune=native -march=native`" OFF)
option(KAHYPAR_ENABLE_THREAD_PINNING "Enables thread pinning in Mt-KaHyPar." OFF)

# algorithm features for CLI build (note: the library always contains all non-experimental features)
option(KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES "Enables graph partitioning features. Can be turned off for faster compilation." OFF)
option(KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES "Enables highest quality preset features. Can be turned off for faster compilation." OFF)
option(KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES "Enables large k partitioning features. Can be turned off for faster compilation." OFF)
option(KAHYPAR_ENABLE_SOED_METRIC "Enables the sum-of-external-degree metric. Can be turned off for faster compilation." OFF)
option(KAHYPAR_ENABLE_STEINER_TREE_METRIC "Enables the Steiner tree metric. Can be turned off for faster compilation." OFF)
option(KAHYPAR_ENABLE_EXPERIMENTAL_FEATURES "Enables some experimental features. Can be turned off for faster compilation." OFF)

# assertions
option(KAHYPAR_DISABLE_ASSERTIONS "Disable all internal assertions." OFF)
option(KAHYPAR_USE_STANDARD_ASSERTIONS "Use standard C++ asserts instead of custom assertions." OFF)
option(KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS "Enable heavy assertions in preprocessing phase." OFF)
option(KAHYPAR_ENABLE_HEAVY_COARSENING_ASSERTIONS "Enable heavy assertions in coarsening phase." OFF)
option(KAHYPAR_ENABLE_HEAVY_INITIAL_PARTITIONING_ASSERTIONS "Enable heavy assertions in initial partitioning phase." OFF)
option(KAHYPAR_ENABLE_HEAVY_REFINEMENT_ASSERTIONS "Enable heavy assertions in refinement phase." OFF)

# developer build options
option(KAHYPAR_CI_BUILD "Indicate that this build is executed on GitHub Actions." OFF)
option(KAHYPAR_USE_GCOV "Compile and run tests with gcov for coverage analysis." OFF)


if(KAHYPAR_ENABLE_STEINER_TREE_METRIC AND NOT KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES)
  message(FATAL_ERROR "Steiner tree metric requires graph features. Add -DKAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES=On to your cmake command")
endif()


#################################################################
## Print header with most important infos                      ##
#################################################################

if(CMAKE_BUILD_TYPE)
  # normalize so that exactly the first letter is uppercase
  string(TOUPPER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_UPPER)
  string(TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_LOWER)
  string(SUBSTRING "${BUILD_TYPE_UPPER}" 0 1 FIRST_LETTER)
  string(SUBSTRING "${BUILD_TYPE_LOWER}" 1 -1 REMAINDER)
  set(KAHYPAR_BUILD_DESCRIPTION "${FIRST_LETTER}${REMAINDER} Build")
  if (KAHYPAR_ENABLE_ARCH_COMPILE_OPTIMIZATIONS AND NOT BUILD_TYPE_UPPER MATCHES "DEBUG")
    set(KAHYPAR_BUILD_DESCRIPTION "${KAHYPAR_BUILD_DESCRIPTION} + march=native")
  endif()
  if(KAHYPAR_ENABLE_TESTING)
    set(KAHYPAR_TESTING_ENABLED "Enabled")
  else()
    set(KAHYPAR_TESTING_ENABLED "Disabled")
  endif()
  if(BUILD_SHARED_LIBS)
    set(KAHYPAR_LINKING_TYPE "Dynamic")
  else()
    set(KAHYPAR_LINKING_TYPE "Static")
  endif()
  message(STATUS "Mt-KaHyPar: ${KAHYPAR_BUILD_DESCRIPTION}, Tests ${KAHYPAR_TESTING_ENABLED}, ${KAHYPAR_LINKING_TYPE} Linking")
endif()


#################################################################
## Compatibility checks which interact with options            ##
#################################################################

if(KAHYPAR_ENABLE_THREAD_PINNING AND KAHYPAR_DISABLE_HWLOC)
  message(WARNING "Thread pinning disabled since hwloc library is disabled.")
  set(KAHYPAR_ENABLE_THREAD_PINNING FALSE CACHE STRING "" FORCE)
endif()

if(UNIX)
  include(CheckIncludeFiles)
  set(CMAKE_REQUIRED_QUIET TRUE)
  set(KAHYPAR_REQUIRED_HEADERS "unistd.h" "sys/mman.h" "sys/ioctl.h")
  check_include_files("${KAHYPAR_REQUIRED_HEADERS}" KAHYPAR_HEADERS_AVAILABLE LANGUAGE CXX)
  if(NOT KAHYPAR_HEADERS_AVAILABLE)
    message(FATAL_ERROR "The following system headers are required, but some are not available: ${KAHYPAR_REQUIRED_HEADERS}")
  endif()

  include(CheckThreadPinning)
  if(KAHYPAR_ENABLE_THREAD_PINNING AND NOT THREAD_PINNING_WORKS)
    message(WARNING "Thread pinning disabled since required system APIs are not available.")
    set(KAHYPAR_ENABLE_THREAD_PINNING FALSE CACHE STRING "" FORCE)
  endif()
elseif(NOT WIN32)
  message(FATAL_ERROR "Only unix-based operating systems and windows are supported.")
endif()


#################################################################
## Options that translate into compile flags                   ##
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

if(KAHYPAR_DISABLE_HWLOC)
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE KAHYPAR_DISABLE_HWLOC)
endif()

if(KAHYPAR_DISABLE_PARLAY)
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE KAHYPAR_DISABLE_PARLAY)
endif()


#################################################################
## Compile options and specialized builds                      ##
#################################################################

if (CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
    set(KAHYPAR_X86 TRUE)
else()
    set(KAHYPAR_X86 FALSE)
endif()

# Add compile flags that enable warnings
target_compile_options(MtKaHyPar-BuildFlags INTERFACE
  -W -Wall -Wextra -Wunused -Wuninitialized -Wfatal-errors -Wcast-qual -Woverloaded-virtual
  -Wredundant-decls -Wno-unused-function -Winit-self -pedantic -DPARANOID -Wno-unused-function)

if(KAHYPAR_USE_GCOV)
  # add coverage anaylsis compile and link flags
  target_compile_options(MtKaHyPar-BuildFlags INTERFACE -fprofile-arcs -ftest-coverage -fprofile-update=atomic)
  target_link_options(MtKaHyPar-BuildFlags INTERFACE -lgcov --coverage)
endif(KAHYPAR_USE_GCOV)

# lld is way faster than ld. If you have it, use it!
find_program(LLD_BIN lld)
if (LLD_BIN AND KAHYPAR_X86)
  message(STATUS "Found and will use LLVM linker " ${LLD_BIN})
  target_compile_options(MtKaHyPar-BuildFlags INTERFACE -fuse-ld=lld)
endif()

target_compile_options(MtKaHyPar-BuildFlags INTERFACE $<$<CONFIG:RELEASE>:-O3>
                                                      $<$<CONFIG:RELWITHDEBINFO>:-g3 -UNDEBUG>  # keep assertions activated
                                                      $<$<CONFIG:DEBUG>:-g3 -fno-omit-frame-pointer>)
target_link_options(MtKaHyPar-BuildFlags INTERFACE $<$<CONFIG:DEBUG>:-fno-omit-frame-pointer>)

if(UNIX AND NOT WIN32)
  target_compile_options(MtKaHyPar-BuildFlags INTERFACE $<$<CONFIG:DEBUG>:-fsanitize=undefined>)
  target_link_options(MtKaHyPar-BuildFlags INTERFACE $<$<CONFIG:DEBUG>:-fsanitize=undefined>)
endif()

if (KAHYPAR_X86 AND KAHYPAR_ENABLE_EXTENDED_INSTRUCTIONS)
  include(CheckCXXCompilerFlag)
  include(CheckSSE4_2)
  block()
    set(CMAKE_REQUIRED_QUIET TRUE)
    check_cxx_compiler_flag(-mcrc32 KAHYPAR_HAS_CRC32)
  endblock()

  target_compile_options(MtKaHyPar-BuildFlags INTERFACE -mcx16)
  if(KAHYPAR_HAS_CRC32)
    target_compile_options(MtKaHyPar-BuildFlags INTERFACE -mcrc32)
  endif()
  if(BUILTIN_POPCNT)
    target_compile_options(MtKaHyPar-BuildFlags INTERFACE -msse4.2)
  endif()
endif()

# check this after include(CheckSSE4_2) so BUILTIN_POPCNT is available
if (NOT KAHYPAR_DISABLE_GROWT AND KAHYPAR_ENABLE_STEINER_TREE_METRIC
    AND KAHYPAR_X86 AND KAHYPAR_ENABLE_EXTENDED_INSTRUCTIONS AND BUILTIN_POPCNT)
    # growt requires SSE instructions
    set(KAHYPAR_USE_GROWT TRUE)
else()
    set(KAHYPAR_USE_GROWT FALSE)
endif()

if(KAHYPAR_ENABLE_ARCH_COMPILE_OPTIMIZATIONS)
  target_compile_options(MtKaHyPar-BuildFlags INTERFACE $<$<CONFIG:RELEASE>:-mtune=native -march=native>
                                                        $<$<CONFIG:RELWITHDEBINFO>:-mtune=native -march=native>)
endif()

if(KAHYPAR_USE_ADDRESS_SANITIZER AND (NOT KAHYPAR_CI_BUILD))
  target_compile_options(MtKaHyPar-BuildFlags INTERFACE $<$<CONFIG:DEBUG>:-fsanitize=address>)
  target_link_options(MtKaHyPar-BuildFlags INTERFACE $<$<CONFIG:DEBUG>:-fsanitize=address>)
endif()


#################################################################
## Setup of dependencies                                       ##
#################################################################

include(FetchContent)
set(KAHYPAR_SHARED_RESOURCES_TAG 6d5c8e2444e4310667ec1925e995f26179d7ee88)
set(KAHYPAR_WHFC_TAG             30b0eeb0e49577d06c3deb09a44b035d81c529d2)
set(KAHYPAR_GROWT_TAG            0c1148ebcdfd4c04803be79706533ad09cc81d37)
set(KAHYPAR_PARLAY_TAG           e1f1dc0ccf930492a2723f7fbef8510d35bf57f5)
set(KAHYPAR_TBB_VERSION          v2022.0.0)
set(KAHYPAR_GOOGLETEST_VERSION   v1.15.2)
set(KAHYPAR_PYBIND11_VERSION     v2.13.6)

message(STATUS "Fetching dependencies...")

# Include header-only dependencies
FetchContent_Populate(
  kahypar-shared-resources QUIET EXLUDE_FROM_ALL
  GIT_REPOSITORY https://github.com/kahypar/kahypar-shared-resources.git
  GIT_TAG        ${KAHYPAR_SHARED_RESOURCES_TAG}
  SOURCE_DIR     external_tools/kahypar-shared-resources
)
FetchContent_Populate(
  WHFC QUIET EXLUDE_FROM_ALL
  GIT_REPOSITORY https://github.com/larsgottesbueren/WHFC.git
  GIT_TAG        ${KAHYPAR_WHFC_TAG}
  SOURCE_DIR     external_tools/WHFC
)
target_include_directories(MtKaHyPar-Include INTERFACE
                           ${CMAKE_CURRENT_BINARY_DIR}/external_tools/kahypar-shared-resources
                           ${CMAKE_CURRENT_BINARY_DIR}/external_tools)

if(KAHYPAR_USE_GROWT)
  FetchContent_Populate(
    growt QUIET EXLUDE_FROM_ALL
    GIT_REPOSITORY https://github.com/TooBiased/growt.git
    GIT_TAG        ${KAHYPAR_GROWT_TAG}
    SOURCE_DIR     external_tools/growt
  )
  target_compile_definitions(MtKaHyPar-BuildFlags INTERFACE KAHYPAR_USE_GROWT)
  target_include_directories(MtKaHyPar-Include INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/external_tools/growt)
endif()

if(NOT KAHYPAR_DISABLE_PARLAY)
  FetchContent_Populate(
    parlay QUIET EXLUDE_FROM_ALL
    GIT_REPOSITORY https://github.com/cmuparlay/parlaylib.git
    GIT_TAG        ${KAHYPAR_PARLAY_TAG}
    SOURCE_DIR     external_tools/parlay
  )
  target_include_directories(MtKaHyPar-Include INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/external_tools/parlay/include)
  target_compile_options(MtKaHyPar-BuildFlags INTERFACE -DPARLAY_TBB)
endif()


if (KAHYPAR_ENABLE_TESTING)
  FetchContent_Declare(
    googletest EXCLUDE_FROM_ALL SYSTEM
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        ${KAHYPAR_GOOGLETEST_VERSION}
  )
  FetchContent_MakeAvailable(googletest)
  include(gmock)
  enable_testing()

  add_library(MtKaHyPar-Test INTERFACE)
  target_link_libraries(MtKaHyPar-Test INTERFACE gmock gtest gtest_main)
endif()

# check for linking problems
if(BUILD_SHARED_LIBS AND KAHYPAR_STATIC_LINK_DEPENDENCIES)
  if(NOT KAHYPAR_DOWNLOAD_BOOST OR NOT CMAKE_POSITION_INDEPENDENT_CODE)
    message(WARNING "Building a shared lib with static linking of transitive dependencies requires building Boost from source with -DCMAKE_POSITION_INDEPENDENT_CODE=On")
  endif()
endif()

block()
  # define required boost libraries
  set(BOOST_INCLUDE_LIBRARIES program_options range dynamic_bitset lexical_cast)

  if(KAHYPAR_STATIC_LINK_DEPENDENCIES)
    set(Boost_USE_STATIC_LIBS ON)
    set(BUILD_SHARED_LIBS OFF)
  endif()
  if(KAHYPAR_DOWNLOAD_BOOST)
    # Download Boost
    set(BOOST_ENABLE_CMAKE ON)
    FetchContent_Declare(
      Boost EXCLUDE_FROM_ALL SYSTEM
      URL https://github.com/boostorg/boost/releases/download/boost-1.86.0/boost-1.86.0-cmake.tar.gz
    )
    FetchContent_MakeAvailable(Boost)
  else()
    # Find system Boost
    find_package(Boost 1.69 REQUIRED COMPONENTS program_options)
    if(NOT Boost_FOUND)
      message(FATAL_ERROR "
        Boost not found. Install Boost on your system or
        add -DKAHYPAR_DOWNLOAD_BOOST=On to the cmake build command.")
    endif()
    message(STATUS "Boost Include: ${Boost_INCLUDE_DIRS}, Boost Library: ${Boost_LIBRARY_DIRS}")
  endif()

  # newer boost versions require to add all targets explicitly
  foreach(BOOST_TARGET ${BOOST_INCLUDE_LIBRARIES})
    if(TARGET Boost::${BOOST_TARGET})
      target_link_libraries(MtKaHyPar-Include INTERFACE Boost::${BOOST_TARGET})
    endif()
  endforeach()
endblock()


if(KAHYPAR_DOWNLOAD_TBB)
  # Download TBB library
  FetchContent_Declare(
    TBB EXCLUDE_FROM_ALL SYSTEM
    GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
    GIT_TAG        ${KAHYPAR_TBB_VERSION}
    GIT_SHALLOW    FALSE  # TBB seems to assume that a git repo is present
  )
  block()
    if(KAHYPAR_STATIC_LINK_TBB)
      set(BUILD_SHARED_LIBS OFF)
    else()
      set(BUILD_SHARED_LIBS ON)
    endif()
    FetchContent_MakeAvailable(TBB)
  endblock()
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
if(NOT KAHYPAR_DISABLE_HWLOC)
  include(DetectHwloc)
  if(NOT TARGET HWLOC::hwloc)
    message(FATAL_ERROR "Hwloc library not found. Install hwloc on your system.")
  endif()
  target_link_libraries(MtKaHyPar-Include INTERFACE HWLOC::hwloc)
endif()


#################################################################
## Include the source code and targets via subdirectories      ##
#################################################################

if(KAHYPAR_PYTHON)
  if (NOT BUILD_SHARED_LIBS)
    message(SEND_ERROR "Python interface must be built as shared library. To do so, add -DBUILD_SHARED_LIBS=On to your cmake command.")
  endif()
  add_subdirectory(python)
endif()
if (KAHYPAR_ENABLE_TESTING)
  add_subdirectory(tests)
endif()

add_subdirectory(mt-kahypar/application)
add_subdirectory(tools)
add_subdirectory(lib)
add_subdirectory(mt-kahypar)

# export target for C interface
add_library(MtKaHyPar::mtkahypar ALIAS mtkahypar)


#################################################################
## Installation                                                ##
#################################################################

# only create installation targets if this is the top level project
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  if(UNIX AND CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    # set default installation prefix if not user-defined
    set(CMAKE_INSTALL_PREFIX "/usr/" CACHE STRING "" FORCE)
  endif()

  if(KAHYPAR_INSTALL_CLI)
    # CLI installation target
    include(GNUInstallDirs)
    set_target_properties(MtKaHyPar-CLI PROPERTIES OUTPUT_NAME mtkahypar)
    install(TARGETS MtKaHyPar-CLI
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      RENAME mtkahypar
      COMPONENT MtKaHyPar_CLI)
    add_custom_target(install-mtkahypar-cli
      ${CMAKE_COMMAND}
      -DBUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_INSTALL_COMPONENT=MtKaHyPar_CLI
      -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
      DEPENDS MtKaHyPar-CLI)
  endif()

  if(NOT WIN32 AND BUILD_SHARED_LIBS)
    # library installation target
    include(SetupInstallation)
  endif()

  if(KAHYPAR_BUILD_DEBIAN_PACKAGE)
    # packaging via CPack
    include(SetupCPack)
  endif()

  if(KAHYPAR_SETUP_PYTHON AND DEFINED KAHYPAR_TBB_DIR)
    # preparation for python wheel
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml.in ${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml @ONLY)
  endif()
  if(KAHYPAR_PYTHON AND DEFINED SKBUILD_PROJECT_NAME)
    # build python wheel
    if(NOT KAHYPAR_DOWNLOAD_TBB OR NOT DEFINED KAHYPAR_TBB_DIR)
      message(FATAL_ERROR "TBB must be built from source for python wheel and KAHYPAR_TBB_DIR must be specified.")
    endif()
    message(STATUS "Python wheel build enabled, installing TBB to ${KAHYPAR_TBB_DIR}")

    install(TARGETS mtkahypar_python
      DESTINATION .
      COMPONENT MtKaHyPar_Python)
    install(TARGETS tbb tbbmalloc
      DESTINATION ${KAHYPAR_TBB_DIR}
      COMPONENT MtKaHyPar_Python)
  endif()
endif()

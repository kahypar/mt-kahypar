cmake_minimum_required(VERSION 3.16)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

project(MtKaHyPar
        DESCRIPTION "Mt-KaHyPar: Multi-Threaded Karlsruhe Hypergraph Partitioning"
        LANGUAGES CXX C)
set(PROJECT_VENDOR "Tobias Heuer")
set(PROJECT_CONTACT "tobias.heuer@kit.edu")
set(PROJECT_URL "https://github.com/kahypar/mt-kahypar")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# Sanity checks
if(NOT CMAKE_BUILD_TYPE)
  message(FATAL_ERROR "CMAKE_BUILD_TYPE must be set. Options are: Debug, Release, RelWithDebInfo")
endif()
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "Build directory must be different from source directory.")
endif()
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  message(FATAL_ERROR "MSVC not supported at the moment")
endif()


#################################################################
## Meta targets                                               ##
#################################################################

# meta targets for flags
add_library(MtKaHyPar-BuildFlags INTERFACE)
set_target_properties(MtKaHyPar-BuildFlags PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)

add_library(MtKaHyPar-WarningFlags INTERFACE)

# meta targets for includes and source files
add_library(MtKaHyPar-Include INTERFACE)
target_include_directories(MtKaHyPar-Include INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(MtKaHyPar-Sources INTERFACE)

# meta target for CLI/test build
add_library(MtKaHyPar-BuildSources OBJECT "")
target_link_libraries(MtKaHyPar-BuildSources PRIVATE MtKaHyPar-Sources)
target_link_libraries(MtKaHyPar-BuildSources PUBLIC MtKaHyPar-Include MtKaHyPar-BuildFlags MtKaHyPar-WarningFlags)

# meta target for library build, which must be built with all features enabled
add_library(MtKaHyPar-LibraryBuildSources OBJECT "")
target_link_libraries(MtKaHyPar-LibraryBuildSources PRIVATE MtKaHyPar-Sources)
target_link_libraries(MtKaHyPar-LibraryBuildSources PUBLIC MtKaHyPar-Include MtKaHyPar-BuildFlags MtKaHyPar-WarningFlags)
target_compile_definitions(MtKaHyPar-LibraryBuildSources PUBLIC MT_KAHYPAR_LIBRARY_MODE
                           KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
                           KAHYPAR_ENABLE_SOED_METRIC KAHYPAR_ENABLE_STEINER_TREE_METRIC)


#################################################################
## Options                                                     ##
#################################################################

option(KAHYPAR_DOWNLOAD_BOOST
  "Download boost automatically and compile required libraries." OFF)

option(KAHYPAR_DOWNLOAD_TBB
  "Download TBB automatically." OFF)

option(KAHYPAR_ENFORCE_MINIMUM_TBB_VERSION
  "Enforces the minimum required TBB version." ON)

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

target_include_directories(MtKaHyPar-Include INTERFACE ${PROJECT_SOURCE_DIR})

find_package(Threads REQUIRED)
message(STATUS "Found Threads: ${CMAKE_THREAD_LIBS_INIT}")  # TODO: do we need this?

# TODO: guard behind option (?)
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

# Include Submodules
target_include_directories(MtKaHyPar-Include INTERFACE external_tools/kahypar-shared-resources)
target_include_directories(MtKaHyPar-Include INTERFACE external_tools/growt)
target_include_directories(MtKaHyPar-Include INTERFACE external_tools/WHFC)

if(NOT MT_KAHYPAR_DISABLE_BOOST)
  if(KAHYPAR_DOWNLOAD_BOOST)
    # Download Boost
    execute_process(COMMAND cmake -P ${CMAKE_CURRENT_SOURCE_DIR}/scripts/download_boost.cmake)
    file(GLOB MTKAHYPAR_MINI_BOOST_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/external_tools/boost/libs/program_options/src/*.cpp)

    add_library(MtKaHyPar-Boost STATIC ${MTKAHYPAR_MINI_BOOST_SOURCES})
    set_target_properties(MtKaHyPar-Boost PROPERTIES LINKER_LANGUAGE CXX)
    target_include_directories(MtKaHyPar-Boost PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/external_tools/boost/boost/)
    target_link_libraries(MtKaHyPar-Boost PRIVATE MtKaHyPar-BuildFlags)
    set(Boost_LIBRARIES MtKaHyPar-Boost)
  else()
    # Find Boost Program Options Library
    find_package(Boost 1.69 REQUIRED COMPONENTS program_options)
    if(Boost_FOUND)
      message(STATUS "Boost Include: ${Boost_INCLUDE_DIRS}")
      message(STATUS "Boost Library Dirs: ${Boost_LIBRARY_DIRS}")
      message(STATUS "Boost Libraries: ${Boost_LIBRARIES}")
      if(WIN32)
        # TODO
        add_definitions(-DBOOST_ALL_NO_LIB)
        add_definitions(-DBOOST_PROGRAM_OPTIONS_DYN_LINK=1)
      endif()
    else()
      MESSAGE(FATAL_ERROR "
        Boost not found. Install Boost on your system or
        add -DKAHYPAR_DOWNLOAD_BOOST=On to the cmake build command.")
    endif()
  endif()
  target_link_libraries(MtKaHyPar-Sources INTERFACE ${Boost_LIBRARIES})
endif()

# Download TBB
if(KAHYPAR_DOWNLOAD_TBB)
  if ( UNIX AND NOT WIN32 )
    execute_process(COMMAND cmake -P ${CMAKE_CURRENT_SOURCE_DIR}/scripts/download_tbb_linux.cmake)
    set(TBB_ROOT ${CMAKE_CURRENT_BINARY_DIR}/external_tools/tbb)
  elseif ( MSVC )
    execute_process(COMMAND cmake -G ${CMAKE_GENERATOR} -P ${CMAKE_CURRENT_SOURCE_DIR}/scripts/download_tbb_windows.cmake)
    set(TBB_ROOT ${CMAKE_CURRENT_BINARY_DIR}/external_tools/tbb)
  else()
    MESSAGE(WARNING "
      We did not find a TBB version for your platform online.
      The build uses the TBB version installed on your system.")
  endif()
endif()

# Find TBB Library
find_package(TBB REQUIRED)
if (TBB_FOUND)
  if ((TBB_INTERFACE_VERSION GREATER_EQUAL 12050) OR (NOT KAHYPAR_ENFORCE_MINIMUM_TBB_VERSION) OR KAHYPAR_CI_BUILD)
    target_include_directories(MtKaHyPar-Include INTERFACE ${TBB_INCLUDE_DIRS} ${TBB_MALLOC_INCLUDE_DIRS})
    target_link_libraries(MtKaHyPar-Sources INTERFACE ${TBB_LIBRARIES} ${TBB_MALLOC_LIBRARIES})
    MESSAGE(STATUS "TBB Version: ${TBB_VERSION_MAJOR}.${TBB_VERSION_MINOR} (Interface Version=${TBB_INTERFACE_VERSION})")
    MESSAGE(STATUS "Found TBB library: inc=${TBB_INCLUDE_DIRS}, lib=${TBB_LIBRARIES}")
    MESSAGE(STATUS "Found TBB Malloc Proxy library: inc=${TBB_MALLOC_PROXY_INCLUDE_DIRS}, lib=${TBB_MALLOC_PROXY_LIBRARIES}")
  ELSE ()
    MESSAGE(FATAL_ERROR "
      Your current TBB version is too old. Your interface version is ${TBB_INTERFACE_VERSION}, but the minimum
      required version is 12050 (OneTBB 2021.5). You can download a suitable version by adding
      -DKAHYPAR_DOWNLOAD_TBB=ON to the cmake build command (may require to clear the cmake cache).
      Note that Mt-KaHyPar also compiles with older versions of TBB. However, we found some unexpected
      behavior of a TBB function that causes on our side segmentation faults in really rare cases.
      If you want to ignore these warnings, you can add -DKAHYPAR_ENFORCE_MINIMUM_TBB_VERSION=OFF to the
      cmake build command which will build Mt-KaHyPar with your current TBB version.")
  ENDIF()
ELSE ()
  MESSAGE(FATAL_ERROR "
    TBB library not found. Install TBB on your system or
    add -DKAHYPAR_DOWNLOAD_TBB=On to the cmake build command
    (only available for Linux and MSVC).")
ENDIF()

# Find HWLOC Library
FIND_PATH(HWLOCK_INCLUDE_DIR NAME hwloc.h
  HINTS $ENV{HOME}/local/include /opt/local/include /usr/local/include /usr/include)
FIND_LIBRARY(HWLOCK_LIBRARY NAME hwloc
  HINTS $ENV{HOME}/local/lib64 $ENV{HOME}/local/lib /usr/local/lib64 /usr/local/lib /opt/local/lib64 /opt/local/lib /usr/lib64 /usr/lib
)

IF (HWLOCK_INCLUDE_DIR AND HWLOCK_LIBRARY)
  MESSAGE(STATUS "Found hwlock library: inc=${HWLOCK_INCLUDE_DIR}, lib=${HWLOCK_LIBRARY}")
  target_include_directories(MtKaHyPar-Include INTERFACE ${HWLOCK_INCLUDE_DIR})
  target_link_libraries(MtKaHyPar-Sources INTERFACE ${HWLOCK_LIBRARY})
ELSE ()
  MESSAGE(FATAL_ERROR "
    HwLoc library not found. Install HwLoc on your system.")
ENDIF ()



#################################################################
## Compile options and specialized builds                      ##
#################################################################

# Add targets for code coverage analysis
if(KAHYPAR_USE_GCOV)
  include(CodeCoverage)
  setup_target_for_coverage(${PROJECT_NAME}_coverage tests coverage)

  # find programs
  find_program(GENHTML genhtml)
  find_program(LCOV lcov)

  if(NOT LCOV OR NOT GENHTML)
    message(SEND_ERROR "Coverage analysis requires lcov and genhtml.")
  endif()

  # add coverage anaylsis compile and link flags
  target_compile_options(MtKaHyPar-BuildFlags INTERFACE -fprofile-arcs -ftest-coverage)
  target_link_options(MtKaHyPar-BuildFlags INTERFACE -lgcov --coverage)

  # add cached variable containing parameters for lcov/genhtml
  set(LCOV_FLAGS "" CACHE STRING "parameters for lcov")
  set(GENHTML_FLAGS --legend --no-branch-coverage
    CACHE STRING "parameters for genhtml")
  set(COVERAGE_LCOV_EXCLUDES "external_tools/*" )

  # custom target to run before tests
  add_custom_target(lcov-reset
    COMMAND ${LCOV} -q --directory ${CMAKE_BINARY_DIR} --zerocounters
    COMMENT "Resetting code coverage counters")

  # custom lcov target to run tests
  add_custom_target(lcov-runtests
    COMMAND make
    DEPENDS lcov-reset
    COMMENT "Running all unit tests")

  # get git version description
  execute_process(COMMAND git describe --tags
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GITDESC
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # command sequence to gather, clean and generate HTML coverage report
  add_custom_target(lcov-html
    COMMAND ${LCOV} -q --directory . --capture --output-file lcov.info
    COMMAND ${LCOV} -q --remove lcov.info '/usr/*' '*/extlib/*' ${LCOV_FLAGS} --output-file lcov-clean.info
    COMMAND ${GENHTML} -q -o coverage --title "KaHyPar ${GITDESC}" --prefix ${PROJECT_SOURCE_DIR} ${GENHTML_FLAGS} lcov-clean.info
    DEPENDS lcov-runtests
    COMMENT "Capturing code coverage counters and create HTML coverage report"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  # top-level target to run tests and generate coverage report
  add_custom_target(test-coverage
    COMMENT "Generate HTML coverage report"
    DEPENDS lcov-html)

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
    set (X86 TRUE)
else()
    set (X86 FALSE)
endif()

# Add compile flags
target_compile_options(MtKaHyPar-WarningFlags INTERFACE
  -W -Wall -Wextra -Wunused -Wuninitialized -Wfatal-errors -Wcast-qual -Woverloaded-virtual
  -Wredundant-decls -Wno-unused-function -Winit-self -pedantic -DPARANOID -Wno-unused-function)

# lld is way faster than ld. If you have it, use it!
find_program(LLD_BIN lld)
if (LLD_BIN AND X86)
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
  include(CheckSSE4_2)

  if(X86)
    target_compile_options(MtKaHyPar-BuildFlags INTERFACE -mcx16)
  endif()
  if( BUILTIN_POPCNT )
    target_compile_options(MtKaHyPar-BuildFlags INTERFACE -msse4.2)
  endif()
  check_cxx_compiler_flag(-mcrc32 KAHYPAR_HAS_CRC32)
  if(KAHYPAR_HAS_CRC32 AND X86)
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

add_custom_target(AnalyzeModifiedSources  perl "${PROJECT_SOURCE_DIR}/codestyle/analyze-source.pl" "-w")
add_custom_target(AnalyzeAllSources  perl "${PROJECT_SOURCE_DIR}/codestyle/analyze-source.pl" "-aw")

message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CMAKE_CXX_FLAGS_RELWITHDEBINFO: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")


if(NOT MT_KAHYPAR_DISABLE_BOOST)
  add_subdirectory(mt-kahypar/application)
endif()

if(KAHYPAR_PYTHON)
  add_subdirectory(python)
endif()

add_subdirectory(tools)
add_subdirectory(lib)
add_subdirectory(tests)
add_subdirectory(mt-kahypar)

cmake_minimum_required(VERSION 3.16)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

project(KaHyPar CXX C)
set(PROJECT_VENDOR "Tobias Heuer")
set(PROJECT_CONTACT "tobias.heuer@kit.edu")
set(PROJECT_URL "https://github.com/kittobi1992/mt-kahypar")
set(PROJECT_DESCRIPTION "Mt-KaHyPar: Multi-Threaded Karlsruhe Hypergraph Partitioning")
set(PROJECT_VERSION "1.0.0")

include_directories(${PROJECT_SOURCE_DIR})

find_package(Threads REQUIRED)
message(STATUS "Found Threads: ${CMAKE_THREAD_LIBS_INIT}")

#ensure that gmock is built before tests are executed
add_subdirectory(external_tools/googletest EXCLUDE_FROM_ALL)
include_directories(SYSTEM ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
include_directories(SYSTEM ${gtest_SOURCE_DIR}/../googlemock/include ${gtest_SOURCE_DIR}/../googlemock/)

#include KaHyPar
include_directories(external_tools/kahypar)
#include WHFC
include_directories(external_tools/WHFC)

option(KAHYPAR_USE_MINIMAL_BOOST
  "Download boost automatically and compile required libraries." OFF)

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

option(KAHYPAR_ENABLE_EXPERIMENTAL_FEATURES
  "Enables experimental and non-evaluated features." OFF)

option(KAHYPAR_TRAVIS_BUILD
  "Indicate that this build is executed on Travis CI." OFF)

if(KAHYPAR_DISABLE_ASSERTIONS)
  add_compile_definitions(KAHYPAR_DISABLE_ASSERTIONS)
endif(KAHYPAR_DISABLE_ASSERTIONS)

if(KAHYPAR_USE_STANDARD_ASSERTIONS)
  add_compile_definitions(KAHYPAR_USE_STANDARD_ASSERTIONS)
endif(KAHYPAR_USE_STANDARD_ASSERTIONS)

if(KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS)
  add_compile_definitions(KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS)
endif(KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS)

if(KAHYPAR_ENABLE_HEAVY_COARSENING_ASSERTIONS)
  add_compile_definitions(KAHYPAR_ENABLE_HEAVY_COARSENING_ASSERTIONS)
endif(KAHYPAR_ENABLE_HEAVY_COARSENING_ASSERTIONS)

if(KAHYPAR_ENABLE_HEAVY_INITIAL_PARTITIONING_ASSERTIONS)
  add_compile_definitions(KAHYPAR_ENABLE_HEAVY_INITIAL_PARTITIONING_ASSERTIONS)
endif(KAHYPAR_ENABLE_HEAVY_INITIAL_PARTITIONING_ASSERTIONS)

if(KAHYPAR_ENABLE_HEAVY_REFINEMENT_ASSERTIONS)
  add_compile_definitions(KAHYPAR_ENABLE_HEAVY_REFINEMENT_ASSERTIONS)
endif(KAHYPAR_ENABLE_HEAVY_REFINEMENT_ASSERTIONS)

if(KAHYPAR_USE_64_BIT_IDS)
  add_compile_definitions(KAHYPAR_USE_64_BIT_IDS)
endif(KAHYPAR_USE_64_BIT_IDS)

if(KAHYPAR_ENABLE_EXPERIMENTAL_FEATURES)
  add_compile_definitions(KAHYPAR_ENABLE_EXPERIMENTAL_FEATURES)
endif(KAHYPAR_ENABLE_EXPERIMENTAL_FEATURES)

if(KAHYPAR_TRAVIS_BUILD)
  add_compile_definitions(KAHYPAR_TRAVIS_BUILD)
endif(KAHYPAR_TRAVIS_BUILD)

if(KAHYPAR_USE_MINIMAL_BOOST)
  execute_process(COMMAND cmake -P ${CMAKE_CURRENT_SOURCE_DIR}/scripts/download_boost.cmake)

  include_directories(SYSTEM ${CMAKE_CURRENT_BINARY_DIR}/external_tools/boost/)

  # glob boost sources
  file(GLOB MINI_BOOST_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/external_tools/boost/libs/program_options/src/*.cpp)

  add_library(mini_boost STATIC ${MINI_BOOST_SOURCES})
  set_target_properties(mini_boost PROPERTIES LINKER_LANGUAGE CXX)
  set(Boost_LIBRARIES mini_boost)
  set(Boost_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/external_tools/boost/boost/)
  include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
else()
  find_package(Boost 1.69 REQUIRED COMPONENTS program_options)
  if(Boost_FOUND)
    include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
    set(KAHYPAR_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} ${KAHYPAR_INCLUDE_DIRS})
    set(KAHYPAR_LINK_LIBRARIES ${Boost_LIBRARIES} ${KAHYPAR_LINK_LIBRARIES})
    message(STATUS "Boost Include: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost Library Dirs: ${Boost_LIBRARY_DIRS}")
    message(STATUS "Boost Libraries: ${Boost_LIBRARIES}")
    if(WIN32)
      add_definitions(-DBOOST_ALL_NO_LIB)
      add_definitions(-DBOOST_PROGRAM_OPTIONS_DYN_LINK=1)
    endif()

  endif()
endif()

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
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgcov")

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

if(NOT CMAKE_BUILD_TYPE)
  set( CMAKE_BUILD_TYPE Debug CACHE STRING
       "Choose the type of build, options are: Debug Release, RelWithDebInfo"
       FORCE )
endif()

#find TBB
find_package(TBB REQUIRED)
include_directories(SYSTEM ${TBB_INCLUDE_DIRS})
link_libraries(${TBB_LIBRARIES})
MESSAGE(STATUS "Found TBB library: inc=${TBB_INCLUDE_DIRS}, lib=${TBB_LIBRARIES}")

# find LIBNUMA
FIND_PATH(NUMA_INCLUDE_DIR NAME numa.h
  HINTS $ENV{HOME}/local/include /opt/local/include /usr/local/include /usr/include)

FIND_LIBRARY(NUMA_LIBRARY NAME numa
  HINTS $ENV{HOME}/local/lib64 $ENV{HOME}/local/lib /usr/local/lib64 /usr/local/lib /opt/local/lib64 /opt/local/lib /usr/lib64 /usr/lib
)

IF (NUMA_INCLUDE_DIR AND NUMA_LIBRARY)
    SET(NUMA_FOUND TRUE)
    MESSAGE(STATUS "Found numa library: inc=${NUMA_INCLUDE_DIR}, lib=${NUMA_LIBRARY}")
    include_directories(${NUMA_INCLUDE_DIR})
    link_libraries(${NUMA_LIBRARY})
ELSE ()
    SET(NUMA_FOUND FALSE)
    MESSAGE(STATUS "WARNING: Numa library not found.")
    MESSAGE(STATUS "Try: 'sudo yum install numactl numactl-devel' (or sudo apt-get install libnuma libnuma-dev)")
ENDIF ()

# find HWLOC

FIND_PATH(HWLOCK_INCLUDE_DIR NAME hwloc.h
  HINTS $ENV{HOME}/local/include /opt/local/include /usr/local/include /usr/include)

FIND_LIBRARY(HWLOCK_LIBRARY NAME hwloc
  HINTS $ENV{HOME}/local/lib64 $ENV{HOME}/local/lib /usr/local/lib64 /usr/local/lib /opt/local/lib64 /opt/local/lib /usr/lib64 /usr/lib
)

IF (HWLOCK_INCLUDE_DIR AND HWLOCK_LIBRARY)
    SET(HWLOCK_FOUND TRUE)
    MESSAGE(STATUS "Found hwlock library: inc=${HWLOCK_INCLUDE_DIR}, lib=${HWLOCK_LIBRARY}")
    include_directories(${HWLOCK_INCLUDE_DIR})
    link_libraries(${HWLOCK_LIBRARY})
ELSE ()
    SET(HWLOCK_FOUND FALSE)
    MESSAGE(STATUS "WARNING: hwlock library not found.")
ENDIF ()

find_package(Git)
include(GetGitRevisionDescription)
get_git_head_revision(MT_KAHYPAR_VERSION_GIT_REFSPEC MT_KAHYPAR_VERSION_GIT_SHA1)
if(MT_KAHYPAR_VERSION_GIT_REFSPEC)
  message(STATUS "Detected git refspec ${MT_KAHYPAR_VERSION_GIT_REFSPEC} sha ${MT_KAHYPAR_VERSION_GIT_SHA1}")
  configure_file(${PROJECT_SOURCE_DIR}/mt-kahypar/application/git_revision.txt.in ${PROJECT_BINARY_DIR}/mt-kahypar/application/git_head.txt)
endif(MT_KAHYPAR_VERSION_GIT_REFSPEC)

if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra ")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused ")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wuninitialized")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfatal-errors")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wcast-qual")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Woverloaded-virtual")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wredundant-decls")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Winit-self")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPARANOID ")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ltbb")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ltbbmalloc_proxy")

  # lld is way faster than ld. If you have it, use it!
  find_program(LLD_BIN lld)
  if (LLD_BIN)
    message(STATUS "Found and will use LLVM linker " ${LLD_BIN})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fuse-ld=lld")
  else()
    message(STATUS "Default linker")
  endif()

  include(CheckSSE4_2)
  if( BUILTIN_POPCNT )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2")
  endif()

  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-mcrc32 KAHYPAR_HAS_CRC32)
  if(KAHYPAR_HAS_CRC32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcrc32")
  endif()


  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=undefined -fno-omit-frame-pointer")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -mtune=native -march=native")
  set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -mtune=native -march=native -g3 ")

  if(NOT KAHYPAR_CI_BUILD)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
    set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=address")
  endif()

  if(ENABLE_PROFILE MATCHES ON)
    message(STATUS "Profiling activated")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DENABLE_PROFILE")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g3 -DENABLE_PROFILE -fno-omit-frame-pointer")
    set(PROFILE_FLAGS "-lprofiler")
  endif()
elseif(MSVC)
  message(STATUS "MSVC not supported")
endif()

add_custom_target(AnalyzeModifiedSources  perl "${PROJECT_SOURCE_DIR}/codestyle/analyze-source.pl" "-w")
add_custom_target(AnalyzeAllSources  perl "${PROJECT_SOURCE_DIR}/codestyle/analyze-source.pl" "-aw")

message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")

include(gmock)
enable_testing()

add_subdirectory(mt-kahypar/application)
add_subdirectory(tools)
add_subdirectory(lib)
add_subdirectory(tests)

# This adds the source files. It's important that this happens after the compile targets have been added
add_subdirectory(mt-kahypar)

add_custom_target(MtKaHyPar
  COMMAND $(MAKE) MtKaHyParDefault
  COMMAND $(MAKE) MtKaHyParQuality
  COMMAND $(MAKE) MtKaHyParGraph
  COMMAND $(MAKE) MtKaHyParGraphQuality
  COMMAND $(MAKE) MtKaHyParWrapper)

add_custom_target(mt_kahypar_tests
  COMMAND $(MAKE) mt_kahypar_fast_tests
  COMMAND $(MAKE) mt_kahypar_strong_tests
  COMMAND $(MAKE) mt_kahypar_graph_tests)

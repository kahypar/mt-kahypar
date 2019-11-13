cmake_minimum_required(VERSION 3.10)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

project(KaHyPar CXX C)
set(PROJECT_VENDOR "Tobias Heuer")
set(PROJECT_CONTACT "tobias.heuer@kit.edu")
set(PROJECT_URL "http://www.kahypar.org")
set(PROJECT_DESCRIPTION "Mt-KaHyPar: Multi-Threaded Karlsruhe Hypergraph Partitioning")

include_directories(${PROJECT_SOURCE_DIR})

find_package(Threads REQUIRED)
message(STATUS "Found Threads: ${CMAKE_THREAD_LIBS_INIT}")

#ensure that gmock is built before tests are executed
add_subdirectory(external_tools/googletest EXCLUDE_FROM_ALL)
include_directories(SYSTEM ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
include_directories(SYSTEM ${gtest_SOURCE_DIR}/../googlemock/include ${gtest_SOURCE_DIR}/../googlemock/)

#include KaHyPar
include_directories(external_tools/kahypar)

option(KAHYPAR_USE_MINIMAL_BOOST
  "Download boost automatically and compile required libraries." OFF)

option(KAHYPAR_DISABLE_ASSERTIONS
  "Disable KaHyPar's internal assertions." OFF)

option(KAHYPAR_USE_STANDARD_ASSERTIONS
  "Use standard C++ asserts instead of custom assertions." OFF)

if(KAHYPAR_DISABLE_ASSERTIONS)
  add_compile_definitions(KAHYPAR_DISABLE_ASSERTIONS)
endif(KAHYPAR_DISABLE_ASSERTIONS)

if(KAHYPAR_USE_STANDARD_ASSERTIONS)
  add_compile_definitions(KAHYPAR_USE_STANDARD_ASSERTIONS)
endif(KAHYPAR_USE_STANDARD_ASSERTIONS)

if(NOT KAHPYAR_HEAVY_ASSERTION_TYPE)
  set( KAHPYAR_HEAVY_ASSERTION_TYPE "NONE" CACHE STRING
       "Choose one of the following heavy assertion types: NONE, PREPROCESSING,
        COARSENING, INITIAL_PARTITIONING, REFINEMENT"
       FORCE )
endif()

string( TOUPPER "${KAHPYAR_HEAVY_ASSERTION_TYPE}" UPPER_CASE_HEAVY_ASSERTION_TYPE )
if ( UPPER_CASE_HEAVY_ASSERTION_TYPE STREQUAL "PREPROCESSING" )
  add_compile_definitions(KAHYPAR_USE_HEAVY_PREPROCESSING_ASSERTIONS)
elseif( UPPER_CASE_HEAVY_ASSERTION_TYPE STREQUAL "COARSENING" )
  add_compile_definitions(KAHYPAR_USE_HEAVY_COARSENING_ASSERTIONS)
elseif( UPPER_CASE_HEAVY_ASSERTION_TYPE STREQUAL "INITIAL_PARTITIONING" )
  add_compile_definitions(KAHYPAR_USE_HEAVY_INITIAL_PARTITIONING_ASSERTIONS)
elseif( UPPER_CASE_HEAVY_ASSERTION_TYPE STREQUAL "REFINEMENT" )
  add_compile_definitions(KAHYPAR_USE_HEAVY_REFINEMENT_ASSERTIONS)
endif()

if(KAHYPAR_USE_MINIMAL_BOOST)
  include(ExternalProject)
  ExternalProject_Add(boost
    URL
    http://downloads.sourceforge.net/project/boost/boost/1.56.0/boost_1_56_0.tar.bz2
    #URL_MD5 6095876341956f65f9d35939ccea1a9f
    DOWNLOAD_DIR /tmp/
    PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/external_tools/boost
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

  include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/external_tools/boost/src/boost/)

  # glob boost sources
  file(GLOB MINI_BOOST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/external_tools/boost/src/boost/libs/program_options/src/*.cpp)

  add_library(mini_boost STATIC ${MINI_BOOST_SOURCES})
  set_target_properties(mini_boost PROPERTIES LINKER_LANGUAGE CXX)
  set(Boost_LIBRARIES mini_boost)
else()
  set(BOOST_MIN_VERSION "1.48.0")
  find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS program_options REQUIRED)
  if(Boost_FOUND)
    include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
    message(STATUS "Boost Include: ${Boost_INCLUDE_DIR}")
    message(STATUS "Boost Library Dirs: ${Boost_LIBRARY_DIRS}")
    message(STATUS "Boost Libraries: ${Boost_LIBRARIES}")
    if(WIN32)
      add_definitions(-DBOOST_ALL_NO_LIB)
      add_definitions(-DBOOST_PROGRAM_OPTIONS_DYN_LINK=1)
    endif()

  endif()
endif()

if(NOT CMAKE_BUILD_TYPE)
  set( CMAKE_BUILD_TYPE Debug CACHE STRING
       "Choose the type of build, options are: Debug Release, RelWithDebInfo"
       FORCE )
endif()

#find TBB
find_package(TBB REQUIRED)
include_directories(SYSTEM ${TBB_INCLUDE_DIRS})
link_libraries(${TBB_LIBRARIES})


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

# find libkahypar

FIND_PATH(KAHYPAR_INCLUDE_DIR NAME libkahypar.h
  HINTS "/usr/local/include" "${CMAKE_CURRENT_SOURCE_DIR}/external_tools/kahypar/include" )

FIND_LIBRARY(KAHYPAR_LIBRARY NAME libkahypar.so
  HINTS "/usr/local/lib" "${CMAKE_CURRENT_SOURCE_DIR}/external_tools/kahypar/build/lib" )

IF (KAHYPAR_INCLUDE_DIR AND KAHYPAR_LIBRARY)
    SET(KAHYPAR_FOUND TRUE)
    MESSAGE(STATUS "Found kahypar library: inc=${KAHYPAR_INCLUDE_DIR}, lib=${KAHYPAR_LIBRARY}")
    include_directories(${KAHYPAR_INCLUDE_DIR})
    link_libraries(${KAHYPAR_LIBRARY})
ELSE ()
    SET(KAHYPAR_FOUND FALSE)
    MESSAGE(STATUS "WARNING: kahypar library not found.")
ENDIF ()

# TODO: Add git revision

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
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ltbb")

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
  set(CMAKE_CXX_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -mtune=native -march=native")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -mtune=native -march=native -g3 ")

  if(ENABLE_PROFILE MATCHES ON)
    message(STATUS "Profiling activated")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DENABLE_PROFILE")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g3 -DENABLE_PROFILE -fno-omit-frame-pointer")
    set(PROFILE_FLAGS "-lprofiler")
  endif()
elseif(MSVC)
  message(STATUS "MSVC not supported")
endif()

message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")

include(gmock)
enable_testing()
add_subdirectory(mt-kahypar/application)
add_subdirectory(tools)
add_subdirectory(tests)
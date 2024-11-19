# custom script for finding the hwloc library, creates a HWLOC::hwloc target
set(HWLOC_TARGET_NAME HWLOC::hwloc)

if(NOT TARGET ${HWLOC_TARGET_NAME})
    if(NOT BUILD_SHARED_LIBS AND KAHYPAR_STATIC_LINK_DEPENDENCIES)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
    endif()

    add_library(hwloc INTERFACE)

    # search predefined paths
    find_path(HWLOCK_INCLUDE_DIR NAME hwloc.h)
    find_library(HWLOCK_LIBRARY NAME hwloc)
    if(HWLOCK_INCLUDE_DIR AND HWLOCK_LIBRARY)
        message(STATUS "Found hwlock library: inc=${HWLOCK_INCLUDE_DIR}, lib=${HWLOCK_LIBRARY}")
        add_library(hwloc INTERFACE IMPORTED)
        target_include_directories(hwloc INTERFACE ${HWLOCK_INCLUDE_DIR})
        target_link_libraries(hwloc INTERFACE ${HWLOCK_LIBRARY})
    elseif(UNIX)
        # search for hwloc via pkg_config
        find_package(PkgConfig QUIET)
        if (PKG_CONFIG_FOUND)
            pkg_search_module(HWLOC hwloc IMPORTED_TARGET)
            if (TARGET PkgConfig::HWLOC)
                target_link_libraries(hwloc INTERFACE PkgConfig::HWLOC)
            endif()
        endif()
    endif()

    if(NOT BUILD_SHARED_LIBS AND KAHYPAR_STATIC_LINK_DEPENDENCIES)
        # for static linking we also need udev as transitive dependency of hwloc
        find_package(PkgConfig QUIET REQUIRED)
        if (PKG_CONFIG_FOUND)
            pkg_search_module(UDEV udev IMPORTED_TARGET)
            if (TARGET PkgConfig::UDEV)
                # -ludev is an annoying workaround since propagating the static dependencies is not properly supported by cmake
                target_link_libraries(hwloc INTERFACE PkgConfig::UDEV -ludev)
            else()
                message(FATAL_ERROR "libudev is required for statically linking hwloc. Install libudev or use -DKAHYPAR_STATIC_LINK_DEPENDENCIES=Off")
            endif()
        endif()
    endif()

    # actual target for outside use
    add_library(${HWLOC_TARGET_NAME} ALIAS hwloc)
endif()

# custom script for finding the hwloc library, creates a HWLOC::hwloc target
set(HWLOC_TARGET_NAME HWLOC::hwloc)

if(NOT TARGET ${HWLOC_TARGET_NAME})
    # search predefined paths
    find_path(HWLOCK_INCLUDE_DIR NAME hwloc.h)
    find_library(HWLOCK_LIBRARY NAME hwloc)
    if(HWLOCK_INCLUDE_DIR AND HWLOCK_LIBRARY)
        message(STATUS "Found hwlock library: inc=${HWLOCK_INCLUDE_DIR}, lib=${HWLOCK_LIBRARY}")
        add_library(hwloc INTERFACE IMPORTED)
        target_include_directories(hwloc INTERFACE ${HWLOCK_INCLUDE_DIR})
        target_link_libraries(hwloc INTERFACE ${HWLOCK_LIBRARY})
        add_library(${HWLOC_TARGET_NAME} ALIAS hwloc)
    elseif(UNIX)
        # search for hwloc via pkg_config
        find_package(PkgConfig QUIET)
        if (PKG_CONFIG_FOUND)
            pkg_search_module(HWLOC hwloc IMPORTED_TARGET)
            if (TARGET PkgConfig::HWLOC)
                add_library(${HWLOC_TARGET_NAME} ALIAS PkgConfig::HWLOC)
            endif()
        endif()
    endif()
endif()

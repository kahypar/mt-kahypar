# Custom script for finding the hwloc library, creates a HWLOC::hwloc target.
# Some complications are necessary in order to also support static linking.
set(HWLOC_TARGET_NAME HWLOC::hwloc)
if(NOT BUILD_SHARED_LIBS AND KAHYPAR_STATIC_LINK_DEPENDENCIES)
    set(HWLOC_LINK_STATICALLY TRUE)
endif()

function(find_hwloc_path)
    # search predefined paths
    find_path(HWLOC_INCLUDE_DIR NAMES hwloc.h)
    find_library(HWLOC_LIBRARY NAMES hwloc)
    if(HWLOC_INCLUDE_DIR AND HWLOC_LIBRARY)
        message(STATUS "Found hwloc library: inc=${HWLOC_INCLUDE_DIR}, lib=${HWLOC_LIBRARY}")
        target_include_directories(hwloc INTERFACE ${HWLOC_INCLUDE_DIR})
        target_link_libraries(hwloc INTERFACE ${HWLOC_LIBRARY})
        set(HWLOC_FOUND TRUE PARENT_SCOPE)
    endif()
endfunction()

if(NOT TARGET ${HWLOC_TARGET_NAME})
    if(HWLOC_LINK_STATICALLY)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
    endif()
    add_library(hwloc INTERFACE IMPORTED)

    find_hwloc_path()
    if(NOT HWLOC_FOUND AND HWLOC_LINK_STATICALLY)
        message(WARNING "Could not find static version of hwloc library. Proceeding to link dynamically.")
        set(HWLOC_LINK_STATICALLY FALSE)
        unset(CMAKE_FIND_LIBRARY_SUFFIXES)
        find_hwloc_path()
    endif()
    if(NOT HWLOC_FOUND)
        # search for hwloc via pkg_config
        find_package(PkgConfig QUIET)
        if (PKG_CONFIG_FOUND)
            pkg_search_module(HWLOC hwloc IMPORTED_TARGET)
            if (TARGET PkgConfig::HWLOC)
                target_link_libraries(hwloc INTERFACE PkgConfig::HWLOC)
                set(HWLOC_FOUND TRUE)
            endif()
        endif()
    endif()
    if(NOT HWLOC_FOUND)
        message(FATAL_ERROR "hwloc library not found. Install hwloc on your system.")
    endif()

    if(HWLOC_LINK_STATICALLY)
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

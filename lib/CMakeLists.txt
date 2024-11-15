include(GNUInstallDirs)

if(BUILD_SHARED_LIBS)
  add_library(mtkahypar SHARED libmtkahypar.cpp)
else()
  add_library(mtkahypar STATIC libmtkahypar.cpp)
endif()
target_link_libraries(mtkahypar PRIVATE MtKaHyPar-LibraryBuildSources)
target_include_directories(mtkahypar INTERFACE "../include")
set_target_properties(mtkahypar PROPERTIES
    PUBLIC_HEADER "../include/libmtkahypar.h;../include/libmtkahypartypes.h")

configure_file(libmtkahypar.pc.in libmtkahypar.pc @ONLY)

install(TARGETS mtkahypar
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${CMAKE_BINARY_DIR}/lib/libmtkahypar.pc
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

configure_file(cmake_uninstall.cmake.in cmake_uninstall.cmake IMMEDIATE @ONLY)

add_custom_target(uninstall-mtkahypar "${CMAKE_COMMAND}" -P cmake_uninstall.cmake)

add_custom_target(install.mtkahypar
    ${CMAKE_COMMAND}
    -DBUILD_TYPE=${CMAKE_BUILD_TYPE}
    -P ${CMAKE_BINARY_DIR}/cmake_install.cmake)
ADD_DEPENDENCIES(install.mtkahypar mtkahypar)

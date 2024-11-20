# Script that creates the installation targets

# this must be set before including GNUInstallDirs
set(CMAKE_INSTALL_PREFIX "/usr/")

include(GNUInstallDirs)

set_target_properties(mtkahypar PROPERTIES
    VERSION ${MT_KAHYPAR_VERSION} SOVERSION ${MT_KAHYPAR_SO_VERSION})
set_target_properties(mtkahypar PROPERTIES
    PUBLIC_HEADER "${MTKAHYPAR_INCLUDE_DIR}/mtkahypar.h;${MTKAHYPAR_INCLUDE_DIR}/mtkahypartypes.h")

# TODO: rename to mtkahypar!
configure_file(cmake/mtkahypar.pc.in mtkahypar.pc @ONLY)
configure_file(cmake/MtKaHyParConfig.cmake.in MtKaHyParConfig.cmake @ONLY)
configure_file(cmake/MtKaHyParConfigVersion.cmake.in MtKaHyParConfigVersion.cmake @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/mtkahypar.pc
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

install(FILES ${CMAKE_BINARY_DIR}/MtKaHyParConfig.cmake
              ${CMAKE_BINARY_DIR}/MtKaHyParConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MtKaHyPar)

install(TARGETS mtkahypar
  EXPORT MtKaHyPar
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT MtKaHyPar
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MtKaHyPar
  FILE MtKaHyPar.cmake
  NAMESPACE MtKaHyPar::)

configure_file(cmake/cmake_uninstall.cmake.in cmake_uninstall.cmake IMMEDIATE @ONLY)
add_custom_target(uninstall-mtkahypar "${CMAKE_COMMAND}" -P cmake_uninstall.cmake)

add_custom_target(install-mtkahypar
  ${CMAKE_COMMAND}
  -DBUILD_TYPE=${CMAKE_BUILD_TYPE}
  -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
  DEPENDS mtkahypar)

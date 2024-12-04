# Script that creates the installation targets for a shared library

include(GNUInstallDirs)

set_target_properties(mtkahypar PROPERTIES
    VERSION ${MT_KAHYPAR_VERSION} SOVERSION ${MT_KAHYPAR_SO_VERSION})
set_target_properties(mtkahypar PROPERTIES
    PUBLIC_HEADER "${MTKAHYPAR_INCLUDE_DIR}/mtkahypar.h;${MTKAHYPAR_INCLUDE_DIR}/mtkahypartypes.h")

configure_file(cmake/mtkahypar.pc.in mtkahypar.pc @ONLY)
configure_file(cmake/MtKaHyParConfig.cmake.in MtKaHyParConfig.cmake @ONLY)
configure_file(cmake/MtKaHyParConfigVersion.cmake.in MtKaHyParConfigVersion.cmake @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/mtkahypar.pc
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
  COMPONENT MtKaHyPar_Lib)

install(FILES ${CMAKE_BINARY_DIR}/MtKaHyParConfig.cmake
              ${CMAKE_BINARY_DIR}/MtKaHyParConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MtKaHyPar
  COMPONENT MtKaHyPar_Lib)

install(TARGETS mtkahypar
  EXPORT MtKaHyPar
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  COMPONENT MtKaHyPar_Lib
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT MtKaHyPar_Lib)

install(EXPORT MtKaHyPar
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MtKaHyPar
  FILE MtKaHyPar.cmake
  NAMESPACE MtKaHyPar::
  COMPONENT MtKaHyPar_Lib)

# custom targets for installing/uninstalling the library
add_custom_target(install-mtkahypar
  ${CMAKE_COMMAND}
  -DBUILD_TYPE=${CMAKE_BUILD_TYPE}
  -DCMAKE_INSTALL_COMPONENT=MtKaHyPar_Lib
  -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
  DEPENDS mtkahypar)

configure_file(cmake/cmake_uninstall.cmake.in cmake_uninstall.cmake IMMEDIATE @ONLY)
add_custom_target(uninstall-mtkahypar
  "${CMAKE_COMMAND}"
  -DMANIFEST_NAME=install_manifest_MtKaHyPar_Lib.txt
  -P ${CMAKE_BINARY_DIR}/cmake_uninstall.cmake)

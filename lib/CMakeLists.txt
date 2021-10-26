include(GNUInstallDirs)

# Library MT-KaHyPar Fast
add_library(mtkahyparfast SHARED libmtkahypar.cc)
target_link_libraries(mtkahyparfast ${Boost_LIBRARIES})
target_compile_definitions(mtkahyparfast PUBLIC MT_KAHYPAR_LIBRARY_MODE)

# Library MT-KaHyPar Strong
add_library(mtkahyparstrong SHARED libmtkahypar.cc)
target_link_libraries(mtkahyparstrong ${Boost_LIBRARIES})
target_compile_definitions(mtkahyparstrong PUBLIC MT_KAHYPAR_LIBRARY_MODE)
target_compile_definitions(mtkahyparstrong PUBLIC USE_STRONG_PARTITIONER)

# Library MT-KaHyPar Fast Graph
add_library(mtkahypargraph SHARED libmtkahypar.cc)
target_link_libraries(mtkahypargraph ${Boost_LIBRARIES})
target_compile_definitions(mtkahypargraph PUBLIC MT_KAHYPAR_LIBRARY_MODE)
target_compile_definitions(mtkahypargraph PUBLIC USE_GRAPH_PARTITIONER)

set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} mtkahyparfast mtkahyparstrong mtkahypargraph PARENT_SCOPE)

set_target_properties(mtkahyparfast PROPERTIES
    PUBLIC_HEADER ../include/libmtkahypar.h)
set_target_properties(mtkahyparstrong PROPERTIES
    PUBLIC_HEADER ../include/libmtkahypar.h)
set_target_properties(mtkahypargraph PROPERTIES
    PUBLIC_HEADER ../include/libmtkahypar.h)

target_include_directories(mtkahyparfast PRIVATE ../include)
target_include_directories(mtkahyparstrong PRIVATE ../include)
target_include_directories(mtkahypargraph PRIVATE ../include)

configure_file(libmtkahyparfast.pc.in libmtkahyparfast.pc @ONLY)
configure_file(libmtkahyparstrong.pc.in libmtkahyparstrong.pc @ONLY)
configure_file(libmtkahypargraph.pc.in libmtkahypargraph.pc @ONLY)

install(TARGETS mtkahyparfast mtkahyparstrong mtkahypargraph
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${CMAKE_BINARY_DIR}/lib/libmtkahyparfast.pc
              ${CMAKE_BINARY_DIR}/lib/libmtkahyparstrong.pc
              ${CMAKE_BINARY_DIR}/lib/libmtkahypargraph.pc
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

configure_file(cmake_uninstall.cmake.in cmake_uninstall.cmake IMMEDIATE @ONLY)

add_custom_target(uninstall-mtkahypar "${CMAKE_COMMAND}" -P cmake_uninstall.cmake)

add_custom_target(install.mtkahypar
    ${CMAKE_COMMAND}
    -DBUILD_TYPE=${CMAKE_BUILD_TYPE}
    -P ${CMAKE_BINARY_DIR}/cmake_install.cmake)
ADD_DEPENDENCIES(install.mtkahypar mtkahyparfast mtkahyparstrong mtkahypargraph)
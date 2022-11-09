include(GNUInstallDirs)

# Library for hypergraph partitioning
add_library(mtkahyparhgp SHARED libmtkahyparhgp.cpp)
target_link_libraries(mtkahyparhgp ${Boost_LIBRARIES})
target_compile_definitions(mtkahyparhgp PUBLIC MT_KAHYPAR_LIBRARY_MODE)
SET_TARGET_PROPERTIES(mtkahyparhgp PROPERTIES COMPILE_FLAGS "-fvisibility=hidden")

# Library for graph partitioning
add_library(mtkahypargp SHARED libmtkahypargp.cpp)
target_link_libraries(mtkahypargp ${Boost_LIBRARIES})
target_compile_definitions(mtkahypargp PUBLIC MT_KAHYPAR_LIBRARY_MODE)
target_compile_definitions(mtkahypargp PUBLIC USE_GRAPH_PARTITIONER)
SET_TARGET_PROPERTIES(mtkahypargp PROPERTIES COMPILE_FLAGS "-fvisibility=hidden")

set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} mtkahyparhgp mtkahypargp PARENT_SCOPE)

# Library for Mt-KaHyPar
add_library(mtkahypar SHARED libmtkahypar.cc)
target_link_libraries(mtkahypar ${Boost_LIBRARIES} mtkahyparhgp mtkahypargp)
target_compile_definitions(mtkahypar PUBLIC MT_KAHYPAR_LIBRARY_MODE)

set(LibMtKaHyParSources
      ../mt-kahypar/partition/context.cpp
      ../mt-kahypar/partition/context_enum_classes.cpp
      ../mt-kahypar/io/command_line_options.cpp)
target_sources(mtkahypar PRIVATE ${LibMtKaHyParSources})

set_target_properties(mtkahyparhgp PROPERTIES
    PUBLIC_HEADER "../include/libmtkahyparhgp.h")
set_target_properties(mtkahypargp PROPERTIES
    PUBLIC_HEADER "../include/libmtkahypargp.h")
set_target_properties(mtkahypar PROPERTIES
    PUBLIC_HEADER "../include/libmtkahypar.h;../include/libmtkahypartypes.h")

target_include_directories(mtkahyparhgp SYSTEM PUBLIC ../include)
target_include_directories(mtkahypargp SYSTEM PUBLIC ../include)
target_include_directories(mtkahypar SYSTEM PUBLIC ../include)

configure_file(libmtkahypar.pc.in libmtkahypar.pc @ONLY)
configure_file(libmtkahyparhgp.pc.in libmtkahyparhgp.pc @ONLY)
configure_file(libmtkahypargp.pc.in libmtkahypargp.pc @ONLY)

install(TARGETS mtkahypar mtkahyparhgp mtkahypargp
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${CMAKE_BINARY_DIR}/lib/libmtkahypar.pc
              ${CMAKE_BINARY_DIR}/lib/libmtkahyparhgp.pc
              ${CMAKE_BINARY_DIR}/lib/libmtkahypargp.pc
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

configure_file(cmake_uninstall.cmake.in cmake_uninstall.cmake IMMEDIATE @ONLY)

add_custom_target(uninstall-mtkahypar "${CMAKE_COMMAND}" -P cmake_uninstall.cmake)

add_custom_target(install.mtkahypar
    ${CMAKE_COMMAND}
    -DBUILD_TYPE=${CMAKE_BUILD_TYPE}
    -P ${CMAKE_BINARY_DIR}/cmake_install.cmake)
ADD_DEPENDENCIES(install.mtkahypar mtkahypar)

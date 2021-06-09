include(GNUInstallDirs)

add_library(mtkahypar SHARED libmtkahypar.cc)
target_link_libraries(mtkahypar ${Boost_LIBRARIES})
target_compile_definitions(mtkahypar PUBLIC MT_KAHYPAR_LIBRARY_MODE)

set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} mtkahypar PARENT_SCOPE)

set_target_properties(mtkahypar PROPERTIES
    PUBLIC_HEADER ../include/libmtkahypar.h)

target_include_directories(mtkahypar PRIVATE ../include)

configure_file(libmtkahypar.pc.in libmtkahypar.pc @ONLY)

if(WIN32)
	install(TARGETS mtkahypar
		RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
		PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
else()
	install(TARGETS mtkahypar
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
		PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

install(FILES ${CMAKE_BINARY_DIR}/lib/libmtkahypar.pc
       DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

configure_file(cmake_uninstall.cmake.in cmake_uninstall.cmake IMMEDIATE @ONLY)

add_custom_target(uninstall-mtkahypar "${CMAKE_COMMAND}" -P cmake_uninstall.cmake)

add_custom_target(install.mtkahypar
    ${CMAKE_COMMAND}
    -DBUILD_TYPE=${CMAKE_BUILD_TYPE}
    -P ${CMAKE_BINARY_DIR}/cmake_install.cmake)
ADD_DEPENDENCIES(install.mtkahypar mtkahypar)
if(BUILD_SHARED_LIBS)
  add_library(mtkahypar SHARED mtkahypar.cpp)
else()
  add_library(mtkahypar STATIC mtkahypar.cpp)
endif()
target_link_libraries(mtkahypar PRIVATE MtKaHyPar-LibraryBuildSources)
target_include_directories(mtkahypar INTERFACE
  $<BUILD_INTERFACE:${MTKAHYPAR_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:include>
)

if(BUILD_SHARED_LIBS)
  target_compile_definitions(mtkahypar PRIVATE MT_KAHYPAR_BUILD_SHARED)
else()
  target_compile_definitions(mtkahypar PUBLIC MT_KAHYPAR_BUILD_STATIC)
endif()

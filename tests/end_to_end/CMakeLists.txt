add_executable(integration_tests mt_kahypar_integration_tests.cc)
target_link_libraries(integration_tests gtest gtest_main ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
set_property(TARGET integration_tests PROPERTY CXX_STANDARD 17)
set_property(TARGET integration_tests PROPERTY CXX_STANDARD_REQUIRED ON)

# exclude integration tests from default build target
set_target_properties(integration_tests PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)


set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} integration_tests PARENT_SCOPE)


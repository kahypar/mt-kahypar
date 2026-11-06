# taken from http://johnlamp.net/cmake-tutorial-5-functionally-improved-testing.html
#
# Note: the tests are not executed as part of the build. Run them via ctest
# (e.g., `ctest --test-dir <build-dir> -R <target> --output-on-failure`) or by
# executing the binary directly.
function(add_gmock_test target)
    add_mtkahypar_executable(${target} ${ARGN})

    target_link_libraries(${target} gtest gtest_main ${CMAKE_THREAD_LIBS_INIT})

    set_property(TARGET ${target} PROPERTY CXX_STANDARD 20)
    set_property(TARGET ${target} PROPERTY CXX_STANDARD_REQUIRED ON)

    add_test(NAME ${target} COMMAND ${target})
endfunction()

find_package(Python3 COMPONENTS Interpreter REQUIRED)

add_test(
    NAME mtkahypar_python_test
    COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/test_mtkahypar.py"
)

set_tests_properties(mtkahypar_python_test PROPERTIES
    ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:mtkahypar_python>;$ENV{PYTHONPATH}"
)

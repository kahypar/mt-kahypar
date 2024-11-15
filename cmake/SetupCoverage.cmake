#################################################################
## Creates a code coverage target with name test-coverage      ##
#################################################################

# find programs
find_program(GENHTML genhtml)
find_program(LCOV lcov)
if(NOT LCOV OR NOT GENHTML)
  message(SEND_ERROR "Coverage analysis requires lcov and genhtml.")
endif()
if(NOT KAHYPAR_ENABLE_TESTING)
  message(SEND_ERROR "Coverage analysis requires that tests are enabled.")
endif()

# add coverage anaylsis compile and link flags
target_compile_options(MtKaHyPar-BuildFlags INTERFACE -fprofile-arcs -ftest-coverage -fprofile-update=atomic)
target_link_options(MtKaHyPar-BuildFlags INTERFACE -lgcov --coverage)

# add cached variable containing parameters for lcov/genhtml
set(GENHTML_FLAGS --legend --no-branch-coverage
  CACHE STRING "parameters for genhtml")

# custom target to run before tests
add_custom_target(lcov-reset
  COMMAND ${LCOV} -q --directory ${CMAKE_BINARY_DIR} --zerocounters
  COMMENT "Resetting code coverage counters")

# custom lcov target to run tests
add_custom_target(lcov-runtests
  # TODO: better handle interface_test and include unit tests
  COMMAND make mtkahypar_interface_test
  DEPENDS lcov-reset MtKaHyPar-LibraryBuildSources MtKaHyPar-Test mtkahypar
  COMMENT "Running all unit tests")

# get git version description
execute_process(COMMAND git describe --tags
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE GITDESC
  OUTPUT_STRIP_TRAILING_WHITESPACE)

# command sequence to gather, clean and generate HTML coverage report
set(COVERAGE_OUTPUT_DIR "coverage")
set(LCOV_EXCLUSIONS "*/external_tools/*" "*/tests/*" "/usr/*" "*/_deps/*")
set(LCOV_COMMAND ${LCOV} -q --directory . --capture --output-file lcov.info)
foreach(EXCLUSION IN LISTS LCOV_EXCLUSIONS)
	list(APPEND LCOV_COMMAND --exclude "${EXCLUSION}")
endforeach()
set(GENHTML_COMMAND ${GENHTML} -q -o ${COVERAGE_OUTPUT_DIR} --title "Mt-KaHyPar ${GITDESC}"
					--prefix ${PROJECT_SOURCE_DIR} ${GENHTML_FLAGS} lcov.info)

add_custom_target(lcov-html
  COMMAND echo ${LCOV_COMMAND} COMMAND ${LCOV_COMMAND}
  COMMAND echo ${GENHTML_COMMAND} COMMAND ${GENHTML_COMMAND}
  DEPENDS lcov-runtests
  COMMENT "Capturing code coverage counters and create HTML coverage report"
  VERBATIM
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# top-level target to run tests and generate coverage report
add_custom_target(test-coverage
  COMMENT "Generate HTML coverage report"
  DEPENDS lcov-html)

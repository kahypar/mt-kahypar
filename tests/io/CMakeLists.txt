target_sources(mt_kahypar_tests PRIVATE
        hypergraph_io_test.cc
        sql_plottools_serializer_test.cc
        )

if(NOT MT_KAHYPAR_DISABLE_BOOST)
        configure_file(context_test.cc.in ${PROJECT_BINARY_DIR}/tests/io/context_test.cc)
        target_sources(mt_kahypar_tests PRIVATE ${PROJECT_BINARY_DIR}/tests/io/context_test.cc)
endif()

target_sources(mtkahypar_tests PRIVATE
        hypergraph_io_test.cc
        sql_plottools_serializer_test.cc
        )

configure_file(context_test.cc.in ${PROJECT_BINARY_DIR}/tests/io/context_test.cc)
target_sources(mtkahypar_tests PRIVATE ${PROJECT_BINARY_DIR}/tests/io/context_test.cc)

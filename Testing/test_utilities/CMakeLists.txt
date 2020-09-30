target_sources(xstell_unit_test_runner
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/unit_test_runner.f>
)

target_sources(xstell_file_opts_runner
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/file_opts_tester.f>
)

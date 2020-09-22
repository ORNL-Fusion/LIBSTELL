target_sources(xstell_unit_test_runner
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/unit_test_runner.f>
)

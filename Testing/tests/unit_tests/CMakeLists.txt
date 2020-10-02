#  Run the built in unit tests of LIBSTELL.
add_test (NAME    libstell_unit_tests 
          COMMAND $<TARGET_PROPERTY:xstell_unit_test_runner,BINARY_DIR>/xstell_unit_test_runner)

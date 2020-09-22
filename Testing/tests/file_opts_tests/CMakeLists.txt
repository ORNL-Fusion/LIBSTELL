#  Run the built in unit tests of LIBSTELL.
add_test (NAME    libstell_file_opts_tests
          COMMAND $<TARGET_PROPERTY:xstell_file_opts_runner,BINARY_DIR>/xstell_file_opts_runner)



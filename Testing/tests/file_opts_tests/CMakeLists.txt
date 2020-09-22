#  Copy input file from the source to the build directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_file
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test_file)

#  Run the built in unit tests of LIBSTELL.
add_test (NAME    libstell_file_opts_tests
          COMMAND $<TARGET_PROPERTY:xstell_file_opts_runner,BINARY_DIR>/xstell_file_opts_runner)

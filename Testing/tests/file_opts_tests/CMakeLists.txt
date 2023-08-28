#  Clean test directory.
add_test (NAME    libstell_file_opts_tests_clean_directory
          COMMAND rm -rf test_file test_file_copy test_directory)

#  Copy input file from the source to the build directory.
add_test (NAME    libstell_file_opts_tests_copy_test_file
          COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/test_file ${CMAKE_CURRENT_BINARY_DIR}/test_file)
set_tests_properties (libstell_file_opts_tests_copy_test_file
                      PROPERTIES
                      DEPENDS libstell_file_opts_tests_clean_directory)

#  Run the built in unit tests of LIBSTELL.
add_test (NAME    libstell_file_opts_tests
          COMMAND $<TARGET_PROPERTY:xstell_file_opts_runner,BINARY_DIR>/xstell_file_opts_runner)
set_tests_properties (libstell_file_opts_tests
                      PROPERTIES
                      DEPENDS libstell_file_opts_tests_copy_test_file)

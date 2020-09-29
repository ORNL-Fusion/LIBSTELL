cmake_minimum_required (VERSION 3.14)

project (libstell C Fortran)

add_library (stell STATIC)

add_subdirectory (Sources)

# FIXME: Should use standard preprocessor names to detect the current system.
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_compile_options (stell PUBLIC -DDARWIN)
endif ()

OPTION (USE_PROFILER "Activate code profiling." OFF)
if (USE_PROFILER)
    target_compile_options (stell PUBLIC -DPROFILE_ON)
endif ()

# Activate or Deactivate fast copies. The fast copies maybe activated by setting -DUSE_FAST_COPY=ON on the cmake command line.
OPTION (USE_FAST_COPY "Activate code fast copies." OFF)
if (USE_FAST_COPY)
    target_compile_options (stell PUBLIC -DFAST_COPY)
endif ()

#  Fortran specific settings. The first setting tells the compiler to use the C
#  preprocessor.
target_compile_options (stell PUBLIC -cpp)

################################################################################
#  Testing                                                                     #
################################################################################

#  Build test utilities.
add_executable (xstell_unit_test_runner)
add_dependencies (xstell_unit_test_runner stell)

target_link_libraries (xstell_unit_test_runner stell)
target_include_directories (xstell_unit_test_runner PUBLIC $<TARGET_PROPERTY:stell,BINARY_DIR>)

#  Build file_opts testers
add_executable (xstell_file_opts_runner)
add_dependencies (xstell_file_opts_runner stell)

target_link_libraries (xstell_file_opts_runner stell)
target_include_directories (xstell_file_opts_runner PUBLIC $<TARGET_PROPERTY:stell,BINARY_DIR>)

add_subdirectory (Testing)

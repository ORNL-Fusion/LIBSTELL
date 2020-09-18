cmake_minimum_required (VERSION 3.14)

project (libstell C Fortran)

#  Fortran specific settings. The first setting tells the compiler to use the C
#  preprocessor.
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")

OPTION (USE_PROFILER "Activate code profiling." OFF)
if (USE_PROFILER)
    add_definitions (-DPROFILE_ON)
endif ()

# Activate or Deactivate fast copies. The fast copies maybe activated by setting -DUSE_FAST_COPY=ON on the cmake command line.
OPTION (USE_FAST_COPY "Activate code fast copies." OFF)
if (USE_FAST_COPY)
    add_definitions (-DFAST_COPY)
endif ()

add_library (stell STATIC)

add_subdirectory (Sources)

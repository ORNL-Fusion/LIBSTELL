#[==[
Provides the following variables:

  * `SCALAPACK_FOUND`: Whether SCALAPACK was found or not.
  * `SCALAPACK_LIBRARIES`: Libraries necessary to use SCALAPACK.
  * `SCALAPACK::SCALAPACK`: A target to use with `target_link_libraries`.
#]==]

# Find it manually.
find_library (SCALAPACK_LIBRARY
              NAMES scalapack scalapack-openmpi mkl_scalapack_ilp64
              DOC "scalapack library")

if (${SCALAPACK_LIBRARY} MATCHES mkl_scalapack_ilp64)
    find_library (BLACS_LIBRARY mkl_blacs_intelmpi_lp64)

endif ()

mark_as_advanced (SCALAPACK_LIBRARY)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (SCALAPACK
                                   REQUIRED_VARS SCALAPACK_LIBRARY)

if (SCALAPACK_FOUND)
    set (SCALAPACK_LIBRARIES "${SCALAPACK_LIBRARY} ${BLACS_LIBRARY}")

    if (NOT TARGET SCALAPACK::SCALAPACK_BASE)
       add_library (SCALAPACK::SCALAPACK_BASE UNKNOWN IMPORTED)
       set_target_properties (SCALAPACK::SCALAPACK_BASE PROPERTIES
                              IMPORTED_LOCATION ${SCALAPACK_LIBRARY})
    endif ()

    if (SCABLACS_LIBRARY)
        if (NOT TARGET SCALAPACK::BLACS)
            add_library (SCALAPACK::BLACS UNKNOWN IMPORTED)
            set_target_properties (SCALAPACK::BLACS PROPERTIES
                                   IMPORTED_LOCATION ${BLACS_LIBRARY})
        endif ()
    endif ()

    if (NOT TARGET  SCALAPACK::SCALAPACK)
        add_library (SCALAPACK::SCALAPACK INTERFACE IMPORTED)
        target_link_libraries (SCALAPACK::SCALAPACK INTERFACE
                               SCALAPACK::SCALAPACK
                               $<$<BOOL:${BLACS_LIBRARY}>:SCALAPACK::BLACS>)
    endif ()
endif ()

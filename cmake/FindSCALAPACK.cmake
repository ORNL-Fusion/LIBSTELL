#[==[
Provides the following variables:

  * `SCALAPACK_FOUND`: Whether SCALAPACK was found or not.
  * `SCALAPACK_LIBRARIES`: Libraries necessary to use SCALAPACK.
  * `SCALAPACK::SCALAPACK`: A target to use with `target_link_libraries`.
#]==]

# Find it manually.
find_library (SCALAPACK_LIBRARY
              NAMES scalapack scalapack-openmpi mkl_scalapack_ilp64 ${SCALAPACK_NAME}
              DOC "scalapack library")

find_library (BLACS_LIBRARY
              NAMES blacs-openmpi mkl_blacs_intelmpi_lp64
              DOC "blacs library")

find_library (BLACS_INIT_LIBRARY
              NAMES blacsF77init-openmpi
              DOC "blacs init library")

mark_as_advanced (SCALAPACK_LIBRARY)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (SCALAPACK
                                   REQUIRED_VARS SCALAPACK_LIBRARY)

if (${SCALAPACK_FOUND})
    set (SCALAPACK_LIBRARIES "${SCALAPACK_LIBRARY} ${BLACS_LIBRARY} ${BLACS_INIT_LIBRARY}")

    if (NOT TARGET SCALAPACK::SCALAPACK_BASE)
       add_library (SCALAPACK::SCALAPACK_BASE UNKNOWN IMPORTED)
       set_target_properties (SCALAPACK::SCALAPACK_BASE PROPERTIES
                              IMPORTED_LOCATION ${SCALAPACK_LIBRARY})
    endif ()

    if (EXISTS ${BLACS_LIBRARY})
        if (NOT TARGET SCALAPACK::BLACS)
            add_library (SCALAPACK::BLACS UNKNOWN IMPORTED)
            set_target_properties (SCALAPACK::BLACS PROPERTIES
                                   IMPORTED_LOCATION ${BLACS_LIBRARY})
        endif ()
    endif ()

    if (EXISTS ${BLACS_INIT_LIBRARY})
        if (NOT TARGET SCALAPACK::BLACS_INIT)
            add_library (SCALAPACK::BLACS_INIT UNKNOWN IMPORTED)
            set_target_properties (SCALAPACK::BLACS_INIT PROPERTIES
                                   IMPORTED_LOCATION ${BLACS_INIT_LIBRARY})
        endif ()
    endif ()

    if (NOT TARGET  SCALAPACK::SCALAPACK)
        add_library (SCALAPACK::SCALAPACK INTERFACE IMPORTED)
        target_link_libraries (SCALAPACK::SCALAPACK INTERFACE
                               SCALAPACK::SCALAPACK_BASE
                               $<$<BOOL:${BLACS_LIBRARY}>:SCALAPACK::BLACS>
                               $<$<BOOL:${BLACS_INIT_LIBRARY}>:SCALAPACK::BLACS_INIT>)
    endif ()
endif ()


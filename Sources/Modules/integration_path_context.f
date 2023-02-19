!*******************************************************************************
!>  @file integration_path_context.f
!>  @brief Contains module @ref integration_path_context
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Module is part of the LIBSTELL. This modules contains code to define and
!>  integrate along an arbitray path.
!*******************************************************************************
      MODULE integration_path_context

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) integration_path_context
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Integration context object interface.
!-------------------------------------------------------------------------------
      TYPE, ABSTRACT :: integration_path_context_class
      CONTAINS
         PROCEDURE (integration_function), DEFERRED :: run
      END TYPE

!-------------------------------------------------------------------------------
!>  Search context object interface.
!-------------------------------------------------------------------------------
      TYPE, ABSTRACT :: search_path_context_class
      CONTAINS
         PROCEDURE (search_function), DEFERRED :: run
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the integration function.
!-------------------------------------------------------------------------------
      INTERFACE
         FUNCTION integration_function(context, xcart, dxcart,                 &
     &                                 length, dx)
         USE stel_kinds
         IMPORT
         REAL (rprec)                                       ::                 &
     &      integration_function
         CLASS (integration_path_context_class), INTENT(in) :: context
         REAL (rprec), DIMENSION(3), INTENT(in)             :: xcart
         REAL (rprec), DIMENSION(3), INTENT(in)             :: dxcart
         REAL (rprec), INTENT(in)                           :: length
         REAL (rprec), INTENT(in)                           :: dx
         END FUNCTION
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for the search function.
!-------------------------------------------------------------------------------
      INTERFACE
         FUNCTION search_function(context, xcart1, xcart2)
         USE stel_kinds
         IMPORT
         LOGICAL                                       ::                      &
     &      search_function
         CLASS (search_path_context_class), INTENT(in) :: context
         REAL (rprec), DIMENSION(3), INTENT(in)        :: xcart1
         REAL (rprec), DIMENSION(3), INTENT(in)        :: xcart2
         END FUNCTION
      END INTERFACE

      END MODULE

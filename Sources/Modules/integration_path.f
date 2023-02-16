!*******************************************************************************
!>  @file integration_path.f
!>  @brief Contains module @ref integration_path
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Module is part of the LIBSTELL. This modules contains code to define and
!>  integrate along an arbitray path.
!*******************************************************************************

      MODULE integration_path
      USE stel_kinds
      USE stel_constants, ONLY : pi
      USE mpi_inc
      USE profiler
      USE integration_path_context

      IMPLICIT NONE

!*******************************************************************************
!  integration_path module parameters
!*******************************************************************************
!>  No qaudrature type.
      INTEGER, PARAMETER :: path_none_type    = -1
!>  Integrate the length in addition.
      INTEGER, PARAMETER :: path_add_type     = 0
!>  Integrate the path using Gauss Legendre Quadrature
      INTEGER, PARAMETER :: path_gleg_type    = 1
!>  Integrate the path using a hp approch.
      INTEGER, PARAMETER :: path_hp_glep_type = 2

!>  Default step size of the integration.
      REAL (rprec), PARAMETER :: path_default_dx = 0.0025

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) vertex
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Class containing the parameters of the integration method to use.
!-------------------------------------------------------------------------------
      TYPE path_int_class
!>  Integration method to use.
         INTEGER      :: method = path_none_type

!>  Step size to use.
         REAL (rprec) :: dx
!>  Quadrature interval length.
         REAL (rprec) :: length

!>  Quadrature weights.
         REAL (rprec), DIMENSION(:), POINTER :: weights
!>  Quadrature abscissas.
         REAL (rprec), DIMENSION(:), POINTER :: absc
      END TYPE

!-------------------------------------------------------------------------------
!>  Base class containing the parameters of the integration method to use.
!-------------------------------------------------------------------------------
      TYPE :: integration_path_class
!>  Step size to use.
         REAL (rprec) :: dx
      CONTAINS
         PROCEDURE    ::                                                       &
     &      integrate_paths => integration_path_integrate_paths
         PROCEDURE    ::                                                       &
     &      integrate_path => integration_path_integrate_path
         GENERIC      :: integrate => integrate_paths,                         &
     &                                integrate_path
      END TYPE

!-------------------------------------------------------------------------------
!>  Subclass to use Gauss Legendre Quadrature.
!-------------------------------------------------------------------------------
      TYPE, EXTENDS(integration_path_class) ::                                 &
     &   integration_path_gleg_class
!>  Quadrature weights.
         REAL (rprec), DIMENSION(:), POINTER :: weights
!>  Quadrature abscissas.
         REAL (rprec), DIMENSION(:), POINTER :: absc
      CONTAINS
         PROCEDURE                           ::                                &
     &      integrate_path => integration_path_gleg_integrate_path
         FINAL                               ::                                &
     &      integration_path_gleg_destruct
      END TYPE

!-------------------------------------------------------------------------------
!>  Subclass to use hp approch.
!-------------------------------------------------------------------------------
      TYPE, EXTENDS(integration_path_gleg_class) ::                            &
     &   integration_path_hp_glep_class
!>  Quadrature interval length.
         REAL (rprec) :: length
      CONTAINS
         PROCEDURE                           ::                                &
     &      integrate_path => integration_path_hp_gleg_integrate_path
      END TYPE

!-------------------------------------------------------------------------------
!>  A single point in space defined by an z, y, z coordinate. A vertex is
!>  structured as a singly linked list.
!-------------------------------------------------------------------------------
      TYPE vertex
!>  Position in cartiesian coordinates.
         REAL (rprec), DIMENSION(3) :: position
!>  Reference to the next vertex.
         TYPE (vertex), POINTER     :: next => null()
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Construction interface for integration_path_class constructor
!-------------------------------------------------------------------------------
      INTERFACE integration_path_class
         MODULE PROCEDURE integration_path_construct
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Construction interface for integration_path_gleg_class construct0r
!-------------------------------------------------------------------------------
      INTERFACE integration_path_gleg_class
         MODULE PROCEDURE integration_path_gleg_construct
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Construction interface for integration_path_gleg_class construct0r
!-------------------------------------------------------------------------------
      INTERFACE integration_path_hp_glep_class
         MODULE PROCEDURE integration_path_hp_glep_construct
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Construction interface using either @ref path_construct_int or
!>  @ref path_construct_vertex
!-------------------------------------------------------------------------------
      INTERFACE path_construct
         MODULE PROCEDURE path_construct_int,                                  &
     &                    path_construct_vertex
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Destruct interface using either @ref path_destruct_int or
!>  @ref path_destruct_vertex
!-------------------------------------------------------------------------------
      INTERFACE path_destruct
         MODULE PROCEDURE path_destruct_int,                                   &
     &                    path_destruct_vertex
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for checking the results of the unit tests
!-------------------------------------------------------------------------------
      INTERFACE check
         MODULE PROCEDURE check_log, check_real, check_int
      END INTERFACE

      PRIVATE :: check, check_real, check_log, check_int,                      &
     &           integrate, integrate_add, integrate_gleg,                     &
     &           integrate_hp_gleg, search, test_function

      CONTAINS

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Factory method to construct a path integrator using a method.
!>
!>  @param[in] method  Integartion method to use.
!>  @param[in] npoints Number of quadrature points to use.
!>  @param[in] length  Length of the interval.
!>  @returns A pointer to a constructed @ref integration_path_class object.
!-------------------------------------------------------------------------------
      FUNCTION make_integrator(method, npoints, length)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (integration_path_class), POINTER :: make_integrator
      CHARACTER (len=*), INTENT(in)           :: method
      INTEGER, INTENT(in)                     :: npoints
      REAL (rprec), INTENT(in)                :: length

!  local variables
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE (method)

         CASE ('add')
            make_integrator => integration_path_class()

         CASE ('gleg')
            make_integrator => integration_path_gleg_class(npoints)

         CASE ('hp_gleg')
            make_integrator => integration_path_hp_glep_class(npoints,         &
     &                                                        length)

      END SELECT

      CALL profiler_set_stop_time('make_integrator', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct an @ref integration_path_class object.
!>
!>  @returns A pointer to a constructed @ref integration_path_class object.
!-------------------------------------------------------------------------------
      FUNCTION integration_path_construct()

      IMPLICIT NONE

!  Declare Arguments
      CLASS (integration_path_class), POINTER ::                               &
     &   integration_path_construct

!  local variables
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(integration_path_construct)
      integration_path_construct%dx = path_default_dx

      CALL profiler_set_stop_time('integration_path_construct',                &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct an @ref integration_path_class object.
!>
!>  @param[in] npoints Number of quadrature points to use.
!>  @returns A pointer to a constructed @ref integration_path_class object.
!-------------------------------------------------------------------------------
      FUNCTION integration_path_gleg_construct(npoints)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (integration_path_gleg_class), POINTER ::                          &
     &   integration_path_gleg_construct
      INTEGER, INTENT(in)                          :: npoints

!  local variables
      REAL (rprec)                                 :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(integration_path_gleg_construct)
      integration_path_gleg_construct%dx = path_default_dx

      ALLOCATE(integration_path_gleg_construct%weights(npoints))
      ALLOCATE(integration_path_gleg_construct%absc(npoints))
      CALL path_get_gaussqad_weights(                                          &
     &        0.0_rprec, 1.0_rprec,                                            &
     &        integration_path_gleg_construct%absc,                            &
     &        integration_path_gleg_construct%weights)

      CALL profiler_set_stop_time('integration_path_gleg_construct',           &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct an @ref integration_path_class object.
!>
!>  @param[in] npoints Number of quadrature points to use.
!>  @param[in] length  Length of the interval.
!>  @returns A pointer to a constructed @ref integration_path_class object.
!-------------------------------------------------------------------------------
      FUNCTION integration_path_hp_glep_construct(npoints, length)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (integration_path_hp_glep_class), POINTER ::                       &
     &   integration_path_hp_glep_construct
      INTEGER, INTENT(in)                             :: npoints
      REAL (rprec), INTENT(in)                        :: length

!  local variables
      REAL (rprec)                                    :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(integration_path_hp_glep_construct)
      integration_path_hp_glep_construct%dx = path_default_dx
      integration_path_hp_glep_construct%length = length

      ALLOCATE(integration_path_hp_glep_construct%weights(npoints))
      ALLOCATE(integration_path_hp_glep_construct%absc(npoints))
      CALL path_get_gaussqad_weights(                                          &
     &        0.0_rprec, 1.0_rprec,                                            &
     &        integration_path_hp_glep_construct%absc,                         &
     &        integration_path_hp_glep_construct%weights)

      CALL profiler_set_stop_time('integration_path_hp_glep_construct',        &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a single @ref path_int_class.
!>
!>  Allocates memory and initializes a @ref path_int_class object.
!>
!>  @param[in] method  Integartion method to use.
!>  @param[in] npoints Number of quadrature points to use.
!>  @param[in] length  Length of the interval.
!>  @returns A pointer to a constructed @ref path_int_class object.
!-------------------------------------------------------------------------------
      FUNCTION path_construct_int(method, npoints, length)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (path_int_class), POINTER :: path_construct_int
      CHARACTER (len=*), INTENT(in)  :: method
      INTEGER, INTENT(in)            :: npoints
      REAL (rprec), INTENT(in)       :: length

!  Start of executable code
      ALLOCATE(path_construct_int)

      path_construct_int%dx = path_default_dx

      SELECT CASE (method)

         CASE ('add')
            path_construct_int%method = path_add_type
            path_construct_int%weights => null()
            path_construct_int%absc => null()
            RETURN

         CASE ('gleg')
            path_construct_int%method = path_gleg_type

         CASE ('hp_gleg')
            path_construct_int%method = path_hp_glep_type
            path_construct_int%length = length

         CASE DEFAULT
            path_construct_int%method = path_none_type

      END SELECT

      ALLOCATE(path_construct_int%weights(npoints))
      ALLOCATE(path_construct_int%absc(npoints))
      CALL path_get_gaussqad_weights(0.0_rprec, 1.0_rprec,                     &
     &                               path_construct_int%absc,                  &
     &                               path_construct_int%weights)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a single @ref vertex.
!>
!>  Allocates memory and initializes a @ref vertex object.
!>
!>  @param[in] position Cartesian position of the vertex object.
!>  @returns A pointer to a constructed @ref vertex object.
!-------------------------------------------------------------------------------
      FUNCTION path_construct_vertex(position)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vertex), POINTER                 :: path_construct_vertex
      REAL (rprec), DIMENSION(3), INTENT(in) :: position

!  Start of executable code
      ALLOCATE(path_construct_vertex)

      path_construct_vertex%position = position

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref integration_path_gleg_class object.
!>
!>  Deallocates memory and uninitializes a @ref integration_path_gleg_class
!>  object.
!>
!>  @param[inout] this A @ref integration_path_gleg_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE integration_path_gleg_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (integration_path_gleg_class), INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%weights)) THEN
         DEALLOCATE(this%weights)
         this%weights => null()
      END IF

      IF (ASSOCIATED(this%absc)) THEN
         DEALLOCATE(this%absc)
         this%absc => null()
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref path_int_class object.
!>
!>  Deallocates memory and uninitializes a @ref path_int_class object.
!>
!>  @param[inout] this A @ref path_int_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE path_destruct_int(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (path_int_class), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%weights)) THEN
         DEALLOCATE(this%weights)
         this%weights => null()
      END IF

      IF (ASSOCIATED(this%absc)) THEN
         DEALLOCATE(this%absc)
         this%absc => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref vertex object.
!>
!>  Deallocates memory and uninitializes a @ref vertex object. This recursively
!>  deconstructed the next vertex until the last in the linked list is found.
!>
!>  @param[inout] this A @ref vertex instance.
!-------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE path_destruct_vertex(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vertex), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%next)) THEN
         CALL path_destruct(this%next)
         this%next => null()
      END IF

      DEALLOCATE(this)
      this => null()

      END SUBROUTINE

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Append a @ref vertex to a path.
!>
!>  Recursively runs through the next vertex to find the last vertex. Once the
!>  last vertex is found, a new vertex is allocated and appended to the path.
!>  This allow works as a constructer and allocates the first vertex if needed.
!>
!>  @param[inout] this     Vertex to append path to.
!>  @param[in]    position Cartesian position of the vertex object.
!-------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE path_append_vertex(this, position)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vertex), POINTER                 :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: position

!  Start of executable code
      IF (ASSOCIATED(this)) THEN
         IF (ASSOCIATED(this%next)) THEN
            CALL path_append_vertex(this%next, position)
         ELSE
            this%next => path_construct(position)
         END IF
      ELSE
         this => path_construct(position)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Integrate along the paths.
!>
!>  Recursively runs through the next vertex to find the last vertex. Once the
!>  last vertex is found, integrate alone that path. The integrand is proveded
!>  by means of call back function.
!>
!>  @param[in] this                 In instance of a @ref integration_path_class
!>                                  instance.
!>  @param[in] path                 Starting @ref vertex to of the path to
!>                                  integrate.
!>  @param     integration_function Function pointer that defines the
!>                                  integrand.
!>  @param[in] context              Generic object that contains data for the
!>                                  integration function.
!>  @returns The total integrated path to the end.
!-------------------------------------------------------------------------------
      RECURSIVE FUNCTION integration_path_integrate_paths(this, path,          &
     &                      integration_function, context)                     &
     &  RESULT(total)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                                       :: total
      CLASS (integration_path_class), INTENT(in)         :: this
      TYPE (vertex), INTENT(in)                          :: path
      CLASS (integration_path_context_class), INTENT(in) :: context
      INTERFACE
         FUNCTION integration_function(context, xcart, dxcart,                 &
     &                                 length, dx)
         USE stel_kinds
         USE integration_path_context
         REAL (rprec)                                       ::                 &
     &      integration_function
         CLASS (integration_path_context_class), INTENT(in) :: context
         REAL (rprec), DIMENSION(3), INTENT(in)             :: xcart
         REAL (rprec), DIMENSION(3), INTENT(in)             :: dxcart
         REAL (rprec), INTENT(in)                           :: length
         REAL (rprec), INTENT(in)                           :: dx
         END FUNCTION
      END INTERFACE

!  local variables
      REAL (rprec)                                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (ASSOCIATED(path%next)) THEN
          total = this%integrate(path%next, integration_function,              &
     &                           context)
          total = this%integrate(context, path, path%next,                     &
     &                           integration_function)
      ELSE
         total = 0.0
      END IF

      CALL profiler_set_stop_time('integration_path_integrate_paths',          &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Integrate along the path.
!>
!>  Recursively runs through the next vertex to find the last vertex. Once the
!>  last vertex is found, integrate alone that path. The integrand is proveded
!>  by means of call back function.
!>
!>  @param[in] this                 In instance of a @ref path_int_class
!>                                  instance.
!>  @param[in] path                 Starting @ref vertex to of the path to
!>                                  integrate.
!>  @param     integration_function Function pointer that defines the
!>                                  integrand.
!>  @param[in] context              Generic object that contains data for the
!>                                  integration function.
!>  @returns The total integrated path to the end.
!-------------------------------------------------------------------------------
      RECURSIVE FUNCTION path_integrate(this, path,                            &
     &                                  integration_function, context)         &
     &  RESULT(total)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                      :: total
      TYPE (path_int_class), INTENT(in) :: this
      TYPE (vertex), INTENT(in)         :: path
      CHARACTER (len=1), INTENT(in)     :: context(:)
      INTERFACE
         FUNCTION integration_function(context, xcart, dxcart,                 &
     &                                 length, dx)
         USE stel_kinds
         REAL (rprec)                           :: integration_function
         CHARACTER (len=1), INTENT(in)          :: context(:)
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
         REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
         REAL (rprec), INTENT(in)               :: length
         REAL (rprec), INTENT(in)               :: dx
         END FUNCTION
      END INTERFACE

!  Start of executable code
      If (ASSOCIATED(path%next)) THEN
         total = path_integrate(this, path%next, integration_function,         &
     &                          context)
         total = total + integrate(this, context, path, path%next,             &
     &                             integration_function)
      ELSE
         total = 0.0
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Search along the path.
!>
!>  Recursively runs through the next vertex to find the last vertex. Once the
!>  last vertex is found, a new vertex is allocated and appended to the path.
!>  The integrand is proveded by means of call back function.
!>
!>  @param[in]  path            Starting vertex to begin search.
!>  @param      search_function Function pointer that defines the search
!>                              criteria.
!>  @param[in]  context         Generic object that contains data for the
!>                              integration function.
!>  @param[out] found           Signals if the condition was met.
!>  @returns The vertex position along the path where the search condition was
!>           found.
!-------------------------------------------------------------------------------
      RECURSIVE FUNCTION path_search(path, search_function, context,           &
     &                               found) RESULT(xcart)

      IMPLICIT NONE

      REAL (rprec), DIMENSION(3)    :: xcart
      TYPE (vertex), INTENT(in)     :: path
      CHARACTER (len=1), INTENT(in) :: context(:)
      LOGICAL, INTENT(out)          :: found
      INTERFACE
         FUNCTION search_function(context, xcart1, xcart2)
         USE stel_kinds
         LOGICAL                                :: search_function
         CHARACTER (len=1), INTENT(in)          :: context(:)
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart1
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart2
         END FUNCTION
      END INTERFACE

!  Start of executable code
      found = .false.

      IF (ASSOCIATED(path%next)) THEN
         found = search(context, path, path%next, search_function,             &
     &                  xcart)
         IF (.not.found) THEN
            xcart = path_search(path%next, search_function, context,           &
     &                          found)
         END IF
      END IF

      END FUNCTION

!*******************************************************************************
!  PRIVATE
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Line integrate between to points.
!>
!>  This chooses the specific integration method.
!>
!>  @param[in] this                 In instance of a @ref integration_path_class
!>                                  instance.
!>  @param[in] context              Generic object that contains data for the
!>                                  integration function.
!>  @param[in] vertex1              Starting point.
!>  @param[in] vertex2              Ending point.
!>  @param[in] dxinput              optional input allowing user to define dx
!>  @param     integration_function Function pointer that defines the
!>                                  integrand.
!>  @returns The path integrated value between the vertex1 and vertex2.
!-------------------------------------------------------------------------------
      FUNCTION integration_path_integrate_path(this, context, vertex1,         &
     &                                         vertex2,                        &
     &                                         integration_function)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                                       ::                    &
     &   integration_path_integrate_path
      CLASS (integration_path_class), INTENT(in)         :: this
      CLASS (integration_path_context_class), INTENT(in) :: context
      TYPE (vertex), INTENT(in)                          :: vertex1
      TYPE (vertex), INTENT(in)                          :: vertex2
      INTERFACE
         FUNCTION integration_function(context, xcart, dxcart,                 &
     &                                 length, dx)
         USE stel_kinds
         USE integration_path_context
         REAL (rprec) :: integration_function
         CLASS (integration_path_context_class), INTENT(in) :: context
         REAL (rprec), DIMENSION(3), INTENT(in)             :: xcart
         REAL (rprec), DIMENSION(3), INTENT(in)             :: dxcart
         REAL (rprec), INTENT(in)                           :: length
         REAL (rprec), INTENT(in)                           :: dx
         END FUNCTION
      END INTERFACE

!  local variables
      REAL (rprec), DIMENSION(3)                         :: xcart
      REAL (rprec), DIMENSION(3)                         :: dxcart
      REAL (rprec)                                       :: length
      REAL (rprec)                                       :: dx
      INTEGER                                            :: i
      INTEGER                                            :: nsteps
      REAL (rprec)                                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Determine the number of integration steps to take by dividing the path length
!  by the step length and rounding to the nearest integer.
      dxcart = vertex2%position - vertex1%position
      length = SQRT(DOT_PRODUCT(dxcart, dxcart))
      nsteps = INT(length/this%dx)
!  Choose the actual step size.
      dxcart = dxcart/nsteps
      dx = SQRT(DOT_PRODUCT(dxcart, dxcart))
      
!  Integrate the length in addition.
      integration_path_integrate_path = 0.0
      length = 0.0
      xcart = vertex1%position
      DO i = 1, nsteps
         xcart = xcart + dxcart
         length = length + dx
         integration_path_integrate_path =                                     &
     &      integration_path_integrate_path +                                  &
     &      integration_function(context, xcart, dxcart, length, dx)
      END DO

      CALL profiler_set_stop_time('integration_path_integrate_path',           &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Line integrate between to points.
!>
!>  This chooses the specific integration method.
!>
!>  @param[in] this                 In instance of a @ref integration_path_class
!>                                  instance.
!>  @param[in] context              Generic object that contains data for the
!>                                  integration function.
!>  @param[in] vertex1              Starting point.
!>  @param[in] vertex2              Ending point.
!>  @param[in] dxinput              optional input allowing user to define dx
!>  @param     integration_function Function pointer that defines the
!>                                  integrand.
!>  @returns The path integrated value between the vertex1 and vertex2.
!-------------------------------------------------------------------------------
      FUNCTION integration_path_gleg_integrate_path(this, context,             &
     &            vertex1, vertex2, integration_function)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                                       ::                    &
     &   integration_path_gleg_integrate_path
      CLASS (integration_path_gleg_class), INTENT(in)    :: this
      CLASS (integration_path_context_class), INTENT(in) :: context
      TYPE (vertex), INTENT(in)                          :: vertex1
      TYPE (vertex), INTENT(in)                          :: vertex2
      INTERFACE
         FUNCTION integration_function(context, xcart, dxcart,                 &
     &                                 length, dx)
         USE stel_kinds
         USE integration_path_context
         REAL (rprec) :: integration_function
         CLASS (integration_path_context_class), INTENT(in) :: context
         REAL (rprec), DIMENSION(3), INTENT(in)             :: xcart
         REAL (rprec), DIMENSION(3), INTENT(in)             :: dxcart
         REAL (rprec), INTENT(in)                           :: length
         REAL (rprec), INTENT(in)                           :: dx
         END FUNCTION
      END INTERFACE

!  local variables
      REAL (rprec), DIMENSION(3)                         :: xcart
      REAL (rprec), DIMENSION(3)                         :: dxcart
      REAL (rprec)                                       :: length
      REAL (rprec)                                       :: temp_length
      REAL (rprec)                                       :: dx
      INTEGER                                            :: i
      REAL (rprec)                                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Determine the number of integration steps to take by dividing the path length
!  by the step length and rounding to the nearest integer.
      dxcart = vertex2%position - vertex1%position
      length = SQRT(DOT_PRODUCT(dxcart, dxcart))
      dxcart(:) = dxcart(:)/length
      
!  Integrate the length using Gauss-Legendre quadrature
      integration_path_gleg_integrate_path = 0.0
      xcart = vertex1%position
      DO i = 1, SIZE(this%weights)
         xcart = vertex1%position + this%absc(i)*dxcart*length
         temp_length = this%absc(i)*length

         integration_path_gleg_integrate_path =                                &
     &      integration_path_gleg_integrate_path +                             &
     &      this%weights(i)*integration_function(context, xcart,               &
     &                                           dxcart, temp_length,          &
     &                                           1.0_rprec)
      END DO
      integration_path_gleg_integrate_path =                                   &
     &   integration_path_gleg_integrate_path*length

      CALL profiler_set_stop_time(                                             &
     &        'integration_path_gleg_integrate_path', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Line integrate between to points.
!>
!>  This chooses the specific integration method.
!>
!>  @param[in] this                 In instance of a @ref integration_path_class
!>                                  instance.
!>  @param[in] context              Generic object that contains data for the
!>                                  integration function.
!>  @param[in] vertex1              Starting point.
!>  @param[in] vertex2              Ending point.
!>  @param[in] dxinput              optional input allowing user to define dx
!>  @param     integration_function Function pointer that defines the
!>                                  integrand.
!>  @returns The path integrated value between the vertex1 and vertex2.
!-------------------------------------------------------------------------------
      FUNCTION integration_path_hp_gleg_integrate_path(this, context,          &
     &            vertex1, vertex2, integration_function)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                                       ::                    &
     &   integration_path_hp_gleg_integrate_path
      CLASS (integration_path_hp_glep_class), INTENT(in) :: this
      CLASS (integration_path_context_class), INTENT(in) :: context
      TYPE (vertex), INTENT(in)                          :: vertex1
      TYPE (vertex), INTENT(in)                          :: vertex2
      INTERFACE
         FUNCTION integration_function(context, xcart, dxcart,                 &
     &                                 length, dx)
         USE stel_kinds
         USE integration_path_context
         REAL (rprec) :: integration_function
         CLASS (integration_path_context_class), INTENT(in) :: context
         REAL (rprec), DIMENSION(3), INTENT(in)             :: xcart
         REAL (rprec), DIMENSION(3), INTENT(in)             :: dxcart
         REAL (rprec), INTENT(in)                           :: length
         REAL (rprec), INTENT(in)                           :: dx
         END FUNCTION
      END INTERFACE

!  local variables
      REAL (rprec), DIMENSION(3)                         :: xcart
      REAL (rprec), DIMENSION(3)                         :: xcart1
      REAL (rprec), DIMENSION(3)                         :: dxcart
      REAL (rprec)                                       :: length
      REAL (rprec)                                       :: lengthp
      REAL (rprec)                                       :: temp_length
      INTEGER                                            :: i
      INTEGER                                            :: j
      INTEGER                                            :: nsteps
      REAL (rprec)                                       :: int_length
      REAL (rprec)                                       :: integratej
      REAL (rprec)                                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Determine the number of integration steps to take by dividing the path length
!  by the step length and rounding to the nearest integer.
      dxcart = vertex2%position - vertex1%position
      length = SQRT(DOT_PRODUCT(dxcart, dxcart))
      dxcart(:) = dxcart(:)/length

      nsteps =  INT(length/this%length)
      int_length = length/REAL(nsteps, rprec)

!  Integrate the length using Gauss-Legendre quadrature
      integration_path_hp_gleg_integrate_path = 0.0
      xcart = vertex1%position
      lengthp = 0.0
      DO j = 1, nsteps
         integratej = 0.0
         DO i = 1, SIZE(this%weights)
            xcart1 = xcart + this%absc(i)*dxcart*int_length
            temp_length = lengthp + this%absc(i)*int_length
            integratej = integratej +                                          &
     &         this%weights(i)*integration_function(context, xcart1,           &
     &                                              dxcart, temp_length,       &
     &                                              1.0_rprec)
         END DO

         integration_path_hp_gleg_integrate_path =                             &
     &      integration_path_hp_gleg_integrate_path +                          &
     &      integratej*int_length
         lengthp = lengthp + int_length
         xcart = xcart + dxcart*int_length
      END DO

      CALL profiler_set_stop_time(                                             &
     &        'integration_path_hp_gleg_integrate_path', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Line integrate between to points.
!>
!>  This chooses the specific integration method.
!>
!>  @param[in] this                 In instance of a @ref path_int_class
!>                                  instance.
!>  @param[in] context              Generic object that contains data for the
!>                                  integration function.
!>  @param[in] vertex1              Starting point.
!>  @param[in] vertex2              Ending point.
!>  @param[in] dxinput              optional input allowing user to define dx
!>  @param     integration_function Function pointer that defines the
!>                                  integrand.
!>  @returns The path integrated value between the vertex1 and vertex2.
!-------------------------------------------------------------------------------
      FUNCTION integrate(this, context, vertex1, vertex2,                      &
     &                   integration_function)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                      :: integrate
      TYPE (path_int_class), INTENT(in) :: this
      CHARACTER(len=1), INTENT(in)      :: context(:)
      TYPE (vertex), INTENT(in)         :: vertex1
      TYPE (vertex), INTENT(in)         :: vertex2
      INTERFACE
         FUNCTION integration_function(context, xcart, dxcart,                 &
     &                                 length, dx)
         USE stel_kinds
         REAL (rprec) :: integration_function
         CHARACTER (len=1), INTENT(in)          :: context(:)
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
         REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
         REAL (rprec), INTENT(in)               :: length
         REAL (rprec), INTENT(in)               :: dx
         END FUNCTION
      END INTERFACE

!  Start of executable code
      SELECT CASE (this%method)

         CASE (path_add_type)
            integrate = integrate_add(this, context, vertex1, vertex2,         &
     &                                integration_function)

         CASE (path_gleg_type)
            integrate = integrate_gleg(this, context, vertex1,                 &
     &                                 vertex2, integration_function)

         CASE (path_hp_glep_type)
            integrate =                                                        &
     &         integrate_hp_gleg(this, context, vertex1, vertex2,              &
     &                           integration_function)

         CASE DEFAULT
            integrate = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Line integrate between to points.
!>
!>  This divids the straight line path defined by two vertices and line
!>  integrates the function. The integrand is proveded by means of call back
!>  function.
!>
!>  @param[in] this                 In instance of a @ref path_int_class
!>                                  instance.
!>  @param[in] context              Generic object that contains data for the
!>                                  integration function.
!>  @param[in] vertex1              Starting point.
!>  @param[in] vertex2              Ending point.
!>  @param     integration_function Function pointer that defines the
!>                                  integrand.
!>  @returns The path integrated value between the vertex1 and vertex2.
!-------------------------------------------------------------------------------
      FUNCTION integrate_add(this, context, vertex1, vertex2,                  &
     &                       integration_function)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                      :: integrate_add
      TYPE (path_int_class), INTENT(in) :: this
      CHARACTER(len=1), INTENT(in)      :: context(:)
      TYPE (vertex), INTENT(in)         :: vertex1
      TYPE (vertex), INTENT(in)         :: vertex2
      INTERFACE
         FUNCTION integration_function(context, xcart, dxcart,                 &
     &                                 length, dx)
         USE stel_kinds
         REAL (rprec) :: integration_function
         CHARACTER (len=1), INTENT(in)          :: context(:)
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
         REAL (rprec), INTENT(in)               :: length, dx
         END FUNCTION
      END INTERFACE

!  local variables
      REAL (rprec), DIMENSION(3)        :: xcart
      REAL (rprec), DIMENSION(3)        :: dxcart
      REAL (rprec)                      :: length
      REAL (rprec)                      :: dx
      INTEGER                           :: i, nsteps

!  Start of executable code
!  Determine the number of integration steps to take by dividing the path length
!  by the step length and rounding to the nearest integer.
      dxcart = vertex2%position - vertex1%position
      length = SQRT(DOT_PRODUCT(dxcart, dxcart))
      nsteps = INT(length/this%dx)
!  Choose the actual step size.
      dxcart = dxcart/nsteps
      dx = SQRT(DOT_PRODUCT(dxcart, dxcart))

      
!  Integrate the length in addition.
      integrate_add = 0.0
      length = 0.0
      xcart = vertex1%position
      DO i = 1, nsteps
         xcart = xcart + dxcart
         length = length + dx
         integrate_add = integrate_add                                         &
     &                 + integration_function(context, xcart, dxcart,          &
     &                                        length, dx)
      END DO

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Integrate between to points using Guass-Legendre quadrature.
!>
!>  This divids the straight line path defined by two vertices and line
!>  integrates the function. The integrand is proveded by means of call back
!>  function.
!>
!>  @param[in] this                 In instance of a @ref path_int_class
!>                                  instance.
!>  @param[in] context              Generic object that contains data for the
!>                                  integration function.
!>  @param[in] vertex1              Starting point.
!>  @param[in] vertex2              Ending point.
!>  @param     integration_function Function pointer that defines the
!>                                  integrand.
!>  @returns The path integrated value between the vertex1 and vertex2.
!-------------------------------------------------------------------------------
      FUNCTION integrate_gleg(this, context, vertex1, vertex2,                 &
     &                        integration_function)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                      :: integrate_gleg
      TYPE (path_int_class), INTENT(in) :: this
      CHARACTER(len=1), INTENT(in)      :: context(:)
      TYPE (vertex), INTENT(in)         :: vertex1
      TYPE (vertex), INTENT(in)         :: vertex2
      INTERFACE
         FUNCTION integration_function(context, xcart, dxcart,                 &
     &                                 length, dx)
         USE stel_kinds
         REAL (rprec) :: integration_function
         CHARACTER (len=1), INTENT(in)          :: context(:)
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
         REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
         REAL (rprec), INTENT(in)               :: length
         REAL (rprec), INTENT(in)               :: dx
         END FUNCTION
      END INTERFACE

!  local variables
      REAL (rprec), DIMENSION(3)        :: xcart
      REAL (rprec), DIMENSION(3)        :: xcart1
      REAL (rprec), DIMENSION(3)        :: dxcart
      REAL (rprec)                      :: length
      REAL (rprec)                      :: temp_length
      REAL (rprec)                      :: dx
      INTEGER                           :: i

!  Start of executable code

!  Determine the number of integration steps to take by dividing the path length
!  by the step length and rounding to the nearest integer.
      dxcart = vertex2%position - vertex1%position
      length = SQRT(DOT_PRODUCT(dxcart, dxcart))
      dxcart(:) = dxcart(:)/length

!  Integrate the length using Gauss-Legendre quadrature
      integrate_gleg = 0.0
      xcart = vertex1%position
      DO i = 1, SIZE(this%weights)
         xcart = vertex1%position + this%absc(i)*dxcart*length
         temp_length = this%absc(i)*length

         integrate_gleg = integrate_gleg +                                     &
     &      this%weights(i)*integration_function(context, xcart,               &
     &                                           dxcart, temp_length,          &
     &                                           1.0_rprec)
      END DO
      integrate_gleg = integrate_gleg*length

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Integrate between to points using hp Guass-Legendre quadrature
!>  method.
!>
!>  This divids the straight line path defined by two vertices into a set of
!>  intergval of length len. Then Gauss-Legendre quadrature is used to
!>  integrate the function across each interval. The integrand is provided by
!>  means of call back function.
!>
!>  @param[in] this                 In instance of a @ref path_int_class
!>                                  instance.
!>  @param[in] context              Generic object that contains data for the
!>                                  integration function.
!>  @param[in] vertex1              Starting point.
!>  @param[in] vertex2              Ending point.
!>  @param     integration_function Function pointer that defines the
!>                                  integrand.
!>  @returns The path integrated value between the vertex1 and vertex2.
!-------------------------------------------------------------------------------
      FUNCTION integrate_hp_gleg(this, context, vertex1, vertex2,              &
     &                           integration_function)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                      :: integrate_hp_gleg
      TYPE (path_int_class), INTENT(in) :: this
      CHARACTER(len=1), INTENT(in)      :: context(:)
      TYPE (vertex), INTENT(in)         :: vertex1
      TYPE (vertex), INTENT(in)         :: vertex2
      INTERFACE
         FUNCTION integration_function(context, xcart, dxcart,                 &
     &                                 length, dx)
         USE stel_kinds
         REAL (rprec) :: integration_function
         CHARACTER (len=1), INTENT(in)          :: context(:)
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
         REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
         REAL (rprec), INTENT(in)               :: length
         REAL (rprec), INTENT(in)               :: dx
         END FUNCTION
      END INTERFACE

!  local variables
      REAL (rprec), DIMENSION(3)        :: xcart
      REAL (rprec), DIMENSION(3)        :: xcart1
      REAL (rprec), DIMENSION(3)        :: dxcart
      REAL (rprec)                      :: length
      REAL (rprec)                      :: lengthp
      REAL (rprec)                      :: temp_length
      INTEGER                           :: i
      INTEGER                           :: j
      INTEGER                           :: nsteps
      REAL (rprec)                      :: int_length
      REAL (rprec)                      :: integratej

!  Start of executable code

!  Determine the number of integration steps to take by dividing the path length
!  by the step length and rounding to the nearest integer.
      dxcart = vertex2%position - vertex1%position
      length = SQRT(DOT_PRODUCT(dxcart, dxcart))
      dxcart(:) = dxcart(:)/length

      nsteps =  INT(length/this%length)
      int_length = length/REAL(nsteps, rprec)

!  Integrate the length using Gauss-Legendre quadrature
      integrate_hp_gleg = 0.0
      xcart = vertex1%position
      lengthp = 0.0
      DO j = 1, nsteps
         integratej = 0.0
         DO i = 1, SIZE(this%weights)
            xcart1 = xcart + this%absc(i)*dxcart*int_length
            temp_length = lengthp + this%absc(i)*int_length
            integratej = integratej +                                          &
     &         this%weights(i)*integration_function(context, xcart1,           &
     &                                              dxcart, temp_length,       &
     &                                              1.0_rprec)
         END DO

         integrate_hp_gleg = integrate_hp_gleg + integratej*int_length
         lengthp = lengthp + int_length
         xcart = xcart + dxcart*int_length
      END DO

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Search line between to points.
!>
!>  This divides the straight line path defined by two vertices and searched for
!>  a condition. The search criteria is proveded by means of call back function.
!>
!>  @param[in] context          Generic object that contains data for the search
!>                              function.
!>  @param[in]  vertex1         Starting point.
!>  @param[in]  vertex2         Ending point.
!>  @param      search_function Function pointer that defines the search
!>                              criteria.
!>  @param[out] xcart           Point where the search criteria was found.
!>  @returns True if the criteria was met between vertex1 and vertex2.
!-------------------------------------------------------------------------------
      FUNCTION search(context, vertex1, vertex2, search_function, xcart)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                                 :: search
      CHARACTER(len=1), INTENT(in)            :: context(:)
      TYPE (vertex), INTENT(in)               :: vertex1
      TYPE (vertex), INTENT(in)               :: vertex2
      REAL (rprec), DIMENSION(3), INTENT(out) :: xcart
      INTERFACE
         FUNCTION search_function(context, xcart1, xcart2)
         USE stel_kinds
         LOGICAL                                :: search_function
         CHARACTER (len=1), INTENT(in)          :: context(:)
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart1
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart2
         END FUNCTION
      END INTERFACE

!  local variables
      REAL (rprec), DIMENSION(3)              :: dxcart
      REAL (rprec)                            :: dx
      INTEGER                                 :: i
      INTEGER                                 :: nsteps

!  local parameters
      REAL (rprec), PARAMETER                 :: dxCourse = 0.01
      REAL (rprec), PARAMETER                 :: dxFine = 1.0E-20

!  Start of executable code
!  Determine the number of integration steps to take by dividing the path length
!  by the step length and rounding to the nearest integer.
      dxcart = vertex2%position - vertex1%position
      nsteps = INT(SQRT(DOT_PRODUCT(dxcart, dxcart))/dxCourse)

!  Choose the actual step size.
      dxcart = dxcart/nsteps

      search = .false.
      xcart = vertex1%position

!  Linearly search the line until an interval containing the point is detected.
      DO i = 1, nsteps
         search = search_function(context, xcart, xcart + dxcart)
         IF (search) THEN

!  Found an interval. Bisect the interval until the length is machine precision.
            DO WHILE (SQRT(DOT_PRODUCT(dxcart, dxcart)) .gt. dxFine)
               dxcart = dxcart/2.0
               IF (.not.search_function(context, xcart,                        &
     &                                  xcart + dxcart)) THEN
                  xcart = xcart + dxcart
               END IF
            END DO

            xcart = xcart + dxcart/2.0
            RETURN

         END IF

         xcart = xcart + dxcart
      END DO

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @breif Calculate the weights and abscissas for Gauss-Legendre integration.
!>
!>
!>  This suprogram calculated the weights and abscissas for Gauss-Legendre
!>  integration. The subroutine is adapted from NIMROD which was adapted
!>  from Numerical Recipies, 2ed.,  Cambridge Press. 
!>
!>  @param[in]  a         Start point for integration.
!>  @param[in]  b         End point for integration.
!>  @param[out] abscissas Array of abscissas.
!>  @param[out] weights   Array of weights.
!-------------------------------------------------------------------------------
      SUBROUTINE path_get_gaussqad_weights(a, b, abscissas, weights)
      
      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), INTENT(in)               :: a
      REAL (rprec), INTENT(in)               :: b
      REAL(rprec), DIMENSION(:), INTENT(out) :: abscissas
      REAL(rprec), DIMENSION(:), INTENT(out) :: weights

!  local variables
      INTEGER                                :: i
      INTEGER                                :: j
      INTEGER                                :: m
      REAL (rprec)                           :: p1
      REAL (rprec)                           :: p2
      REAL (rprec)                           :: p3
      REAL (rprec)                           :: pp
      REAL (rprec)                           :: xl
      REAL (rprec)                           :: xm
      REAL (rprec)                           :: current
      REAL (rprec)                           :: last

!  local parameters
      REAL (rprec), PARAMETER                :: eps = 1.0E-14

!  Start of executable code
      m = (SIZE(weights) + 1)/2

!  Change basis from [-1,1] to [a.b]
      xm = 0.5*(b + a)
      xl = 0.5*(b - a)

      DO i = 1, m
!  current guess for i-th root of P_n(x)
         current = COS(pi*(REAL(i, rprec) - 0.25_rprec)/                       &
     &                 (REAL(SIZE(weights), rprec) + 0.5_rprec))

!  Newtons methods to find the roots  of P_n(x)
         DO

!  Calculate P_n(x)
            p1 = 1.0
            p2 = 0.0
            DO j = 1, SIZE(weights)
               p3 = p2
               p2 = p1
               p1 = (REAL(2.0*j - 1, rprec)*current*p2 -                       &
     &               REAL(j - 1, rprec)*p3)/REAL(j, rprec)
            END DO

!  P'_n = n(xP_n - P_(n-1)/(x^2-1)
            pp = REAL(SIZE(weights), rprec)*(current*p1 - p2)                  &
     &         / (current**2 - 1._rprec)

!  Previous guess for i-th root
            last = current
            current = last - p1/pp
            IF (ABS(current - last) .le. eps) EXIT
         END DO

         abscissas(i) = xm - xl*current
         abscissas(SIZE(weights) + 1 - i) = xm + xl*current

!  wi = 2/[(1 - x^2)P'n(x)^2]
         weights(i) = 2.0*xl/((1.0 - current*current)*pp*pp)
         weights(SIZE(weights) + 1 - i) = weights(i)
      END DO

      END SUBROUTINE

!*******************************************************************************
!  UNIT TESTS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Path unit test function.
!>
!>  This runs the associated unit tests and returns the result.
!>
!>  @returns True if the tests pass and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION path_test()

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                                         :: path_test

!  local variables
      REAL (rprec)                                    :: result
      TYPE (vertex), POINTER                          ::                       &
     &   test_path => null()
      CLASS (integration_path_class), POINTER         :: int_params
      CLASS (integration_path_context_class), POINTER :: context
      REAL (rprec), DIMENSION(3)                      :: absc
      REAL (rprec), DIMENSION(3)                      :: wgts
      REAL (rprec)                                    :: work

!  Start of executable code
!  Test to make sure the vertices begin in an unallocated state.
      path_test = check(.false., ASSOCIATED(test_path), 1, 'ASSOCIATED')
      IF (.not.path_test) RETURN

!  Test to make sure first vertex is created.
      CALL path_append_vertex(test_path,                                       &
     &                        (/ 1.0_rprec, 2.0_rprec, 3.0_rprec /))
      path_test = check(.true., ASSOCIATED(test_path), 1,                      &
     &                  'path_append_vertex')
      IF (.not.path_test) RETURN
      path_test = check(1.0_rprec, test_path%position(1), 2,                   &
     &                  'path_append_vertex')
      IF (.not.path_test) RETURN
      path_test = check(2.0_rprec, test_path%position(2), 3,                   &
     &                  'path_append_vertex')
      IF (.not.path_test) RETURN
      path_test = check(3.0_rprec, test_path%position(3), 4,                   &
     &                  'path_append_vertex')
      IF (.not.path_test) RETURN

!  Test to make sure second vertex is appended.
      CALL path_append_vertex(test_path,                                       &
     &                        (/ 4.0_rprec, 5.0_rprec, 6.0_rprec /))
      path_test = check(.true., ASSOCIATED(test_path%next), 5,                 &
     &                  'path_append_vertex')
      IF (.not.path_test) RETURN
      path_test = check(1.0_rprec, test_path%position(1), 6,                   &
     &                  'path_append_vertex')
      IF (.not.path_test) RETURN
      path_test = check(2.0_rprec, test_path%position(2), 7,                   &
     &                  'path_append_vertex')
      IF (.not.path_test) RETURN
      path_test = check(3.0_rprec, test_path%position(3), 8,                   &
     &                  'path_append_vertex')
      IF (.not.path_test) RETURN
      path_test = check(4.0_rprec, test_path%next%position(1), 9,              &
     &                  'path_append_vertex')
      IF (.not.path_test) RETURN
      path_test = check(5.0_rprec, test_path%next%position(2), 10,             &
     &                  'path_append_vertex')
      IF (.not.path_test) RETURN
      path_test = check(6.0_rprec, test_path%next%position(3), 11,             &
     &                  'path_append_vertex')
      IF (.not.path_test) RETURN

!  Test to make sure third vertex is appended.
      CALL path_append_vertex(test_path,                                       &
     &                        (/ 7.0_rprec, 8.0_rprec, 9.0_rprec /))
      path_test = check(.true., ASSOCIATED(test_path%next%next),               &
     &                  12, 'path_append_vertex')
      IF (.not.path_test) RETURN
      path_test = check(7.0_rprec, test_path%next%next%position(1), 12,        &
     &                  'path_append_vertex')
      IF (.not.path_test) RETURN
      path_test = check(8.0_rprec, test_path%next%next%position(2), 13,        &
     &                  'path_append_vertex')
      IF (.not.path_test) RETURN
      path_test = check(9.0_rprec, test_path%next%next%position(3), 14,        &
     &                  'path_append_vertex')
      IF (.not.path_test) RETURN

!  Test to make sure path object is destroyed.
      CALL path_destruct(test_path)
      path_test = check(.false., ASSOCIATED(test_path), 1,                     &
     &                  'path_destruct')
      IF (.not.path_test) RETURN

!  Test to make sure a path_in_class object was allocated.
      int_params => make_integrator('add', 0, 0.0_rprec)
      path_test = check(.true., ASSOCIATED(int_params), 1,                     &
     &                  'path_construct_int')
      IF (.not.path_test) RETURN

!  Test path integration.
      CALL path_append_vertex(test_path,                                       &
     &                        (/ 2.0_rprec, 0.0_rprec, 0.0_rprec /))
      CALL path_append_vertex(test_path,                                       &
     &                        (/ 0.0_rprec, 0.0_rprec, 0.0_rprec /))
      result = int_params%integrate(test_path, test_function, context)
      path_test = check(2.0_rprec, result, 1, 'path_integrate')
      DEALLOCATE(int_params)

!  Test gll integrator weights and abscissas
      CALL path_get_gaussqad_weights(-1.0_rprec, 1.0_rprec, absc, wgts)
      work = SQRT(0.6_rprec)
      path_test = check(-1.0_rprec*work, absc(1), 1, 'get_gaussqad_a')
      path_test = check( 0.0_rprec,      absc(2), 2, 'get_gaussqad_a')
      path_test = check(           work, absc(3), 3, 'get_gaussqad_a')
      path_test = check(5._rprec/9._rprec, wgts(1), 4, 'get_gaussqad_w')
      path_test = check(8._rprec/9._rprec, wgts(2), 5, 'get_gaussqad_w')
      path_test = check(5._rprec/9._rprec, wgts(3), 6, 'get_gaussqad_w')
      IF (.not.path_test) RETURN

!  Test gleg integrator.
      int_params => make_integrator('gleg', 100, 0.0_rprec)
      result = int_params%integrate(test_path, test_function, context)
      path_test = check(2.0_rprec, result, 1, 'path_integrate_gleg')
      IF (.not.path_test) RETURN
      DEALLOCATE(int_params)

!  Test hp_gleg integrator.
      int_params => make_integrator('hp_gleg', 100, 0.01_rprec)
      result = int_params%integrate(test_path, test_function, context)
      path_test = check(2.0_rprec, result, 1,                                &
     &                  'path_integrate_hp_gleg')
      IF (.not.path_test) RETURN
      CALL path_destruct(test_path)
      DEALLOCATE(int_params)

      END FUNCTION

!*******************************************************************************
!  CHECK FUNCTIONS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Check a logical value.
!>
!>  Checks that the expected value matches the recieved. Otherwise report an
!>  error.
!>
!>  @param[in] expected The known value.
!>  @param[in] received The known test value.
!>  @param[in] testNum  The number of the test.
!>  @param[in] name     The name of the test.
!>  @returns True if the check passes and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION check_log(expected, received, testNum, name)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: check_log
      LOGICAL, INTENT(in)           :: expected
      LOGICAL, INTENT(in)           :: received
      INTEGER, INTENT(in)           :: testNum
      CHARACTER (LEN=*), INTENT(in) :: name

!  Start of executable code
      check_log = expected .eqv. received
      IF (.not.check_log) THEN
         WRITE(*,*) "integration_path.f: ", name, " test", testNum,            &
     &              "failed."
         WRITE(*,*) "Expected", expected, "Received", received
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Check a real value.
!>
!>  Checks that the expected value matches the recieved. Otherwise report an
!>  error.
!>
!>  @param[in] expected The known value.
!>  @param[in] received The known test value.
!>  @param[in] testNum  The number of the test.
!>  @param[in] name     The name of the test.
!>  @returns True if the check passes and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION check_real(expected, received, testNum, name)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: check_real
      REAL(rprec), INTENT(in)       :: expected
      REAL(rprec), INTENT(in)       :: received
      INTEGER, INTENT(in)           :: testNum
      CHARACTER (LEN=*), INTENT(in) :: name

! local variables
      REAL(rprec)                   :: eps = 1.0e-8

!  Start of executable code
      check_real = (ABS(expected - received) .lt. eps)
      IF (.not.check_real) THEN
         WRITE(*,*) "integration_path.f: ", name, " test", testNum,            &
     &              "failed."
         WRITE(*,*) "Expected", expected, "Received", received
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Check an integer value.
!>
!>  Checks that the expected value matches the recieved. Otherwise report an
!>  error.
!>
!>  @param[in] expected The known value.
!>  @param[in] received The known test value.
!>  @param[in] testNum  The number of the test.
!>  @param[in] name     The name of the test.
!>  @returns True if the check passes and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION check_int(expected, received, testNum, name)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: check_int
      INTEGER, INTENT(in)           :: expected
      INTEGER, INTENT(in)           :: received
      INTEGER, INTENT(in)           :: testNum
      CHARACTER (LEN=*), INTENT(in) :: name

!  Start of executable code
      check_int = expected .eq. received
      IF (.not.check_int) THEN
         WRITE(*,*) "integration_path.f: ", name, " test", testNum,            &
     &              "failed."
         WRITE(*,*) "Expected", expected, "Received", received
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Call back function to test the integration.
!>
!>  This defined a function of f(x) = 1.
!>
!>  @param[in] context Generic object that contains data for the integration
!>                     function.
!>  @param[in] xcart   Current point in the integration.
!>  @param[in] dxcart  Vector direction of the current integration path.
!>  @param[in] length  Length along the current integration.
!>  @param[in] dx      Scaler length of the current integration path.
!>  @returns One
!-------------------------------------------------------------------------------
      FUNCTION test_function(context, xcart, dxcart, length, dx)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                                       ::                    &
     &   test_function
      CLASS (integration_path_context_class), INTENT(in) :: context
      REAL (rprec), DIMENSION(3), INTENT(in)             :: xcart
      REAL (rprec), DIMENSION(3), INTENT(in)             :: dxcart
      REAL (rprec), INTENT(in)                           :: length
      REAL (rprec), INTENT(in)                           :: dx

!  Start of executable code
      test_function = dx

      END FUNCTION

      END MODULE

!*******************************************************************************
!>  @file gradient_descent.f
!>  @brief Contains module @ref gradient_descent
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Module is part of the LIBSTELL. This modules contains code to define and
!>  integrate along an arbitray path.
!*******************************************************************************

      MODULE gradient_descent
      USE stel_kinds
      USE profiler

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) context
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class containing the parameters for the gradient_descent.
!-------------------------------------------------------------------------------
      TYPE :: gradient_descent_class
!>  Step size to use.
         REAL (rprec)                        :: step
!>  Minimization tolarance.
         REAL (rprec)                        :: tolarance
!>  Variables.
         REAL (rprec), DIMENSION(:), POINTER :: x_var => null()
      CONTAINS
         PROCEDURE :: chi2 => gradient_descent_chi2
         PROCEDURE :: minimize => gradient_descent_minimize
         FINAL     :: gradient_descent_destruct
      END TYPE

!-------------------------------------------------------------------------------
!>  Base class containing the parameters for the gradient_descent.
!-------------------------------------------------------------------------------
      TYPE, EXTENDS(gradient_descent_class) ::                                  &
     &   gradient_descent_test_class
      CONTAINS
         PROCEDURE :: chi2 => gradient_descent_test_chi2
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      INTERFACE gradient_descent_test_class
         MODULE PROCEDURE gradient_descent_test_construct
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for checking the results of the unit tests
!-------------------------------------------------------------------------------
      INTERFACE check
         MODULE PROCEDURE check_log, check_real
      END INTERFACE

      PRIVATE :: check, check_log, check_real

      CONTAINS

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref gradient_descent_test_class object.
!>
!>  Deallocates memory and uninitializes a @ref gradient_descent_test_class
!>  object.
!>
!>  @param[inout] this A @ref gradient_descent_test_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE gradient_descent_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (gradient_descent_class), INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%x_var)) THEN
         DEALLOCATE(this%x_var)
         this%x_var => null()
      END IF

      END SUBROUTINE

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Compute the χ^2 value and the gradients.
!>
!>  This method should not be called directly and needs to be overloaded in a
!>  subclass.
!>
!>  @param[in]  this A @ref gradient_descent_context instance.
!>  @param[out] chi2 The χ^2 value and the gradients.
!-------------------------------------------------------------------------------
      SUBROUTINE gradient_descent_chi2(this, chi2)
      USE v3_utilities

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) ::  gradient_descent_chi
      CLASS (gradient_descent_class), INTENT(in) :: this
      REAL (rprec), DIMENSION(:), INTENT(out)    :: chi2

!  local variables
      REAL (rprec)                               :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      CALL assert(.false., 'This method needs to be ' //                       &
     &                     'overloaded in a subclass.')

      CALL profiler_set_stop_time('gradient_descent_chi2',                     &
     &                            start_time)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Gradient descent to minimize the χ^2 function.
!>
!>  @param[inout] this A @ref gradient_descent_context instance.
!>  @returns The χ^2 residule.
!-------------------------------------------------------------------------------
      FUNCTION gradient_descent_minimize(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: gradient_descent_minimize
      CLASS (gradient_descent_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                               :: start_time
      REAL (rprec), DIMENSION(:), ALLOCATABLE    :: chi2_temp

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE (chi2_temp(SIZE(this%x_var) + 1))

      CALL this%chi2(chi2_temp)
      DO WHILE (chi2_temp(1) .gt. this%tolarance)
         this%x_var = this%x_var - this%step*chi2_temp(2:)
         CALL this%chi2(chi2_temp)
      END DO

      gradient_descent_minimize = chi2_temp(1)

      DEALLOCATE (chi2_temp)

      CALL profiler_set_stop_time('gradient_descent_minimize',                 &
     &                            start_time)

      END FUNCTION

!*******************************************************************************
!  UNIT TESTS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref vertex.
!>
!>  Allocates memory and initializes a @ref vertex object.
!>
!>  @param[in] position Cartesian position of the vertex object.
!>  @returns A pointer to a constructed @ref vertex object.
!-------------------------------------------------------------------------------
      FUNCTION gradient_descent_test_construct()

      IMPLICIT NONE

!  Declare Arguments
      CLASS (gradient_descent_test_class), POINTER ::                          &
     &   gradient_descent_test_construct

!  local variables
      REAL (rprec)                                 :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(gradient_descent_test_construct)

      ALLOCATE(gradient_descent_test_construct%x_var(1))
      gradient_descent_test_construct%step = 1.0
      gradient_descent_test_construct%tolarance = 1.0E-20

      CALL profiler_set_stop_time('gradient_descent_test_construct',           &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Compute the χ^2 value and the gradients.
!>
!>  This method should not be called directly and needs to be overloaded in a
!>  subclass.
!>
!>  @param[in]  this A @ref gradient_descent_context instance.
!>  @param[out] chi2 The χ^2 value and the gradients.
!-------------------------------------------------------------------------------
      SUBROUTINE gradient_descent_test_chi2(this, chi2)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) ::  gradient_descent_chi
      CLASS (gradient_descent_test_class), INTENT(in) :: this
      REAL (rprec), DIMENSION(:), INTENT(out)         :: chi2

!  local variables
      REAL (rprec)                                    :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      chi2(1) = (10.0*this%x_var(1) - 2.0 - 5.0)**2
      chi2(2) = 200.0*this%x_var(1) - 140.0

      CALL profiler_set_stop_time('gradient_descent_test_chi2',                &
     &                            start_time)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Gradient descent unit test function.
!>
!>  This runs the associated unit tests and returns the result.
!>
!>  @returns True if the tests pass and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION gradient_descent_test()

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: gradient_descent_test

!  local variables
      REAL (rprec)                                 :: start_time
      REAL (rprec)                                 :: residule
      CLASS (gradient_descent_test_class), POINTER :: context

!  Start of executable code
      start_time = profiler_get_start_time()

!  Test to make sure the vertices begin in an unallocated state.
      gradient_descent_test = check(.false., ASSOCIATED(context),              &
     &                              1, 'ASSOCIATED')
      IF (.not.gradient_descent_test) THEN
         RETURN
      END IF

      context => gradient_descent_test_class()
!  Test to make sure the vertices begin in an unallocated state.
      gradient_descent_test = check(.true., ASSOCIATED(context),              &
     &                              2, 'ASSOCIATED')
      IF (.not.gradient_descent_test) THEN
         RETURN
      END IF

      context%x_var(1) = 0.0
      context%step = 0.001
      residule = context%minimize()
!  Test to make sure the vertices begin in an unallocated state.
      gradient_descent_test = check(.true., residule .le. 1.0E-20,             &
     &                              1, 'residule')
      IF (.not.gradient_descent_test) THEN
         RETURN
      END IF

      gradient_descent_test = check(context%x_var(1), 0.7_rprec,               &
     &                              1, 'value')
      IF (.not.gradient_descent_test) THEN
         RETURN
      END IF

      DEALLOCATE(context)

      CALL profiler_set_stop_time('gradient_descent_test',                     &
     &                            start_time)

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
         WRITE(*,*) 'gradient_descent.f: ', name, ' test', testNum,            &
     &              'failed.'
         WRITE(*,*) 'Expected', expected, 'Received', received
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
         WRITE(*,*) 'gradient_descent.f: ', name, ' test', testNum,            &
     &              'failed.'
         WRITE(*,*) 'Expected', expected, 'Received', received
      END IF

      END FUNCTION

      END MODULE

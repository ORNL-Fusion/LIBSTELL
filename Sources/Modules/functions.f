!*******************************************************************************
!>  @file File functions.f
!>  @brief Contains module @ref functions.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  This module containes functions used by the profiles.
!*******************************************************************************
      MODULE functions
      USE stel_kinds

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      PUBLIC :: two_power, two_power_gs, function_test
      PRIVATE :: check

      CONTAINS

!*******************************************************************************
!  PUBLIC FUNCTIONS AND SUBROUTINES
!*******************************************************************************
!*******************************************************************************
!>  @brief A two power profile.
!>
!>  A two power profiles is defined using 3 parameters.
!>
!>    b(0) * (1.0 - x ** b(1)) ** b(2)                                       (1)
!>
!>  b(0) controls the profile amplitude.
!>  b(1) controls the profile width.
!>  b(2) controls the profile tail.
!>
!>  @param[in] x Point to evaluate function from.
!>  @param[in] b Array of profile parameters.
!>  @returns The function evaluated at x.
!*******************************************************************************
      REAL(rprec) PURE FUNCTION two_power(x, b)

      IMPLICIT NONE

!  Declare Arguments
      REAL(rprec), INTENT(in)                  :: x
      REAL(rprec), DIMENSION(0:20), INTENT(in) :: b

!  Start of executable code
      two_power = b(0)*((1.0 - x**b(1))**b(2))

      END FUNCTION

!*******************************************************************************
!>  @brief A two power profile with added guassian peaks.
!>
!>  The power profile is defined as @ref functions::two_power. Then up to 9
!>  additional guassians peaks can be added.
!>
!>    two_power(x)*(1 + Sum[b(i)*Exp(-(x - b(i+1))/b(i+2)) ** 2])            (1)
!>
!>  Parameters beyond the first three control the additional guassians.
!>  b(i) controls the amplitude.
!>  b(i+1) controls the mean.
!>  b(i+2) controls the width.
!>
!>  @param[in] x Point to evaluate function from.
!>  @param[in] b Array of profile parameters.
!>  @returns The function evaluated at x.
!*******************************************************************************
      REAL(rprec) PURE FUNCTION two_power_gs(x, b)

      IMPLICIT NONE

!  Declare Arguments
      REAL(rprec), INTENT(in)                  :: x
      REAL(rprec), DIMENSION(0:20), INTENT(in) :: b

!  local variables
      INTEGER                                  :: i

!  Start of executable code
      two_power_gs = 1.0
      DO i = 3, 18, 3
         two_power_gs = two_power_gs                                           &
     &                + b(i)*exp(-((x - b(i+1))/b(i+2))**2.0)
      END DO
      two_power_gs = two_power_gs*two_power(x,b)

      END FUNCTION

!*******************************************************************************
!>  @brief A profile with a single guassian peak.
!>
!>  This defines a profile with a single guassian peak.
!>
!>    b(0)*Exp(-(x - b(1))/b(2)) ** 2                                        (1)
!>
!>  b(0) controls the amplitude.
!>  b(1) controls the mean.
!>  b(2) controls the width.
!>
!>  @param[in] x Point to evaluate function from.
!>  @param[in] b Array of profile parameters.
!>  @returns The function evaluated at x.
!*******************************************************************************
      REAL(rprec) PURE FUNCTION guassian(x, b)

      IMPLICIT NONE

!  Declare Arguments
      REAL(rprec), INTENT(in)                  :: x
      REAL(rprec), DIMENSION(0:20), INTENT(in) :: b

!  Start of executable code
      guassian = b(0)*exp(-((x - b(1))/b(2))**2.0)

      END FUNCTION

!*******************************************************************************
!>  @brief A profile with a integrated single guassian peak.
!>
!>  This defines a profile with a with an integarted single guassian peak.
!>
!>    Sqrt(2Pi)/2*b(0)*b(1)*(Erf(b(1)/b(2)) - Erf((b(1) - x)/b(2)))          (1)
!>
!>  b(0) controls the amplitude.
!>  b(1) controls the mean.
!>  b(2) controls the width.
!>
!>  @note Fortran 90 does not have a built in Error function must use the C
!>        function instead.
!>
!>  @param[in] x Point to evaluate function from.
!>  @param[in] b Array of profile parameters.
!>  @returns The function evaluated at x.
!*******************************************************************************
      REAL(rprec) PURE FUNCTION guassian_int(x, b)
      USE stel_constants, only : pi

      IMPLICIT NONE

!  Declare Arguments
      REAL(rprec), INTENT(in)                  :: x
      REAL(rprec), DIMENSION(0:20), INTENT(in) :: b

!  local parameters
      REAL(rprec), PARAMETER                   :: sqrpi = 0.5*SQRT(pi)

!  Start of executable code

      guassian_int = sqrpi*b(0)*b(1)*(ERF(b(1)/b(2)) -                         &
     &                                ERF((b(1) - x)/b(2)))

      END FUNCTION

!*******************************************************************************
!  UNIT TESTS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Profile unit test function.
!>
!>  This runs the associated unit tests and returns the result.
!>
!>  @returns True if the tests pass and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION function_test()

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                      :: function_test

!  local variables
      REAL(rprec)                  :: result
      REAL(rprec), DIMENSION(0:20) :: b(0:20) = 0

!  Start of executable code
! Test two_power function for x = 0, b = {1,10,2} is 1
      b(0:2) = (/ 1.0d+0, 10.0d+0, 2.0d+0 /)
      result = two_power(0.0d+0, b)
      function_test = check(1.0d+0,result,1,"two_power()")
      IF (.not.function_test) RETURN

! Test two_power function for x = 1, b = {1,10,2} is 0
      result = two_power(1.0d+0, b)
      function_test = check(0.0d+0,result,2,"two_power()")
      IF (.not.function_test) RETURN

! Test two_power function for x = 0.5, b = {1,1,1} is 0.5
      b(0:2) = (/ 1.0d+0, 1.0d+0, 1.0d+0 /)
      result = two_power(0.5d+0, b)
      function_test = check(0.5d+0,result,3,"two_power()")
      IF (.not.function_test) RETURN

! Test two_power function for x = 0.5, b = {1,1,2} is 0.25
      b(0:2) = (/ 1.0d+0, 1.0d+0, 2.0d+0 /)
      result = two_power(0.5d+0, b)
      function_test = check(0.25d+0,result,4,"two_power()")
      IF (.not.function_test) RETURN

! Test two_power_gs function for x = 0.4, b = {1,1,1,0,0,1} is twopower(x,b)
      b(0:5) = (/ 1.0d+0, 1.0d+0, 1.0d+0, 0.0d+0, 0.0d+0, 1.0d+0 /)
      result = two_power_gs(0.4d+0, b)
      function_test = check(two_power(0.4d+0, b),                              &
     &                      result,1,"two_power_gs")
      IF (.not.function_test) RETURN

! Test two_power_gs function for x = 0.8, b = {1,1,0,1,0.8,0.1} is 2
      b(0:5) = (/ 1.0d+0, 1.0d+0, 0.0d+0, 1.0d+0, 0.8d+0, 0.1d+0 /)
      result = two_power_gs(0.8d+0, b)
      function_test = check(2.0d+0,result,1,"two_power_gs")
      IF (.not.function_test) RETURN

      END FUNCTION

!*******************************************************************************
!  CHECK FUNCTIONS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Check a value.
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
      FUNCTION check(expected, received, testNum, name)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: check
      REAL(rprec), INTENT(in)       :: expected, received
      INTEGER, INTENT(in)           :: testNum
      CHARACTER (LEN=*), INTENT(in) :: name

!  Start of executable code
      check = expected .eq. received
      IF (.not.check) THEN
         write(*,*) "functions.f: ", name, " test", testNum,                   &
     &              "failed."
         write(*,*) "Expected", expected, "Received", received
      END IF

      END FUNCTION

      END MODULE

!*******************************************************************************
!>  @file unit_test_runner.f
!>  @brief Utility to run the LIBSTELL unit tests.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  These tests run module tests embeded within specific modules.
!*******************************************************************************
      PROGRAM LIBSTELL_UNIT_TEST_RUNNER
      USE line_segment
      USE coordinate_utilities
      USE functions
      USE integration_path

      IMPLICIT NONE

!  local variables
      LOGICAL :: tests_failed

!  Start of executable code
      tests_failed = .false.

! Line segment tests
      IF (line_seg_test()) THEN
         WRITE (*,*) 'Line Segment Test Passed'
      ELSE
         tests_failed = .true.
      END IF

! Coordinate tests
      IF (cood_utils_test()) THEN
         WRITE (*,*) 'Coord Utilites Test Passed'
      ELSE
         tests_failed = .true.
      END IF

! Profile function tests
      IF (function_test()) THEN
         WRITE (*,*) 'Functions Test Passed'
      ELSE
         tests_failed = .true.
      END IF

! Integration path tests
      IF (path_test()) THEN
         WRITE (*,*) 'Integration Path Test Passed'
      ELSE
         tests_failed = .true.
      END IF

!  The CTest unit test frame work chooses a pass fail based on the exit code of
!  the program.
      IF (tests_failed) THEN
         CALL EXIT(1)
      END IF

      END PROGRAM

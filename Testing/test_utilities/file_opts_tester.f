!*******************************************************************************
!>  @file file_opts_tester.f
!>  @brief Utility to run the file opts operations.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  These tests run subroutine form the @ref file_opts module.
!*******************************************************************************
     PROGRAM file_opts_tester
     USE file_opts

     IMPLICIT NONE

!  Local Variables
     INTEGER :: local_error
     INTEGER :: test_error

!  Start of executable code.
     local_error = 0
     test_error = 0

#  A second python utility will check that the copy was made.
     CALL copy_file('test_file', 'test_file_copy', local_error)
     IF (local_error .ne. 0) THEN
        WRITE (*,*) 'Copy file test failed.'
        test_error = local_error
     END IF

     CALL create_directory('test_directory', local_error)
     IF (local_error .ne. 0) THEN
        WRITE (*,*) 'Create directory test failed.'
        test_error = local_error
     END IF

#  The delete directory should check if the file was successfully moved.
     CALL move_file('test_file', 'test_directory/test_file_move', local_error)
     IF (local_error .ne. 0) THEN
        WRITE (*,*) 'Move file test failed.'
        test_error = local_error
     END IF

     CALL change_directory('test_directory', local_error)
     IF (local_error .ne. 0) THEN
        WRITE (*,*) 'Change directory test failed.'
        test_error = local_error
     END IF

     CALL delete_file('test_file_move', local_error)
     IF (local_error .ne. 0) THEN
        WRITE (*,*) 'Delete file test failed.'
        test_error = local_error
     END IF

     CALL EXIT(test_error)

     END PROGRAM

      SUBROUTINE multiprocess(numprocess, maxprocess, wrapperfcn, fcn)
      USE system_mod
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER numprocess, maxprocess
      REAL, EXTERNAL :: fcn
      EXTERNAL wrapperfcn
#if !defined(MPI_OPT)
!C-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, status, iretpid, ierror
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      EXTERNAL myfork
!-----------------------------------------------
 
 
!   This Program was written by S. P. Hirshman (7/1/99), under contract
!   the U.S. DOE, and should not be used without explicit permission.
!
!   Input Variable Names
!   maxprocess:   user defined constant, maximum number of processes that this routine will
!                 request. It should be about MAX[NCPU/(1+NSPAWN), 1], WHERE NCPU is the number of
!                 machine cpus, and NSPAWN is the maximum number of processes spawned by the
!                 CALL to loc_function (for stellopt code, NSPAWN = 1)
!   NumProcess:   the TOTAL number of processes to be launched in parallel (IF possible). If
!                 this exceeds the number of available processors (max_process, THEN the next
!                 available processor scheduled by the operating system
!                 returns only after all processes are completed.
!   WrapperFcn:   Fortran SUBROUTINE that performs the desired task.
!                 It takes for arguments
!                    (1) the INDEX of the process to execute, WHERE 1 <= index <= NumProcess
!                    (2) the SUBROUTINE (Fcn) which is called ito evaluate specific information
!                        (such as FUNCTIONal minima, etc.).
!
!   Calling convention from a Fortran PROGRAM:
!
!   CALL MultiProcess(nprocess, maxprocess, Wrapper_Subroutine, Worker_Subroutine)
 
      iretpid = 0; ierror = 0

      IF (maxprocess .GT. 1) THEN

         WRITE(6,*)
         WRITE(6,*) ' Begin multi-processing: request ', numprocess,    &
            ' processes distributed among ', maxprocess, ' processors'
         CALL flush(6)

         DO i = 1, numprocess
            CALL myfork (i, maxprocess, wrapperfcn, fcn)
         END DO
 
!        Wait for ALL processes to finish...
         DO i = 1, numprocess
            CALL pxfwait (status, iretpid, ierror)
         END DO
 
      ELSE

         DO i = 1, numprocess
            CALL wrapperfcn (i, fcn)
         END DO
      END IF
#endif
      END SUBROUTINE multiprocess

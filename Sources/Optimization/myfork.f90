      SUBROUTINE myfork(i, maxprocess, wrapper, fcn)
      USE system_mod
      IMPLICIT NONE      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER i, maxprocess
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      EXTERNAL wrapper, fcn
#if !defined(MPI_OPT)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: pid, status, iretpid, ierror, werror
      INTEGER, SAVE :: nprocess = 0
!-----------------------------------------------
 
      IF (i .eq. 1) nprocess = 0
      ierror = -1
 
!     Child process: limit number to max_process to avoid potential system hang-up
 
      DO WHILE(ierror .ne. 0)
 
         IF (nprocess .lt. maxprocess) CALL pxffork (pid, ierror)
 
         IF (ierror .ne. 0) THEN
!           wait for next available processor
            CALL pxfwait (status, iretpid, werror)
!           IF (status.gt.0 .and. nprocess.ge.1) THEN
            IF (nprocess .ge. 1) nprocess = nprocess - 1
!           ELSE
!              nprocess = 0
!           END IF
         END IF
      END DO
 
      IF (pid .eq. 0) THEN
         CALL wrapper (i, fcn)
#if defined(CRAY)
         CALL EXIT(1)
#else
         STOP
#endif
      END IF
 
      nprocess = nprocess + 1
#endif
      END SUBROUTINE myfork

      SUBROUTINE de_mpi(np, fcn, funcval)
      USE de_mod
      USE mpi_params
      USE mpi_inc
      IMPLICIT NONE
      INTEGER :: np
      REAL(rprec) :: funcval(np)
      EXTERNAL fcn

#if defined(MPI_OPT)
      INTEGER :: status(MPI_STATUS_size)                     !mpi stuff
      INTEGER :: i,  j, iflag                                !mpi stuff
      INTEGER :: numsent, sender, ierr                       !mpi stuff
      INTEGER :: anstype, column                             !mpi stuff
      REAL(rprec), DIMENSION(n_free) :: x
      REAL(rprec), DIMENSION(nopt) :: fvec

!******************************************
!
!  mpi : set barrier so ALL processors get here before starting
!
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)                 !mpi stuff

!******************************************
!
!     ****Master node controls this portion of the code****
!
      IF (myid .eq. master) THEN
         numsent = 0    !numsent is a counter used to track how many
                        !jobs have been sent to workers
c
c     SEND forward difference displacements from master to each
c           worker process. Tag with these with the column number.
c
         DO j = 1,MIN(numprocs-1,np)
            x(:) = ui_XC(j,:)
            CALL MPI_SEND(x, n_free, MPI_REAL8, j,
     1                  j, MPI_COMM_WORLD, ierr)
            IF (ierr .ne. 0) STOP 'MPI_SEND error(1) in de_mpi'
            numsent = numsent+1
         END DO          !j = 1,MIN(numprocs-1,n)
c
c      Looping through the columns, collect answers from the workers.
c      As answers are received, new uncalculated columns are sent
c      out to these same workers.
c
         DO j = 1,np
            CALL MPI_RECV(fvec, nopt, MPI_REAL8,
     1           MPI_any_SOURCE, MPI_any_TAG,
     2           MPI_COMM_WORLD, status, ierr)
            IF (ierr .ne. 0) STOP 'MPI_RECV error(1) in de_mpi'
            sender     = status(MPI_SOURCE)
            anstype    = status(MPI_TAG)       ! column is tag value
            IF (anstype .gt. np) STOP 'ANSTYPE > NP IN de_mpi'

            funcval(anstype) = SUM(fvec(:nopt)**2)
c           WRITE(6,'(a,1pe10.3,a,i3,a,i3)')' FUNCVAL = ',
c    1       funcval(anstype),
c    2      ' for iteration ', anstype+nfev,' processor = ', sender

c
c           If more columns are left, then send another column to the worker(sender)
c           that just sent in an answer
c
            IF (numsent .lt. np) THEN
               numsent = numsent+1
               x(:) = ui_XC(numsent,:)

               CALL MPI_SEND(x, n_free, MPI_REAL8,
     1                       sender, numsent, MPI_COMM_WORLD, ierr)
               IF (ierr .ne. 0) STOP 'MPI_SEND error(2) in de_mpi'

            ELSE                ! Tell worker that there is no more work to DO

               CALL MPI_SEND(MPI_BOTTOM, 0, MPI_REAL8,
     1                       sender, 0, MPI_COMM_WORLD, ierr)
               IF (ierr .ne. 0) STOP 'MPI_end error(3) in de_mpi'
            ENDIF      ! IF( myid .eq. master ) THEN
         END DO     ! DO j = 1,n
c
c     ****Worker portion of the code****
c        Skip this when processor id exceeds work to be done
c
      ELSE IF (myid .le. np) THEN               !!IF( myid .ne. master )
c
c        Otherwise accept the next available column, check the tag,
c        and IF the tag is non-zero CALL SUBROUTINE fcn.
c        If the tag is zero, there are no more columns
c        and worker skips to the END.
c
 90      CALL MPI_RECV(x, n_free, MPI_REAL8, master,
     1                 MPI_any_TAG, MPI_COMM_WORLD, status, ierr)
         IF (ierr .ne. 0) STOP 'MPI_RECV error(2) in de_mpi'

         column = status(MPI_TAG)                !!ID of pseudo-processor issuing this message
         IF (column .eq. 0) THEN
            GOTO 200
         ELSE
            iflag = column
c           CALL the chisq fcn for the portion of displacement vector which
c           was just received. Note that WA stores the local fvec_min array

            CALL fcn(nopt, n_free, x, fvec, iflag, nfev)
            IF (iflag.ne.0) GOTO 300
c
c           Send this function evaluation back to the master process tagged
c           with the column number so the master knows where to put it
c
            CALL MPI_SEND(fvec, nopt, MPI_REAL8, master,
     1                    column, MPI_COMM_WORLD, ierr)
            IF (ierr .ne. 0) STOP 'MPI_SEND error(4) in de_mpi'
            GOTO 90    !Return to 90 and check IF master process has sent any more jobs
         END IF
 200     CONTINUE
      ENDIF       ! IF( myid .ne. master )

!
!     Broadcast the funcval array to all processors FROM master
!
      CALL MPI_BCAST(funcval, np, MPI_REAL8, master,
     1     MPI_COMM_WORLD, ierr)
      IF (ierr .ne. 0) GOTO 100

      RETURN

 100  CONTINUE
      PRINT *,' MPI_BCAST error in de_mpi: IERR=', ierr

      RETURN

 300  CONTINUE
      PRINT *,' IFLAG = ', iflag, ' in de_mpi CALL to fcn'
      STOP
#endif
      END SUBROUTINE de_mpi

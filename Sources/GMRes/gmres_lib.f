      MODULE GMRES_LIB
      USE stel_kinds, ONLY: dp
      USE stel_constants, ONLY: one, zero
      USE mpi_inc
      IMPLICIT NONE

      TYPE GMRES_INFO
         INTEGER  :: m, mblk_size, icntl(9), info(3)
         INTEGER  :: endglobrow, startglobrow, iam, nprocs,
     1               ngmres_type=2
#if defined(MPI_OPT)
         INTEGER  :: my_comm=MPI_COMM_WORLD, 
     1               my_comm_world=MPI_COMM_WORLD
#else
         INTEGER  :: my_comm=0,
     1               my_comm_world=0
#endif
         INTEGER, POINTER  :: rcounts(:), disp(:)
         LOGICAL  :: lactive = .TRUE.
         REAL(dp) :: cntl(5), ftol
         LOGICAL  :: lverbose = .TRUE.
         LOGICAL  :: l_nonlinear = .TRUE.
      END TYPE GMRES_INFO

      INTEGER, PARAMETER :: finish=0, matveci=1, precondLeft=2,
     1                      precondRight=3, dotProd=4, peek=5

      CONTAINS

      SUBROUTINE gmres_ser (n, gi, yAx, apply_precond,
     &                      getnlforce, x0, b)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: n
      TYPE(GMRES_INFO)        :: gi
      REAL(dp), INTENT(IN)    :: b(n)
      REAL(dp), INTENT(INOUT) :: x0(n)
!      REAL(dp) :: cntl(5)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: revcom, colx, coly, colz, nbscal, m, iam
      INTEGER :: irc(5)
      INTEGER :: nout, lwork, icount, jcount, kout, kprint
      REAL(dp)  :: rinfo(2), fsq_nl=-1, fsq_lin, fsq_min,
     1             delfsq, fsq_last, bnorm, xmod, xmax,
     2             fsq_min_lin, gsum
      REAL(dp), TARGET, ALLOCATABLE :: work(:)
      REAL(dp), ALLOCATABLE         :: xmin(:)
      REAL(dp), POINTER :: sx(:), sy(:), sz(:)
	LOGICAL   :: lprint, lstore, l_nonlinear
!-----------------------------------------------
      EXTERNAL yAx, apply_precond, getnlforce
!-----------------------------------------------
!
!     EASY-TO-USE WRAPPER FOR GMRES DRIVER CALL
!
!     X0: on input, initial guess if icntl(6) == 1
!         on output, solution of Ax = b
!         NOTE: it is not overwritten until the END of this routine
!
      bnorm = SQRT(SUM(b*b))
      IF (bnorm .EQ. 0) RETURN

!      PRINT *,'SERIAL BNORM: ', bnorm

      fsq_min = 1; fsq_min_lin = 1
      delfsq = 1
      jcount = 0; icount = 0
      kout = 0; kprint = 0
	l_nonlinear = gi%l_nonlinear

      m = gi%m
      iam = gi%iam
	  lprint = (iam.EQ.0 .AND. gi%lverbose)
	
!      lwork = m**2 + m*(n+6) + 5*n + 1
      lwork = m**2 + m*(n+6) + 6*n + 1    !Additional space for peek revcom (5 -> 6)

      ALLOCATE (work(lwork), stat=nout)
      IF (nout .NE. 0) STOP 'Allocation error in gmres!'
      work = 0
      IF (gi%icntl(6) .EQ. 1) work(1:n) = x0/bnorm
      work(n+1:2*n) = b(1:n)/bnorm

      IF (lprint) PRINT 900
 900  FORMAT(/,1x,'GMRES (SERIAL) CONVERGENCE SUMMARY',/,1x,15('-'))     

      !****************************************
      !* Reverse communication implementation
      !****************************************

 10   CONTINUE
      CALL drive_dgmres(n, n, m, lwork, work, irc,
     &                  gi%icntl, gi%cntl, gi%info, rinfo)
      revcom = irc(1)
      colx   = irc(2)
      coly   = irc(3)
      colz   = irc(4)
      nbscal = irc(5)
      sx => work(colx:);  sy => work(coly:);  sz => work(colz:)

      IF (revcom .EQ. matveci) THEN
! perform the matrix vector product work(colz) <-- A * work(colx)
         CALL yAx (sx, sz, n) 
!Debug
!         gsum = SUM(work(colz:colz+n-1)**2)
!         IF (lprint) PRINT *,'CALL MATVECI, |Ax| = ', SQRT(gsum)
!End Debug
         GOTO 10

      ELSE IF (revcom.eq.precondLeft) THEN
! perform the left preconditioning
!          IF (lprint) PRINT *,'CALL PRECONDL'        
!         work(colz) <-- M^{-1} * work(colx)
!         WRITE(10000,*) "left_dcopy"; CALL FLUSH(10000)
         CALL dcopy(n,sx,1,sz,1)
!         WRITE(10000,*) "precondLeft"; CALL FLUSH(10000)
         CALL apply_precond(sz)
         GOTO 10

      ELSE IF (revcom .EQ. precondRight) THEN
! perform the right preconditioning
!          IF (lprint) PRINT *,'CALL PRECONDR'        
!         WRITE(10000,*) "right_dcopy"; CALL FLUSH(10000)
         CALL dcopy(n,sx,1,sz,1)
!         WRITE(10000,*) "precondRight"; CALL FLUSH(10000)
         CALL apply_precond(sz)
         GOTO 10

      ELSE IF (revcom .EQ. dotProd) THEN
!      perform the scalar product
!      work(colz) <-- work(colx) work(coly)
!
!         CALL dgemv('C',n,nbscal,ONE,sx,n,sy,1,ZERO,sz,1)
!        WRITE(10000,*) "dotProd/truncate"; CALL FLUSH(10000)
        DO nout = 0, nbscal-1
          work(colz+nout) = SUM(work(colx:colx+n-1)*work(coly:coly+n-1))
          CALL Truncate(work(colz+nout), 15)
          colx = colx+n
        END DO

! Debug
!        IF (lprint) THEN
!           gsum = SUM(work(colz:colz+nbscal-1)**2)
!           PRINT 200,' CALL DOTPROD: nbscal = ', nbscal, 
!     1               ' |WORK|: ', SQRT(gsum)
! 200    FORMAT(a,i4,a,1p,e14.6)
!        END IF
! End Debug
        GOTO 10

      ELSE IF (revcom .EQ. peek) THEN
        fsq_last = fsq_nl                                                !need for delfsq criteria
        IF (l_nonlinear) CALL GetNLForce(work(colx), fsq_nl, bnorm)      !get nonlinear force

        fsq_lin = (rinfo(1)*bnorm)**2
        icount = icount+1
        lStore = (fsq_nl.LT.fsq_min .OR. icount.EQ.1)
        fsq_min = MIN(fsq_min, fsq_nl)
        fsq_min_lin = MIN(fsq_min_lin, fsq_lin)
        IF (fsq_lin .LT. 1.E-30_dp) gi%icntl(7)=gi%info(1)

        IF (lprint .AND. icount.EQ.1) THEN
	     IF (l_nonlinear) THEN
		    PRINT 905
	     ELSE
              PRINT 906
	     END IF
	  END IF

        kout = gi%info(1)
        IF (icount.EQ.1 .OR. MOD(icount,10).EQ.0) THEN
           kprint = kout
           xmod = bnorm*SQRT(SUM(work(colx:colx+n-1)**2))
           xmax = bnorm*MAXVAL(ABS(work(colx:colx+n-1)))
           IF (icount .EQ. 1) THEN
              fsq_min = fsq_nl; fsq_min_lin = fsq_lin
           END IF
           IF (lprint) THEN
		    IF (l_nonlinear) THEN
			   PRINT 910, kout, fsq_nl, xmod, xmax, fsq_lin
              ELSE
			   PRINT 911, kout, xmod, xmax, fsq_lin
              END IF
	     END IF
        END IF

        IF (gi%ngmres_type .EQ. 1) GOTO 10                               !OLDSTYLE:IGNORE LOGIC BELOW

        delfsq = (fsq_last-fsq_nl)/fsq_min                              
        IF (delfsq .LT. 0.05_dp) THEN
          jcount = jcount+1
        ELSE
          jcount = 0
        END IF

! STOPPING CRITERIA (reset max iterations to current iteration)         
        IF (fsq_nl.GT.(3*fsq_min) .OR. jcount.GT.3 
     &            .OR. fsq_nl.LE.gi%ftol) gi%icntl(7) = gi%info(1)
        IF (lStore) THEN
          IF (.NOT.ALLOCATED(xmin)) ALLOCATE (xmin(n))
          xmin = work(colx:colx+n-1)
        END IF
        GOTO 10

      ENDIF
 905  FORMAT(1x,'ITER',7x,'FSQ_NL',10x,'||X||',9x,'MAX|X|',9x,'FSQ_ARN')
 906  FORMAT(1x,'ITER',7x,'||X||',11x,'MAX|X|',9x,'FSQ_ARN')
 910  FORMAT(i5, 4(3x,1pe12.3))
 911  FORMAT(i5, 3(3x,1pe12.3))
  
 100   CONTINUE

      IF (ALLOCATED (xmin)) THEN
        x0(1:n) = xmin(1:n)
        DEALLOCATE (xmin)
      ELSE
        x0(1:n) = work(1:n)
        fsq_min = fsq_nl
	  IF (.NOT.l_Nonlinear) fsq_min = fsq_lin
      END IF

      gi%ftol = fsq_min
      x0 = bnorm*x0
      IF (lprint .AND. kprint.NE.kout) THEN
         xmod = SQRT(SUM(x0*x0));  xmax = MAXVAL(ABS(x0))
	   IF (l_NonLinear) THEN
            PRINT 910, kout, fsq_min, xmod, xmax, fsq_min_lin
         ELSE
	      PRINT 911, kout, xmod, xmax, fsq_lin
	   END IF
      END IF

      DEALLOCATE (work)

!Sanity check
!      CALL GetNLForce(x0, fsq_nl, one)
!      IF (iam .EQ. 0) PRINT *,'FSQ_MIN: ', fsq_min,' FSQ_NL: ', fsq_nl

      END SUBROUTINE gmres_ser

      SUBROUTINE gmres_par (n, gi, yAx, apply_precond, GetNLForce, 
     &                      x0, b)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)     :: n
      TYPE(GMRES_INFO)        :: gi
      REAL(dp), INTENT(IN)    :: b(n)
      REAL(dp), INTENT(INOUT) :: x0(n)
      EXTERNAL apply_precond, yAx, GetNLForce
#if defined(MPI_OPT)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: revcom, colx, coly, colz, nbscal, m, iam
      INTEGER  :: irc(5), jcount, icount
      INTEGER  :: nout, lwork, kout, kprint, MPI_ERR, 
     &            MY_COMM_WORLD, MY_COMM
      INTEGER, ALLOCATABLE :: itest(:)
      REAL(dp) :: rinfo(2), fsq_nl, fsq_min, fsq_last, delfsq, fsq_lin, 
     &            fsq_min_lin
      REAL(dp), ALLOCATABLE :: work(:), xmin(:)
      INTEGER  :: nloc, myrowstart, myrowend
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: aux, tmpbuf
      REAL(dp) :: skston, skstoff, xmax, xmod, bnorm, lsum, gsum
      LOGICAL, PARAMETER :: lfast = .true.
      LOGICAL  :: lactive = .TRUE., lprint, lStore, l_NonLinear
!-----------------------------------------------
!
!     EASY-TO-USE WRAPPER FOR GMRES PARALLEL DRIVER CALL
!
!     X0: on input, initial guess if icntl(6) == 1
!         on output, solution of Ax = b
!         NOTE: it is not overwritten UNTIL the end of this routine
      fsq_nl = -1
      lactive = gi%lactive
      MY_COMM_WORLD = gi%MY_COMM_WORLD
      MY_COMM       = gi%MY_COMM
      l_nonlinear = gi%l_nonlinear

      fsq_min = 1
      fsq_min_lin = 1
      delfsq = 1
      jcount = 0
      icount = 0
      kout = 0
      kprint = 0

      nloc = (gi%endglobrow - gi%startglobrow + 1)*gi%mblk_size
      myrowstart = (gi%startglobrow - 1)*gi%mblk_size + 1
      myrowend = myrowstart + nloc - 1

      IF (.NOT.lactive) THEN
         nloc = 1
         myrowstart = 1
         myrowend = myrowstart
      END IF

      ALLOCATE(tmpbuf(n), aux(n), itest(gi%nprocs), stat=nout)
      IF (nout .NE. 0) THEN
         STOP 'Allocation error in gmres_fun!'
      END IF

      m = gi%m
      iam = gi%iam
      lprint = (iam .EQ. 0 .AND. gi%lverbose .AND. lactive)
      lwork = m**2 + m*(nloc + 6) + 6*nloc + 1    !Additional space for peek revcom (5 -> 6)
      ALLOCATE(work(lwork), stat=nout)
      IF (nout .NE. 0) THEN
         STOP 'Allocation error in gmres_fun!'
      END IF
      work = 0

!      PRINT *,'In gmres_par, iam: ',iam,' lactive: ',lactive

      LACTIVE0: IF (lactive) THEN
      lsum = SUM(b(myrowstart:myrowend)**2)
      CALL MPI_ALLREDUCE(lsum, gsum, 1, MPI_REAL8, MPI_SUM, MY_COMM,
     &                   MPI_ERR)
      bnorm = SQRT(gsum)
      IF (bnorm .EQ. 0) THEN
         RETURN
      END IF

      IF (gi%icntl(6) .EQ. 1) THEN
         work(1:nloc) = x0(myrowstart:myrowend)/bnorm
      END IF

      work(nloc + 1:2*nloc) = b(myrowstart:myrowend)/bnorm

      IF (lprint) THEN
         PRINT 900
      END IF
900   FORMAT(/,1x,'GMRES (PARALLEL) CONVERGENCE SUMMARY',/,1x,15('-'))

      END IF LACTIVE0

      CALL MPI_BCAST(bnorm, 1, MPI_REAL8, 0, MY_COMM_WORLD,
     &               MPI_ERR)
 
!****************************************
!* Reverse communication implementation
!****************************************

 10   CONTINUE

      IF (lactive) THEN
         CALL drive_dgmres(n, nloc, m, lwork, work, irc,
     &                     gi%icntl, gi%cntl, gi%info, rinfo)
      ELSE
         irc = 1
      END IF

!BCast to all processors in the local world group
      CALL MPI_BCAST(irc(1), 1, MPI_INTEGER, 0, MY_COMM_WORLD,
     &               MPI_ERR)
      CALL MPI_BCAST(rinfo(1), 1, MPI_REAL8, 0, MY_COMM_WORLD,
     &               MPI_ERR)

      revcom = irc(1)
      colx   = irc(2)
      coly   = irc(3)
      colz   = irc(4)
      nbscal = irc(5)

      REVCOM_LOOP: IF (revcom .EQ. matveci) THEN
! perform the matrix vector product work(colz) <-- A * work(colx)
!        WRITE(10000+iam,*) "matvec"; CALL FLUSH(10000+iam)
         CALL yAx (work(colx), work(colz), nloc)
         GOTO 10

      ELSE IF (revcom .EQ. precondLeft) THEN
!        IF (lprint) PRINT *,'CALL PRECONL'
! perform the left preconditioning work(colz) <-- M^{-1} * work(colx)
!        WRITE(10000+iam,*) "left_dcopy"; CALL FLUSH(10000+iam)
         IF (lactive) THEN
            CALL dcopy(nloc, work(colx), 1, work(colz), 1)
         END IF
         GOTO 10

      ELSE IF (revcom .EQ. precondRight) THEN
         CALL dcopy(nloc,work(colx),1,work(colz),1)

         IF (lactive) THEN
            aux(myrowstart:myrowend) = work(colz:colz + nloc - 1)
!        WRITE(10000+iam,*) "precondRight"; CALL FLUSH(10000+iam)
            CALL apply_precond(aux)
            work(colz:colz + nloc - 1) = aux(myrowstart:myrowend)
         END IF
         GOTO 10

      ELSE IF (revcom .EQ. dotProd) THEN
! perform the scalar product (uses nbscal columns of A starting at work(colx))
        ! work(colz) <-- work(colx) work(coly)

         LACTIVE_DP: IF (lactive) THEN
        !MAKE SURE nbscal is the same on all processors -
            CALL MPI_ALLGATHER(nbscal, 1, MPI_INTEGER, itest, 1,
     &                         MPI_INTEGER, MY_COMM, MPI_ERR)
            IF (lprint .AND. ANY(itest(1:gi%nprocs) .NE. nbscal)) THEN
               PRINT *,'itest: ',itest(1:gi%nprocs)
               STOP 'nbscal not same on all procs!'
            END IF

!       THIS IS FASTER THAN ALLGATHERV LOOP, BUT MAY LEAD TO SLIGHTLY DIFFERENT CONVERGENCE
!       SEQUENCES FOR DIFFERENT # PROCESSORS. THE DGEMV CALL IS SLIGHTLY SLOWER THAN THE LOOP
            IF (lFast) THEN

!!        CALL dgemv('C',nloc,nbscal,one,work(colx),nloc,work(coly),1,zero,aux,1)

               DO nout = 1, nbscal
                  aux(nout) = SUM(work(colx:colx + nloc - 1) *
     &                            work(coly:coly + nloc - 1))
                  colx = colx + nloc
               END DO

               CALL MPI_ALLREDUCE(aux, work(colz), nbscal, MPI_REAL8,
     &                            MPI_SUM, MY_COMM, MPI_ERR)

            ELSE

!       THE TIMING FOR THIS DOES NOT SCALE AS WELL WITH # PROCESSORS, BUT
!       LEADS TO A CONVERGENCE SEQUENCE THAT IS INDEPENDENT OF # PROCESSORS
               CALL MPI_GATHERV(work(coly), nloc, MPI_REAL8, aux,
     &                          gi%rcounts, gi%disp, MPI_REAL8, 0,
     &                          MY_COMM, MPI_ERR)

               DO nout = 0, nbscal - 1
                  CALL MPI_GATHERV(work(colx), nloc, MPI_REAL8, tmpbuf,
     &                             gi%rcounts, gi%disp, MPI_REAL8, 0,
     &                             MY_COMM, MPI_ERR)
                  IF (iam .EQ. 0) THEN
                     work(colz + nout) = SUM(tmpbuf*aux)              !DO FOR ALL PROCS IF ALLGATHER FORM USED
                  END IF
                  colx = colx+nloc
               END DO

               CALL MPI_BCAST(work(colz), nbscal, MPI_REAL8, 0, MY_COMM,
     &                        MPI_ERR)
            END IF

            DO nout = 0,nbscal-1
               CALL Truncate(work(colz+nout), 15)
            END DO
         END IF LACTIVE_DP

         GOTO 10

      ELSE IF (revcom .EQ. peek) THEN
         IF (lactive) THEN
            fsq_last = fsq_nl                                              !need for delfsq criteria
            aux(myrowstart:myrowend) = work(colx:colx + nloc - 1)
         END IF
!get nonlinear force
         IF (l_NonLinear) THEN
            CALL GetNLForce(aux, fsq_nl, bnorm)
         END IF

         fsq_lin = (rinfo(1)*bnorm)**2
         icount = icount + 1
         lStore = (fsq_nl .LT. fsq_min .OR. icount .EQ. 1)
         fsq_min = MIN(fsq_min, fsq_nl)
         fsq_min_lin = MIN(fsq_min_lin, fsq_lin)
           
         IF (fsq_lin .LT. 1.E-30_dp) THEN
            gi%icntl(7) = gi%info(1)
         END IF

         IF (lprint .AND. icount .EQ. 1) THEN
            IF (l_nonlinear) THEN
		       PRINT 905
            ELSE
               PRINT 906
	        END IF
         END IF
         kout = gi%info(1)

         LACTIVE1: IF (lactive) THEN
        
            IF (icount .EQ. 1 .OR. MOD(icount,10) .EQ. 0) THEN
              lsum = SUM(aux(myrowstart:myrowend)**2)
              CALL MPI_REDUCE(lsum, gsum, 1, MPI_REAL8, MPI_SUM, 0,
     &                        MY_COMM, MPI_ERR)
              xmod = bnorm*SQRT(gsum)
              lsum = MAXVAL(ABS(aux(myrowstart:myrowend)))
              CALL MPI_REDUCE(lsum, xmax, 1, MPI_REAL8, MPI_MAX, 0,
     &                        MY_COMM, MPI_ERR)
              xmax = bnorm*xmax
              kprint = kout
              IF (lprint) THEN
                 IF (l_nonlinear) THEN
                    PRINT 910, kout, fsq_nl, xmod, xmax, fsq_lin
                 ELSE
                    PRINT 911, kout, xmod, xmax, fsq_lin
                 END IF
              END IF
           END IF

        END IF LACTIVE1

        IF (gi%ngmres_type .EQ. 1) THEN
           GOTO 10  !OLDSTYLE:IGNORE LOGIC BELOW
        END IF

!       LACTIVE1: IF (lactive) THEN
        delfsq = (fsq_last - fsq_nl)/fsq_min
        IF (delfsq .LT. 0.05_dp) THEN                                     !PROGRESS COUNTER
          jcount = jcount + 1
        ELSE
          jcount = 0
        END IF

! STOPPING CRITERIA (reset max iterations to current iteration)         
         IF (fsq_nl .GT. (3*fsq_min) .OR.
     &       jcount .GT. 3           .OR.
     &       fsq_nl .LE. gi%ftol) THEN
            gi%icntl(7) = gi%info(1)
         END IF
         IF (lStore) THEN
            IF (.NOT.ALLOCATED(xmin)) THEN
               ALLOCATE (xmin(nloc))
            END IF
            xmin = work(colx:colx+nloc-1)
         END IF
         GOTO 10

      END IF REVCOM_LOOP

 905  FORMAT(1x,'ITER',7x,'FSQ_NL',10x,'||X||',9x,'MAX|X|',9x,'FSQ_ARN')
 906  FORMAT(1x,'ITER',7x,'||X||',11x,'MAX|X|',9x,'FSQ_ARN')
 910  FORMAT(i5, 4(3x,1pe12.3))
 911  FORMAT(i5, 3(3x,1pe12.3))

!*******************************
! end reverse loop: dump the solution to a file for debugging
!******************************
      GOTO 100

      LACTIVE2: IF (lprint .AND. lactive) THEN

         nout = 11
         OPEN(nout,FILE='sol_dTestgmres', STATUS='unknown')
         IF (gi%icntl(5) .EQ. 0) THEN
            WRITE(nout,*) 'Orthogonalisation : MGS'
         ELSE IF (gi%icntl(5) .eq. 1) THEN
            WRITE (nout,*) 'Orthogonalisation : IMGS'
         ELSE IF (gi%icntl(5) .eq. 2) THEN
            WRITE (nout,*) 'Orthogonalisation : CGS'
         ELSE IF (gi%icntl(5) .eq. 3) THEN
            WRITE (nout,*) 'Orthogonalisation : ICGS'
         END IF
         WRITE (nout,*) 'Restart : ', m
         WRITE (nout,*) 'info(1) = ',gi%info(1),
     &                  '  info(2) = ',gi%info(2)
         WRITE (nout,*) 'rinfo(1) = ',rinfo(1),'  rinfo(2) = ',rinfo(2)
         WRITE (nout,*) 'Optimal workspace = ', gi%info(3)
         WRITE (nout,*) 'Solution : '
         DO jcount = 1, n
            WRITE (nout,*) work(jcount)
         END DO
         WRITE (nout,*) '   '

      END IF LACTIVE2

 100  CONTINUE

      LACTIVE3: IF (lactive) THEN
         IF (ALLOCATED(xmin)) THEN ! otherwise, work(1:nloc) holds final state
            work(1:nloc) = xmin(1:nloc)
            DEALLOCATE(xmin)
         ELSE
            fsq_min = fsq_nl
         END IF

         CALL MPI_ALLGATHERV(work, nloc, MPI_REAL8, tmpbuf, gi%rcounts,
     &                       gi%disp, MPI_REAL8, MY_COMM, MPI_ERR)
         x0 = tmpbuf*bnorm

         IF (lprint .AND. kout.NE.kprint) THEN
            xmod = SQRT(SUM(x0*x0)); xmax = MAXVAL(ABS(x0))
            IF (l_NonLinear) THEN
               PRINT 910, kout, fsq_min, xmod, xmax, fsq_min_lin
            ELSE
	           PRINT 911, kout, xmod, xmax, fsq_lin
	        END IF
         END IF

         gi%ftol = fsq_min

!Sanity check
!         CALL GetNLForce(tmpbuf, fsq_nl, bnorm)
!         IF (iam .EQ. 0) PRINT *,'FSQ_MIN: ',fsq_min,' FSQ_NL: ', fsq_nl

      END IF LACTIVE3

      DEALLOCATE (work, tmpbuf, aux, itest)
#endif
      END SUBROUTINE gmres_par

      SUBROUTINE Truncate(num, iprec0)
      USE stel_kinds, ONLY: dp
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)     :: iprec0
      REAL(dp), INTENT(INOUT) :: num
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER*24 :: chnum, tchnum
!-----------------------------------------------
!
!     TRUNCATES double-precision to precision iprec0 digits, keeping exponent range of double
!     WRITE TO INTERNAL FILE TO DO TRUNCATION
!
!      RETURN

      WRITE (chnum, '(a,i2,a,i2,a)') '(1p,e',iprec0+7,'.',iprec0,')'
      WRITE (tchnum, chnum) num

      READ (tchnum, chnum) num

      END SUBROUTINE Truncate


      END MODULE GMRES_LIB

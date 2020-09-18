      MODULE lmpar_mod
      USE stel_kinds
      INTEGER :: nscan, ldfjac
      INTEGER, DIMENSION(:), POINTER :: ipvt
      REAL(rprec) :: pnorm, fnorm1, delta, par, spread_ratio
      REAL(rprec), DIMENSION(:), POINTER :: wa2p, wa3p, wa4p
      REAL(rprec), DIMENSION(:), POINTER :: diag, qtf
      REAL(rprec), DIMENSION(:,:), POINTER :: fjac
      LOGICAL :: lfirst_lm

      CONTAINS

      SUBROUTINE levmarq_param_mp(x, wa1, wa2, wa3, wa4,
     1     nfev, m, n, iflag, fcn)
      USE fdjac_mod, ONLY: flag_cleanup
      USE mpi_params
      USE mpi_inc
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: n, m
      INTEGER :: iflag, nfev
      REAL(rprec), INTENT(in) :: x(n)
      REAL(rprec) :: wa1(n), wa2(n), wa3(n), wa4(m)
      EXTERNAL fcn
#if defined(MPI_OPT)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
      REAL(rprec), DIMENSION(11), PARAMETER :: factors =
     1  (/ 1.0_dp, 0.5_dp, 0.25_dp, 0.128_dp, 0.75_dp,
     2      1.25_dp, 1.5_dp, 0.9_dp, 1.1_dp, 1.75_dp, 2.1_dp /)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iproc, iproc_min, nfact, num_lev, istat, ierr, j,
     1           iflag_min(1)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: iflag_array
      REAL(rprec) :: scale_factor
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: fnorm_array
      CHARACTER(LEN=1) :: ext, low_mark
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec), EXTERNAL :: enorm
C-----------------------------------------------
!     Perform Numprocs different function calls (in parallel) to find the minimum norm (chi-sq).
!     MPI calls are used to determine which processor has the minimum norm and then the
!     information is sent to all processors using MPI Broadcasts (modifications made by D. A. Spong 8/27/2000).
!
!     Uses processors 0,...,numprocs-1 (which belong to the user-defined MPI_COMM_WORKERS communicator)

      nfact = SIZE(factors)
      num_lev = numprocs
      ALLOCATE (iflag_array(numprocs), fnorm_array(numprocs),
     1          stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in levmarq_param_mp'

      iproc = myid
      IF (lfirst_lm .and. num_lev > 2) THEN
!
!       Do an exponential spread the first call (lfirst_lm=.true.) to see where we are
!
         scale_factor = 10._dp**(-iproc)
!        scale_factor = 10._dp**(one-iproc)
         lfirst_lm = .false.
      ELSE IF (num_lev > 2*nfact) THEN
         scale_factor = ((iproc+1)*MAXVAL(factors))/num_lev
      ELSE IF (iproc .lt. nfact) THEN
         scale_factor = factors(iproc+1)
      ELSE
         scale_factor = ((iproc-nfact)*MINVAL(factors))/
     1                   (num_lev-nfact)
      END IF

      delta = scale_factor*delta

      CALL lmpar (n, fjac, ldfjac, ipvt, diag, qtf,
     1            delta, par, wa1, wa2, wa3, wa4)
!
!     store the direction p and x + p. calculate the norm of p.
!
      IF (par .eq. zero) wa1 = wa1*scale_factor
      wa1 = -wa1
      wa2 = x + wa1
      wa3 = diag*wa1
      pnorm = enorm(n,wa3)
!
!     evaluate the function at x + p and calculate its norm.
!     Only do for 0 < myid < n processors (the MPI_COMM_WORKERS group),
!     which the v3post routines are expecting. To use other processors, must
!     pass the communicator id (maybe through hi bits of iflag, using IAND...)
!
      MPI_COMM_WORKERS = MPI_COMM_WORLD
      iflag = iproc
      CALL fcn (m, n, wa2, wa4, iflag, nfev)

!
!     Gather iflag information to all processors and check for iflag < 0
!
      CALL MPI_ALLGATHER(iflag, 1, MPI_INTEGER, iflag_array, 1,
     1     MPI_INTEGER, MPI_COMM_WORLD, ierr)
      IF (ierr .ne. 0) STOP 'MPI_ALLGATHER failed in levmarq_param_mp'

      fnorm1 = enorm(m,wa4)

!
!     Find processor with minimum fnorm1 value
!
      CALL MPI_ALLGATHER(fnorm1, 1, MPI_REAL8, fnorm_array, 1,
     1     MPI_REAL8, MPI_COMM_WORLD, ierr)
      IF (ierr .ne. 0) STOP 'MPI_ALLGATHER failed in LMDIF'
      iflag_min = MINLOC(fnorm_array)
      iproc_min = iflag_min(1) - 1

!      iflag = MINVAL(iflag_array)
!      IF (iflag .lt. 0) GOTO 100

      ext = ' '
      low_mark = ' '
      IF (myid .eq. master) ext = '*'
      IF (iproc .eq. iproc_min) low_mark = '*'
      iproc = iproc+1
      WRITE(6, '(2x,i6,8x,i3,4x,2(3x,1es12.4,a),3x,1es12.4)')
     1         iproc+nfev, iproc, fnorm1**2, low_mark, par, ext,
     2         delta

      CALL flush(6)

      fnorm1 = MINVAL(fnorm_array)

 100  CONTINUE

      DEALLOCATE (iflag_array, fnorm_array)

!
!     Broadcast all relevant scalars and arrays from the
!     processor with minimum fnorm1 to the other processors,
!     overwriting their data. Note: diag, ipvt are same already on
!     ALL processors. wa3 is overwritten...
!
      CALL MPI_BCAST(pnorm,1,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
      IF (ierr .ne. 0) GOTO 3000
      CALL MPI_BCAST(par,1,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
      IF (ierr .ne. 0) GOTO 3000
      CALL MPI_BCAST(delta,1,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
      IF (ierr .ne. 0) GOTO 3000
      CALL MPI_BCAST(wa1,n,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
      IF (ierr .ne. 0) GOTO 3000
      CALL MPI_BCAST(wa2,n,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
      IF (ierr .ne. 0) GOTO 3000
      CALL MPI_BCAST(wa4,m,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
      IF (ierr .ne. 0) GOTO 3000

!
!     BROADCAST JACOBIAN fjac(:n,j), j=1,n ROW BY ROW (CHANGED IN LMPAR) TO OTHER PROCESSES
!
      DO j = 1, n
         IF (myid .eq. iproc_min) wa3(:n) = fjac(:n,j)
         CALL MPI_BCAST(wa3, n, MPI_REAL8, iproc_min,
     1        MPI_COMM_WORLD, ierr)
         IF (ierr .ne. 0) GOTO 3000
         IF (myid .ne. iproc_min) fjac(:n,j) = wa3(:n)
      END DO

!
!     CLEANUP AFTER LEVENBERG-MARQUARDT LOOP AS NEEDED (WA4 IS NOT CHANGED)
!
      iflag = flag_cleanup
      CALL fcn (m, n, x, wa4, iflag, nfev)             !Contains Bcast Barrier

      nfev = nfev + num_lev

      RETURN

 3000 CONTINUE

      WRITE (6, *) 'MPI_BCAST error in LEVMARQ_PARAM_MP, ierr = ', ierr
      STOP
#endif
      END SUBROUTINE levmarq_param_mp

      SUBROUTINE levmarq_param(x, wa1, wa2, wa3, wa4,
     1     time, nfev, m, n, iflag, fcn)
      USE fdjac_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nfev, n, m, iflag
      REAL(rprec) :: time
      REAL(rprec), TARGET :: x(n), wa1(n), wa2(n), wa3(n), wa4(m)
      EXTERNAL fcn
#if !defined(MPI_OPT)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero=0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: j, istat, iread, ic1, ic2, irate, count_max,
     1     jmin
      REAL(rprec), DIMENSION(num_lm_params) ::
     1      fnorm_min, pnorm_min, delta_min, par_min
      CHARACTER(LEN=1) :: ext, low_mark
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL multiprocess, lmpar_parallel
C-----------------------------------------------
!
!     Initialize module variables for use in lmpar_parallel
!
      xp => x
      wap => wa1
      wa2p => wa2
      wa3p => wa3
      wa4p => wa4
      np = n
      mp = m
      ncntp = nfev

      CALL SYSTEM_CLOCK(ic1, irate)

      CALL multiprocess(num_lm_params, max_processors,
     1    lmpar_parallel, fcn)

      CALL SYSTEM_CLOCK(ic2, irate, count_max)
      IF (ic2 .lt. ic1) ic2 = ic2 + count_max

      nfev = nfev + num_lm_params

!
!     Read in optimal wa1, wa2, wa4, par, delta, pnorm, fnorm1 value from file
!
      DO j = 1, num_lm_params

        READ (j+1000, iostat=iread) istat, iflag, pnorm_min(j),
     1        fnorm_min(j), par_min(j), delta_min(j)
        IF (iread .ne. 0) THEN
           WRITE (6, *) 'Error reading from file fort.', j+1000,
     1       ' in levmarq_param', ' IOSTAT = ', iread
           iflag = -15
        ELSE IF (j .ne. istat) THEN
           WRITE (6, *)
     1        'Incorrect value READ in for INDEX j in levmarq_param'
           iflag = -15
        END IF

        IF (iflag .NE. 0) RETURN

        IF (j .EQ. 1) fnorm1 = fnorm_min(j)
        IF (fnorm_min(j) .LE. fnorm1) THEN
           jmin = j
           fnorm1 = fnorm_min(jmin)
           pnorm  = pnorm_min(jmin)
           par    = par_min(jmin)
           delta  = delta_min(jmin)
#if defined(CRAY)
           DO k = 1, n
              READ (j+1000) wa1(k), wa2(k)
              DO istat = 1, n
                 READ (j+1000) fjac(k, istat)
              END DO
           END DO
           DO k = 1, m
              READ (j+1000) wa4(k)
           END DO
#else
           READ (j+1000) wa1, wa2, wa4, fjac(1:n, 1:n)
#endif

        END IF

        CLOSE (j+1000, status='delete')                        !!Needed to run correctly in multi-tasking...

      END DO

      DO j = 1, num_lm_params
         ext = ' '
         low_mark = ' '
         IF (j .eq. 1) ext = '*'
         IF (j .eq. jmin) low_mark = '*'
         WRITE (6, '(2x,i6,4x,2(3x,1es12.4,a),3x,1es12.4)') j+ncntp,
     1         fnorm_min(j)**2, low_mark, par_min(j), ext, delta_min(j)
      END DO

!
!     Do any special cleanup now for IFLAG = flag_cleanup. WA4 LEFT UNCHANGED
!
      iflag = flag_cleanup
      CALL fcn(m, n, x, wa4, iflag, ncntp)

      time = time + REAL(ic2 - ic1)/REAL(irate)                !!Time in multi-process CALL
#endif
      END SUBROUTINE levmarq_param

      END MODULE lmpar_mod

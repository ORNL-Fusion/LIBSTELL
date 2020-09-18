      SUBROUTINE lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev,    &
         epsfcn, diag, mode, factor, nprint, info, nfev, fjac,          &
         ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4)
      USE stel_kinds
      USE lmpar_mod, fjac_mod=>fjac, ldfjac_mod=>ldfjac,                &
         ipvt_mod=>ipvt, qtf_mod=>qtf, diag_mod=>diag
#if defined(MPI_OPT)
      USE fdjac_mod, ONLY: flip,flag_singletask,flag_cleanup,fdjac2_mp
#else
      USE fdjac_mod, ONLY: max_processors, flip, flag_singletask,       &
                           flag_cleanup, fdjac2
#endif
      USE mpi_params
      USE mpi_inc
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: m, n, maxfev, mode, nprint, info, nfev, ldfjac
      REAL(rprec), INTENT(in) ::  ftol, xtol, gtol, epsfcn, factor
      REAL(rprec), DIMENSION(n) :: x, wa1, wa2, wa3
      REAL(rprec), DIMENSION(m) :: fvec, wa4
      INTEGER, DIMENSION(n), TARGET :: ipvt
      REAL(rprec), DIMENSION(n), TARGET :: diag, qtf
      REAL(rprec), DIMENSION(ldfjac,n), TARGET :: fjac
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1,                      &
         p1=0.1_dp, p5=0.5_dp, p25=0.25_dp, p75=0.75_dp, p0001=1.e-4_dp
      CHARACTER(LEN=130), DIMENSION(0:11) :: info_array 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: iflag, iter, j, l, istat, ikey
      INTEGER :: cycle_count, subcycle
      REAL(rprec) :: actred, dirder, epsmch, fnorm,                     &
         gnorm, prered, ratio, sum0, temp,                              &
         temp1, temp2, xnorm
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: fjac_save
#if !defined(MPI_OPT)
      REAL(rprec) :: wall_time, wall_time_lev
#endif
      REAL(rprec) :: fnorm_min
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: x_min, fvec_min
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      EXTERNAL fcn
      REAL(rprec), EXTERNAL :: dpmpar, enorm
!-----------------------------------------------
!
!     SUBROUTINE lmdif
!
!     the purpose of lmdif is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the levenberg-marquardt algorithm. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.
!
!     the subroutine statement is
!
!       subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
!                        diag,mode,factor,nprint,info,nfev,fjac,
!                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling. see lmdif1 for
!         documentation and should be written as follows.
!
!         subroutine fcn(m, n, x, fvec, iflag, ncnt)
!         integer m,n,iflag
!         real(rprec) x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!
!       ftol is a nonnegative input variable. termination
!         occurs when both the actual and predicted relative
!         reductions in the sum of squares are at most ftol.
!         therefore, ftol measures the relative error desired
!         in the sum of squares.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol. therefore, xtol measures the
!         relative error desired in the approximate solution.
!
!       gtol is a nonnegative input variable. termination
!         occurs when the cosine of the angle between fvec and
!         any column of the jacobian is at most gtol in absolute
!         value. therefore, gtol measures the orthogonality
!         desired between the function vector and the columns
!         of the jacobian.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn is at least
!         maxfev by the end of an iteration.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.). 100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x and fvec available
!         for printing. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set
!       nfev is an integer output variable set to the number of
!         calls to fcn.
!
!       fjac is an output m by n array. the upper n by n submatrix
!         of fjac contains an upper triangular matrix r with
!         diagonal elements of nonincreasing magnitude such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower trapezoidal
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular
!         with diagonal elements of nonincreasing magnitude.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       qtf is an output array of length n which contains
!         the first n elements of the vector (q transpose)*fvec.
!
!       wa1, wa2, and wa3 are work arrays of length n.
!
!       wa4 is a work array of length m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac
!
!       fortran-supplied ... ABS,max,min,sqrt,mod
!
!     argonne national laboratory. MINpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********

#if defined(MPI_OPT)
!     Get mpi parameters
      CALL MPI_COMM_RANK (MPI_COMM_WORLD, myid, ierr_mpi)       !mpi stuff
      IF (ierr_mpi .ne. 0) STOP 'MPI_COMM_RANK error in LMDIF'
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD, numprocs, ierr_mpi)   !mpi stuff
      IF (ierr_mpi .ne. 0) STOP 'MPI_COMM_SIZE error in LMDIF'
#endif
      info = 0;      iflag = 0;      nfev = 0;      cycle_count = 0

      info_array(9) =                                                   &
      "improper input parameters (number constraints MUST " //          &
      "be greater than number variables (>0)."

      info_array(1) =                                                   &
      "algorithm estimates that the relative error in the " //          &
      "sum of squares is at most ftol."

      info_array(2) =                                                   &
      "algorithm estimates that the relative error between" //          &
      " x and the solution is at most xtol."

      info_array(3) =                                                   &
      "algorithm estimates that the relative error in the " //          &
      "sum of squares and between x and the solution is at" //          &
      " most ftol and xtol."

      info_array(4) =                                                   &
      "the cosine of the angle between fvec and any column" //          &
      " of the jacobian is at most gtol in absolute value."

      info_array(5) =                                                   &
      "number of calls to fcn has reached or exceeded maxfev."

      info_array(6) =                                                   &
      "ftol is too small. no further reduction in the sum " //          &
      "of squares is possible."

      info_array(7) =                                                   &
      "xtol is too small. no further improvement in the " //            &
      "approximate solution x is possible."

      info_array(8) =                                                   &
      "gtol is too small. fvec is orthogonal to the columns" //         &
      " of the jacobian to machine precision."

      info_array(0) =                                                   &
      "levenberg-marquardt optimizer terminated properly."              

      info_array(10) =                                                  &
      "Levenberg-Marquardt terminated prematurely due to " //           &
      "function evaluation error." 

      info_array(11) =                                                  &
      "Levenberg-Marquardt terminated prematurely: " //                 &
      "must request more than one (1) processor" 

#if defined(MPI_OPT)
      IF (numprocs > n) THEN
         IF (myid .eq. master) THEN
            WRITE (6, *)'Warning: more processors have been requested', &
            ' than the maximum (nvar) required = ',n
         END IF
      ELSE IF (numprocs <= 1) THEN
         info = 11
         GOTO 400
      END IF
#endif
!
!     check the input parameters for errors.
!

      ALLOCATE (x_min(n), fvec_min(m), flip(n), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation error in lmdif'

!
!     epsmch is the machine precision.
!     flip is control for direction flipping in fdjac2!

      epsmch = dpmpar(1)
      flip = .false.
#if !defined(MPI_OPT)
      myid = 0;      wall_time = 0;      wall_time_lev = 0
#endif
!
!     ASSIGN VALUES TO MODULE VARIABLES (FACILITATES PASSING TO SUBROUTINES)
!
      ldfjac_mod = ldfjac
      ipvt_mod => ipvt
      fjac_mod => fjac
      diag_mod => diag
      qtf_mod => qtf
!
!     check the input parameters for errors.
!
      IF (n.le.0 .or. m.lt.n .or. ldfjac.lt.m .or. ftol.lt.zero         &
         .or. xtol.lt. zero .or. gtol.lt.zero .or. maxfev.le.0          &
         .or. factor.le.zero) THEN
         info = 9
         GOTO 400
      END IF

      IF (mode .EQ. 2) THEN
         DO j = 1, n
            IF (diag(j) .le. zero) GOTO 300
         END DO
      END IF

!
!     evaluate the function at the starting point
!     and calculate its norm.
!
!     Set up workers communicator (only master processor here) for initial run
!
#if defined(MPI_OPT)
      IF (myid .ne. master) THEN
         ikey = MPI_UNDEFINED
      ELSE
         ikey = WORKER_SPLIT_KEY+1
      END IF
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, ikey, worker_id,              &
                          MPI_COMM_WORKERS, ierr_mpi)
      IF (ierr_mpi .ne. 0) STOP 'MPI_COMM_SPLIT error in lsfun1'
#endif
      IF (myid .EQ. master) THEN
         iflag = flag_singletask
         CALL fcn (m, n, x, fvec, iflag, nfev)
#if defined(MPI_OPT)
         CALL MPI_COMM_FREE(MPI_COMM_WORKERS, ierr_mpi)      
#endif
      END IF


#if defined(MPI_OPT)
      CALL MPI_BCAST(iflag,1, MPI_INTEGER, master,                      &
                     MPI_COMM_WORLD, ierr_mpi)
      IF (ierr_mpi .NE. 0) GOTO 3000
      IF (iflag .ge. 0) CALL                                            &
           MPI_BCAST(fvec, m, MPI_REAL8, master,                        &
                     MPI_COMM_WORLD, ierr_mpi)
      IF (ierr_mpi .NE. 0) GOTO 3000
#endif

      IF (iflag .LT. 0) GOTO 300
      fnorm = enorm(m,fvec)
      IF (nfev.GE.maxfev .OR. maxfev.EQ.1) info = 5
      IF (info .NE. 0) GOTO 300
!
!     initialize levenberg-marquardt parameter (par) and iteration counter.
!
      par = 0
      iter = 1
      lfirst_lm = .true.
#if defined(MPI_OPT)
      IF (myid .EQ. master) WRITE (6, 1000) numprocs
 1000 FORMAT (/,' Beginning Levenberg-Marquardt Iterations',/,          &
              ' Number of Processors: ',i4,//,                          &
!     1        ' Number of Processors: ',i4,' (1 controller proc)',//,
              70('='),/,2x,'Iteration',3x,'Processor',7x,'Chi-Sq',7x,   &
             'LM Parameter',6x,'Delta Tol'/,70('='))
#else
      WRITE (6, 1000) max_processors
 1000 FORMAT (/,' Beginning Levenberg-Marquardt Iterations',/,          &
              ' Number processors requested: ', i4,//,                  &
              59('='),/,2x,'Iteration',8x,'Chi-Sq',7x,                  &
             'LM Parameter',6x,'Delta Tol',/,59('='))
#endif
!
!     beginning of the outer loop.
!
      outerloop: DO WHILE (nfev .lt. maxfev)
!
!        calculate the jacobian matrix.
!
#if defined(MPI_OPT)
         CALL fdjac2_mp(fcn, m, n, x, fvec, fjac, ldfjac,               &
                 iflag, nfev, epsfcn, fnorm_min, x_min, fvec_min)
#else
         iflag = 2
         CALL fdjac2(fcn, m, n, x, fvec, fjac, ldfjac, iflag, nfev,     &
                     epsfcn, wa4, wall_time, fnorm_min, x_min, fvec_min)
#endif
         nfev = nfev + n
         IF (iflag .LT. 0) EXIT outerloop

!
!        compute the qr factorization of the jacobian.
!
         CALL qrfac(m, n, fjac, ldfjac, .true., ipvt,                   &
                    n, wa1, wa2, wa3)

!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
         IF (iter .EQ. 1) THEN
            IF (mode .ne. 2) THEN
               diag = wa2
               WHERE (wa2 .eq. zero) diag = one
            END IF

!
!        also on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
            wa3 = diag*x
            xnorm = enorm(n,wa3)
            delta = factor*xnorm
            IF (delta .eq. zero) delta = factor
         END IF

!
!        form (q TRANSPOSE)*fvec and store the first n components in qtf.
!
         wa4 = fvec
         DO j = 1, n
            IF (fjac(j,j) .ne. zero) THEN
               sum0 = SUM(fjac(j:m,j)*wa4(j:m))
               temp = -sum0/fjac(j,j)
               wa4(j:m) = wa4(j:m) + fjac(j:m,j)*temp
            END IF
            fjac(j,j) = wa1(j)
            qtf(j) = wa4(j)
         END DO

!
!        compute the norm of the scaled gradient.
!
         gnorm = zero
         IF (fnorm .ne. zero) THEN
            DO j = 1, n
               l = ipvt(j)
               IF (wa2(l) .ne. zero) THEN
                  sum0 = SUM(fjac(1:j,j)*(qtf(1:j)/fnorm))
                  gnorm = MAX(gnorm,ABS(sum0/wa2(l)))
               END IF
            END DO
         END IF

!
!        test for convergence of the gradient norm.
!
         IF (gnorm .LE. gtol) info = 4
         IF (info .NE. 0) EXIT outerloop

!
!        rescale if necessary.
!
         IF (mode .ne. 2) diag = MAX(diag,wa2)

!
!        set up for inner loop (levmarqloop) to determine x update.
!
         subcycle = 0
         ratio = 0

         ALLOCATE (fjac_save(n,n), stat=istat)
         IF (istat .NE. 0) STOP 'Fjac_save allocation error'

         fjac_save(:n,:n) = fjac(:n,:n)

         levmarqloop: DO WHILE (ratio.lt.p0001 .and. subcycle.lt.2)

           subcycle = subcycle + 1
           fjac(:n,:n) = fjac_save(:n,:n)
           spread_ratio = epsfcn/delta/10

!
!        Determine the levenberg-marquardt parameter.
!
#if defined(MPI_OPT)
!        for parallel processing, scan a range of values for this
!        parameter to find the 'optimal' choice
!
           CALL levmarq_param_mp (x, wa1, wa2, wa3, wa4,                &
                                  nfev, m, n, iflag, fcn)
#else
           CALL levmarq_param(x, wa1, wa2, wa3, wa4,                    &
                wall_time_lev, nfev, m, n, iflag, fcn)
#endif
           IF (iflag .LT. 0) EXIT

!
!        on the first iteration, adjust the initial step bound.
!
           IF (iter .EQ. 1) delta = MIN(delta, pnorm)

!
!        compute the scaled actual reduction.
!
           actred = -one
           IF (p1*fnorm1 .LT. fnorm) actred = one - (fnorm1/fnorm)**2

!
!        compute the scaled predicted reduction (prered) and
!        the scaled directional derivative (dirder).
!
           DO j = 1, n
              wa3(j) = zero
              l = ipvt(j)
              temp = wa1(l)
              wa3(1:j) = wa3(1:j) + fjac(1:j,j)*temp
           END DO

           temp1 = enorm(n,wa3)/fnorm
           temp2 = (SQRT(par)*pnorm)/fnorm
           prered = temp1**2 + temp2**2/p5
           dirder = -(temp1**2 + temp2**2)

!
!        compute the ratio of the actual to the predicted reduction.
!
           ratio = zero
           IF (prered .NE. zero) ratio = actred/prered

!
!        test if Jacobian calculation gave best improvement
!        (added by MZ and SH)
!
           IF (fnorm_min .LT. MIN(fnorm,fnorm1)) THEN
              IF (myid .EQ. master) WRITE (6, *)                        &
                 ' Using minimum state from Jacobian evaluation'
              wa2 = x_min
              wa4 = fvec_min
              fnorm1 = fnorm_min
              IF (p1*fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
              IF (prered .ne. zero) ratio = actred/prered
           END IF

!
!        update the step bound.
!
           IF (ratio .le. p25) THEN
              IF (actred .ge. zero) THEN
                 temp = p5
              ELSE
                 temp = p5*dirder/(dirder + p5*actred)
              END IF
              IF (p1*fnorm1.ge.fnorm .or. temp.lt.p1) temp = p1
              delta = temp*MIN(delta,pnorm/p1)
              par = par/temp
           ELSE IF (par.eq.zero .or. ratio.ge.p75) THEN
              delta = pnorm/p5
              par = p5*par
           END IF

!
!        test for successful iteration.
!        update x, fvec, and their norms.
!
           IF (ratio .GE. p0001) THEN
              x = wa2
              wa2 = diag*x
              fvec = wa4
              xnorm = enorm(n,wa2)
              fnorm = fnorm1
              iter = iter + 1
              IF (myid .EQ. master) WRITE(6,'(/,3(2x,a,1es10.3)/)')     &
                'new minimum = ', fnorm**2,'lm-par = ', par,            &
                'delta-tol = ', delta
           END IF

           cycle_count = cycle_count + 1

!
!        tests for convergence.
!
           IF (ABS(actred).LE.ftol .AND. prered.LE.ftol                 &
              .AND. p5*ratio.LE.one) info = 1
!
!        next test made more stringent to avoid premature declaration of
!        completion observed on some problems.
!
           IF (delta.LE.xtol*xnorm .AND. cycle_count.GT.2) info = 2
!          IF (delta .le. xtol*xnorm) info = 2

           IF (ABS(actred).LE.ftol .AND. prered.LE.ftol                 &
              .AND. p5*ratio.LE.one .AND. info.EQ.2) info = 3
           IF (info .NE. 0) EXIT levmarqloop

!
!        tests for termination and stringent tolerances.
!
           IF (nfev .ge. maxfev) info = 5
           IF (ABS(actred).le.epsmch .and. prered.le.epsmch             &
              .and. p5*ratio.le.one) info = 6
           IF (delta .le. epsmch*xnorm) info = 7
           IF (gnorm .le. epsmch) info = 8
           IF (info .ne. 0) EXIT levmarqloop
!
!        END of the inner loop. repeat if iteration unsuccessful (ratio < p0001)
!
         END DO levmarqloop

         DEALLOCATE (fjac_save)
         IF (info.ne.0 .or. iflag.ne.0) EXIT outerloop

      END DO outerloop

!     termination, either normal or user imposed.

 300  CONTINUE

      IF (info.EQ.0 .AND. iflag.NE.0) info = 10

 400  CONTINUE

      IF (myid .EQ. master) THEN
         IF (info.ge.LBOUND(info_array,1) .and.                         &
             info.le.UBOUND(info_array,1)) THEN
            WRITE (6, '(2(/,1x,a))')                                    &
        'Levenberg-Marquardt optimizer status: ',TRIM(info_array(info))
         ELSE
            WRITE (6, '(a,i5)')' LM info out of bounds: ', info
         END IF
      ENDIF                                                              ! MPI

      IF (ALLOCATED(x_min)) DEALLOCATE (x_min, fvec_min, flip)

      IF (iflag .LT. 0) info = iflag

      IF (nfev.LE.1 .and. info.NE.11) THEN
         nfev = 1
         iflag = flag_cleanup                        !!Clean-up last time through for master
         CALL fcn (m, n, x, fvec, iflag, nfev)
      END IF

#if defined(MPI_OPT)
      RETURN

 3000 CONTINUE
      WRITE (6, *) 'MPI_BCAST error in LMDIF, ierr = ', ierr_mpi
#else
      WRITE(*, '(2(/,a, f10.2))')                                       &
           ' Total wall clock time in jacobian multi-process call  = ', &
           wall_time,                                                   &
           ' Total wall clock time in lev param multi-process call = ', &
           wall_time_lev
#endif
      END SUBROUTINE lmdif

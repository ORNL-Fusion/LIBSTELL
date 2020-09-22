      SUBROUTINE lmdif1(fcn, m, n, x, fvec, tol, epsfcn, nfev_end,      &
                 diag, mode, info, lwa, max_processors, num_lm_params)
!!    ADDED EPSFCN TO ARG LIST: BY SPH (2/97)
!!    ADDED NFEV_end == MAXFEV TO ARG LIST (6/31/99)
!!    ADDED MAX_PROCESSORS, NUM_LM_PARAMS TO ARG LIST (11/23/99)
!!    (NEED MAX_PROCESSORS, NUM_LM_PARAMS FOR MULTI-PROCESSOR APPLICATIONS)

      USE fdjac_mod, ONLY: maxj_processors=>max_processors,             &
                           numj_lm_params=>num_lm_params
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: m, n, lwa, nfev_end, mode
      INTEGER, INTENT(out) :: info
      REAL(rprec), INTENT(in) :: tol, epsfcn
      REAL(rprec), DIMENSION(n), INTENT(inout) :: x, diag
      REAL(rprec), DIMENSION(m), INTENT(out) :: fvec
      INTEGER, INTENT(in) :: max_processors, num_lm_params
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: zero=0, factor=10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(:), ALLOCATABLE :: iwa
      INTEGER :: maxfev, mp5n, nfev, nprint
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: wa
      REAL(rprec) :: ftol, gtol, xtol
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      EXTERNAL fcn
!-----------------------------------------------
      INTERFACE
         SUBROUTINE lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, &
                  epsfcn, diag, mode, factor, nprint, info, nfev, fjac, &
                  ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4)
         USE stel_kinds
         INTEGER :: m, n, maxfev, mode, nprint, info, nfev, ldfjac
         REAL(rprec), INTENT(in) ::  ftol, xtol, gtol, epsfcn, factor
         REAL(rprec), DIMENSION(n) :: x, wa1, wa2, wa3
         REAL(rprec), DIMENSION(m) :: fvec, wa4
         INTEGER, DIMENSION(n), TARGET :: ipvt
         REAL(rprec), DIMENSION(n), TARGET :: diag, qtf
         REAL(rprec), DIMENSION(ldfjac,n), TARGET :: fjac
         EXTERNAL fcn
         END SUBROUTINE lmdif
      END INTERFACE

!
!     SUBROUTINE lmdif1
!
!     the purpose of lmdif1 is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of the
!     levenberg-marquardt algorithm. this is done by using the more
!     general least-squares solver lmdif. the user must provide a
!     SUBROUTINE which calculates the functions. the jacobian is
!     THEN calculated by a forward-difference approximation.
!
!     the SUBROUTINE statement is
!
!       SUBROUTINE lmdif1(fcn,m,n,x,fvec,tol,info,lwa)
!
!     WHERE
!
!       fcn is the name of the user-supplied SUBROUTINE which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         SUBROUTINE fcn(m, n, x, fvec, iflag, ncnt)
!         INTEGER m,n,iflag
!         REAL(rprec) x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         RETURN this vector in fvec.
!         ----------
!         RETURN
!         END
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmdif1.
!         in this CASE set iflag to a negative INTEGER. On a multi-processor
!         machine, iflag will be initialized to the particular processor id.
!
!
!       m is a positive INTEGER input variable set to the number
!         of functions.
!
!       n is a positive INTEGER input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which CONTAINS
!         the functions evaluated at the output x.
!
!       ncnt is a positive INTEGER input variable set to the current
!         iteration count (added by SPH - 7/99)
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates either that the relative
!         error in the sum of squares is at most tol or that
!         the relative error between x and the solution is at
!         most tol.
!
!       info is an INTEGER output variable. IF the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn.
!
!       lwa is a positive INTEGER input variable not less than
!         m*n+5*n+m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... lmdif
!
!     argonne national laboratory. MINpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     modified to accept improvements from jacobian calc. MZarnstorff Oct 2001
!     modified to flip sign of jacobian offsets for more efficient search
!     and start with an exponential levenberg spread to settle on scale
!             M. Zarnstorff                          Jan 2002
!
!     **********
!
!     CALL lmdif.
!
      IF (lwa .lt. m*n+5*n+m) STOP 'lwa too small'

      ALLOCATE (wa(lwa), iwa(n), stat=info)
      IF (info .NE. 0) STOP 'Allocation error in lmdif1!'

#if !defined(MPI_OPT)
!
!     Load fdjac module values
!
      maxj_processors = MAX(max_processors,1)
      numj_lm_params  = MAX(num_lm_params,1)
#endif
      maxfev = 200*(n + 1)
      maxfev = MIN (maxfev, nfev_end)             !!SPH-Added 7/99
      ftol = tol
      xtol = tol
      gtol = zero

!!    ADDED BY SPH -- PASSED IN ARG LIST(2/97)
!!    epsfcn = zero
!!    mode = 1       (DAS, passed through arg list 9/13/00)

      nprint = 0
      mp5n = m + 5*n
      fvec = 0
      CALL lmdif (fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev,         &
                  epsfcn, diag, mode, factor, nprint, info, nfev,       &
                  wa(mp5n+1), m, iwa, wa(n+1), wa(2*n+1), wa(3*n+1),    &
                  wa(4*n+1), wa(5*n+1))

#if !defined(MPI_OPT)
      IF (info .eq. 8) info = 4
#endif
      DEALLOCATE(wa, iwa)

      END SUBROUTINE lmdif1

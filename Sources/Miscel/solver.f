      SUBROUTINE solver(amat, b, m, nrhs, info)
      USE stel_kinds
#ifdef ACCELERATE_NEW_LAPACK
      USE accelerate_lapack_lp64
#endif
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in)  :: m, nrhs
      INTEGER, INTENT(out) :: info
      REAL(rprec), INTENT(inout) :: amat(m,m), b(m,nrhs)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, ALLOCATABLE :: ipiv(:)
C-----------------------------------------------
      info = 0
      ALLOCATE (ipiv(m))

!     Compute the solution to a REAL system of linear equations
!       AMAT * X = B,
!     WHERE AMAT is an M-by-M matrix and X and B are N-by-NRHS matrices.
!
!     FACTOR AMATRIX INTO LU FORM
!     AND SOLVE BY GAUSSIAN ELIMINATION
!
#ifdef ACCELERATE_NEW_LAPACK
      CALL dgesv (m, nrhs, amat(1,1), m, ipiv(1), b(1,1), m, info)
#else
      CALL dgesv (m, nrhs, amat, m, ipiv, b, m, info)
#endif
!     IF (info .ne. 0) PRINT *, ' Condition No. = 0 in Solver'

      DEALLOCATE (ipiv)

      END SUBROUTINE solver

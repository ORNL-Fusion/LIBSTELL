      SUBROUTINE lmpar_parallel(j, fcn)
      USE fdjac_mod, ONLY: wa1p => wap, m=>mp, n=>np,                   &
         num_lm_params, xp, ncnt=>ncntp
      USE lmpar_mod
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: j
      EXTERNAL fcn
#if !defined(MPI_OPT)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
      REAL(rprec), DIMENSION(11), PARAMETER :: factors =                &
        (/ 1.0_dp, 0.5_dp, 0.25_dp, 0.128_dp, 0.75_dp,                  &
            1.25_dp, 1.5_dp, 0.9_dp, 1.1_dp, 1.75_dp, 2.1_dp /)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: iflag, nfact
#if defined(CRAY)
      INTEGER :: istat, k
#endif
      REAL(rprec) :: deltain, parin, fnorm_in, pnorm_in, scale_factor
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      REAL(rprec), EXTERNAL :: enorm
!-----------------------------------------------
!
!     THIS ROUTINE IS PASSED TO THE MULTI-PROCESSOR HANDLING
!     ROUTINE


! ***************************************************
!  stepping algorithm similar to that used in the original parallel optimizer
!  by M.Zarnstorff and S. Ethier,  Feb. 1999
!
!  Re-implemented,  MCZ  July 2000
! ***************************************************
      nfact = SIZE(factors)

      IF (lfirst_lm .and. num_lm_params > 2) THEN
!
!       do an exponential spread the first time to see where we are
!
!SPH        scale_factor = EXP((j-1)*log(spread_ratio)/num_lm_params)
         scale_factor = 10._dp**(1._dp - j)
      ELSE IF (num_lm_params > 2*nfact) THEN
         scale_factor = (j*MAXVAL(factors))/num_lm_params
      ELSE IF (j .le. nfact) THEN
         scale_factor = factors(j)
      ELSE
         scale_factor =((j-nfact)*MINVAL(factors))/(num_lm_params-nfact)
      ENDif

      deltain = delta * scale_factor

!
!     Compute perturbation vector (wa1p) and Lev/Marq PARAMETER (par)
!     for different tolerances, delta
!
      parin = par

      CALL lmpar (n, fjac, ldfjac, ipvt, diag, qtf, deltain, parin,    &
                  wa1p, wa2p, wa3p, wa4p)

!
!     store the direction p and x + p. calculate the norm of p.
!
      IF (parin.eq.zero .and. j.ne.1) wa1p = wa1p*scale_factor
      wa1p = -wa1p
      wa2p = xp + wa1p
      wa3p = diag*wa1p
      pnorm_in = enorm(n, wa3p)

!
!     evaluate the function at x + p and calculate its norm.
!c
      iflag = j
      CALL fcn (m, n, wa2p, wa4p, iflag, ncnt)

      fnorm_in = enorm(m, wa4p)

!
!     OPEN A UNIQUE FILE FOR I/O IN MULTI-PROCESSOR SYSTEM
!
      WRITE (j+1000) j, iflag, pnorm_in, fnorm_in, parin, deltain
#if defined(CRAY)
      DO k = 1, n
         WRITE (j+1000) wa1p(k), wa2p(k)
         DO istat = 1, n
            WRITE (j+1000) fjac(k, istat)
         END DO
      END DO
      DO k = 1, m
         WRITE (j+1000) wa4p(k)
      END DO
#else
      WRITE (j+1000) wa1p, wa2p, wa4p, fjac(1:n, 1:n)
#endif
      CLOSE (j+1000)                      !!Needed to run correctly in multi-tasking...
#endif
      END SUBROUTINE lmpar_parallel

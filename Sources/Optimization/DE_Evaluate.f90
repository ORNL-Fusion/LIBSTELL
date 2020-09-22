      SUBROUTINE DE_Evaluate (num_proc, fcn, val, NP, Dim_XC, nfeval)
      USE stel_kinds
      USE DE_mod
      IMPLICIT NONE

      INTEGER, PARAMETER :: iflag_cleanup = -100
      INTEGER :: num_proc, NP, Dim_XC, nfeval
      REAL(rprec), DIMENSION(NP) :: val

      REAL(rprec), DIMENSION(nopt) :: fvec
      REAL(rprec) :: funcval

      INTEGER :: iflag, j, istat, jstat
      EXTERNAL fcn
#if !defined(MPI_OPT)
      EXTERNAL de_parallel
#endif

      n_pop = NP
      n_free = Dim_XC
      nfev = nfeval

#if defined(MPI_OPT)
      CALL de_mpi(np, fcn, val)
#else
      CALL multiprocess(NP, num_proc, de_parallel, fcn)

!
!  gather results here from ALL processors
!
      DO  j=1, NP

!
!  Read in the results from the individual evaluations
!
         READ(j+1000, iostat=istat) jstat, iflag, funcval
         IF( istat .ne. 0) WRITE(6,*) 'Iostat =',istat,' for CASE ',j

         IF( jstat .ne. j ) THEN
            WRITE(6,*) "wrong INDEX READ in de_evaluate"
            iflag=-14
            EXIT
         END IF

         val(j)=funcval

         CLOSE(j+1000, status='delete')

      ENDDO
#endif
      nfeval = nfeval + np
      iflag = iflag_cleanup
      CALL fcn (nopt, n_free, ui_XC(1,:), fvec, iflag, nfeval)

      END SUBROUTINE DE_Evaluate

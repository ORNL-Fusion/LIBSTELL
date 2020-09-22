      SUBROUTINE DE_DRIVER(fcn, n_opt, n_var, x, fvec, tol, eps,        &
           num_iter_opt, max_processors, filename, info, lwa, lrestart )
! ******************************************************************
      USE stel_kinds
      USE de_mod
      USE gade_mod, ONLY: ga_de
      USE system_mod
      USE mpi_params, ONLY: master, myid
      USE safe_open_mod
      USE mpi_inc
      IMPLICIT NONE

      INTEGER :: ierr
      INTEGER, PARAMETER :: iflag_cleanup = -100
      INTEGER :: n_opt, n_var, info, lwa, num_iter_opt, max_processors
      REAL(rprec), DIMENSION(n_opt), TARGET :: fvec
      REAL(rprec), DIMENSION(n_var) :: x
      REAL(rprec) :: tol, eps, chi_sq
      EXTERNAL fcn
      CHARACTER(LEN=*), INTENT(in) :: filename
      CHARACTER(LEN=LEN_TRIM(filename)+25) :: temp
      LOGICAL :: lrestart
!
!     local variables
!
      INTEGER :: iunit
      INTEGER :: i, istat, iflag, irestart = 25
      REAL(rprec), DIMENSION(n_var) :: par_max, par_min, partemp
      REAL(rprec) ::  tgt

! ******************************************************************
!  entries for the 'de' NAMELIST
!
! npopsiz    -  population SIZE
! ngen       -  number of generations
! idum       -  IF < 0, THEN |idum| is used as seed for random-number gen.
! strategy   - The strategy of the mutation operations is used in HDE,
!              except that cross-over strategy is separated into separate
!              specifier CR_strategy below
! CR_strategy - cross-over strategy 0=exponential, 1=binomial
! f_cross    - crossover scaling factor
! pcross     -  crossover probability (CR)
! parmin     - array specifying minimum value for each free-parameter,
! parmax     - array specifying maximum value for each free-parameter,
! ibound     - =1 then interpret parmin and parmax as scale-factors to be
!              multiplied by the initial guess values for each parameter
!              =0 interpret parmin and parmax as absolute values
! out_iter   - The intermediate output will be produced after "out_iter"
!              iterations. no intermediate output will be produced if
!              "out_iter < 1".
!
!     IMPORTANT: IT IS ASSUMED THAT THE MPI_PARAMS MODULE
!                HAS BEEN LOADED EXTERNALLY WITH MPI_INIT, MPI_COMM_RANK CALLS
! ******************************************************************
      info = 0

      n_free = n_var
      nopt = n_opt
      n_pop = npopsiz
      IF( n_pop <= 0) n_pop = 10*n_free
      IF( (pcross < 0._dp) .or. (pcross > 1.0_dp) ) THEN
         IF (myid .eq. master)                                          &
             WRITE(6,*) 'Illegal value of pcross for de:', pcross
         STOP
      ENDif
!     out_iter = MAX(1, out_iter)


      IF( ibound .eq. 1 ) THEN
         par_max(:n_var) = x(:n_var)*parmax(:n_var)
         par_min(:n_var) = x(:n_var)*parmin(:n_var)

         WHERE (par_max(:n_var) < par_min(:n_var) )
            partemp(:n_var) = par_max(:n_var)
            par_max(:n_var) = par_min(:n_var)
            par_min(:n_var) = partemp(:n_var)
         END WHERE

      ELSE
         par_max(:n_var) = parmax(:n_var)
         par_min(:n_var) = parmin(:n_var)
      ENDif

      nfev=1

      IF (myid .eq. master) THEN
!         WRITE(6, nml = ga_de)

! ------setup restart file -------------

         IF( lrestart ) THEN
            temp = "cp ../de_restart." // TRIM(filename) // " ."
            CALL system(TRIM(temp))
         END IF

         temp = "de_restart." // filename
         CALL safe_open(irestart, istat, TRIM(temp),                    &
                        'unknown', 'formatted')
         IF (istat .ne. 0) STOP 'Error opening de_restart file'
      END IF

      iflag=1
      tgt = 0
      CALL DE_Evolve(fcn,n_free,par_min,par_max, tgt, n_pop, ngen,      &
              f_cross,pcross,strategy,CR_strategy,out_iter,iunit,       &
              irestart, max_processors,lrestart, x, chi_sq,nfev)

      IF (myid .eq. master) THEN
         WRITE(6,*) "final solution: "
!        WRITE(6,*) "x ", x(:n_var)
!        WRITE(6,*) "fvec ",(fvec(i),i=1,n_opt)
         WRITE(6,*) "x "
         WRITE(6,990) (x(i), i=1,n_var)
 990  FORMAT(5(1pe22.14))
         WRITE(6,*) "y ", chi_sq
      END IF

      iflag = iflag_cleanup
      CALL fcn(n_opt, npopsiz, x, fvec, iflag, nfev)

      END SUBROUTINE DE_DRIVER

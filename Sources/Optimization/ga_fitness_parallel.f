      SUBROUTINE ga_fitness_parallel(j,fcn)
      USE ga_mod
      IMPLICIT NONE
      INTEGER :: j
      EXTERNAL fcn
#if !defined(MPI_OPT)
      INTEGER :: iflag
      REAL(rprec) :: funcval, ga_evaluate

      iflag=j
      funcval = ga_evaluate(fcn, num_obj, f_obj, nparam, parent(1,j),
     >                   iflag, nfit_eval)
      WRITE (j+1000) j, iflag
      WRITE (j+1000) funcval
      CLOSE (j+1000)
#endif
      END SUBROUTINE ga_fitness_parallel

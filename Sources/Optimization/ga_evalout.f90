      SUBROUTINE ga_evalout(fbar, best, fcn, nopt, fvec, nfev,          &
                            num_processors, iflag, myid)
!#######################################################################
!
!  this subroutine evaluates the population, assigns fitness,
!c  establishes the best individual, and outputs information.
      USE ga_mod
      USE mpi_params, ONLY: master
      IMPLICIT NONE
      EXTERNAL fcn
#if !defined(MPI_OPT)
      EXTERNAL ga_fitness_parallel
#endif
      INTEGER :: nopt, n, j, k, kk, iflag, myid
      REAL(rprec), DIMENSION(nopt) :: fvec
      REAL(rprec), DIMENSION(nparmax) :: paramsm,paramav
      INTEGER :: nfev, num_processors
      REAL(rprec) :: fitsum, funcval, fbar, best
      INTEGER :: jstart, jend, istat, jstat
      LOGICAL :: ldiag_opt

      SAVE

      fitSUM = 0
      best=-1.0e30_dp

      ldiag_opt = .false.

      IF (myid .eq. master) WRITE(6,*) 'in ga_evalout',num_processors
         WRITE(6,*) fbar,best,nopt,nfev,iflag
!  ,iflag
!                   fbar,best,nopt,nfev,num_processors,iflag

      DO 29 n=1,nparam
         paramsm(n)=0
 29   CONTINUE
      jstart=1
      jend=npopsiz
      IF(iskip.ne.0) jstart=iskip
      IF(iend.ne.0) jend=iend

      DO  j=jstart,jend

         CALL ga_decode(j,parent,iparent)

!        IF(iskip.ne.0 .and. iend.ne.0 .and. iskip.eq.iend) THEN
!        IF(lscreen) THEN
!        IF(nchrome .le. 120) THEN
!        WRITE(6,1075) j,(iparent(k,j),k=1,nchrome)
!        ELSE
!        WRITE(6,1075) j,(iparent(k,j),k=1,120)
!        WRITE(6,1077) (iparent(k,j),k=121,nchrome)
!        END IF
!        WRITE(6,1076)   (parent(kk,j),kk=1,nparam),0.0
!        END IF
!        END IF
         IF(LDIAG_OPT .and. myid.eq.master) THEN
            IF(nchrome .le. 120) THEN
               WRITE(iunit_ga_out,1075) j,(iparent(k,j),k=1,nchrome)
            ELSE
               WRITE(iunit_ga_out,1075) j,(iparent(k,j),k=1,120)
               WRITE(iunit_ga_out,1077) (iparent(k,j),k=121,nchrome)
            END IF
            WRITE(iunit_ga_out,1076) (parent(kk,j),kk=1,nparam)
         END IF

      END DO
#if defined(MPI_OPT)
         CALL ga_fitness_mpi (jend-jstart+1, f_obj, num_obj,            &
                              fcn, nfev, fitness)
#else
         IF (myid .eq. master) WRITE(6,'(1x,i4,a,i4,a)') jend-jstart+1, &
               ' processes started on ',num_processors, ' processors'

!        flush out buffer before multiprocessing
         CALL flush(6)
         CALL flush(iunit_ga_out)

         CALL multiprocess(jend-jstart+1, num_processors,               &
                           ga_fitness_parallel, fcn )
#endif
        nfev=nfev+jend-jstart+1
        nfit_eval=nfev

!       Clean up...
        iflag=-100
        CALL fcn(nopt, npopsiz, parent(1,jbest), fvec, iflag, nfev)
#if !defined(MPI_OPT)
        DO j=jstart, jend
!
!  Call function evaluator, Write out individual and fitness, and add
!  to the summation for later averaging.
!        iflag=j
!        funcval = ga_evaluate(fcn, nopt, fvec, nparam, parent(1,j),
!    1                      iflag, nfev)

         READ(j+1000, iostat=istat) jstat, iflag
         IF( jstat .ne. j ) THEN
         WRITE(6,*) "wrong INDEX READ in evalout"
         iflag=-14
         EXIT
         END IF

         READ(j+1000, iostat=istat) funcval
         fitness(j)=funcval
         CLOSE(j+1000, status='delete')
        END DO
#endif
        DO 30 j = jstart, jend
         fitsum=fitsum+fitness(j)
         DO 22 n=1,nparam
            paramsm(n)=paramsm(n)+parent(n,j)
 22      CONTINUE

!  Check to see IF fitness of individual j is the best fitness.
         IF (fitness(j).gt.best) THEN
            best=fitness(j)
            jbest=j
            DO 24 k=1,nchrome
               ibest(k)=iparent(k,j)
 24         CONTINUE
         END IF
 30   CONTINUE

!  compute parameter and fitness averages.
      fbar=fitsum/npopsiz
      DO 23 n=1,nparam
         paramav(n)=paramsm(n)/npopsiz
 23   CONTINUE


!  write output information
      IF (myid.eq.master) THEN
      IF (ldiag_opt) THEN
         IF (npopsiz.eq.1) THEN
            IF(nchrome .le. 120) THEN
               WRITE(iunit_ga_out,1075) 1,(iparent(k,1),k=1,nchrome)
            ELSE
            WRITE(iunit_ga_out,1075) 1,(iparent(k,1),k=1,120)
            WRITE(iunit_ga_out,1077) (iparent(k,j),k=121,nchrome)
            END IF
            WRITE(iunit_ga_out,1076)   (parent(k,1),k=1,nparam)
            WRITE(iunit_ga_out,1078)  fitness(1)
            WRITE(iunit_ga_out,*) ' Average Values:'
            WRITE(iunit_ga_out,1275) (parent(k,1),k=1,nparam)
            WRITE(iunit_ga_out,1276) fbar
         ELSE
            WRITE(iunit_ga_out,1275) (paramav(k),k=1,nparam)
            WRITE(iunit_ga_out,1276) (fitness(j),j=1,npopsiz)
         END IF
      END IF
      WRITE(6,1100) fbar
      WRITE(iunit_ga_out,1100) fbar
      WRITE(6,1200) best
      WRITE(iunit_ga_out,1200) best
      END IF

 1075 FORMAT(i3,1x,(120i1))
 1077 FORMAT(3x,1x,(120i1))
 1076 FORMAT(3x,1x,10(1x,e10.4))
 1078 FORMAT(10x,e12.5)
 1100 FORMAT(1x,'Average Function Value of Generation=',e12.5)
 1200 FORMAT(1x,'Maximum Function Value              =',e12.5/)
 1275 FORMAT(/' Average Values:',18x,10(1x,e10.4))
 1276 FORMAT(10x,10e12.5)

      END SUBROUTINE ga_evalout

      MODULE mgrid_mod
      USE v3_utilities
      USE stel_kinds
      USE vmec_input, ONLY: lfreeb
      USE vparams, ONLY: nigroup
      IMPLICIT NONE

      INTEGER, PARAMETER :: nlimset = 2       !number of different limiters
      CHARACTER(LEN=*), PARAMETER ::
     1   vn_br0 = 'br', vn_bp0 = 'bp', vn_bz0 = 'bz',
     2   vn_ar0 = 'ar', vn_ap0 = 'ap', vn_az0 = 'az',
     3   vn_ir = 'ir', vn_jz = 'jz',
     4   vn_kp = 'kp', vn_nfp = 'nfp',
     5   vn_rmin='rmin', vn_rmax='rmax', vn_zmin='zmin',
     6   vn_zmax='zmax', vn_coilgrp='coil_group'
      CHARACTER(LEN=*), PARAMETER ::
     1  vn_nextcur = 'nextcur',  vn_mgmode='mgrid_mode',
     2  vn_coilcur = 'raw_coil_cur',
     3  vn_flp = 'nobser', vn_nobd = 'nobd', vn_nbset = 'nbsets',
     4  vn_nbfld = 'nbfld',
     2  ln_flp = 'flux loops', ln_nobd = 'Connected flux loops',
     3  ln_nbset = 'B-coil loops', ln_next = 'External currents',
     4  ln_nbfld = 'B-coil measurements'

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
!     nr0b, np0b, nz0b:
!              :  grid dimensions of magnetic field in mgrid file
!     nbvac    :  total number of grid points (nr0b*np0b*nz0b) in mgrid file
!     bvac(:,1):  br (radial component of external magnetic field)
!     bvac(:,2):  bp (toroidal component)
!     bvac(:,3) = bz (z-component)
!     rminb, rmaxb : min (max) radial dimension of grid in mgrid
!     zminb, zmaxb : min (max) vertical dimension of grid in mgrid
!
!     nextcur:         no. of EXTERNAL current groups (eg., TF, PF, helical)
!     raw_coil_current  array of raw currents for each coil group
!     mgrid_mode     = 'S', scaled mode; = 'R', raw mode
!     curlabel:   array of labels describing each current group
!                     included in green''s FUNCTION BFIELD response
!
      INTEGER :: nr0b, np0b, nfper0, nz0b
      INTEGER :: nobd, nobser, nextcur, nbfldn, nbsets, nbcoilsn
      INTEGER :: nbvac, nlim, nsets, nrgrid, nzgrid
      INTEGER, DIMENSION(:), ALLOCATABLE :: needflx, nbcoils
      INTEGER, DIMENSION(:), ALLOCATABLE :: limitr, nsetsn
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: iconnect, needbfld
      REAL(rprec) :: rminb, zminb, rmaxb, zmaxb, delrb, delzb
      REAL(rprec) ::rx1, rx2, zy1, zy2, condif
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE, TARGET :: bvac
      REAL(rprec), DIMENSION(:,:,:), POINTER :: brvac, bzvac, bpvac
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: unpsiext,
     1   plbfld, rbcoil, zbcoil, abcoil, bcoil, rbcoilsqr
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: raw_coil_current
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xobser, zobser,
     1   xobsqr, dsiext, psiext, plflux, b_chi
      CHARACTER(LEN=300) :: mgrid_path
      CHARACTER(LEN=300) :: mgrid_path_old = " "
      CHARACTER(LEN=30), DIMENSION(:), ALLOCATABLE :: curlabel
      CHARACTER(LEN=15), DIMENSION(:), ALLOCATABLE ::
     1                                           dsilabel, bloopnames
      CHARACTER(LEN=30) :: tokid
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: dbcoil, pfcspec
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1    rlim, zlim, reslim, seplim
      CHARACTER(LEN=1) :: mgrid_mode

      PRIVATE :: read_mgrid_nc

      CONTAINS

      SUBROUTINE read_mgrid (mgrid_file, extcur, nv, nfp, lscreen,
     1                       ier_flag, comm)
      USE system_mod
      USE mpi_inc
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   V a r i a b l e s
C-----------------------------------------------
!
!     mgrid_file:     full path to mgrid file
!     lscreen   :     logical controlling output to screen
!     ier_flag  ;     error flag returned to caller
!     extcur(n)    :  external current multiplier for bfield(n) components
!     comm      : Optional mpi communicator.
!
      INTEGER, INTENT(out)          :: ier_flag
      INTEGER, INTENT(in)           :: nv, nfp
      LOGICAL, INTENT(in)           :: lscreen
      REAL(rprec), INTENT(in)       :: extcur(:)
      CHARACTER(len=*), INTENT(in)  :: mgrid_file
      INTEGER, INTENT(in), OPTIONAL :: comm
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
#if defined(VMS)
      CHARACTER(LEN=*), PARAMETER :: mgrid_defarea='vmec$:[makegrid]'
#else
      CHARACTER(LEN=*), PARAMETER :: mgrid_defarea='$HOME/vmec/MAKEGRID'
#endif
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
!
!     lgrid_exist  :   logical set if mgrid file is found in given path
!
      INTEGER :: istat, ii
      CHARACTER(LEN=200) :: home_dir
      LOGICAL :: lgrid_exist, lfind
      INTEGER :: mpi_comm
C-----------------------------------------------

      mpi_comm = MPI_COMM_WORLD
      IF (PRESENT(comm)) mpi_comm = comm

      mgrid_path = TRIM(mgrid_file)

      IF ((mgrid_path .eq. TRIM(mgrid_path_old)) .and.
     1     ALLOCATED(curlabel)) THEN
         PRINT *,' mgrid file previously parsed!'
         RETURN
      END IF

      INQUIRE (file=mgrid_path,exist=lgrid_exist,iostat=istat)
      IF (istat.ne.0 .or. .not.lgrid_exist) THEN
          IF (lscreen) PRINT *,' MGRID FILE NOT FOUND IN SPECIFIED ',
     1       'PATH: SEARCHING DEFAULT AREA'
          ii = INDEX(mgrid_file,'/',back=.true.)
          istat = INDEX(mgrid_defarea, '$HOME')
          IF (istat .ne. 0) THEN
             CALL getenv('HOME', home_dir)
             IF (istat .gt. 1) THEN
                home_dir = mgrid_defarea(1:istat-1) // TRIM(home_dir)
     1                   // mgrid_defarea(istat+5:)
             ELSE
                home_dir = TRIM(home_dir) // mgrid_defarea(istat+5:)
             END IF
          ELSE
             home_dir = mgrid_defarea
          END IF
          mgrid_path = TRIM(home_dir) // mgrid_file(ii+1:)
          INQUIRE (file=mgrid_path,exist=lgrid_exist,iostat=istat)
      END IF

      mgrid_path_old = mgrid_path

      ier_flag = 0

      IF (lgrid_exist) THEN
         IF (lscreen) PRINT '(2x,2a)',
     1     'Opening vacuum field file: ', TRIM(mgrid_file)
!
!        Parse mgrid file name, look for .nc extension (netcdf format)
!
         ii = LEN_TRIM(mgrid_path) - 2
         lfind = (mgrid_path(ii:ii+2) == '.nc')
         IF (lfind) THEN
            CALL read_mgrid_nc (mgrid_path, extcur, nv, nfp,
     1                          ier_flag, lscreen, comm)
         END IF

!SPH060517         IF (np0b .ne. nv) THEN
         IF (nv.EQ.0 .OR. MOD(np0b, nv).NE.0) THEN
            PRINT *,' NZETA=',nv,
     1      ' DOES NOT DIVIDE EVENLY INTO NP0B=',np0b,' IN MGRID FILE'
            ier_flag = 9
         ELSE IF (nfper0.ne.nfp) THEN
            PRINT *,' NFP(READ in) = ',nfp,' DOES NOT AGREE WITH ',
     1      'NFPER (in vacuum field file) = ',nfper0
            ier_flag = 9
         END IF

      END IF

      IF (ier_flag .ne. 0) RETURN

      IF (.not.lgrid_exist .or. ier_flag.ne.0) THEN
         lfreeb = .false.
         IF (lscreen) THEN
            PRINT *, ' Error opening/reading mgrid file in dir: ',
     1                TRIM(home_dir)
            PRINT *, ' User must supply vacuum bfield in mgrid to ',
     1                'run vmec in free-boundary mode!'
            PRINT *, ' Proceeding to run vmec in',
     1               ' fixed boundary mode'
         END IF
      END IF

      END SUBROUTINE read_mgrid

!
!     PARALLEL MPI MODIFICATIONS ADDED BY Mark R. Cianciosa <cianciosamr@ornl.gov>, 092315
!
      SUBROUTINE read_mgrid_nc (filename, extcur, nv, nfp,
     1                          ier_flag, lscreen, comm)
      USE ezcdf
      USE mpi_inc
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y  A r g u m e n t s
C-----------------------------------------------
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER, INTENT(in)          :: nv, nfp
      REAL(rprec), INTENT(in)      :: extcur(:)
      INTEGER, INTENT(in)          :: comm
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE ::
     1                       brtemp, bztemp, bptemp
      INTEGER :: ier_flag, ngrid
      INTEGER :: istat, ig
      LOGICAL :: lscreen
      INTEGER, DIMENSION(3)   :: dimlens
      CHARACTER(LEN=100) :: temp

      INTEGER :: mpi_rank, mpi_size, lMPIInit, MPI_ERR

      CALL MPI_INITIALIZED(lMPIInit, MPI_ERR)
      IF (lMPIInit.NE.0) THEN
         CALL MPI_COMM_RANK(comm, mpi_rank, istat)
         CALL MPI_COMM_SIZE(comm, mpi_size, istat)
      ELSE
         mpi_rank = 0; mpi_size = 1
      END IF

C-----------------------------------------------
      call cdf_open(ngrid, filename,'r', istat)
      IF (istat .ne. 0) THEN
         ier_flag = 9
         RETURN
      END IF

!
!     READ IN DATA
!
      CALL cdf_read(ngrid, vn_ir, nr0b)
      CALL cdf_read(ngrid, vn_jz, nz0b)
      CALL cdf_read(ngrid, vn_kp, np0b)
      CALL cdf_read(ngrid, vn_nfp, nfper0)

      IF (nfper0.NE.nfp .OR. MOD(np0b, nv).NE.0) RETURN

      CALL cdf_read(ngrid, vn_nextcur, nextcur)

      IF (nextcur .eq. 0) THEN
        PRINT *,' NEXTCUR = 0 IN READING MGRID FILE'
        ier_flag = 9
        RETURN
      ELSE IF (nextcur .gt. nigroup) THEN
        PRINT *,' NEXTCUR > NIGROUP IN MGRID FILE'
        ier_flag = 9
        RETURN
      END IF

      CALL cdf_read(ngrid, vn_rmin, rminb)
      CALL cdf_read(ngrid, vn_zmin, zminb)
      CALL cdf_read(ngrid, vn_rmax, rmaxb)
      CALL cdf_read(ngrid, vn_zmax, zmaxb)

      delrb = (rmaxb-rminb)/(nr0b-1)
      delzb = (zmaxb-zminb)/(nz0b-1)

      CALL cdf_inquire(ngrid, vn_coilgrp, dimlens)
      IF (.NOT. ALLOCATED(curlabel)) THEN
         ALLOCATE (curlabel(nextcur), stat=istat)
      ELSE IF (SIZE(curlabel) .ne. nextcur) THEN
         DEALLOCATE (curlabel)
         ALLOCATE (curlabel(nextcur), stat=istat)
      END IF
!THIS IS A GLITCH WITH cdf_read: must distinguish 1D char array from multi-D
      IF (nextcur .eq. 1) THEN
         IF (istat .eq. 0)
     1     CALL cdf_read(ngrid, vn_coilgrp, curlabel(1))
      ELSE IF (istat .eq. 0) THEN
           CALL cdf_read(ngrid, vn_coilgrp, curlabel(1:nextcur))
      END IF

      IF (istat .ne. 0) STOP 'Error allocating CURLABEL in mgrid_mod'
!
!     READ 3D Br, Bp, Bz ARRAYS FOR EACH COIL GROUP
!
      nbvac = nr0b*nz0b*nv
      IF (.NOT. ALLOCATED(bvac)) THEN
         ALLOCATE (bvac(nbvac,3), stat=istat)
      ELSE IF (SIZE(bvac,1) .ne. nbvac) THEN
         DEALLOCATE (bvac);  ALLOCATE(bvac(nbvac,3), stat=istat)
      END IF
      IF (istat .ne. 0) STOP 'Error allocating bvac in mgrid_mod'

      bvac = 0

      ALLOCATE (brtemp(nr0b, nz0b, np0b),                                      &
     &          bptemp(nr0b, nz0b, np0b),                                      &
     &          bztemp(nr0b, nz0b, np0b), stat=istat)
      IF (istat .ne. 0)STOP 'Error allocating bXtemp in mgrid_mod'
      brtemp = 0
      bptemp = 0
      bztemp = 0


      DO ig = mpi_rank + 1, nextcur, mpi_size

         WRITE (temp, 1000) vn_br0, ig
         CALL cdf_read(ngrid, temp, brtemp)

         WRITE (temp, 1000) vn_bp0, ig
         CALL cdf_read(ngrid, temp, bptemp)

         WRITE (temp, 1000) vn_bz0, ig
         CALL cdf_read(ngrid, temp, bztemp)

!
!        STORE SUMMED BFIELD (OVER COIL GROUPS) IN BVAC
!
         CALL sum_bfield(bvac(1,1), brtemp, extcur(ig), nv)
         CALL sum_bfield(bvac(1,2), bptemp, extcur(ig), nv)
         CALL sum_bfield(bvac(1,3), bztemp, extcur(ig), nv)
      END DO

	  np0b = nv

      IF (lMPIInit.NE.0) THEN
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, bvac, SIZE(bvac), MPI_REAL8,  &
     &                      MPI_SUM, comm, istat)
         CALL assert_eq(istat,0,'MPI_ALLREDUCE failed in read_mgrid_nc')
      END IF

      CALL cdf_inquire(ngrid, vn_mgmode, dimlens, ier=istat)
      IF (istat .eq. 0) THEN
         CALL cdf_read(ngrid, vn_mgmode, mgrid_mode)
      ELSE
         mgrid_mode = 'N'
      END IF

      CALL cdf_inquire(ngrid, vn_coilcur, dimlens, ier=istat)
      IF (istat .eq. 0) THEN
	   IF (ALLOCATED(raw_coil_current)) DEALLOCATE(raw_coil_current)
         ALLOCATE (raw_coil_current(nextcur), stat=istat)
         IF (istat .ne. 0) STOP 'Error allocating RAW_COIL in mgrid_mod'
         CALL cdf_read(ngrid, vn_coilcur, raw_coil_current)
      END IF

      CALL cdf_close(ngrid)

      IF (ALLOCATED(brtemp))
     1    DEALLOCATE (brtemp, bptemp, bztemp)

      CALL assign_bptrs(bvac)

1000  FORMAT(a,'_',i3.3)

      END SUBROUTINE read_mgrid_nc

      SUBROUTINE sum_bfield(bfield, bf_add, cur, nv)
	  INTEGER, INTENT(IN)        :: nv
      REAL(rprec), INTENT(INOUT) :: bfield(nr0b*nz0b,nv)
      REAL(rprec), INTENT(IN)    :: bf_add(nr0b*nz0b,np0b)
	  INTEGER     :: nskip
      REAL(rprec) :: cur

      nskip = np0b/nv
      bfield = bfield + cur*bf_add(:,1:np0b:nskip)

      END SUBROUTINE sum_bfield

      SUBROUTINE assign_bptrs(bptr)
      IMPLICIT NONE
      REAL(rprec), TARGET, INTENT(in) :: bptr(nr0b,nz0b,np0b,3)

      brvac => bptr(:,:,:,1)
      bpvac => bptr(:,:,:,2)
      bzvac => bptr(:,:,:,3)

      END SUBROUTINE assign_bptrs

      SUBROUTINE free_mgrid (istat)
      INTEGER :: istat

      istat = 0

      IF (ALLOCATED(bvac)) DEALLOCATE (bvac,stat=istat)
      IF (ALLOCATED(xobser))
     1   DEALLOCATE (xobser, xobsqr, zobser, unpsiext, dsiext,
     2      psiext,plflux, iconnect, needflx, needbfld, plbfld,
     3      nbcoils, rbcoil, zbcoil, abcoil, bcoil, rbcoilsqr, dbcoil,
     4      pfcspec,dsilabel, bloopnames, curlabel, b_chi, stat=istat)
      IF (ALLOCATED(raw_coil_current)) DEALLOCATE(raw_coil_current)

      IF (ALLOCATED(rlim))
     1   DEALLOCATE (rlim,zlim, reslim,seplim,stat=istat)

!  Reset mgrid_path_old, so that can reread an mgrid file. SL, JDH 2012-07-16
      mgrid_path_old = " "

      END SUBROUTINE free_mgrid

      END MODULE mgrid_mod

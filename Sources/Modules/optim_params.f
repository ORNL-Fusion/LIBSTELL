      MODULE optim_params
      USE vparams, ONLY: rprec, dp, nsd, ntord, ntor1d, mpol1d
      USE vsvd0, ONLY: nigroup
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NumJstard = 10
      INTEGER, PARAMETER :: ini_max=100                                  !! COBRA
      REAL(rprec) :: bigno = 1.e10_dp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: niter_opt, NumJstar, NumJinvariant,
     1   mboz_opt, nboz_opt,  nbmn, nproc, jboot                         !LPK
      INTEGER :: num_processors, num_levmar_params
      INTEGER, DIMENSION(20) :: n_jac, m_jac
      INTEGER, DIMENSION(20) :: n_vac_island, m_vac_island
      REAL(rprec) :: r00_scale, b00_scale, epsfcn, rgrid_min, rgrid_max,
     1               zgrid_min, zgrid_max, r00_opt
      REAL(rprec) :: sigma_jstar(nsd, NumJstard),
     1   sigma_jinvariant(nsd, NumJstard)
      REAL(rprec), DIMENSION(nsd) ::
     1   sigma_iota, sigma_mercier, sigma_vp
      REAL(rprec), DIMENSION(nsd) :: sigma_bmin,
     1   sigma_bmax, sigma_ripple, sigma_bmn, nsurf_mask
      REAL(rprec) :: sigma_aspect, sigma_ellipticity, sigma_maxcurrent,
     1   sigma_coil_complex, sigma_curv, sigma_beta, target_aspectratio,
     2   target_maxcurrent, target_beta, target_rmax, target_rmin,
     3   target_ellipticity, sigma_rmax, sigma_rmin, sigma_iota_max,
     4   sigma_iota_min, target_iota_max, target_iota_min
     5   ,sigma_centering                                                 !!Obsolete, replaced by sigma_rmax, sigma_rmin
      REAL(rprec) :: coil_separation, phiedge_min, phiedge_max
      REAL(rprec), DIMENSION(0:10) :: target_iota, target_well,
     1   at, aseedcur, fboot                                             !!LPK
      REAL(rprec) :: sigma_bal, sigma_boot, zeff_boot,                   !!LPK
     1    sigma_fluxp, target_fluxp, sigma_zmax, target_zmax,
     2    sigma_pseudo, sigma_pseudo2, sigma_rbtor, target_rbtor
      REAL(rprec), DIMENSION(9) :: sigma_kink, target_kink
      REAL(rprec), DIMENSION(20) :: sigma_jac, sigma_vac_island

      COMPLEX :: helicity
      LOGICAL :: lfix_ntor(-ntord:ntord), 
     1           lfix_rhob(-ntord:ntord,0:mpol1d+1)
      LOGICAL, DIMENSION(nsd) :: lsurf_mask, ldkes_mask
      LOGICAL :: lextcur(nigroup), lcoilp_sep
      LOGICAL :: lreset_opt, lbmn, lcoil_complex, ledge_current,
     1   lcur_prof_opt, liota_prof_opt, lbootsj_opt, lj_star,
     2   lj_invariant, lcoil_opt, lcoil_geom, laspect_max, lbeta_min,
     3   lvac_opt, lphiedge,
     4   lbal_opt, lbootstrap, lseedcur, lkink_opt, lnescoil_opt         !!LPK
     5  ,lcurprof_opt, lprof_opt, lpress_opt, lboundary                  !!Obsolete, replaced by lcur_prof_opt, liota_prof_opt, lpres_prof_opt
      LOGICAL :: lreconp, lreconj, lp1zero, lj1zero
      INTEGER :: kpp, kjj
      CHARACTER(len=100) :: seq_ext, opt_ext
      CHARACTER(len=200) :: v3rfun_dir, v3post_in
      LOGICAL :: lv3post
      CHARACTER(len=4) :: sym_type

      REAL(rprec) :: target_coil_complex, target_coil_jmax               !!LPK
      REAL(rprec) :: sigma_coil_jmax, sigma_berr_ave                     !!LPK, SPH

      REAL(rprec), DIMENSION(nsd) :: sigma_neo, nneo_mask                !! NEO
      LOGICAL, DIMENSION(nsd) :: lneo_mask                               !! NEO
      LOGICAL :: lneo_opt                                                !! NEO
      REAL(rprec), DIMENSION(nsd) :: sigma_dsubr
      LOGICAL :: ldsubr_opt
      REAL(rprec) :: sigma_orbit
      LOGICAL :: lorbit_opt
      INTEGER :: nopt_alg, nopt_boundary
      LOGICAL :: ldiag_opt, lkeep_mins                                   !! Diagnostic output

      REAL(rprec) :: sigma_kappa, target_kappa, sigma_oh
      REAL(rprec), DIMENSION(nigroup) :: sigma_extcur, oh_coefs
      REAL(rprec), DIMENSION(0:10) :: target_iota_p
      REAL(rprec), DIMENSION(nsd) :: sigma_iota_pmax, sigma_iota_pmin

      INTEGER :: nini_theta, nini_zeta, nini_tot                         !! COBRA
      REAL(rprec), DIMENSION(nsd) :: sigma_balloon, sigma_pgrad,
     1                               target_balloon                      !! COBRA
      REAL(rprec), DIMENSION(ini_max):: bal_theta0, bal_zeta0            !! COBRA

      REAL(rprec) :: sigma_pedge(1)                                      !! COBRA
      REAL(rprec) :: sigma_bootsj(nsd)
      LOGICAL :: lballoon_mask(nsd)                                      !! COBRA
      LOGICAL :: lballoon_opt, lpres_prof_opt                            !! COBRA

      REAL(rprec) :: nballoon_mask(nsd)                                  !! VMECCOBRA (RS)
      LOGICAL :: l_legendre                                              !! LEGENDRE (RS)

      REAL(rprec), DIMENSION(nsd) :: ndkes_mask, dkes_nu, dkes_efield,
     1     sigma_dkes                                                    !! RHF
      LOGICAL :: ldkes_opt                                               !! RHF

      REAL(rprec) :: sigma_vv, sigma_vv_rms, vv_dist, vv_dist_rms,
     1               target_vv, target_vv_rms                            !! RH & MZ
      INTEGER :: mpol_vv, ntor_vv, nu_vv, nv_vv
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d) ::
     1   rbc_vv, zbs_vv
      LOGICAL :: shapeweight
CEAL    !       deviation weighting defaults
      REAL(rprec) :: theta0_bw(3), phi0_bw,
     1   wtheta_bw, wphi_bw, amplw_bw(3),  planes_bw(3)

      REAL(rprec) :: sigma_diagno(1000), data_diagno(1000)
      CHARACTER(len=30) :: name_diagno(1000)

      NAMELIST /optimum/ epsfcn, niter_opt, num_processors,             & 
     &   num_levmar_params, nsurf_mask, nopt_alg, nopt_boundary,        &
     &   lreset_opt, ldiag_opt, lkeep_mins, lbmn, lj_star, laspect_max, &
     &   lbeta_min, lj_invariant, liota_prof_opt, lcur_prof_opt,        &
     &   ledge_current, lphiedge, lbootsj_opt, lkink_opt, lballoon_opt, &
     &   l_legendre, ldkes_opt, lneo_opt, ldsubr_opt, lorbit_opt,       &
     &   lpres_prof_opt, lnescoil_opt, lcoil_geom, lv3post, lvac_opt,   &
     &   lfix_ntor, lfix_rhob, lextcur, sigma_extcur, oh_coefs,sigma_oh,&
     &   r00_opt, r00_scale, b00_scale, rgrid_min, rgrid_max,           &
     &   zgrid_min, zgrid_max, mboz_opt, nboz_opt, phiedge_min,         &
     &   phiedge_max, coil_separation, target_aspectratio, sigma_aspect,&
     &   target_beta, sigma_beta, target_kink, sigma_kink,              &
     &   target_maxcurrent, sigma_maxcurrent, target_rmax, sigma_rmax,  &
     &   target_rmin, sigma_rmin, target_zmax, sigma_zmax, target_kappa,&
     &   sigma_kappa,target_ellipticity,sigma_ellipticity, target_fluxp,&
     &   sigma_fluxp, target_rbtor, sigma_rbtor, target_coil_complex,   &
     &   sigma_coil_complex, target_coil_jmax, sigma_coil_jmax,         &
     &   target_iota, sigma_iota, target_iota_p, sigma_iota_pmax,       &
     &   sigma_iota_pmin, target_iota_min, sigma_iota_min,              &
     &   target_iota_max, sigma_iota_max, sigma_curv, sigma_berr_ave,   &
     &   sigma_pseudo, sigma_pseudo2, sigma_mercier, sigma_jac, n_jac,  &
     &   m_jac, n_vac_island, m_vac_island, sigma_vac_island,  helicity,&
     &   sigma_bmin, sigma_bmax, sigma_bmn, sigma_ripple, NumJstar,     &
     &   NumJInvariant, sigma_jstar, sigma_jinvariant, nballoon_mask,   &
     &   target_balloon, sigma_balloon, sigma_pgrad, sigma_pedge,       &
     &   bal_theta0, bal_zeta0, fboot, aseedcur,sigma_bootsj,ndkes_mask,&
     &   dkes_nu, dkes_efield, sigma_dkes, nneo_mask, sigma_neo,        &
     &   sigma_dsubr, sigma_orbit, v3rfun_dir, v3post_in, name_diagno,  &
     &   data_diagno, sigma_diagno, target_vv, sigma_vv, target_vv_rms, &
     &   lreconp, lp1zero, kpp, lreconj, lj1zero, kjj,                  &
     &   sigma_vv_rms, mpol_vv, ntor_vv, nu_vv, nv_vv, rbc_vv, zbs_vv,  &
     &   shapeweight, planes_bw, amplw_bw, theta0_bw,phi0_bw, wtheta_bw,&
     &   wphi_bw,  
     &   lcoil_complex, lballoon_mask, ldkes_mask, lseedcur, lneo_mask  &
     &   ,lsurf_mask, lcurprof_opt, lprof_opt, target_well, lcoil_opt   & !!Obsolete, retain for consistency
     &   ,sigma_centering, lpress_opt, lbootstrap, lbal_opt, lboundary  &
     &   ,sigma_vp, sym_type, sigma_boot, sigma_bal, nproc, nbmn        &
     &   ,zeff_boot, jboot, at
!
!        VARIABLE DESCRIPTIONS
!------------------------------------
!        SCALARS/VECTORS
!------------------------------------
!        nopt_alg             determines the algorithm for optimization
!                             0 = Levenberg-Marquardt (default)
!                             1 = Genetic (GA)
!                             2 = Differential Evolution (DE)
!        nopt_boundary        determines internal representation used for optimization 
!                             in fixed-boundary runs (LFREEB=.f)
!                             0 = Hirshman-Breslau (default)
!                             1 = Advanced HB
!                             2 = Garabedian delta_mn
!        niter_opt            Maximum number of optimization iterations
!        num_processors       Maximum number of processors to try and use simultaneously
!                             (obsolescent: nproc)
!        num_levmar_params
!        nproc
!        epsfcn               'Annealing' parameter, 1.E-2 or less. If
!                             lreset_opt=F, epsfcn <= 1.E-3; if lreset_opt=T,
!                             epsfcn <= 1.E-4 is recommended
!        r00_scale 
!        b00_scale
!        r00_opt
!        mboz_opt             User-input number of poloidal boozer modes
!        nboz_opt             User-input number of toroidal boozer modes
!        numJstar
!        numJInvariant
!        sym_type 
!        m_vac_island, n_vac_island 
!        m_jac, n_jac
!        jboot
!        fboot
!        at
!        aseedcur
!        dkes_nu
!        dkes_efield
!        mpol_vv, ntor_vv
!        coil_separation,
!        bal_theta0, bal_zeta0
!        rgrid_min, rgrid_max
!        zgrid_min, zgrid_max
!        nbmn
!        zeff_boot
!        oh_coefs
!        rbc_vv, zbs_vv
!        theta0_bw, phi0_bw, amplw_bw  
!        planes_bw
!------------------------------------
!        CONTROL-MASK ARRAYS
!------------------------------------
!        nsurf_mask           Real array (size=nrad) giving the fractional radial s-values
!                             at which optimizations are to be done for bmin, bmax, Jinvariant, 
!                             DKES, NEO calculations
!        lsurf_mask (obs)     Logical array (size=nrad); use nsurf_mask 
!        nballoon_mask        Real array (size=nrad) containing fractional radial s-values 
!                             where ballooning growth rates are to be computed
!        lballoon_mask (obs)  Logical array (size=nrad) designating surfaces where
!                             ballooning growth rates are to be computed (=T); use nballoon mask
!        lextcur              Logical array (size>=nextcur); if (i-th) element =T and lfreeb=True, 
!                             extcur(i) will be varied as an independent variable.
!        ndkes_mask
!        ldkes_mask
!        nneo_mask
!        lneo_mask
!------------------------------------
!        LOGICAL CONTROL VARIABLES
!------------------------------------
!        lreset_opt           =F, VMEC runs without resetting (after convergence)
!                             to a coarse radial grid
!                             =T (default), VMEC resets after each new iteration
!        lcur_prof_opt        =T, ncurr set to 1 and ac current series expansion coefficients
!                             are varied as independent variables
!                             =F (default), ac coefficients are fixed (if ncurr = 1)
!        liota_prof_opt       =T, ncurr set to 0 and ai iota series expansion coefficients 
!                             are varied as independent variables
!        ledge_current        =T, vary curtor (edge current), provided ncurr=1
!        lphiedge             =T, vary phiedge (edge toroidal flux)
!        lprof_opt (obs)      (see lcur_prof_opt and liota_prof_opt)
!        lcurprof_opt (obs)   (see lcur_prof_opt)
!        lpres_prof_opt       =T  optimize pressure profile shape (vary mass expansion coefficients
!                             as independent variables); used to achieve ballooning stability
!                             =F  keep pressure profile shape fixed
!        lpress_opt (obs)     (see lpres_prof_opt)
!        laspect_max          =T, then sigma_aspect = bigno if aspect ratio <= target_aspectratio.
!                             used to guarantee a no larger than the target.
!        lbeta_min            =T, then sigma_beta = bigno if beta >= Target_Beta. Used to
!                             guarantee a minimum beta
!        lkink_opt            =T, do global kink stability calculation
!        lballoon_opt         =T, do ballooning calculation
!                             =F, do not do ballooning calculation if nballoon_mask prescribed
!        lbal_opt (obs)       (use nballoon_mask)
!                             =F (or lfreeb=False), extcur(i) is fixed during optimization
!        lnescoil_opt         =T, evaluate NESCOIL coil optimization targets
!
!        lfix_ntor(n)         Logical array. =F, r0n's, z0n's are free to vary (i.e., they ARE NOT fixed)
!                             for the toroidal index 'n' of the arrays rbc(n,m=0), zbs(n,m=0)
!                             =T, r0n's, z0n's are fixed (not allowed to varied) during optimization
!                             Useful if certain n-s are to be kept fixed (such as the
!                             axi-symmetric boundary components: lfix_ntor(0) = T)
!                             NO EFFECT IN FREE-BDY OPTIMIZATION
!        lfix_rhob(n,m)       2D logical array. =F, rhobc(n,m) component is varied during optimization (default)
!                             =T, rhobc(n,m) compoents is fixed during optimization
!                             NO EFFECT IN FREE-BDY OPTIMIZATION
!        lbmn                 =T, add Boozer spectra target to chi-sq (for targetting symmetries: QA, QH, QP)
!        lbootsj_opt
!        lbootstrap (obs)     (see lbootsj_opt)
!        lseedcur
!        lj_star
!        lj_invariant
!        l_legendre
!        lkeep_mins
!        lcoil_complex
!        lcoil_opt            (obsolete)
!        lcoil_geom           =T, do coil geometry optimization
!                             (coil harmonics are independent variables in the optimization)
!                             Note that in this case, the initial coil shapes
!                             and the weights for the coil evaluations targets
!                             (e.g. coil curvatures and separations) are
!                             specified in the COILSIN NAMELIST (in coilsnamin.f in LIBSTELL).
!
!                             =F, coil geometry is FIXED and either the
!                             EXTERNAL currents are varied (lfreeb=.true., lextcur)
!                             or the bdy coefficients (rbc, zbs) are varied.
!        lcoilp_sep           a COMPUTED variable, which = F if lcoil_geom = T and
!                             MAY = T when lcoil_geom=F. Used to force computation of
!                             minimum coil-plasma separation even when the coil geometry is not
!                             evolving in the optimizer
!        lvac_opt             =T, use rotations and shifts of coils for independent optimization variables
!        ldkes_opt
!        lneo_opt
!        ldiag_opt
!        ldsubr_opt
!        lorbit_opt
!        phiedge_max,min      max/min allowable range for variation of phiedge, when lphiedge = TRUE
!------------------------------------
!        TARGET/WEIGHTS (SIGMA = 1/WEIGHT)
!------------------------------------
!        Equilibrium and Geometry
!------------------------------------
!        Target_AspectRatio   Aspect Ratio 
!        sigma_aspect         
!        Target_MaxCurrent    Maximum integrated toroidal current
!                             (bounding current, matched at each radial position)
!        sigma_maxcurrent
!        Target_Beta          Volume averaged beta to match
!        sigma beta
!        Target_Iota          Coefficients of power series (in flux s)
!                             for iota profile
!        sigma_iota
!        Target_Iota_P        Coefficients of power series (in s) for iota-prime
!                             = d(iota)/ds
!        sigma_iota_pmax      Array (size nrad) of sigmas for target_iota_p as
!                                   upper bound on d-iota/ds as
!        sigma_iota_pmin      Array (size nrad) of sigmas for target_iota_p as
!                                   lower bound on d-iota/ds as
!        Target_iota_min      minimum iota value (as lower bound)
!        sigma_iota_min
!        Target_iota_max      maximum iota value (as upper bound)
!        sigma_iota_max
!        Target_rmin          minimum major radius over the boundary
!        sigma_rmin           
!        Target_rmax          maximum major radius over the boundary
!        sigma_rmax           
!        Target_zmax          maximum height over the boundary
!        sigma_zmax           
!        Target_ellipticity   desired elongation of phi=0 cross section
!        sigma_ellipticity
!        Target_kappa         desired <kappa> (n=0 component of elongation)
!        sigma_kappa          
!        Target_fluxp         Minimum poloidal flux (Wb)
!        sigma_fluxp          
!        sigma_curv           Forces curvature kurtosis to zero at phi=0,90,180,270
!
!        Stability
!------------------------------------
!        Target_Well          Coefficients Power series (in flux s)
!                             for the magnetic well
!                             [Vp(s) - Vp(0)]/Vp(0) [Replaced with Mercier criterion]
!        sigma_vp             Array (nrad) of sigmas for well
!        sigma_Mercier        Sigma for Mercier stability (not relative)
!        Target_balloon       Real array (size=nrad) of ballooning eigenvalue (=0 for marginal stability)
!        sigma_balloon        Real array (size=nrad) 
!        sigma_pgrad          Real array (size=nrad) forcing local pressure gradient to zero
!        sigma_pedge          Real array (size=1) for forcing pressure value at edge to zero
!        Target_kink          array of kink eigenvalues. Can be used to bias the chi-square
!                             to achieve stability
!        sigma_kink           Array of sigmas for full kink stability calc.,
!                             The first sigma value is used in a call to xtprp/ft5tpr
!                             the second value is used for xtprp2/ft5tpr2,and so-on for successive values.  
!                             This allows targetting of multiple mode families. The sigma values are  *not relative*
!
!        Coils Currents
!------------------------------------
!        sigma_extcur         Array (dim=nextcur) of sigmas for the coil-currents used for regularizing 
!                             (forcing to zero) the coil currents in the extcur array for which the mask array 
!                             lextcur = T
!        Target_rbtor         The effective poloidal coil currents (R*Bt) used to constrain sum of varied external 
!                             coil currents ([T] - [m]). Only needed if lcoil_geom=F (fixed coil geometry); otherwise,
!                             handled by xcoilgeom.
!        sigma_rbtor  
!        oh_coefs             Array (size=nextcur) for imposing a linear sum constraint on the coil currents, 
!                             e.g. for ensuring that a PF set does not generate net poloidal flux.
!                             Constraint is sum_i (oh_coefs(i)*extcur(i))
!        sigma_oh             sigma for the OH flux constraint (weighted linear sum
!                             coil currents, targetted to zero)
!
!        NESCOIL current-sheet (lnescoil_opt=T)
!------------------------------------
!        Target_coil_complex  desired maximum coil complexity measure (<M> for coils on sheet)
!        sigma_coil_complex   
!        Target_coil_jmax     desired maximum coil current density
!        sigma_coil_jmax
!        sigma_berr_ave       forces NESCOIL berr to zero
!
!        COILOPT
!------------------------------------
!
!        Transport (see lfix_ntor, lfix_rhob, lbmn described above)
!------------------------------------
!        helicity             If lbmn=T, used to select helicity of Boozer spectra that
!                             influence symmetry hence transport properties.
!                             For kh = real(helicity), lh = imag(helicity), then
!                             kh = 0         =>   quasi-poloidal symmetry
!                             lh = 0         =>   quasi-axisymmetry
!                             m*kh+n*lh = 0  =>  quasi-helical symmetry (lu + kv)
!        sigma_bmn            Real array (size=nrad) for forcing specific bmns to zero
!                             (relative to largest bmn). Used to obtain quasi-symmetric spectra
!                             in conjunction with helicity specification.
        

!        target_iota_max, sigma_iota_max
!        target_iota_min, sigma_iota_min
!        target_vv, sigma_vv
!        target_vv_rms, sigma_vv_rms
!        sigma_bootsj (sigma_boot)
!        sigma_pseudo
!        sigma_pseudo2
!        sigma_bal
!        sigma_vac_island,
!        sigma_jstar 
!        sigma_jinvariant
!        sigma_jac
!        sigma_oh         
!        sigma_bmin
!        sigma_bmax 
!        sigma_ripple
!        sigma_curv
!        sigma_pgrad 
!        sigma_pedge
!        sigma_neo
!        sigma_dsubr
!        sigma_orbit
!        sigma_dkes
!        shapeweight
!        wtheta_bw
!        wphi_bw
!        sigma_centering
!
!        Magnetic Diagnostics
!------------------------------------
!        lv3post              =T, match diagnostic signals, =F, ignore diagnostic signals
!                             in optimization
!        v3rfun_dir           fully-qualified path to directory containing pre-computed
!                             response tables for the diagnostic set used for matching
!        v3post_in            fully-qualified name of v3post input namelist file

      CONTAINS

      SUBROUTINE read_optimum_namelist (iunit, istat)
      INTEGER, INTENT(in)    :: iunit
      INTEGER, INTENT(inout) :: istat

      istat = 0

!
!     DEFAULT VALUES
!
      r00_opt = -1
      r00_scale = 1
      b00_scale = 1
      kpp = 10 ; kjj = 10
      lreconp = .false.  ;  lreconj = .false.
      lp1zero = .false.  ;  lj1zero = .false. 

      mboz_opt = 0; nboz_opt = 0; nbmn=0; nproc=0
      opt_ext = 'none'
      v3rfun_dir = ' '
      v3post_in = ' '
      lbeta_min = .false.
      laspect_max = .false.
      lcoil_geom = .false.
      lvac_opt = .false.
      lv3post = .false.
      coil_separation = 0
      NumJstar = 0
      NumJinvariant = 0
      target_iota_max = 1
      target_iota_min = 1
      sigma_iota_max = bigno
      sigma_iota_min = bigno
      sigma_iota = bigno
      sigma_vp   = bigno
      sigma_mercier = bigno
      target_kink = 0
      sigma_kink = bigno
      sigma_curv = bigno
      sigma_aspect = bigno
      sigma_coil_complex = bigno
      sigma_MaxCurrent = bigno
      sigma_beta = bigno
      sigma_centering = bigno
      sigma_rmin = bigno
      sigma_rmax = bigno
      sigma_zmax = bigno
      sigma_ellipticity = bigno
      sigma_bmin = bigno
      sigma_bmax = bigno
      sigma_bmn  = bigno
      sigma_jstar = bigno
      sigma_jinvariant = bigno
      sigma_ripple = bigno
      sigma_bootsj = bigno
      sigma_boot   = bigno
      at = 0
      sigma_jac = bigno
      sigma_vac_island = bigno
      sigma_diagno = bigno
      data_diagno = 0
      name_diagno = ""

!
!     INITIALIZATION OBSOLETE VARS, TOO
!
      lcurprof_opt = .false.; lcoil_opt = .false.; lboundary = .false.

      num_levmar_params = 1
      num_processors =    1

!
!     LPK ADDITIONS
!
      sigma_fluxp = bigno
      sigma_pseudo = bigno
      sigma_pseudo2 = bigno
      nproc = 1
      lkink_opt = .false.
      lbal_opt = .false.
      lbootstrap = .false.
      jboot = 0
      fboot = 0
      fboot(0) = 1
      zeff_boot = 1
      lseedcur = .false.
      lpress_opt = .false.
      lcoil_complex = .false.
      lnescoil_opt = .false.
      sigma_coil_jmax = bigno
      sigma_berr_ave = bigno
      target_coil_complex = 1
      target_coil_jmax = 0
      ldkes_opt = .false.
      ldkes_mask = .false.
      sigma_dkes = bigno
      ndkes_mask = 0
      dkes_efield = 0
      dkes_nu = 0.001_dp
      aseedcur = 0
      Target_RBtor = 1
      Sigma_RBtor = bigno
!!
!!    NEO
!!
      sigma_neo = bigno
      lneo_opt = .false.
      nneo_mask = 0
      lneo_mask = .false.
      ldiag_opt = .false.
      lkeep_mins = .false.
      ldsubr_opt = .false.
      sigma_dsubr = bigno
      lorbit_opt = .false.
      sigma_orbit = bigno
!!
!!    END LPK ADDITIONS
!!
      Target_AspectRatio = 3
      Target_MaxCurrent = 0
      Target_Iota = 0
      Target_Well = (/ 0.0_dp, -0.39_dp, 0.19_dp, (0.0_dp,istat=1,8) /)
      niter_opt = 1
      nopt_alg = 0
      nopt_boundary = 0
      lreset_opt = .true.
      lcur_prof_opt = .false.
      liota_prof_opt = .false.
      ledge_current = .false.
      lphiedge = .false.
      lbmn = .false.
      lj_star = .false.
      lj_invariant = .false.
      lprof_opt = .false.
      lextcur = .false.
      lbootsj_opt = .false.
      lsurf_mask = .false.
      lfix_ntor = .false.
      lfix_rhob = .false.
      phiedge_max = 0;  phiedge_min = 0
      nsurf_mask = 0
      epsfcn = -1
      helicity = CMPLX(0._dp, 0._dp)
      sym_type = 'NONE'

!---------------------------------------------------------------------------------
!   CODE added by R.SANCHEZ (01/19/99): ballooning related variables and sigmas.

      target_balloon = 0
      sigma_balloon = bigno
      sigma_bal = bigno                                                  !Old style - LPK
      sigma_pgrad = bigno
      sigma_pedge = bigno
      lballoon_opt = .false.
      lballoon_mask = .false.
      lpres_prof_opt = .false.
      nballoon_mask = 0                                                  !VMECCOBRA
      l_legendre = .false.                                               !LEGENDRE

!------------------------------------------------------------------------------
!  initialize 3D boundary limiter parameters,  RH & MZ  June 2000

      target_vv = 0
      target_vv_rms = 0
      sigma_vv = bigno
      sigma_vv_rms = bigno
      rbc_vv = 0
      zbs_vv = 0
      mpol_vv = 0
      ntor_vv = 0
      nu_vv = 5
      nv_vv = 2
      shapeweight = .FALSE.; amplw_bw = 0
      theta0_bw(1) = 0 ; theta0_bw(2) = 0 ; theta0_bw(3) = 0 
      phi0_bw = 0 ; wtheta_bw = .5 ; wphi_bw = .5
      planes_bw(1) = 0 ; planes_bw(2) = 1.5707 ; planes_bw(3) = 3.14159

!------------------------------------------------------------------------------

      sigma_kappa = bigno
      target_iota_p = 0
      sigma_iota_pmax = bigno
      sigma_iota_pmin = bigno
      sigma_extcur = bigno
      sigma_oh = bigno
      oh_coefs = 0
      bal_theta0 = 0; bal_zeta0 = 0

      READ (iunit, nml=optimum, iostat=istat)

      END SUBROUTINE read_optimum_namelist

      END MODULE optim_params

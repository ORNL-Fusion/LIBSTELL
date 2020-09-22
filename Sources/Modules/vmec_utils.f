      MODULE vmec_utils
      USE stel_kinds
      USE stel_constants, ONLY: twopi, one, zero
      USE cyl_flux, ONLY: flx2cyl, cyl2flx
      IMPLICIT NONE

!-------------------------------------------------------------------------------!
!     THIS MODULE CONTAINS USEFUL UTILITIES FOR PROCESSING VMEC 
!     DATA. MOST FUNCTIONS ARE OVERLOADED TO BE ABLE TO USE EITHER 
!     INTERNALLY DATA (LOCAL FROM WITHIN VMEC) OR DATA FROM WOUT FILE
!
!  CONTAINS:
!     GetBcyl_WOUT - call as GetBcyl
!     GetBcyl_VMEC - call as GetBcyl
!     GetJcyl_WOUT - call as GetJcyl
!     MSE_pitch_WOUT - call as MSE_pitch
!     MSE_pitch_VMEC - call as MSE_pitch
!
!  2011-09-08 JDH. rzl* declaration 2*ntmax) -> 3*ntmax), consistency
!    with rzl* usage in VMEC, even though lambdas aren't used here.
!
!  2011-09-06 JDH
!    Split subroutines flx2cyl and cyl2flx into the module cyl_flux
!    - to clarify interfaces and dependencies
!-------------------------------------------------------------------------------
!

!
!     OVERLOADED FUNCTIONS
!
      INTERFACE GetBcyl
          MODULE PROCEDURE GetBcyl_WOUT, GetBcyl_VMEC
      END INTERFACE

      INTERFACE GetJcyl
          MODULE PROCEDURE GetJcyl_WOUT
      END INTERFACE

      INTERFACE MSE_pitch
          MODULE PROCEDURE MSE_pitch_WOUT
          MODULE PROCEDURE MSE_pitch_VMEC
      END INTERFACE

      CONTAINS

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      SUBROUTINE GetBcyl_WOUT(R1, Phi, Z1, Br, Bphi, Bz, 
     1                        sflx, uflx, info)
      USE read_wout_mod, phi_wout=>phi, ns_w=>ns, ntor_w=>ntor,
     1     mpol_w=>mpol, ntmax_w=>ntmax, lthreed_w=>lthreed,
     2     lasym_w=>lasym
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, OPTIONAL, INTENT(out) :: info
      REAL(rprec), INTENT(in)  :: R1, Z1, Phi
      REAL(rprec), INTENT(out) :: Br, Bphi, Bz
      REAL(rprec), INTENT(out), OPTIONAL :: sflx, uflx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: fmin_acceptable = 1.E-6_dp
      INTEGER     :: nfe, info_loc
      REAL(rprec) :: r_cyl(3), c_flx(3), fmin
      REAL(rprec) :: Ru1, Zu1, Rv1, Zv1
      REAL(rprec) :: bsupu1, bsupv1
C-----------------------------------------------
      IF (.not.lwout_opened) THEN
         WRITE(6, '(2a,/,a)')
     1   ' This form of GetBcyl can only be called if WOUT has been',
     2   ' previously opened!',' Try GetBcyl_VMEC form instead.'
         RETURN
      END IF

      CALL LoadRZL

!     Computes cylindrical components of the magnetic field, Br, Bphi, Bz,
!     at the specified cylindrical coordinate point (R1, Phi, Z1), where
!     Phi is the true geometric toroidal angle (NOT N*Phi)
!
!     INPUT
!     R1, Phi, Z1  : cylindrical coordinates at which evaluation is to take place
!     
!     OUTPUT
!     Br, Bphi, Bz : computed cylindrical components of B at input point
!     sflx, uflx   : computed flux and theta angle at the cylindrical point
!
!     1. Convert to point in flux-coordinates: cflux = (s, u, v=N*phi)
!        and evaluate Ru, Zu, Rv, Zv at that point
!
      r_cyl(1) = R1;  r_cyl(2) = nfp*Phi;  r_cyl(3) = Z1
      c_flx(1) = 0;   c_flx(2) = 0;        c_flx(3) = r_cyl(2)
      IF (PRESENT(sflx)) c_flx(1) = sflx
      IF (PRESENT(uflx)) c_flx(2) = uflx
      CALL cyl2flx(rzl_local, r_cyl, c_flx, ns_w, ntor_w, mpol_w, 
     1     ntmax_w, lthreed_w, lasym_w, info_loc, nfe, fmin, 
     2     RU=Ru1, ZU=Zu1, RV=Rv1, ZV=Zv1)
      Rv1 = nfp*Rv1;  Zv1 = nfp*Zv1

      IF (info_loc.eq.-1 .and. (fmin .le. fmin_acceptable)) info_loc = 0

      IF (PRESENT(info)) info = info_loc
      IF (info_loc .ne. 0) RETURN

      IF (PRESENT(sflx)) sflx = c_flx(1)  
      IF (PRESENT(uflx)) uflx = c_flx(2)

      IF (c_flx(1) .gt. one) THEN
         Br = 0;  Bphi = 0;  Bz = 0
         RETURN
      END IF
!
!     2. Evaluate Bsupu, Bsupv at this point
!
      CALL tosuvspace (c_flx(1), c_flx(2), c_flx(3), 
     1                 BSUPU=bsupu1, BSUPV=bsupv1)
!
!     3. Form Br, Bphi, Bz
!
      Br   = Ru1*bsupu1 + Rv1*bsupv1
      Bphi = R1 *bsupv1
      Bz   = Zu1*bsupu1 + Zv1*bsupv1
      
      END SUBROUTINE GetBcyl_WOUT

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      SUBROUTINE GetBcyl_VMEC(R1, Phi, Z1, Br, Bphi, Bz, sflx, uflx, 
     1     bsupu, bsupv, rzl_array, ns_in, ntor_in, mpol_in, ntmax_in, 
     2     nzeta, ntheta3, nper, mscale, nscale, lthreed_in, lasym_in,  
     3     info)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ns_in, ntor_in, mpol_in, ntmax_in, 
     1                       nzeta, ntheta3, nper
      INTEGER, OPTIONAL, INTENT(out) :: info
      LOGICAL, INTENT(in) :: lthreed_in, lasym_in
      REAL(rprec), INTENT(in)  :: R1, Z1, Phi
      REAL(rprec), INTENT(in)  :: 
     1             rzl_array(ns_in,0:ntor_in,0:mpol_in-1,3*ntmax_in),
     2             mscale(0:mpol_in-1), nscale(0:ntor_in)
      REAL(rprec), DIMENSION(ns_in,nzeta,ntheta3), INTENT(in) 
     1                         :: bsupu, bsupv
      REAL(rprec), INTENT(out) :: Br, Bphi, Bz, sflx, uflx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: c1p5 = 1.5_dp
      REAL(rprec), PARAMETER :: fmin_acceptable = 1.E-6_dp
      INTEGER     :: nfe, info_loc, jslo, jshi, julo, juhi, 
     1               kvlo, kvhi, ntheta1
      REAL(rprec) :: r_cyl(3), c_flx(3), vflx, vflx_norm, 
     1               uflx_norm, fmin
      REAL(rprec) :: wgt_s, wgt_u, wgt_v, hs1, hu1, hv1
      REAL(rprec) :: Ru1, Zu1, Rv1, Zv1
      REAL(rprec) :: bsupu1, bsupv1, bsupu2, bsupv2
C-----------------------------------------------

!     Computes cylindrical components of the magnetic field, Br, Bphi, Bz,
!     at the specified cylindrical coordinate point (R1, Phi, Z1), where
!     Phi is the true geometric toroidal angle (NOT NPER*Phi)
!     Also, sflx, uflx are the computed flux and theta angle at the point
!
!     This routine is callable from within the VMEC code (in contrast to
!     the GetBcyl routine, which requires WOUT output file).
!
!     1. Convert to point in flux-coordinates: cflux = (s, u, v=N*phi)
!        and evaluate Ru, Zu, Rv, Zv at that point
!
      r_cyl(1) = R1;  r_cyl(2) = nper*Phi;  r_cyl(3) = Z1
      c_flx(1) = 0;   c_flx(2) = 0;         c_flx(3) = r_cyl(2)
      CALL cyl2flx(rzl_array, r_cyl, c_flx, ns_in, ntor_in, mpol_in, 
     1     ntmax_in, lthreed_in, lasym_in, info_loc, nfe, fmin, 
     2     mscale, nscale, RU=Ru1, ZU=Zu1, RV=Rv1, ZV=Zv1)
      Rv1 = nper*Rv1;  Zv1 = nper*Zv1

      IF (info_loc.eq.-1 .and. (fmin .le. fmin_acceptable)) info_loc = 0

      IF (PRESENT(info)) info = info_loc
      IF (info_loc .ne. 0) RETURN

      sflx = c_flx(1);  uflx = c_flx(2);  vflx = c_flx(3)
      IF (c_flx(1) .gt. one) THEN
         Br = 0;  Bphi = 0;  Bz = 0
         RETURN
      END IF

!
!     2. Evaluate Bsupu, Bsupv at this flux coordinate point by 2D interpolation in s, u space
!        This is not quite as accurate as the 1D (s) interpolation based on the Fourier coefficients
!        of bsupu, bsupv...done in GetBcyl...
!        Formula 25.2.66 (Bivariate, 4pt Formula) in Abramowitz and Stegun
!
      hs1 = one/(ns_in - 1)
      jslo = INT(c1p5 + sflx/hs1)
      jshi = jslo+1
      wgt_s = (sflx - hs1*(jslo-c1p5))/hs1
      IF (jslo .eq. ns_in) THEN
!        USE Xhalf(ns+1) = 2*Xhalf(ns) - Xhalf(ns-1) FOR "GHOST" POINT VALUE hs/2 OUTSIDE EDGE
!        THEN, X = wlo*Xhalf(ns) + whi*Xhalf(ns+1) == Xhalf(ns) + whi*(Xhalf(ns) - Xhalf(ns-1)) 
!        WHERE wlo = 1 - wgt_s, whi = wgt_s
         jshi = jslo-1
         wgt_s = 1+wgt_s
      ELSE IF (jslo .eq. 1) THEN
         jslo = 2
      END IF

      IF (lasym_in) THEN
         ntheta1 = ntheta3
      ELSE
         ntheta1 = 2*(ntheta3 - 1)
      END IF
      
      uflx = MOD(uflx, twopi)
      DO WHILE (uflx .lt. zero) 
         uflx = uflx+twopi
      END DO

      hu1 = one/ntheta1
      uflx_norm = uflx/twopi
      julo = INT(1 + uflx_norm/hu1)
      IF (julo .gt. ntheta3) THEN
         IF (ABS(uflx_norm - 1) .lt. 1.E-2*hu1) THEN
            julo = 1
            uflx_norm = 0
         ELSE IF (ABS(uflx_norm - .5_dp) .lt. 1.E-2*hu1) THEN
            julo = ntheta3
            uflx_norm = .5_dp
         ELSE
            PRINT *, 'julo=', julo,' > ntheta3=', ntheta3,
     1      ' uflx_norm=', uflx_norm, ' in GetBcyl!'
            IF (PRESENT(info)) info = -10
            RETURN
         END IF
      END IF
      juhi = julo + 1
      IF (julo .eq. ntheta3) juhi = 1         !Periodic point at u = 0
      wgt_u = (uflx_norm - hu1*(julo-1))/hu1

      
      DO WHILE (vflx .lt. zero) 
         vflx = vflx+twopi
      END DO
      vflx = MOD(vflx, twopi)
      hv1 = one/nzeta
      vflx_norm = vflx/twopi
      kvlo = INT(1 + vflx_norm/hv1)
      kvhi = kvlo+1
      IF (kvlo .eq. nzeta) kvhi = 1
      wgt_v = (vflx_norm - hv1*(kvlo-1))/hv1

!
!     BIVARIATE INTERPOLATION IN S, U AT 2 kv PLANES
!
      bsupu1 = (1-wgt_s)*((1-wgt_u)*bsupu(jslo,kvlo,julo)
     2       +               wgt_u *bsupu(jslo,kvlo,juhi))
     1       + wgt_s*    ((1-wgt_u)*bsupu(jshi,kvlo,julo)
     3       +               wgt_u *bsupu(jshi,kvlo,juhi))

      bsupv1 = (1-wgt_s)*((1-wgt_u)*bsupv(jslo,kvlo,julo)
     2       +               wgt_u *bsupv(jslo,kvlo,juhi))
     1       + wgt_s*    ((1-wgt_u)*bsupv(jshi,kvlo,julo)
     3       +               wgt_u *bsupv(jshi,kvlo,juhi))

      bsupu2 = (1-wgt_s)*((1-wgt_u)*bsupu(jslo,kvhi,julo)
     2       +               wgt_u *bsupu(jslo,kvhi,juhi))
     1       + wgt_s*    ((1-wgt_u)*bsupu(jshi,kvhi,julo)
     3       +               wgt_u *bsupu(jshi,kvhi,juhi))

      bsupv2 = (1-wgt_s)*((1-wgt_u)*bsupv(jslo,kvhi,julo)
     2       +               wgt_u *bsupv(jslo,kvhi,juhi))
     1       + wgt_s*    ((1-wgt_u)*bsupv(jshi,kvhi,julo)
     3       +               wgt_u *bsupv(jshi,kvhi,juhi))

!
!     LINEAR INTERPOLATION IN V
!
      bsupu1 = (1-wgt_v)*bsupu1 + wgt_v*bsupu2
      bsupv1 = (1-wgt_v)*bsupv1 + wgt_v*bsupv2

!
!     3. Form Br, Bphi, Bz
!
      Br   = Ru1*bsupu1 + Rv1*bsupv1
      Bphi = R1 *bsupv1
      Bz   = Zu1*bsupu1 + Zv1*bsupv1
      
      END SUBROUTINE GetBcyl_VMEC

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      SUBROUTINE GetJcyl_WOUT(R1, Phi, Z1, JR, JPHI, JZ, 
     1                        sflx, uflx, info)
      USE read_wout_mod, phi_wout1=>phi, ns_w1=>ns, ntor_w1=>ntor,
     1     mpol_w1=>mpol, ntmax_w1=>ntmax, lthreed_w1=>lthreed,
     2     lasym_w1=>lasym
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, OPTIONAL, INTENT(out) :: info
      REAL(rprec), INTENT(in)  :: R1, Z1, Phi
      REAL(rprec), INTENT(out) :: JR, JPHI, JZ
      REAL(rprec), INTENT(out), OPTIONAL :: sflx, uflx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: fmin_acceptable = 1.E-6_dp
      INTEGER     :: nfe, info_loc
      REAL(rprec) :: r_cyl(3), c_flx(3), fmin
      REAL(rprec) :: Ru1, Zu1, Rv1, Zv1
      REAL(rprec) :: jsupu1, jsupv1, gsqrt1
C-----------------------------------------------
      IF (.not.lwout_opened) THEN
         WRITE(6, '(2a,/,a)')
     1   ' This form of GetBcyl can only be called if WOUT has been',
     2   ' previously opened!'
         RETURN
      END IF

      CALL LoadRZL

!     Computes cylindrical components of the current, Jr, Jphi, Jz,
!     at the specified cylindrical coordinate point (R1, Phi, Z1), where
!     Phi is the true geometric toroidal angle (NOT N*Phi)
!
!     INPUT
!     R1, Phi, Z1  : cylindrical coordinates at which evaluation is to take place
!     
!     OUTPUT
!     Br, Bphi, Bz : computed cylindrical components of B at input point
!     sflx, uflx   : computed flux and theta angle at the cylindrical point
!
!     1. Convert to point in flux-coordinates: cflux = (s, u, v=N*phi)
!        and evaluate Ru, Zu, Rv, Zv at that point
!
      r_cyl(1) = R1;  r_cyl(2) = nfp*Phi;  r_cyl(3) = Z1
      c_flx(1) = 0;   c_flx(2) = 0;        c_flx(3) = r_cyl(2)
      CALL cyl2flx(rzl_local, r_cyl, c_flx, ns_w1, ntor_w1, mpol_w1, 
     1     ntmax_w1, lthreed_w1, lasym_w1, info_loc, nfe, fmin, 
     2     RU=Ru1, ZU=Zu1, RV=Rv1, ZV=Zv1)
      Rv1 = nfp*Rv1;  Zv1 = nfp*Zv1

      IF (info_loc.eq.-1 .and. (fmin .le. fmin_acceptable)) info_loc = 0

      IF (PRESENT(info)) info = info_loc
      IF (info_loc .ne. 0) RETURN

      IF (PRESENT(sflx)) sflx = c_flx(1)  
      IF (PRESENT(uflx)) uflx = c_flx(2)

      IF (c_flx(1) .gt. one) THEN
         Jr = 0;  Jphi = 0;  Jz = 0
         RETURN
      END IF

!     3. Evaluate d(Bsubs)/du and d(Bsubs)/dv, d(Bsubu)/ds, d(Bsubv)/ds at this point
      CALL tosuvspace (c_flx(1), c_flx(2), c_flx(3), 
     1                 GSQRT=gsqrt1, JSUPU=jsupu1, JSUPV=jsupv1)

!      WRITE (36, '(1p4e12.4)') R1*jsupv1, dbsubuds1, dbsubsdu1, gsqrt1
!
!     4. Return Jr, Jphi, Jz
!
      Jr   = Ru1*jsupu1 + Rv1*jsupv1
      Jphi =              R1 *jsupv1
      Jz   = Zu1*jsupu1 + Zv1*jsupv1
      
      END SUBROUTINE GetJcyl_WOUT

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      FUNCTION MSE_pitch_WOUT(r1, phi1, z1, acoef, efield, info)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(out) :: info
      REAL(rprec), INTENT(in) :: r1, phi1, z1, acoef(6)
      REAL(rprec), INTENT(in), OPTIONAL :: efield(2)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: MSE_pitch_WOUT, Er, Ez, Br, Bphi, Bz
C-----------------------------------------------
!
!     Computes the Motional Stark Effect tan(pitch angle) == mse_pitch
!     at a given cylindrical point (R, phi, Z) inside the plasma
!
!     INPUT
!     Acoef : array of constants defined by the viewing geometry and beam velocity
!     r1, f1, z1: cylindrical coordinate of measurement point
!     Efield: (optional) electric field components (Er, Ez) in rest frame
!
!     OUTPUT
!     MSE_pitch  pitch angle at the input point
!     info       info = 0, calculation is valid
!
      IF (PRESENT(efield)) THEN
         Er = efield(1);  Ez = efield(2)
      ELSE
         Er = 0; Ez = 0
      END IF
      info = -1
!
!     Compute cylindrical components of B-field at given point R1, phi=f1, Z1
!
      CALL GetBcyl_WOUT(r1, phi1, z1, br, bphi, bz, INFO=info)

      MSE_pitch_WOUT = (acoef(1)*Bz   + acoef(5)*Er)/
     1                 (acoef(2)*Bphi + acoef(3)*Br 
     2               + acoef(4)*Bz   + acoef(6)*Ez)

      END FUNCTION MSE_pitch_WOUT

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      FUNCTION MSE_pitch_VMEC(r1, phi1, z1, acoef, efield, sflx, uflx, 
     1     bsupu, bsupv, rzl_array, ns_in, ntor_in, mpol_in, ntmax_in, 
     2     nzeta, ntheta3, nper, mscale, nscale, lthreed_in, lasym_in,  
     3     info)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), INTENT(in) :: r1, phi1, z1, acoef(6), efield(2)
      INTEGER, INTENT(in) :: ns_in, ntor_in, mpol_in, ntmax_in, 
     1                       nzeta, ntheta3, nper
      LOGICAL, INTENT(in) :: lthreed_in, lasym_in
      REAL(rprec), INTENT(in)  :: 
     1             rzl_array(ns_in,0:ntor_in,0:mpol_in-1,3*ntmax_in),
     2             mscale(0:mpol_in-1), nscale(0:ntor_in)
      REAL(rprec), DIMENSION(ns_in,nzeta,ntheta3), INTENT(in) 
     1                         :: bsupu, bsupv
      REAL(rprec), INTENT(out) :: sflx, uflx
      INTEGER, INTENT(out) :: info
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: MSE_pitch_VMEC, Er, Ez, Br, Bphi, Bz
C-----------------------------------------------
!
!     Computes the Motional Stark Effect tan(pitch angle) == mse_pitch
!     at a given cylindrical point (R, phi, Z) inside the plasma
!
!     INPUT
!     Acoef : array of constants defined by the viewing geometry and beam velocity
!     r1, f1, z1: cylindrical coordinate of measurement point
!     Efield: (optional) electric field components (Er, Ez) in rest frame
!
!     OUTPUT
!     MSE_pitch  pitch angle at the input point
!     info       info = 0, calculation is valid
!
      Er = efield(1);  Ez = efield(2)
      info = -1
!
!     Compute cylindrical components of B-field at given point R1, phi=f1, Z1
!
      CALL GetBcyl_VMEC(r1, phi1, z1, br, bphi, bz, sflx, uflx, 
     1     bsupu, bsupv, rzl_array, ns_in, ntor_in, mpol_in, ntmax_in, 
     2     nzeta, ntheta3, nper, mscale, nscale, lthreed_in, lasym_in,  
     3     info)

      MSE_pitch_VMEC = (acoef(1)*Bz   + acoef(5)*Er)/
     1                 (acoef(2)*Bphi + acoef(3)*Br 
     2               +  acoef(4)*Bz   + acoef(6)*Ez)

      END FUNCTION MSE_pitch_VMEC

      END MODULE vmec_utils

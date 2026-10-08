      SUBROUTINE DFLUX(FLUX,SOL,KSTEP,KINC,TIME,NOEL,NPT,COORDS,
     1                 JLTYP,TEMP,PRESS,SNAME)
C======================================================================
C DFLUX - Conical Gaussian volumetric heat source travelling along a
C         circular path around the global Z axis (circumferential
C         laser / electron-beam weld of two concentric parts).
C
C Units   : length = mm, time = s, energy = mJ, temperature = C
C           -> volumetric flux returned in mJ/(s mm3) = mW/mm3
C Path    : circle of radius RADIUS centred on the Z axis, starting on
C           the +X axis (THETA = 0) and rotating counter-clockwise.
C Power   : optional linear ramp-up, constant segment, linear
C           ramp-down, all defined as angles travelled by the beam.
C Language: fixed-form Fortran 77 (72 columns, no tabs, ASCII only).
C Author  : Victor Maia (vhfm08@gmail.com)
C
C Heat source (single cone, Gaussian in r, linear radius in depth):
C
C   q(r,z) = Q0 * exp( -3 r^2 / R0(z)^2 )
C   R0(z)  = r_e + (r_i - r_e) * z / zi          0 <= z <= zi
C   Q0     = 3 ETA QTOT / ( pi * (zi/3)(r_e^2 + r_e r_i + r_i^2) )
C
C   so that the volume integral of q equals ETA*QTOT exactly.
C======================================================================
      INCLUDE 'ABA_PARAM.INC'

C======================================================================
C ---------------------- Abaqus-mandated arguments --------------------
C======================================================================
      DIMENSION FLUX(2), TIME(2), COORDS(3)
      CHARACTER*80 SNAME

C======================================================================
C -------------------------- Local variables --------------------------
C======================================================================
      REAL*8 X, Y, ZINT, XC, YC, DX, DY
      REAL*8 THETA, THETA_DEG, AMP
      REAL*8 R_LOC, Z_REL, R0Z, R_CUT, EXP_R, QVOL

C======================================================================
C ------------------------- User parameters ---------------------------
C QTOT       = nominal beam power (mJ/s)   [1 W = 1.0D3 mJ/s]
C ETA        = absorption efficiency (-)
C R_E        = beam radius at the top surface (mm)
C R_I        = beam radius at depth ZI (mm)
C ZI         = depth of the heat source (mm)
C RADIUS     = radius of the circular weld path (mm)
C VBEAM      = beam travel speed (mm/min)
C RAMP_U_DEG = ramp-up segment (deg)
C WELD_DEG   = full-power segment (deg)
C RAMP_D_DEG = ramp-down segment (deg)
C Z0         = Z coordinate of the top (irradiated) surface (mm)
C TOL        = relative Gaussian cut-off, exp(-3 r^2/R0^2) < TOL (-)
C
C NOTE: values below are illustrative and not tied to any real process.
C======================================================================
      REAL*8 QTOT, ETA, R_E, R_I, ZI, RADIUS, VBEAM
      REAL*8 RAMP_U_DEG, WELD_DEG, RAMP_D_DEG, Z0, TOL
      PARAMETER (QTOT       = 1.8D6 )
      PARAMETER (ETA        = 0.60D0)
      PARAMETER (R_E        = 1.00D0)
      PARAMETER (R_I        = 0.75D0)
      PARAMETER (ZI         = 3.00D0)
      PARAMETER (RADIUS     = 15.0D0)
      PARAMETER (VBEAM      = 1200.0D0)
      PARAMETER (RAMP_U_DEG = 0.0D0 )
      PARAMETER (WELD_DEG   = 360.0D0)
      PARAMETER (RAMP_D_DEG = 20.0D0)
      PARAMETER (Z0         = 20.0D0)
      PARAMETER (TOL        = 1.0D-8)

C======================================================================
C ------------- Derived constants (DO NOT CHANGE) ---------------------
C OMEGA = angular speed (rad/s)
C DENOM = (zi/3)(r_e^2 + r_e r_i + r_i^2)  -> cone volume / pi
C PREFAC= peak volumetric flux Q0 (mJ/(s mm3))
C TOTDEG= total angle travelled with the beam switched on (deg)
C Computed as compile-time PARAMETERs: no SAVE, thread-safe.
C======================================================================
      REAL*8 PI, OMEGA, DENOM, PREFAC, TOTDEG
      PARAMETER (PI     = 3.141592653589793D0)
      PARAMETER (OMEGA  = (VBEAM/60.0D0)/RADIUS)
      PARAMETER (DENOM  = (ZI/3.0D0)*(R_I*R_I + R_E*R_I + R_E*R_E))
      PARAMETER (PREFAC = (ETA*QTOT)*3.0D0/(PI*DENOM))
      PARAMETER (TOTDEG = RAMP_U_DEG + WELD_DEG + RAMP_D_DEG)

C======================================================================
C --- Input checks -----------------------------------------------------
C======================================================================
      IF (JLTYP .NE. 1) THEN
        WRITE(6,*) '**ERROR: DFLUX must be applied as body flux (BFNU)'
        CALL XIT
      END IF
      IF (RADIUS .LE. 0.0D0) THEN
        WRITE(6,*) '**ERROR: RADIUS must be > 0'
        CALL XIT
      END IF

      FLUX(1) = 0.0D0
      FLUX(2) = 0.0D0

C======================================================================
C --- Integration point coordinates ------------------------------------
C======================================================================
      X    = COORDS(1)
      Y    = COORDS(2)
      ZINT = COORDS(3)

C======================================================================
C --- Current angular position of the source (TIME(1) = step time) ----
C======================================================================
      THETA     = OMEGA*TIME(1)
      THETA_DEG = THETA*180.0D0/PI

C --- Beam switched off after ramp-up + weld + ramp-down --------------
      IF (THETA_DEG .GT. TOTDEG) RETURN

C======================================================================
C --- Depth below top surface (points above Z0 clamped to surface) ----
C======================================================================
      Z_REL = Z0 - ZINT
      IF (Z_REL .LT. 0.0D0) Z_REL = 0.0D0

C --- Deeper than the heat-source depth: no flux -----------------------
      IF (Z_REL .GT. ZI) RETURN

C======================================================================
C --- Horizontal distance to the instantaneous beam centre ------------
C======================================================================
      XC    = RADIUS*COS(THETA)
      YC    = RADIUS*SIN(THETA)
      R0Z   = R_E + (R_I - R_E)*(Z_REL/ZI)
      DX    = X - XC
      DY    = Y - YC
      R_LOC = SQRT(DX*DX + DY*DY)
      R_CUT = SQRT(-LOG(TOL)/3.0D0)*R0Z

C --- Outside the Gaussian cut-off radius: no flux ---------------------
      IF (R_LOC .GT. R_CUT) RETURN

C======================================================================
C --- Power amplitude: ramp-up / constant / ramp-down -----------------
C======================================================================
      IF (THETA_DEG .LT. RAMP_U_DEG) THEN
        AMP = THETA_DEG/RAMP_U_DEG
      ELSE IF (THETA_DEG .LT. RAMP_U_DEG + WELD_DEG) THEN
        AMP = 1.0D0
      ELSE IF (RAMP_D_DEG .GT. 0.0D0) THEN
        AMP = (TOTDEG - THETA_DEG)/RAMP_D_DEG
      ELSE
        AMP = 0.0D0
      END IF

C======================================================================
C --- Gaussian volumetric flux; no temperature dependence -> dq/dT = 0
C======================================================================
      EXP_R   = EXP(-3.0D0*R_LOC*R_LOC/(R0Z*R0Z))
      QVOL    = PREFAC*EXP_R*AMP
      FLUX(1) = QVOL
      FLUX(2) = 0.0D0

      RETURN
      END

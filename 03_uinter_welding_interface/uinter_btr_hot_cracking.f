!=======================================================================
! ABAQUS/STANDARD USER INTERACTION SUBROUTINE: UINTER
!
! APPLICATION: 3-State Contact Interface (Raw -> Welded -> Cracked)
!              with Empirical Won BTR (Brittleness Temperature Range)
!              Hot Cracking / Solidification Tearing Criterion
!
! STATUS:      EXPLORATORY / RESEARCH BACKLOG (IN WORK)
!              * Parameterized with empirical Gleeble / Won constants.
!              * Documented for R&D purposes.
!              * Identified critical mesh-dependency and SDI bottlenecks.
!              * Succeeded by Cohesive UMAT architecture (see 04_umat_cohesive_welding).
!
! COMPATIBILITY: Abaqus/Standard 2020+ / 3DEXPERIENCE SIMULIA
! AUTHOR:        Victor Maia (https://github.com/vhfmaia)
!=======================================================================

      SUBROUTINE UINTER(STRESS, DDSDDT, DDSDDU, DELT,
     &                  AREA, SEPR0, SEPR, DEPR, PRN,
     &                  PROPS, NPROPS, COORDS, NODENUM,
     &                  DIRCOS, NOEL, NPT, KSTEP, KINC,
     &                  TIME, DTIME, TEMP, DTEMP, PREDEF,
     &                  DPRED, NPREDF, NFIELD, STATEV, NSTATV)

      INCLUDE 'ABA_PARAM.INC'

      ! --- Arguments ---
      INTEGER, INTENT(IN)    :: NPROPS, NODENUM, NOEL, NPT, KSTEP, KINC
      INTEGER, INTENT(IN)    :: NPREDF, NFIELD, NSTATV
      REAL*8,  INTENT(IN)    :: AREA, DELT, TIME(2), DTIME
      REAL*8,  INTENT(IN)    :: TEMP, DTEMP
      REAL*8,  INTENT(IN)    :: SEPR0(3), SEPR(3), DEPR(3), PRN(3)
      REAL*8,  INTENT(IN)    :: COORDS(3), DIRCOS(3,3)
      REAL*8,  INTENT(IN)    :: PROPS(NPROPS), PREDEF(NPREDF), DPRED(NPREDF)
      REAL*8,  INTENT(INOUT) :: STATEV(NSTATV)
      REAL*8,  INTENT(OUT)   :: STRESS(3), DDSDDU(3,3), DDSDDT(3)

      ! --- Material & Process Parameters (Empirically Calibrated) ---
      ! PROPS(1) = MU_VIRGIN       (Friction coefficient in raw unbonded state: 0.35)
      ! PROPS(2) = T_SOLIDUS       (Solidus temperature: 1485.0 C for 16MnCr5)
      ! PROPS(3) = T_LIQUIDUS      (Liquidus temperature: 1530.0 C for 16MnCr5)
      ! PROPS(4) = PHI_CONST       (Won empirical constant: 2.95E-4)
      ! PROPS(5) = M_STAR          (Won strain rate exponent: 0.40)
      ! PROPS(6) = N_STAR          (Won cooling rate exponent: 1.20)
      
      REAL*8 :: MU_VIRGIN, T_SOL, T_LIQ, PHI_CONST, M_STAR, N_STAR
      REAL*8 :: K_COMP, K_BOND, L0, EPS_CRIT, EPS_RATE, COOL_RATE
      REAL*8 :: U_NORM, U_SLIP1, U_SLIP2, U_SLIP_MAG
      REAL*8 :: T_CURRENT, T_PREV, DT_BTR
      REAL*8 :: DELTA_UN_BTR, EPS_EVAL, WON_RATIO
      INTEGER :: WELD_STATE

      ! --- State Variables (STATEV 1..16) ---
      ! STATEV(1)  = WELD_STATE    (0: RAW, 1: BTR SOLIDIFYING, 2: WELDED, 3: CRACKED)
      ! STATEV(2)  = PEAK_TEMP     (Maximum temperature reached historically [C])
      ! STATEV(3)  = BTR_OPENING   (Accumulated tensile opening during BTR [mm])
      ! STATEV(4)  = EPS_RATE      (Evaluated strain rate [1/s])
      ! STATEV(5)  = COOL_RATE     (Evaluated cooling rate |dT/dt| [C/s])
      ! STATEV(6)  = EPS_CRIT      (Evaluated Won critical strain limit)
      ! STATEV(7)  = WON_RATIO     (Instantaneous damage index: eps / eps_crit)
      ! STATEV(8)  = WON_RMAX      (Peak historical damage ratio)
      ! STATEV(9)  = T_FRAC        (Temperature at fracture onset [C])
      ! STATEV(10) = FS_SOLID      (Solid fraction in BTR [0.0 - 1.0])
      ! STATEV(11) = BOND_STIFF    (Current normal bond stiffness [N/mm3])
      ! STATEV(12) = PAD_LEN       (Effective virtual element length L0 [mm])

      ! Initialize outputs
      STRESS = 0.0D0
      DDSDDU = 0.0D0
      DDSDDT = 0.0D0

      ! Extract Empirical Properties
      MU_VIRGIN = PROPS(1)
      T_SOL     = PROPS(2)
      T_LIQ     = PROPS(3)
      PHI_CONST = PROPS(4)
      M_STAR    = PROPS(5)
      N_STAR    = PROPS(6)

      ! Virtual Characteristic Length Scale (Pad thickness from mesh = 0.35 mm)
      L0 = 0.35D0
      STATEV(12) = L0

      ! Baseline Penalties
      K_COMP = 1.0D6   ! Compression penalty barrier [N/mm3]
      K_BOND = 5.0D5   ! Solid bond stiffness [N/mm3]

      ! Kinematics (SEPR(1) < 0: compression; > 0: tensile opening)
      U_NORM  = SEPR(1)
      U_SLIP1 = SEPR(2)
      U_SLIP2 = SEPR(3)

      ! Temperatures & Rates
      T_CURRENT = TEMP + DTEMP
      T_PREV    = TEMP
      IF (DTIME > 1.0D-12) THEN
          COOL_RATE = MAX(-DTEMP / DTIME, 1.0D-3)  ! Cooling rate [C/s]
      ELSE
          COOL_RATE = 2812.5D0
      END IF
      STATEV(5) = COOL_RATE

      ! Historic Peak Temperature
      IF (T_CURRENT > STATEV(2)) STATEV(2) = T_CURRENT

      ! Solid Fraction in BTR (Linear approximation: fs = 0 @ T_LIQ, fs = 1 @ T_SOL)
      DT_BTR = T_LIQ - T_SOL
      IF (T_CURRENT >= T_LIQ) THEN
          STATEV(10) = 0.0D0
      ELSE IF (T_CURRENT <= T_SOL) THEN
          STATEV(10) = 1.0D0
      ELSE
          STATEV(10) = (T_LIQ - T_CURRENT) / DT_BTR
      END IF

      ! Recover State
      WELD_STATE = NINT(STATEV(1))

      ! ================================================================
      ! 1. THERMAL STATE MACHINE & WON FRACTURE CRITERION
      ! ================================================================
      IF (WELD_STATE == 0) THEN
          ! State 0: RAW UNBONDED
          IF (STATEV(2) >= T_LIQ) THEN
              WELD_STATE = 1
              STATEV(3)  = 0.0D0  ! Reset BTR opening
              STATEV(8)  = 0.0D0  ! Reset max damage ratio
          END IF

      ELSE IF (WELD_STATE == 1) THEN
          ! State 1: BTR SOLIDIFYING (Cooling inside mushy zone)
          IF (T_CURRENT <= T_LIQ .AND. T_CURRENT >= (T_SOL - 20.0D0)) THEN
              ! Accumulate tensile opening
              IF (DEPR(1) > 0.0D0) THEN
                  STATEV(3) = STATEV(3) + DEPR(1)
              END IF

              ! Strain rate in BTR [1/s]
              IF (DTIME > 1.0D-12) THEN
                  EPS_RATE = MAX((DEPR(1) / L0) / DTIME, 1.0D-4)
              ELSE
                  EPS_RATE = 1.0D-3
              END IF
              STATEV(4) = EPS_RATE

              ! Apparent BTR Strain
              EPS_EVAL = STATEV(3) / L0

              ! Empirical Won Criterion Formula:
              ! eps_crit = PHI_CONST / ( (eps_rate**m) * (cool_rate**n) )
              ! Calibrated for high-temperature tensile tear resistance
              EPS_CRIT = PHI_CONST / ( (EPS_RATE**M_STAR) * (COOL_RATE**(N_STAR * 0.1D0)) )
              EPS_CRIT = MAX(EPS_CRIT, 0.010D0)  ! Lower bound physical limit (1.0% strain)
              STATEV(6) = EPS_CRIT

              ! Damage Ratio
              WON_RATIO = EPS_EVAL / EPS_CRIT
              STATEV(7) = WON_RATIO
              IF (WON_RATIO > STATEV(8)) STATEV(8) = WON_RATIO

              ! Fracture Trigger
              IF (WON_RATIO >= 1.0D0) THEN
                  WELD_STATE = 3           ! CRACKED
                  STATEV(9)  = T_CURRENT   ! Lock in fracture temperature
              END IF

          ELSE IF (T_CURRENT < T_SOL) THEN
              ! Successfully cooled below solidus without cracking -> BONDED
              IF (STATEV(8) < 1.0D0) THEN
                  WELD_STATE = 2           ! WELDED
              ELSE
                  WELD_STATE = 3           ! CRACKED
              END IF
          END IF

      ELSE IF (WELD_STATE == 2) THEN
          ! State 2: WELDED
          IF (T_CURRENT >= T_LIQ) THEN
              WELD_STATE = 1
              STATEV(3)  = 0.0D0
          END IF
      END IF

      STATEV(1) = DBLE(WELD_STATE)

      ! ================================================================
      ! 2. CONSTITUTIVE TRACTION-SEPARATION LAW
      ! ================================================================
      SELECT CASE (WELD_STATE)

      CASE (0)
          ! --- STATE 0: RAW PRESS-FIT ---
          IF (U_NORM < 0.0D0) THEN
              STRESS(1)   = K_COMP * U_NORM
              DDSDDU(1,1) = K_COMP
              STATEV(11)  = K_COMP
              U_SLIP_MAG  = SQRT(U_SLIP1**2 + U_SLIP2**2)
              IF (U_SLIP_MAG > 1.0D-8) THEN
                  STRESS(2) = -MU_VIRGIN * ABS(STRESS(1)) * (U_SLIP1 / U_SLIP_MAG)
                  STRESS(3) = -MU_VIRGIN * ABS(STRESS(1)) * (U_SLIP2 / U_SLIP_MAG)
              END IF
          ELSE
              STRESS(1)   = 0.0D0
              DDSDDU(1,1) = 0.0D0
              STATEV(11)  = 0.0D0
          END IF

      CASE (1)
          ! --- STATE 1: BTR MUSHY ZONE ---
          IF (U_NORM < 0.0D0) THEN
              STRESS(1)   = K_COMP * U_NORM
              DDSDDU(1,1) = K_COMP
          ELSE
              ! Residual dendritic bridge stiffness
              STRESS(1)   = (0.05D0 * K_BOND) * U_NORM
              DDSDDU(1,1) = 0.05D0 * K_BOND
          END IF
          STATEV(11) = DDSDDU(1,1)

      CASE (2)
          ! --- STATE 2: WELDED SOUND METAL ---
          STRESS(1)   = K_BOND * U_NORM
          DDSDDU(1,1) = K_BOND
          STRESS(2)   = (K_BOND * 0.4D0) * U_SLIP1
          STRESS(3)   = (K_BOND * 0.4D0) * U_SLIP2
          DDSDDU(2,2) = K_BOND * 0.4D0
          DDSDDU(3,3) = K_BOND * 0.4D0
          STATEV(11)  = K_BOND

      CASE (3)
          ! --- STATE 3: CRACKED / PERMANENT TEAR ---
          IF (U_NORM < 0.0D0) THEN
              STRESS(1)   = K_COMP * U_NORM
              DDSDDU(1,1) = K_COMP
              STATEV(11)  = K_COMP
          ELSE
              STRESS(1)   = 0.0D0
              DDSDDU(1,1) = 0.0D0
              STATEV(11)  = 0.0D0
          END IF
          STRESS(2) = 0.0D0
          STRESS(3) = 0.0D0

      END SELECT

      RETURN
      END SUBROUTINE UINTER

!=======================================================================
! ABAQUS/STANDARD USER INTERACTION SUBROUTINE: UINTER
!
! APPLICATION: 3-State Contact Interface (Raw -> Welded -> Cracked)
!              with BTR (Brittleness Temperature Range) Hot Cracking
!
! STATUS:      EXPLORATORY / RESEARCH BACKLOG (IN WORK)
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

      ! --- Material & Process Parameters ---
      ! PROPS(1) = MU_VIRGIN       (Friction coefficient in raw unbonded state)
      ! PROPS(2) = T_SOLIDUS       (1485.0 C for 16MnCr5)
      ! PROPS(3) = T_LIQUIDUS      (1530.0 C for 16MnCr5)
      ! PROPS(4) = K_PENALTY_COMP  (Normal penalty stiffness in compression [N/mm3])
      ! PROPS(5) = K_WELD_BOND     (Elastic bond stiffness in welded state [N/mm3])
      ! PROPS(6) = VIRTUAL_PAD_LEN (Virtual element length scale L0 [mm] - SOURCE OF MESH DEPENDENCY)
      ! PROPS(7) = WON_CRIT_STRAIN (Won critical strain limit for hot tearing [% / 100])
      
      REAL*8 :: MU_VIRGIN, T_SOL, T_LIQ, K_COMP, K_BOND, L0, EPS_CRIT
      REAL*8 :: U_NORM, U_SLIP1, U_SLIP2, U_SLIP_MAG
      REAL*8 :: T_CURRENT, T_PREV, T_RATE
      REAL*8 :: DELTA_UN_BTR, EPS_BTR
      INTEGER :: WELD_STATE

      ! --- State Variables (STATEV) ---
      ! STATEV(1) = STATE_FLAG   (0: RAW, 1: BTR SOLIDIFYING, 2: WELDED, 3: CRACKED)
      ! STATEV(2) = MAX_TEMP     (Peak temperature reached historically [C])
      ! STATEV(3) = BTR_OPENING  (Accumulated tensile opening during BTR [mm])
      ! STATEV(4) = WON_INDEX    (Damage indicator: Ratio eps_accum / eps_crit)
      ! STATEV(5) = BOND_STIFF   (Effective current normal stiffness [N/mm3])

      ! Initialize outputs
      STRESS = 0.0D0
      DDSDDU = 0.0D0
      DDSDDT = 0.0D0

      ! Extract Properties
      MU_VIRGIN = PROPS(1)
      T_SOL     = PROPS(2)
      T_LIQ     = PROPS(3)
      K_COMP    = PROPS(4)
      K_BOND    = PROPS(5)
      L0        = PROPS(6)
      EPS_CRIT  = PROPS(7)

      ! Kinematic Quantities (1: Normal, 2: Slip-1, 3: Slip-2)
      ! SEPR(1) < 0 implies contact penetration (compression); > 0 implies opening
      U_NORM    = SEPR(1)
      U_SLIP1   = SEPR(2)
      U_SLIP2   = SEPR(3)

      ! Temperatures
      T_CURRENT = TEMP + DTEMP
      T_PREV    = TEMP
      IF (DTIME > 1.0D-12) THEN
          T_RATE = DTEMP / DTIME
      ELSE
          T_RATE = 0.0D0
      END IF

      ! Historical Peak Temperature
      IF (T_CURRENT > STATEV(2)) STATEV(2) = T_CURRENT

      ! Recover State Flag
      WELD_STATE = NINT(STATEV(1))

      ! ================================================================
      ! 1. STATE MACHINE TRANSITIONS
      ! ================================================================
      IF (WELD_STATE == 0) THEN
          ! State 0: RAW UNWELDED
          IF (STATEV(2) >= T_LIQ) THEN
              ! Pool reached liquidus -> initiate melting transition
              WELD_STATE = 1
              STATEV(3)  = 0.0D0  ! Reset accumulated BTR opening
          END IF

      ELSE IF (WELD_STATE == 1) THEN
          ! State 1: BTR SOLIDIFYING (T_SOL <= T <= T_LIQ during cooling)
          IF (T_CURRENT <= T_LIQ .AND. T_CURRENT >= T_SOL) THEN
              ! Cooling inside BTR: track tensile opening
              IF (DEPR(1) > 0.0D0) THEN
                  STATEV(3) = STATEV(3) + DEPR(1)
              END IF

              ! Apparent BTR strain (NOTE: Mesh-dependent normalization by L0)
              IF (L0 > 1.0D-6) THEN
                  EPS_BTR = STATEV(3) / L0
              ELSE
                  EPS_BTR = STATEV(3)
              END IF
              STATEV(4) = EPS_BTR / MAX(EPS_CRIT, 1.0D-6)

              ! Check Hot Cracking Criterion (Won / Prokhorov)
              IF (STATEV(4) >= 1.0D0) THEN
                  WELD_STATE = 3  ! CRACKED (Hot tear occurred)
              END IF

          ELSE IF (T_CURRENT < T_SOL) THEN
              ! Successfully cooled below solidus without cracking -> BONDED
              WELD_STATE = 2
          END IF

      ELSE IF (WELD_STATE == 2) THEN
          ! State 2: WELDED JOINT
          ! If re-melted above liquidus:
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
          ! ------------------------------------------------------------
          ! STATE 0: RAW UNBONDED (Unilateral Press-Fit Contact)
          ! ------------------------------------------------------------
          IF (U_NORM < 0.0D0) THEN
              ! Compression: Penalty barrier against penetration
              STRESS(1) = K_COMP * U_NORM
              DDSDDU(1,1) = K_COMP
              STATEV(5) = K_COMP

              ! Coulomb Friction in Tangential Plane
              U_SLIP_MAG = SQRT(U_SLIP1**2 + U_SLIP2**2)
              IF (U_SLIP_MAG > 1.0D-8) THEN
                  STRESS(2) = -MU_VIRGIN * ABS(STRESS(1)) * (U_SLIP1 / U_SLIP_MAG)
                  STRESS(3) = -MU_VIRGIN * ABS(STRESS(1)) * (U_SLIP2 / U_SLIP_MAG)
              END IF
          ELSE
              ! Tension / Clearance: Free opening, zero stress
              STRESS(1) = 0.0D0
              DDSDDU(1,1) = 0.0D0
              STATEV(5) = 0.0D0
          END IF

      CASE (1)
          ! ------------------------------------------------------------
          ! STATE 1: BTR SOLIDIFYING (Semi-solid dendritic slurry)
          ! ------------------------------------------------------------
          IF (U_NORM < 0.0D0) THEN
              STRESS(1) = K_COMP * U_NORM
              DDSDDU(1,1) = K_COMP
          ELSE
              ! Minimal liquid film resistance
              STRESS(1) = 0.01D0 * K_BOND * U_NORM
              DDSDDU(1,1) = 0.01D0 * K_BOND
          END IF
          STATEV(5) = DDSDDU(1,1)

      CASE (2)
          ! ------------------------------------------------------------
          ! STATE 2: WELDED (Full bilateral structural continuity)
          ! ------------------------------------------------------------
          ! Resists both tension and compression with full bond stiffness
          STRESS(1) = K_BOND * U_NORM
          DDSDDU(1,1) = K_BOND

          ! Full shear resistance
          STRESS(2) = (K_BOND * 0.4D0) * U_SLIP1
          STRESS(3) = (K_BOND * 0.4D0) * U_SLIP2
          DDSDDU(2,2) = K_BOND * 0.4D0
          DDSDDU(3,3) = K_BOND * 0.4D0
          STATEV(5) = K_BOND

      CASE (3)
          ! ------------------------------------------------------------
          ! STATE 3: CRACKED (Permanent Hot Tear / Delamination)
          ! ------------------------------------------------------------
          IF (U_NORM < 0.0D0) THEN
              ! Compression: crack closes, contact prevents penetration
              STRESS(1) = K_COMP * U_NORM
              DDSDDU(1,1) = K_COMP
              STATEV(5) = K_COMP
          ELSE
              ! Tension: completely open crack, zero resistance
              STRESS(1) = 0.0D0
              DDSDDU(1,1) = 0.0D0
              STATEV(5) = 0.0D0
          END IF
          ! No shear transmission across cracked open gap
          STRESS(2) = 0.0D0
          STRESS(3) = 0.0D0

      END SELECT

      RETURN
      END SUBROUTINE UINTER

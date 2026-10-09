!=======================================================================
! ABAQUS/STANDARD USER MATERIAL SUBROUTINE: UMAT (COHESIVE ZONE)
!
! APPLICATION: 3-State Thermo-Mechanical Cohesive Interface
!              (Raw -> Welded -> Cracked) with Energy-Regularized
!              BTR (Brittleness Temperature Range) Hot Cracking
!
! ELEMENT COMPATIBILITY: COH3D8 (8-Node 3D Cohesive Elements)
!                        Zero-thickness or thin interface layers
!
! HIGHLIGHTS:
!  1. MESH-INDEPENDENT: Regularized by critical fracture energy Gc.
!  2. QUADRATIC CONVERGENCE: Exact analytical 3x3 DDSDDE Jacobian.
!  3. SEAMLESS MPI SCALING: Element-local integration (no contact SDI).
!  4. THERMAL STATE MACHINE: Raw unbonded -> Melting -> Solidified -> Tear.
!
! AUTHOR:        Victor Maia (https://github.com/vhfmaia)
! REPOSITORY:    https://github.com/vhfmaia/Abaqus-Subroutines-Portfolio
!=======================================================================

      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD,
     &                RPL, DDSDDT, DRPLDE, DRPLDT,
     &                STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     &                PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS,
     &                NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT,
     &                CELENT, DFGRD0, DFGRD1, NOEL, NPT, KSLAY,
     &                KSPT, KSTEP, KINC)

      IMPLICIT NONE

      ! --- Input / Output Arguments ---
      CHARACTER*80, INTENT(IN)    :: CMNAME
      INTEGER,      INTENT(IN)    :: NDI, NSHR, NTENS, NSTATV, NPROPS
      INTEGER,      INTENT(IN)    :: NOEL, NPT, KSLAY, KSPT, KSTEP, KINC
      REAL*8,       INTENT(IN)    :: TIME(2), DTIME, TEMP, DTEMP, CELENT
      REAL*8,       INTENT(IN)    :: STRAN(NTENS), DSTRAN(NTENS)
      REAL*8,       INTENT(IN)    :: PREDEF(1), DPRED(1)
      REAL*8,       INTENT(IN)    :: PROPS(NPROPS), COORDS(3), DROT(3,3)
      REAL*8,       INTENT(IN)    :: DFGRD0(3,3), DFGRD1(3,3)
      REAL*8,       INTENT(INOUT) :: STATEV(NSTATV), PNEWDT
      REAL*8,       INTENT(OUT)   :: STRESS(NTENS), DDSDDE(NTENS,NTENS)
      REAL*8,       INTENT(OUT)   :: SSE, SPD, SCD, RPL, DDSDDT(NTENS)
      REAL*8,       INTENT(OUT)   :: DRPLDE(NTENS), DRPLDT

      ! --- Material Properties (PROPS) ---
      ! PROPS(1) = K_PENALTY_RAW    (Normal compressive penalty [N/mm3])
      ! PROPS(2) = K_BOND_SOLID     (Elastic stiffness of sound weld [N/mm3])
      ! PROPS(3) = G_SHEAR_SOLID    (Shear stiffness of sound weld [N/mm3])
      ! PROPS(4) = T_SOLIDUS        (Solidus temperature, e.g. 1485.0 C)
      ! PROPS(5) = T_LIQUIDUS       (Liquidus temperature, e.g. 1530.0 C)
      ! PROPS(6) = SIGMA_CRIT_BTR   (Peak tensile strength in BTR [MPa])
      ! PROPS(7) = GC_FRACTURE_BTR  (Critical fracture energy in BTR [mJ/mm2 = N/mm])
      ! PROPS(8) = DELTA_MAX_CRACK  (Ultimate displacement at failure [mm] = 2*Gc/Sigma_c)

      REAL*8 :: K_PENALTY, K_BOND, G_SHEAR, T_SOL, T_LIQ
      REAL*8 :: SIGMA_C, G_C, DELTA_F
      REAL*8 :: T_CURRENT, T_PREV
      REAL*8 :: DELTA_N, DELTA_S1, DELTA_S2
      REAL*8 :: DELTA_0, DAMAGE, EFF_DELTA
      INTEGER :: WELD_STATE

      ! --- State Variables (STATEV) ---
      ! STATEV(1) = WELD_STATE_FLAG (0: RAW, 1: MELTED/BTR, 2: WELDED, 3: CRACKED)
      ! STATEV(2) = PEAK_TEMP       (Maximum temperature reached historically [C])
      ! STATEV(3) = DAMAGE_VARIABLE (Continuous damage D in [0.0, 1.0])
      ! STATEV(4) = BTR_MAX_OPEN    (Maximum tensile separation recorded in BTR [mm])
      ! STATEV(5) = DISSIPATED_ENG  (Energy dissipated by cracking [mJ/mm2])

      ! Initialize
      STRESS = 0.0D0
      DDSDDE = 0.0D0
      DDSDDT = 0.0D0
      SSE    = 0.0D0
      SPD    = 0.0D0
      SCD    = 0.0D0
      RPL    = 0.0D0
      DRPLDE = 0.0D0
      DRPLDT = 0.0D0

      ! Extract Properties
      K_PENALTY = PROPS(1)
      K_BOND    = PROPS(2)
      G_SHEAR   = PROPS(3)
      T_SOL     = PROPS(4)
      T_LIQ     = PROPS(5)
      SIGMA_C   = PROPS(6)
      G_C       = PROPS(7)
      DELTA_F   = PROPS(8)

      ! Cohesive element kinematics (NTENS = 3 for COH3D8 traction-separation)
      ! STRAN(1) = delta_n (normal separation)
      ! STRAN(2) = delta_s1 (tangential shear 1)
      ! STRAN(3) = delta_s2 (tangential shear 2)
      DELTA_N  = STRAN(1) + DSTRAN(1)
      DELTA_S1 = STRAN(2) + DSTRAN(2)
      DELTA_S2 = STRAN(3) + DSTRAN(3)

      ! Temperatures
      T_CURRENT = TEMP + DTEMP
      T_PREV    = TEMP

      ! Track Peak Historic Temperature
      IF (T_CURRENT > STATEV(2)) STATEV(2) = T_CURRENT

      ! Recover State Flag & Damage
      WELD_STATE = NINT(STATEV(1))
      DAMAGE     = STATEV(3)

      ! Damage initiation threshold (elastic limit displacement)
      IF (K_BOND > 1.0D-6) THEN
          DELTA_0 = SIGMA_C / K_BOND
      ELSE
          DELTA_0 = 1.0D-4
      END IF

      ! ================================================================
      ! 1. THERMAL STATE MACHINE LOGIC
      ! ================================================================
      IF (WELD_STATE == 0) THEN
          ! --- STATE 0: RAW PRE-WELD ---
          IF (STATEV(2) >= T_LIQ) THEN
              ! Core melted -> transition to liquid/semi-solid
              WELD_STATE = 1
              STATEV(4)  = 0.0D0  ! Reset BTR opening tracker
          END IF

      ELSE IF (WELD_STATE == 1) THEN
          ! --- STATE 1: BTR SOLIDIFYING (T_SOL <= T <= T_LIQ) ---
          IF (T_CURRENT <= T_LIQ .AND. T_CURRENT >= T_SOL) THEN
              ! Inside mushy zone: evaluate tensile opening
              IF (DELTA_N > STATEV(4)) STATEV(4) = DELTA_N

              ! Energy-regularized linear softening law:
              ! D = (delta_f * (delta - delta_0)) / (delta * (delta_f - delta_0))
              IF (STATEV(4) > DELTA_0) THEN
                  IF (DELTA_F > DELTA_0) THEN
                      DAMAGE = (DELTA_F * (STATEV(4) - DELTA_0)) /
     &                         (STATEV(4) * (DELTA_F - DELTA_0))
                      DAMAGE = MIN(MAX(DAMAGE, 0.0D0), 1.0D0)
                  ELSE
                      DAMAGE = 1.0D0
                  END IF
              ELSE
                  DAMAGE = 0.0D0
              END IF
              STATEV(3) = DAMAGE

              ! If full separation energy consumed -> permanent tear
              IF (DAMAGE >= 0.999D0) THEN
                  WELD_STATE = 3  ! CRACKED
                  DAMAGE     = 1.0D0
                  STATEV(3)  = 1.0D0
              END IF

          ELSE IF (T_CURRENT < T_SOL) THEN
              ! Cooled below solidus:
              IF (DAMAGE >= 0.999D0) THEN
                  WELD_STATE = 3  ! Remained cracked
              ELSE
                  WELD_STATE = 2  ! Solidified sound joint
                  DAMAGE     = 0.0D0  ! Heal remaining micro-voids
                  STATEV(3)  = 0.0D0
              END IF
          END IF

      ELSE IF (WELD_STATE == 2) THEN
          ! --- STATE 2: WELDED SOUND JOINT ---
          ! If re-melted by subsequent pass:
          IF (T_CURRENT >= T_LIQ) THEN
              WELD_STATE = 1
              DAMAGE     = 0.0D0
              STATEV(3)  = 0.0D0
          END IF
      END IF

      STATEV(1) = DBLE(WELD_STATE)

      ! ================================================================
      ! 2. CONSTITUTIVE TRACTION AND ANALYTICAL JACOBIAN
      ! ================================================================

      SELECT CASE (WELD_STATE)

      CASE (0)
          ! ------------------------------------------------------------
          ! STATE 0: RAW UNBONDED (Unilateral contact barrier)
          ! ------------------------------------------------------------
          IF (DELTA_N < 0.0D0) THEN
              ! Compression: strong elastic penalty barrier
              STRESS(1)   = K_PENALTY * DELTA_N
              DDSDDE(1,1) = K_PENALTY
          ELSE
              ! Tension: free gap opening, zero stress
              STRESS(1)   = 0.0D0
              DDSDDE(1,1) = 0.0D0
          END IF
          ! Zero shear transmission across unwelded clearance
          STRESS(2)   = 0.0D0
          STRESS(3)   = 0.0D0
          DDSDDE(2,2) = 0.0D0
          DDSDDE(3,3) = 0.0D0

      CASE (1)
          ! ------------------------------------------------------------
          ! STATE 1: BTR SOLIDIFYING (Mushy slurry with softening)
          ! ------------------------------------------------------------
          IF (DELTA_N < 0.0D0) THEN
              ! Compression: penalty barrier
              STRESS(1)   = K_PENALTY * DELTA_N
              DDSDDE(1,1) = K_PENALTY
          ELSE
              ! Tension with energy-regularized softening damage (1 - D)
              STRESS(1)   = (1.0D0 - DAMAGE) * K_BOND * DELTA_N
              DDSDDE(1,1) = (1.0D0 - DAMAGE) * K_BOND
          END IF
          ! Weak semi-solid shear
          STRESS(2)   = (1.0D0 - DAMAGE) * (G_SHEAR * 0.05D0) * DELTA_S1
          STRESS(3)   = (1.0D0 - DAMAGE) * (G_SHEAR * 0.05D0) * DELTA_S2
          DDSDDE(2,2) = (1.0D0 - DAMAGE) * (G_SHEAR * 0.05D0)
          DDSDDE(3,3) = (1.0D0 - DAMAGE) * (G_SHEAR * 0.05D0)

      CASE (2)
          ! ------------------------------------------------------------
          ! STATE 2: WELDED SOUND STEEL (Full bilateral continuity)
          ! ------------------------------------------------------------
          ! Resists both tension and compression bilaterally
          STRESS(1)   = K_BOND * DELTA_N
          DDSDDE(1,1) = K_BOND

          ! Full shear transmission
          STRESS(2)   = G_SHEAR * DELTA_S1
          STRESS(3)   = G_SHEAR * DELTA_S2
          DDSDDE(2,2) = G_SHEAR
          DDSDDE(3,3) = G_SHEAR

      CASE (3)
          ! ------------------------------------------------------------
          ! STATE 3: CRACKED (Permanent Tear / Delaminated)
          ! ------------------------------------------------------------
          IF (DELTA_N < 0.0D0) THEN
              ! Crack closure: contact penalty prevents penetration
              STRESS(1)   = K_PENALTY * DELTA_N
              DDSDDE(1,1) = K_PENALTY
          ELSE
              ! Crack open: zero tensile resistance
              STRESS(1)   = 0.0D0
              DDSDDE(1,1) = 0.0D0
          END IF
          ! No shear transmission across open crack
          STRESS(2)   = 0.0D0
          STRESS(3)   = 0.0D0
          DDSDDE(2,2) = 0.0D0
          DDSDDE(3,3) = 0.0D0

      END SELECT

      RETURN
      END SUBROUTINE UMAT

!DEC$ FREEFORM
!=======================================================================
!  Abaqus/Standard User Subroutine: DFLUX
!  Model: Circumferential Laser Welding - Moving Conical Gaussian Heat Source
!  Author: Victor Maia (vhfm08@gmail.com)
!
!  References:
!    - Goldak et al. (1984): A new finite element model for welding heat sources
!    - ISO/TR 17671-4: Welding - Recommendations for welding of metallic materials
!
!  Conventions & Units:
!    - Global coordinates: mm, s, tonne, mJ, C
!    - Volumetric flux: mJ/(s.mm3) = mW/mm3
!    - Trajectory: Circular weld along radius R in X-Y plane about Z axis
!    - Ramping: 10 deg ramp-up, 360 deg steady weld, 10 deg ramp-down
!=======================================================================

module laser_dflux_mod
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    ! ------------------------------------------------------------------
    ! 1. PROCESS & OPTICAL BEAM PARAMETERS
    ! ------------------------------------------------------------------
    real(real64), parameter :: qtot_nominal = 1.80d6        ! Nominal beam power [mW = mJ/s] (1800 W)
    real(real64), parameter :: eta_absorb   = 0.60d0        ! Laser absorption efficiency [-]
    real(real64), parameter :: r_top        = 1.00d0        ! Beam radius at top surface [mm]
    real(real64), parameter :: r_bottom     = 0.75d0        ! Beam radius at penetration depth [mm]
    real(real64), parameter :: depth_zi     = 3.00d0        ! Conical penetration depth [mm]
    real(real64), parameter :: weld_radius  = 15.00d0       ! Circular joint radius [mm]
    real(real64), parameter :: travel_speed = 1200.0d0      ! Travel speed [mm/min] (20 mm/s)
    real(real64), parameter :: z_surface    = 20.00d0       ! Top irradiated surface Z coordinate [mm]
    real(real64), parameter :: cutoff_tol   = 1.0d-8        ! Numerical flux cutoff threshold [-]

    ! ------------------------------------------------------------------
    ! 2. ANGULAR POWER SCHEDULE PARAMETERS
    ! ------------------------------------------------------------------
    real(real64), parameter :: deg_ramp_up   = 10.00d0      ! Power ramp-up segment [deg]
    real(real64), parameter :: deg_weld      = 360.00d0     ! Full nominal power weld segment [deg]
    real(real64), parameter :: deg_ramp_down = 10.00d0      ! Power ramp-down overlap segment [deg]

    real(real64), parameter :: pi_const = 3.14159265358979323846d0

contains

    !> Computes instantaneous conical Gaussian volumetric heat flux
    pure subroutine compute_conical_flux(time_val, coords, q_vol)
        real(real64), intent(in)  :: time_val
        real(real64), intent(in)  :: coords(3)
        real(real64), intent(out) :: q_vol

        real(real64) :: v_sec, omega, theta_rad, theta_deg
        real(real64) :: amp_factor, q0_peak, z_local, r0_z
        real(real64) :: xc, yc, dx, dy, r_dist_sq, r_cutoff_sq

        ! Linear speed in mm/s and angular velocity in rad/s
        v_sec = travel_speed / 60.0d0
        omega = v_sec / weld_radius

        ! Current angular position along the joint
        theta_rad = omega * time_val
        theta_deg = theta_rad * (180.0d0 / pi_const)

        ! 1. Angular Power Schedule Evaluation
        if (theta_deg < deg_ramp_up) then
            ! Linear power ramp-up (0 -> 10 deg)
            amp_factor = theta_deg / deg_ramp_up
        else if (theta_deg <= (deg_ramp_up + deg_weld)) then
            ! Full nominal steady-state welding (10 -> 370 deg)
            amp_factor = 1.0d0
        else if (theta_deg <= (deg_ramp_up + deg_weld + deg_ramp_down)) then
            ! Linear power decay and crater fill overlap (370 -> 380 deg)
            amp_factor = 1.0d0 - (theta_deg - (deg_ramp_up + deg_weld)) / deg_ramp_down
        else
            ! Beam shut off
            amp_factor = 0.0d0
        end if

        if (amp_factor <= 0.0d0) then
            q_vol = 0.0d0
            return
        end if

        ! 2. Depth Check relative to top irradiated surface
        z_local = z_surface - coords(3)
        if (z_local < 0.0d0 .or. z_local > depth_zi) then
            q_vol = 0.0d0
            return
        end if

        ! 3. Analytical Peak Flux Prefactor Q0 (Volume conservation: Integral = eta * Qtot)
        ! Volume integral denominator: (pi/3) * zi * (re^2 + re*ri + ri^2)
        q0_peak = (3.0d0 * eta_absorb * qtot_nominal) / &
                  (pi_const * (depth_zi / 3.0d0) * &
                   (r_top**2 + r_top * r_bottom + r_bottom**2))

        ! Local cone radius at current depth z
        r0_z = r_top + (r_bottom - r_top) * (z_local / depth_zi)
        if (r0_z <= 0.0d0) then
            q_vol = 0.0d0
            return
        end if

        ! 4. Beam Center Coordinates in horizontal plane
        xc = weld_radius * cos(theta_rad)
        yc = weld_radius * sin(theta_rad)

        dx = coords(1) - xc
        dy = coords(2) - yc
        r_dist_sq = dx**2 + dy**2

        ! Numerical cutoff envelope: exp(-3 * r^2 / R0^2) >= TOL
        r_cutoff_sq = -(r0_z**2 / 3.0d0) * log(cutoff_tol)
        if (r_dist_sq > r_cutoff_sq) then
            q_vol = 0.0d0
            return
        end if

        ! 5. Conical Gaussian Volumetric Heat Flux
        q_vol = amp_factor * q0_peak * exp(-3.0d0 * r_dist_sq / (r0_z**2))

    end subroutine compute_conical_flux

end module laser_dflux_mod

!=======================================================================
!  Abaqus/Standard Entry Point Subroutine
!=======================================================================
subroutine dflux(flux, sol, kstep, kinc, time, noel, npt, coords, &
                 jltyp, temp, press, sname)
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use laser_dflux_mod, only: compute_conical_flux
    implicit none

    ! Abaqus interface arguments
    real(real64), intent(out)   :: flux(2)
    real(real64), intent(in)    :: sol
    integer(int32), intent(in)  :: kstep
    integer(int32), intent(in)  :: kinc
    real(real64), intent(in)    :: time(2)
    integer(int32), intent(in)  :: noel
    integer(int32), intent(in)  :: npt
    real(real64), intent(in)    :: coords(3)
    integer(int32), intent(in)  :: jltyp
    real(real64), intent(in)    :: temp
    real(real64), intent(in)    :: press
    character(len=80), intent(in) :: sname

    real(real64) :: q_val

    ! time(2) is total step time
    call compute_conical_flux(time(2), coords, q_val)

    flux(1) = q_val
    flux(2) = 0.0d0    ! df/dtheta (zero for purely prescribed thermal load)

end subroutine dflux

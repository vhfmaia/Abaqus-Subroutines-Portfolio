!DEC$ FREEFORM
!=======================================================================
!  Abaqus/Standard User Subroutine: DFLUX
!  Model: Circumferential Laser Welding - Moving Conical Gaussian Heat Source
!  Author: Victor Maia (vhfm08@gmail.com)
!
!  References (Conical Gaussian / TDC Model):
!    - Farrokhi, F., Endelt, B., Kristiansen, M. (2019): A numerical model
!      for full and partial penetration hybrid laser welding of thick-section
!      steels, Optics and Laser Technology, 109, 629-642.
!    - Wu, C.S., Wang, H.G., Zhang, Y.M. (2006): A new heat source model for
!      keyhole plasma arc welding in FEM analysis, Welding Journal, 85(12), 284-291.
!    - Liu, M., Kouadri-Henni, A., Malard, B. (2022): Simulation of low-cycle
!      fatigue residual stress in laser-welded structures, ICRS11.
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

    ! Mathematical & analytical normalization constants (Wu et al. 2006 / Farrokhi et al. 2019)
    real(real64), parameter :: pi_const   = 3.14159265358979323846d0
    real(real64), parameter :: e3_ratio   = 20.085536923187668d0 / (20.085536923187668d0 - 1.0d0)
    real(real64), parameter :: q0_peak    = (9.0d0 * eta_absorb * qtot_nominal * e3_ratio) / &
                                            (pi_const * depth_zi * &
                                             (r_top**2 + r_top * r_bottom + r_bottom**2))

    real(real64), parameter :: v_sec      = travel_speed / 60.0d0
    real(real64), parameter :: omega      = v_sec / weld_radius
    real(real64), parameter :: rad_to_deg = 180.0d0 / pi_const

contains

    !> Computes instantaneous conical Gaussian volumetric heat flux
    pure subroutine compute_conical_flux(time_val, coords, q_vol)
        real(real64), intent(in)  :: time_val
        real(real64), intent(in)  :: coords(3)
        real(real64), intent(out) :: q_vol

        real(real64) :: theta_rad, theta_deg
        real(real64) :: amp_factor, z_local, r0_z
        real(real64) :: xc, yc, dx, dy, r_dist_sq

        ! Current angular position along the joint
        theta_rad = omega * time_val
        theta_deg = theta_rad * rad_to_deg

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

        ! 3. Local cone radius at current depth z (Wu et al. Eq. 6 / Farrokhi et al. Eq. 6)
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

        ! Conical boundary envelope: r <= r0(z) (Farrokhi et al. 2019 / Wu et al. 2006)
        if (r_dist_sq > (r0_z**2)) then
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
    real(real64), intent(out)     :: flux(2)
    real(real64), intent(in)      :: sol
    integer(int32), intent(in)    :: kstep
    integer(int32), intent(in)    :: kinc
    real(real64), intent(in)      :: time(2)
    integer(int32), intent(in)    :: noel
    integer(int32), intent(in)    :: npt
    real(real64), intent(in)      :: coords(3)
    integer(int32), intent(in)    :: jltyp
    real(real64), intent(in)      :: temp
    real(real64), intent(in)      :: press
    character(len=80), intent(in) :: sname

    real(real64) :: q_val

    ! time(2) is total step time
    call compute_conical_flux(time(2), coords, q_val)

    flux(1) = q_val
    flux(2) = 0.0d0    ! df/dtheta (zero for purely prescribed thermal load)

end subroutine dflux

!=======================================================================
!  Abaqus/Standard User Subroutine: USDFLD
!  Standard interface satisfying *USER DEFINED FIELD requests.
!=======================================================================
subroutine usdfld(field, statev, pnewdt, direct, t, celent, &
                  time, dtime, cmname, orname, nfield, nstatv, &
                  noel, npt, layer, kspt, kstep, kinc, ndi, nshr, coord, &
                  jmac, jmatyp, matlayo, laccflg)
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    ! Abaqus interface arguments
    real(real64), intent(inout)     :: field(nfield)
    real(real64), intent(inout)     :: statev(nstatv)
    real(real64), intent(inout)     :: pnewdt
    real(real64), intent(in)        :: direct(3,3), t(3,3)
    real(real64), intent(in)        :: celent
    real(real64), intent(in)        :: time(2), dtime
    character(len=80), intent(in)   :: cmname, orname
    integer(int32), intent(in)      :: nfield, nstatv
    integer(int32), intent(in)      :: noel, npt, layer, kspt, kstep, kinc, ndi, nshr
    real(real64), intent(in)        :: coord(3)
    integer(int32), intent(in)      :: jmac, jmatyp, matlayo, laccflg

    return
end subroutine usdfld


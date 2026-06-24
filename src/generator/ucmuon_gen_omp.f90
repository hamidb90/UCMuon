!=============================================================================
! ucmuon_gen_omp.f90  – CosmoALEPH Muon Generator (OpenMP version)
! UCLouvain Muography Group
!
! Changes from v10:
!   1. uses rng_parallel  (par_ranlux replaces RANLUX everywhere)
!   2. uses omp_lib
!   3. removed: external :: RLUXGO, external :: RANLUX
!   4. RNG init: par_init_rng(iranlux) instead of RLUXGO(3,iranlux,0,0)
!   5. Main loop: !$OMP PARALLEL + infinite DO/EXIT (guarantees exactly nmuons saved)
!   6. Removed max_tries fixed ceiling — loop runs until i == nmuons regardless of acceptance rate
!=============================================================================
program ucmuon_gen_omp
  use ucmuon_source_module
  use geom_module
  use phits_module
  use rng_parallel        ! <-- OMP: thread-safe RNG
  use omp_lib             ! <-- OMP: thread count, get_thread_num
  use parma_path          ! configurable PARMA data directory
  implicit none

  !---------------------------------------------------------------------------
  ! External PARMA functions from parma_subroutines.f90
  !---------------------------------------------------------------------------
  real(8), external :: getMuonSpec
  real(8), external :: getSpecAngFinal
  real(8), external :: getd
  real(8), external :: getr
  real(8), external :: getFFPfromW
  real(8), external :: getHP

  !---------------------------------------------------------------------------
  ! Generator parameters
  !---------------------------------------------------------------------------
  real(8) :: e_min, e_max, p_min, p_max

  ! --- Source geometry ---
  integer :: source_mode   ! 1=disc  2=rectangle  3=hemisphere
  integer :: source_plane  ! 1=XY  2=XZ  3=YZ  (flat sources only)

  ! Disk source (source_mode == 1): center + radius + tilt
  real(8) :: src_disk_cx_m, src_disk_cy_m, src_disk_r_m   ! center [m], radius [m]
  real(8) :: src_tilt_deg, src_tilt_az_deg                 ! tilt from plane normal [deg], azimuth [deg]
  real(8) :: src_disk_cx_cm, src_disk_cy_cm, src_disk_r_cm ! same in cm
  real(8) :: src_tilt_rad, src_tilt_az_rad                 ! tilt angles [rad]

  ! Rectangle source (source_mode == 2): bounding box in metres
  real(8) :: src_u1_m, src_u2_m   ! first-axis range  [m]
  real(8) :: src_v1_m, src_v2_m   ! second-axis range [m]
  real(8) :: src_w_m              ! fixed coordinate  [m]

  ! Derived working quantities (cm, set in DERIVED QUANTITIES section)
  real(8) :: half_lu_cm, half_lv_cm  ! half-extents [cm]
  real(8) :: centre_u_cm, centre_v_cm, src_w_cm

  ! Hemisphere uses separate radius + centre-z (plane selector not applicable)
  real(8) :: hemi_radius_m, hemi_radius_cm, hemi_cz_m, hemi_cz_cm

  ! Kept for compatibility with generate_muon / sample_position call signatures
  real(8) :: radius_cm, half_lx_cm, half_ly_cm, source_z_cm

  ! Temporary scalars for plane-rotation swap
  real(8) :: tmp_pos, tmp_dir
  real(8) :: theta_max, theta_max_deg
  integer :: nmuons, angular_mode
  integer :: spectrum_mode_in

  !---------------------------------------------------------------------------
  ! PARMA mode parameters
  !---------------------------------------------------------------------------
  real(8)       :: parma_lat, parma_lon, parma_alt_km
  integer       :: parma_year, parma_month, parma_day
  integer       :: parma_charge_mode
  character(200):: parma_datapath
  real(8)       :: parma_d_gcm2
  real(8)       :: parma_rc_GV
  real(8)       :: parma_s_W
  real(8)       :: parma_ffp_MV
  integer       :: parma_ic
  real(8)       :: parma_ratio_plus

  integer, parameter :: NPARMA_E = 500
  real(8) :: parma_E_MeV(NPARMA_E)
  real(8) :: parma_cdf_plus(NPARMA_E)
  real(8) :: parma_cdf_minus(NPARMA_E)

  integer, parameter :: NPARMA_ANG = 200
  real(8) :: parma_cos_arr(NPARMA_ANG)
  real(8) :: parma_cdf_ang_plus(NPARMA_ANG)
  real(8) :: parma_cdf_ang_minus(NPARMA_ANG)

  !---------------------------------------------------------------------------
  ! Detector geometry
  !---------------------------------------------------------------------------
  integer, parameter :: MAX_DET = 10
  integer            :: ndet, i_det
  type(cyl_t)        :: det_cyl_arr(MAX_DET)
  type(aabb_t)       :: det_box_arr(MAX_DET)
  integer            :: det_shape_arr(MAX_DET)
  real(8)            :: safety_margin

  !---------------------------------------------------------------------------
  ! Output flags
  !---------------------------------------------------------------------------
  integer :: use_detector, save_all, save_phits
  integer :: write_surface

  !---------------------------------------------------------------------------
  ! Per-event variables (will be PRIVATE in OMP parallel region)
  !---------------------------------------------------------------------------
  real(8)  :: x, y, z, emu, cx, cy, cz
  integer  :: muon_charge
  real(8)  :: momentum, theta, phi, px, py, pz
  real(8)  :: ekin_MeV
  integer  :: kf, hit_flag, det_mask

  logical :: hit, hit_i
  real(8) :: t_hit, t_enter, t_exit

  !---------------------------------------------------------------------------
  ! Misc — serial section
  !---------------------------------------------------------------------------
  integer  :: j
  integer(8)  :: ntry           ! total attempts (SHARED, ATOMIC CAPTURE)
  integer(8)  :: i              ! accepted count (SHARED, updated in CRITICAL)
  integer     :: tim(8), iranlux
  integer     :: use_defaults
  character(120) :: output_all, output_sel, output_phits
  real(8)        :: depth
  character(1)   :: dummy
  real(8)        :: alt_z_m, alt_p_avg, alt_L_p, alt_corr   ! altitude correction

  ! Working scalars for PARMA sampling (PRIVATE in parallel region)
  real(8) :: rnd, rnd2, frac_plus, E_sampled_MeV, cos_sampled
  real(8) :: phi_p, phi_m
  real(8) :: sum_p, sum_m, dE, dcos
  real(8) :: ang_fac                ! angular shape factor (init only)
  integer :: ibin, k                ! k = inner energy loop index (init only)
  real(4) :: rnd4

  !---------------------------------------------------------------------------
  ! OMP-specific variables
  !---------------------------------------------------------------------------
  integer(8) :: my_ntry            ! per-thread captured event ID
  integer(8) :: progress_interval  ! how often to print progress

  !===========================================================================
  write(*,*)
  write(*,*) ' ============================================================'
  write(*,*) '    UCMuon Generator - UCLouvain Muography Group              '
  write(*,*) '    ** OpenMP parallel version **                            '
  write(*,*) ' ============================================================'
  write(*,*)

  write(*,*) ' Use default parameters? (1=Yes, 0=No)'
  write(*,*) '   Defaults: E=100-2500 GeV, R=800 m, N=100000'
  write(*,*) '             Angular: cos^2(theta), theta_max=85 deg'
  write(*,*) '             No detector filter, no PHITS output'
  read(*,*) use_defaults
  write(*,*)

  !===========================================================================
  ! DEFAULTS
  !===========================================================================
  if (use_defaults == 1) then
    e_min            = 100.0d0
    e_max            = 2500.0d0
    source_mode      = 1       ! disc
    source_plane     = 1       ! XY
    src_disk_cx_m    = 0.0d0;    src_disk_cy_m = 0.0d0
    src_disk_r_m     = 800.0d0
    src_tilt_deg     = 0.0d0;    src_tilt_az_deg = 0.0d0
    src_u1_m         = -800.0d0;  src_u2_m = 800.0d0   ! rectangle defaults
    src_v1_m         = -800.0d0;  src_v2_m = 800.0d0
    src_w_m          = 0.0d0
    hemi_radius_m    = 800.0d0;   hemi_cz_m = 0.0d0
    nmuons           = 100000
    angular_mode     = 2
    theta_max        = 85.0d0 * PI / 180.0d0
    spectrum_mode_in = 1
    use_detector     = 0
    ndet             = 0
    save_all         = 1
    save_phits       = 0
    output_all       = 'muons_surface.dat'
    output_sel       = 'muons_selected.dat'
    output_phits     = 'muons_for_phits.dat'

  !===========================================================================
  ! INTERACTIVE INPUT  (unchanged from v10)
  !===========================================================================
  else

    write(*,*) ' --- [1/7] Energy Range ---'
    write(*,*) ' Min muon energy [GeV]  (e.g. 100):'
    read(*,*) e_min
    write(*,*) ' Max muon energy [GeV]  (e.g. 2500):'
    read(*,*) e_max
    write(*,*)

    write(*,*) ' --- [2/7] Spectrum Model ---'
    write(*,*) '   1 = CosmoALEPH   dN/dp ~ p^-3.1952  (Schmelling 2013)'
    write(*,*) '   2 = Power-law   dN/dE ~ E^-3.7      (legacy MUSIC cross-check)'
    write(*,*) '   3 = PARMA / EXPACS  (location & date-aware, mu+/mu- correct)'
    write(*,*) '   4 = Guan et al. (2015)  Modified Gaisser, arXiv:1509.06176'
    write(*,*) '   5 = Frosin et al. (2025)  J.Phys.G 52, 035002'
    write(*,*) '   6 = Bugaev/Gaisser (1990)  dN/dE ~ E^-2.7*(pi+K)  (pair with cos^2)'
    write(*,*) '   7 = Reyna-Bugaev (2006)  log-poly p^3*F  (pair with cos^3, best surface)'
    write(*,*) '   8 = Cosmic electrons  dN/dE ~ E^-3.0  (surface/shallow, 10 MeV-1 GeV)'
    read(*,*) spectrum_mode_in
    if (spectrum_mode_in < 1 .or. spectrum_mode_in > 8) spectrum_mode_in = 1
    write(*,*)

    if (spectrum_mode_in == 3) then
      write(*,*) ' --- [2b/7] PARMA Location & Date ---'
      write(*,*) ' Latitude  [deg, -90 to +90]:'
      read(*,*) parma_lat
      write(*,*) ' Longitude [deg, -180 to +180]:'
      read(*,*) parma_lon
      write(*,*) ' Altitude  [km, 0 = sea level]:'
      read(*,*) parma_alt_km
      write(*,*) ' Year  (e.g. 2026):'
      read(*,*) parma_year
      write(*,*) ' Month (1-12):'
      read(*,*) parma_month
      write(*,*) ' Day   (1-31):'
      read(*,*) parma_day
      write(*,*) ' Muon charge selection:  0=both   1=mu+ only   -1=mu- only'
      read(*,*) parma_charge_mode
      write(*,*) ' Path to PARMA data directory:'
      read(*,'(A)') parma_datapath
      if (len_trim(parma_datapath) == 0) parma_datapath = '.'
      write(*,*) ' W (Wolf/sunspot) index (0=solar min, ~150=solar max):'
      read(*,*) parma_s_W
      write(*,*)
    end if

    write(*,*) ' --- [3/7] Generation Surface ---'
    write(*,*) '   1 = Circular disk  (center + radius)'
    write(*,*) '   2 = Rectangle      (bounding box)'
    write(*,*) '   3 = Hemisphere     (radius R, XY plane only)'
    read(*,*) source_mode
    if (source_mode < 1 .or. source_mode > 3) source_mode = 1
    write(*,*)

    source_plane    = 1
    src_tilt_deg    = 0.0d0
    src_tilt_az_deg = 0.0d0

    if (source_mode == 1 .or. source_mode == 2) then

      write(*,*) ' Generation plane (orientation of the source surface):'
      write(*,*) '   1 = XY  (horizontal, z=const, muons travel in -Z)'
      write(*,*) '   2 = XZ  (vertical,   y=const, muons travel in -Y)'
      write(*,*) '   3 = YZ  (vertical,   x=const, muons travel in -X)'
      read(*,*) source_plane
      if (source_plane < 1 .or. source_plane > 3) source_plane = 1
      write(*,*)

      if (source_mode == 1) then
        ! Disk: center + radius (no bounding box)
        write(*,*) ' Disk center U [m]  (0 = centred):'
        read(*,*) src_disk_cx_m
        write(*,*) ' Disk center V [m]  (0 = centred):'
        read(*,*) src_disk_cy_m
        write(*,*) ' Disk radius [m]:'
        read(*,*) src_disk_r_m
        if (src_disk_r_m <= 0.0d0) src_disk_r_m = 1.0d0
        write(*,*) ' Fixed W [m]  (0 = surface):'
        read(*,*) src_w_m
        write(*,*) ' Tilt angle from plane normal [deg]  (0 = no tilt):'
        read(*,*) src_tilt_deg
        write(*,*) ' Tilt azimuth [deg]  (0 = toward U axis):'
        read(*,*) src_tilt_az_deg
      else
        ! Rectangle: bounding box
        if (source_plane == 1) then
          write(*,*) ' X range: Xmin [m]:'
          read(*,*) src_u1_m
          write(*,*) ' X range: Xmax [m]:'
          read(*,*) src_u2_m
          write(*,*) ' Y range: Ymin [m]:'
          read(*,*) src_v1_m
          write(*,*) ' Y range: Ymax [m]:'
          read(*,*) src_v2_m
          write(*,*) ' Fixed Z [m]  (0 = surface):'
          read(*,*) src_w_m
        else if (source_plane == 2) then
          write(*,*) ' X range: Xmin [m]:'
          read(*,*) src_u1_m
          write(*,*) ' X range: Xmax [m]:'
          read(*,*) src_u2_m
          write(*,*) ' Z range: Zmin [m]:'
          read(*,*) src_v1_m
          write(*,*) ' Z range: Zmax [m]:'
          read(*,*) src_v2_m
          write(*,*) ' Fixed Y [m]:'
          read(*,*) src_w_m
        else
          write(*,*) ' Y range: Ymin [m]:'
          read(*,*) src_u1_m
          write(*,*) ' Y range: Ymax [m]:'
          read(*,*) src_u2_m
          write(*,*) ' Z range: Zmin [m]:'
          read(*,*) src_v1_m
          write(*,*) ' Z range: Zmax [m]:'
          read(*,*) src_v2_m
          write(*,*) ' Fixed X [m]:'
          read(*,*) src_w_m
        end if
        if (src_u2_m <= src_u1_m) src_u2_m = src_u1_m + 1.0d0
        if (src_v2_m <= src_v1_m) src_v2_m = src_v1_m + 1.0d0
        write(*,*) ' Tilt angle from plane normal [deg]  (0 = no tilt):'
        read(*,*) src_tilt_deg
        write(*,*) ' Tilt azimuth [deg]  (0 = toward U axis):'
        read(*,*) src_tilt_az_deg
      end if

    else
      ! Hemisphere — radius + centre z
      write(*,*) ' Hemisphere radius [meters]:'
      read(*,*) hemi_radius_m
      write(*,*) ' Hemisphere centre Z [meters]  (0 = surface):'
      read(*,*) hemi_cz_m
      src_u1_m = 0.0d0;  src_u2_m = 0.0d0
      src_v1_m = 0.0d0;  src_v2_m = 0.0d0
      src_w_m  = hemi_cz_m
    end if
    write(*,*)

    write(*,*) ' --- [4/7] Angular Distribution ---'
    write(*,*) '   1 = Vertical only'
    write(*,*) '   2 = cos^2(theta)'
    write(*,*) '   3 = Uniform cone'
    write(*,*) '   4 = Guan/Frosin self-consistent P(theta|E)'
    write(*,*) '   5 = cos^3(theta)  Reyna-Bugaev / cosmic electrons'
    read(*,*) angular_mode
    if (angular_mode < 1 .or. angular_mode > 5) angular_mode = 2
    theta_max = 0d0
    if (angular_mode == 2 .or. angular_mode == 3 .or. angular_mode == 4 .or. angular_mode == 5) then
      write(*,*) ' Max zenith angle [degrees]  (e.g. 85):'
      read(*,*) theta_max_deg
      theta_max = theta_max_deg * PI / 180d0
    end if
    write(*,*)

    write(*,*) ' --- [5/7] Number of Muons ---'
    read(*,*) nmuons
    write(*,*)

    write(*,*) ' --- [6/7] Detector Filter ---'
    write(*,*) ' Enable detector filter? (1=Yes, 0=No)'
    read(*,*) use_detector
    write(*,*)

    ndet          = 0
    safety_margin = 0.0d0
    det_shape_arr = 1

    if (use_detector == 1) then
      write(*,'(A,I2,A)') ' How many detectors? (1 to ', MAX_DET, ')'
      read(*,*) ndet
      if (ndet < 1)       ndet = 1
      if (ndet > MAX_DET) ndet = MAX_DET
      do i_det = 1, ndet
        write(*,'(A,I2,A,I2,A)') ' --- Detector [', i_det, '/', ndet, '] ---'
        write(*,*) ' Shape:  1=Cylinder   2=Box (AABB)'
        read(*,*) det_shape_arr(i_det)
        write(*,*) ' Safety margin [cm]  (0=exact):'
        read(*,*) safety_margin
        if (det_shape_arr(i_det) == 1) then
          write(*,*) ' Bottom axis Ax Ay Az [cm]:'
          read(*,*) det_cyl_arr(i_det)%ax, det_cyl_arr(i_det)%ay, det_cyl_arr(i_det)%az
          write(*,*) ' Top axis Bx By Bz [cm]:'
          read(*,*) det_cyl_arr(i_det)%bx, det_cyl_arr(i_det)%by, det_cyl_arr(i_det)%bz
          write(*,*) ' Radius [cm]:'
          read(*,*) det_cyl_arr(i_det)%r
          det_cyl_arr(i_det)%margin = safety_margin
          det_cyl_arr(i_det)%caps   = .true.
        else
          write(*,*) ' Xmin Xmax [cm]:'
          read(*,*) det_box_arr(i_det)%xmin, det_box_arr(i_det)%xmax
          write(*,*) ' Ymin Ymax [cm]:'
          read(*,*) det_box_arr(i_det)%ymin, det_box_arr(i_det)%ymax
          write(*,*) ' Zmin Zmax [cm]:'
          read(*,*) det_box_arr(i_det)%zmin, det_box_arr(i_det)%zmax
          det_box_arr(i_det)%margin = safety_margin
        end if
      end do
    end if

    write(*,*) ' --- [7/7] Output ---'
    write(*,*) ' Save ALL muons? (1=Yes, 0=No)'
    read(*,*) save_all
    write(*,*) ' Write PHITS file? (1=Yes, 0=No)'
    read(*,*) save_phits
    write(*,*) ' Output filename for surface muons (Enter = muons_surface.dat):'
    read(*,'(A)') output_all
    if (len_trim(output_all) == 0) output_all = 'muons_surface.dat'
    if (use_detector == 1) then
      write(*,*) ' Selected muons filename (Enter = muons_selected.dat):'
      read(*,'(A)') output_sel
      if (len_trim(output_sel) == 0) output_sel = 'muons_selected.dat'
    end if
    if (save_phits == 1) then
      write(*,*) ' PHITS filename (Enter = muons_for_phits.dat):'
      read(*,'(A)') output_phits
      if (len_trim(output_phits) == 0) output_phits = 'muons_for_phits.dat'
    end if
    write(*,*)

  end if  ! use_defaults

  !===========================================================================
  ! PARMA INITIALISATION  (single-threaded, before parallel region)
  !===========================================================================
  if (spectrum_mode_in == 3) then

    call parma_set_datadir(trim(parma_datapath))

    parma_d_gcm2 = getd(parma_alt_km, parma_lat)
    parma_rc_GV  = getr(parma_lat, parma_lon)
    parma_ffp_MV = getHP(parma_year, parma_month, parma_day, parma_ic)

    parma_s_W = max(parma_s_W, -135.4d0)

    write(*,'(A,F8.2,A)') '  PARMA atm. depth:   ', parma_d_gcm2, ' g/cm2'
    write(*,'(A,F8.3,A)') '  PARMA cutoff rigid: ', parma_rc_GV,  ' GV'
    write(*,'(A,F8.1,A)') '  PARMA FFP:          ', parma_ffp_MV, ' MV'
    write(*,'(A,F8.1)')   '  PARMA W index:      ', parma_s_W
    write(*,*)

    do j = 1, NPARMA_E
      parma_E_MeV(j) = exp( log(e_min*1.0d3) + &
           dble(j-1)/dble(NPARMA_E-1) * log(e_max/e_min) )
    end do

    parma_cdf_plus(1)  = 0.0d0
    parma_cdf_minus(1) = 0.0d0
    do j = 2, NPARMA_E
      dE = parma_E_MeV(j) - parma_E_MeV(j-1)
      phi_p = getMuonSpec(1, parma_s_W, parma_rc_GV, parma_d_gcm2, parma_E_MeV(j-1))
      phi_m = getMuonSpec(2, parma_s_W, parma_rc_GV, parma_d_gcm2, parma_E_MeV(j-1))
      parma_cdf_plus(j)  = parma_cdf_plus(j-1)  + max(0.0d0, phi_p) * dE
      parma_cdf_minus(j) = parma_cdf_minus(j-1) + max(0.0d0, phi_m) * dE
    end do

    sum_p = parma_cdf_plus(NPARMA_E)
    sum_m = parma_cdf_minus(NPARMA_E)

    if (sum_p + sum_m < 1.0d-30) then
      write(*,*) ' ERROR: PARMA muon spectrum integral is zero.'
      stop
    end if

    if (sum_p > 0.0d0) parma_cdf_plus  = parma_cdf_plus  / sum_p
    if (sum_m > 0.0d0) parma_cdf_minus = parma_cdf_minus / sum_m

    if (parma_charge_mode == 1) then
      parma_ratio_plus = 1.0d0
    else if (parma_charge_mode == -1) then
      parma_ratio_plus = 0.0d0
    else
      frac_plus = sum_p / max(sum_p + sum_m, 1.0d-30)
      parma_ratio_plus = frac_plus
    end if

    do j = 1, NPARMA_ANG
      parma_cos_arr(j) = cos(theta_max) + &
           dble(j-1)/dble(NPARMA_ANG-1) * (1.0d0 - cos(theta_max))
    end do

    ! Build energy-averaged angular CDF:
    !   P(cos θ) ∝ ∫ Φ(E) · F_ang(E, cos θ) dE
    ! where Φ(E) = getMuonSpec  (mu+/mu- energy spectrum)
    !   and F_ang = getSpecAngFinal(ip=4)  (PARMA muon zenith shape factor).
    ! Integrating over the full energy range avoids the bias of using a single
    ! representative energy; mu+ and mu- get separate CDFs because their spectra
    ! have different shapes (rigidity / solar modulation corrections differ).
    parma_cdf_ang_plus(1)  = 0.0d0
    parma_cdf_ang_minus(1) = 0.0d0
    do j = 2, NPARMA_ANG
      dcos  = parma_cos_arr(j) - parma_cos_arr(j-1)
      sum_p = 0.0d0
      sum_m = 0.0d0
      do k = 1, NPARMA_E - 1
        dE      = parma_E_MeV(k+1) - parma_E_MeV(k)
        phi_p   = max(0.0d0, getMuonSpec(1, parma_s_W, parma_rc_GV, &
                             parma_d_gcm2, parma_E_MeV(k)))
        phi_m   = max(0.0d0, getMuonSpec(2, parma_s_W, parma_rc_GV, &
                             parma_d_gcm2, parma_E_MeV(k)))
        ang_fac = max(0.0d0, getSpecAngFinal(4, parma_s_W, parma_rc_GV, &
                             parma_d_gcm2, parma_E_MeV(k), 0.0d0, parma_cos_arr(j-1)))
        sum_p   = sum_p + phi_p * ang_fac * dE
        sum_m   = sum_m + phi_m * ang_fac * dE
      end do
      parma_cdf_ang_plus(j)  = parma_cdf_ang_plus(j-1)  + sum_p * dcos
      parma_cdf_ang_minus(j) = parma_cdf_ang_minus(j-1) + sum_m * dcos
    end do
    if (parma_cdf_ang_plus(NPARMA_ANG)  > 0.0d0) &
      parma_cdf_ang_plus  = parma_cdf_ang_plus  / parma_cdf_ang_plus(NPARMA_ANG)
    if (parma_cdf_ang_minus(NPARMA_ANG) > 0.0d0) &
      parma_cdf_ang_minus = parma_cdf_ang_minus / parma_cdf_ang_minus(NPARMA_ANG)

  end if  ! spectrum_mode_in == 3

  !===========================================================================
  ! DERIVED QUANTITIES
  !===========================================================================
  write_surface = 1
  if (use_detector == 1 .and. save_all == 0) write_surface = 0

  if (spectrum_mode_in == 8) then
    if (e_min <= ELECTRON_MASS) e_min = ELECTRON_MASS * 1.01d0
    p_min = sqrt(max(e_min**2 - ELECTRON_MASS**2, 0.0d0))
    p_max = sqrt(max(e_max**2 - ELECTRON_MASS**2, 0.0d0))
  else
    if (e_min <= MUON_MASS) then
      write(*,*) ' ERROR: E_min must be > 0.106 GeV. Stopping.'
      stop
    end if
    p_min = sqrt(e_min**2 - MUON_MASS**2)
    p_max = sqrt(e_max**2 - MUON_MASS**2)
  end if

  ! --- Derived quantities per source mode ---
  src_w_cm       = src_w_m * 100.0d0
  hemi_radius_cm = hemi_radius_m * 100.0d0
  hemi_cz_cm     = hemi_cz_m    * 100.0d0

  if (source_mode == 1) then
    ! Disk: use center + radius directly
    src_disk_cx_cm  = src_disk_cx_m  * 100.0d0
    src_disk_cy_cm  = src_disk_cy_m  * 100.0d0
    src_disk_r_cm   = src_disk_r_m   * 100.0d0
    src_tilt_rad    = src_tilt_deg    * PI / 180.0d0
    src_tilt_az_rad = src_tilt_az_deg * PI / 180.0d0
    centre_u_cm     = src_disk_cx_cm
    centre_v_cm     = src_disk_cy_cm
    radius_cm       = src_disk_r_cm
    half_lu_cm      = src_disk_r_cm   ! used in weight formula
    half_lv_cm      = src_disk_r_cm
  else if (source_mode == 2) then
    ! Rectangle: derive from bounding box
    half_lu_cm  = (src_u2_m - src_u1_m) * 50.0d0
    half_lv_cm  = (src_v2_m - src_v1_m) * 50.0d0
    centre_u_cm = (src_u2_m + src_u1_m) * 50.0d0
    centre_v_cm = (src_v2_m + src_v1_m) * 50.0d0
    radius_cm   = sqrt(half_lu_cm**2 + half_lv_cm**2)
    src_tilt_rad    = src_tilt_deg    * PI / 180.0d0
    src_tilt_az_rad = src_tilt_az_deg * PI / 180.0d0
  else
    ! Hemisphere
    half_lu_cm  = hemi_radius_cm
    half_lv_cm  = hemi_radius_cm
    centre_u_cm = 0.0d0
    centre_v_cm = 0.0d0
    radius_cm   = hemi_radius_cm
    src_tilt_rad    = 0.0d0
    src_tilt_az_rad = 0.0d0
  end if
  half_lx_cm  = half_lu_cm
  half_ly_cm  = half_lv_cm
  source_z_cm = 0.0d0
  if (source_mode == 3) source_z_cm = hemi_cz_cm

  depth       = p_min * 2.4d0

  !===========================================================================
  ! CONFIGURATION SUMMARY
  !===========================================================================
  write(*,*)
  write(*,*) ' ============================================================'
  write(*,*) '                GENERATION CONFIGURATION'
  write(*,*) ' ============================================================'
  write(*,'(A,F12.2,A)') '  Min Energy:     ', e_min, ' GeV'
  write(*,'(A,F12.2,A)') '  Max Energy:     ', e_max, ' GeV'
  write(*,'(A,I12)')     '  Muons to save:  ', nmuons
  write(*,'(A,I4,A)')    '  OMP threads:    ', omp_get_max_threads(), &
                         '  (set OMP_NUM_THREADS to change)'
  if (source_mode == 1) write(*,'(A)') '  Source shape:   Disc (center + radius)'
  if (source_mode == 2) write(*,'(A)') '  Source shape:   Rectangle'
  if (source_mode == 3) write(*,'(A)') '  Source shape:   Hemisphere'
  if (source_mode == 1 .or. source_mode == 2) then
    if (source_plane == 1) write(*,'(A)') '  Source plane:   XY  (muons travel in -Z)'
    if (source_plane == 2) write(*,'(A)') '  Source plane:   XZ  (muons travel in -Y)'
    if (source_plane == 3) write(*,'(A)') '  Source plane:   YZ  (muons travel in -X)'
    if (source_mode == 1) then
      write(*,'(A,F9.2,A,F9.2,A)') '  Disk center:    (', src_disk_cx_m, ',', src_disk_cy_m, ') m'
      write(*,'(A,F9.2,A)')         '  Disk radius:    ', src_disk_r_m, ' m'
      write(*,'(A,F9.2,A)')         '  W fixed:        ', src_w_m, ' m'
      if (src_tilt_deg > 0.01d0) then
        write(*,'(A,F7.2,A)')       '  Tilt angle:     ', src_tilt_deg, ' deg'
        write(*,'(A,F7.2,A)')       '  Tilt azimuth:   ', src_tilt_az_deg, ' deg'
      end if
    else
      if (source_plane == 1) then
        write(*,'(A,F9.2,A,F9.2,A)') '  X range:        ', src_u1_m, ' ..', src_u2_m, ' m'
        write(*,'(A,F9.2,A,F9.2,A)') '  Y range:        ', src_v1_m, ' ..', src_v2_m, ' m'
        write(*,'(A,F9.2,A)')         '  Z fixed:        ', src_w_m,  ' m'
      else if (source_plane == 2) then
        write(*,'(A,F9.2,A,F9.2,A)') '  X range:        ', src_u1_m, ' ..', src_u2_m, ' m'
        write(*,'(A,F9.2,A,F9.2,A)') '  Z range:        ', src_v1_m, ' ..', src_v2_m, ' m'
        write(*,'(A,F9.2,A)')         '  Y fixed:        ', src_w_m,  ' m'
      else
        write(*,'(A,F9.2,A,F9.2,A)') '  Y range:        ', src_u1_m, ' ..', src_u2_m, ' m'
        write(*,'(A,F9.2,A,F9.2,A)') '  Z range:        ', src_v1_m, ' ..', src_v2_m, ' m'
        write(*,'(A,F9.2,A)')         '  X fixed:        ', src_w_m,  ' m'
      end if
      if (src_tilt_deg > 0.01d0) then
        write(*,'(A,F7.2,A)')         '  Tilt angle:     ', src_tilt_deg,    ' deg'
        write(*,'(A,F7.2,A)')         '  Tilt azimuth:   ', src_tilt_az_deg, ' deg'
      end if
    end if
  else
    write(*,'(A,F9.2,A)') '  Hemi radius:    ', hemi_radius_m, ' m'
    write(*,'(A,F9.2,A)') '  Hemi centre Z:  ', hemi_cz_m,     ' m'
  end if
  write(*,'(A,A)')       '  Output file:    ', trim(output_all)

  ! Altitude correction — Cecchini & Spurio (2012) §3
  ! I_μ(H) = I_μ(0)·exp(+H/L(p)),  L(p) = 4900 + 750·p [m],  p in GeV/c
  ! Valid for p > 10 GeV/c, 0 < H < 1000 m above sea level.
  if (source_mode == 3) then
    alt_z_m = hemi_cz_m
  else
    alt_z_m = src_w_m
  end if
  if (spectrum_mode_in == 8) then
    alt_p_avg = p_min
  else
    alt_p_avg = sqrt(max(e_min, 0.106d0)**2 - 0.10566d0**2)
  end if
  ! p at E_min [GeV/c] — altitude correction (Cecchini & Spurio 2012) not valid for electrons
  alt_L_p   = 4900.0d0 + 750.0d0 * alt_p_avg                ! scale length [m]
  if (alt_z_m > 0.0d0 .and. alt_z_m < 1000.0d0) then
    alt_corr = exp(alt_z_m / alt_L_p)
    write(*,'(A,F9.2,A)') '  Surface alt.:   ', alt_z_m,  ' m above sea level'
    write(*,'(A,F8.4,A,F8.1,A)') &
         '  Altitude corr.: x', alt_corr, &
         '  (L(p_min)=', alt_L_p/1000.0d0, ' km)'
    write(*,'(A)') &
         '  NOTE: multiply integrated flux by this factor for eq. measurement time.'
  else if (alt_z_m < 0.0d0) then
    write(*,'(A,F9.2,A)') '  Surface depth:  ', -alt_z_m, ' m underground'
    write(*,'(A)') '  NOTE: set E_min >= threshold from Groom formula for this depth.'
  end if
  write(*,*) ' ============================================================'
  write(*,*)
  write(*,*) ' Press Enter to start...'
  read(*,'(A)') dummy

  !===========================================================================
  ! RNG INITIALISATION  (OMP change: par_init_rng instead of RLUXGO)
  !===========================================================================
  write(*,*) ' Initializing RNG streams...'
  call DATE_AND_TIME(VALUES=tim)
  iranlux = tim(6) + tim(5)*60 + tim(4)*3600
  call par_init_rng(iranlux)      ! <-- OMP: initialises one stream per thread

  if (spectrum_mode_in /= 3) then
    call build_cosmoaleph_cdf(p_min, p_max, spectrum_mode_in)
  end if

  !===========================================================================
  ! OPEN OUTPUT FILES
  !===========================================================================
  if (write_surface == 1) then
    open(unit=20, file=trim(output_all), status='unknown')
    write(20,'(A)') &
      '# EventID  x_cm  y_cm  z_cm  p_GeV  px_GeV  py_GeV  pz_GeV' // &
      '  theta_rad  phi_rad  E_GeV  charge  hit_flag  det_mask'
  end if
  if (use_detector == 1) then
    open(unit=21, file=trim(output_sel), status='unknown')
    write(21,'(A)') &
      '# EventID  x_cm  y_cm  z_cm  p_GeV  px_GeV  py_GeV  pz_GeV' // &
      '  theta_rad  phi_rad  E_GeV  charge  det_mask'
  end if
  if (save_phits == 1) open(unit=22, file=trim(output_phits), status='unknown')

  i    = 0_8
  ntry = 0_8
  write(*,*) ' Generating...'

  progress_interval = max(int(nmuons, 8) / 20_8, 1_8)   ! print ~20 progress lines

  !===========================================================================
  ! MAIN GENERATION LOOP  — OpenMP parallel
  !
  ! Pattern: !$OMP PARALLEL + infinite DO loop + EXIT condition.
  ! This guarantees exactly nmuons accepted muons regardless of acceptance rate.
  !
  ! DO NOT use !$OMP PARALLEL DO with a fixed iteration count — that pattern
  ! cannot guarantee reaching nmuons when acceptance << 1 (detector filter on).
  !
  ! Each thread loops independently until the shared counter i reaches nmuons.
  ! The exit check at the loop top uses a non-atomic read of i: at worst one
  ! thread does one extra iteration after target is reached, which is harmless
  ! (the CRITICAL block's  if (i < nmuons)  is the authoritative gate).
  !
  ! Shared state protected by:
  !   !$OMP ATOMIC CAPTURE  → ntry  (unique monotonic event ID per attempt)
  !   !$OMP CRITICAL(c_surf) → write to surface file (unit 20)
  !   !$OMP CRITICAL(c_acc)  → accepted counter i + selected/PHITS file writes
  !
  ! Output row order is non-deterministic (threads interleave).
  ! EventIDs are unique but not necessarily sequential — this is normal OMP MC.
  !===========================================================================

  !$OMP PARALLEL                                                               &
  !$OMP   DEFAULT(SHARED)                                                      &
  !$OMP   PRIVATE(x, y, z, emu, cx, cy, cz, muon_charge,                      &
  !$OMP           momentum, theta, phi, px, py, pz, ekin_MeV, kf,             &
  !$OMP           hit, hit_i, hit_flag, det_mask,                              &
  !$OMP           t_hit, t_enter, t_exit, i_det,                               &
  !$OMP           rnd, rnd2, rnd4, frac_plus, E_sampled_MeV,                   &
  !$OMP           cos_sampled, ibin, j, my_ntry, tmp_pos, tmp_dir)

  do   ! ← infinite loop: each thread keeps generating until i == nmuons

    ! Non-atomic early exit: if target already reached, stop this thread.
    ! At worst one extra muon is generated per thread; CRITICAL below discards it.
    if (i >= int(nmuons, 8)) exit

    !==========================================================================
    ! GENERATE ONE MUON
    !==========================================================================

    !--- (A)  PARMA mode: sample energy and charge from PARMA CDFs -----------
    if (spectrum_mode_in == 3) then

      call par_ranlux(rnd4);  rnd2 = dble(rnd4)
      if (rnd2 < parma_ratio_plus) then
        muon_charge = 1
      else
        muon_charge = -1
      end if

      call par_ranlux(rnd4);  rnd = dble(rnd4)
      if (muon_charge == 1) then
        ibin = 1
        do j = 2, NPARMA_E
          if (parma_cdf_plus(j) >= rnd) then
            ibin = j-1;  exit
          end if
        end do
        if (parma_cdf_plus(ibin+1) > parma_cdf_plus(ibin)) then
          frac_plus = (rnd - parma_cdf_plus(ibin)) / &
                      (parma_cdf_plus(ibin+1) - parma_cdf_plus(ibin))
        else
          frac_plus = 0.0d0
        end if
      else
        ibin = 1
        do j = 2, NPARMA_E
          if (parma_cdf_minus(j) >= rnd) then
            ibin = j-1;  exit
          end if
        end do
        if (parma_cdf_minus(ibin+1) > parma_cdf_minus(ibin)) then
          frac_plus = (rnd - parma_cdf_minus(ibin)) / &
                      (parma_cdf_minus(ibin+1) - parma_cdf_minus(ibin))
        else
          frac_plus = 0.0d0
        end if
      end if
      ibin = min(ibin, NPARMA_E-1)
      E_sampled_MeV = parma_E_MeV(ibin) + &
                      frac_plus * (parma_E_MeV(ibin+1) - parma_E_MeV(ibin))
      emu = E_sampled_MeV * 1.0d-3

      call par_ranlux(rnd4);  rnd = dble(rnd4)
      if (angular_mode == 1) then
        cos_sampled = 1.0d0
      else if (angular_mode == 3) then
        cos_sampled = cos(theta_max) + rnd * (1.0d0 - cos(theta_max))
      else
        if (muon_charge == 1) then
          ibin = 1
          do j = 2, NPARMA_ANG
            if (parma_cdf_ang_plus(j) >= rnd) then
              ibin = j-1;  exit
            end if
          end do
          if (parma_cdf_ang_plus(ibin+1) > parma_cdf_ang_plus(ibin)) then
            frac_plus = (rnd - parma_cdf_ang_plus(ibin)) / &
                        (parma_cdf_ang_plus(ibin+1) - parma_cdf_ang_plus(ibin))
          else
            frac_plus = 0.0d0
          end if
        else
          ibin = 1
          do j = 2, NPARMA_ANG
            if (parma_cdf_ang_minus(j) >= rnd) then
              ibin = j-1;  exit
            end if
          end do
          if (parma_cdf_ang_minus(ibin+1) > parma_cdf_ang_minus(ibin)) then
            frac_plus = (rnd - parma_cdf_ang_minus(ibin)) / &
                        (parma_cdf_ang_minus(ibin+1) - parma_cdf_ang_minus(ibin))
          else
            frac_plus = 0.0d0
          end if
        end if
        ibin = min(ibin, NPARMA_ANG-1)
        cos_sampled = parma_cos_arr(ibin) + &
                      frac_plus * (parma_cos_arr(ibin+1) - parma_cos_arr(ibin))
      end if

      theta = acos(max(-1.0d0, min(1.0d0, cos_sampled)))

      call par_ranlux(rnd4);  rnd = dble(rnd4)
      phi = rnd * 2.0d0 * PI

      cx =  sin(theta) * cos(phi)
      cy =  sin(theta) * sin(phi)
      cz = -cos(theta)

      call sample_position(source_mode, radius_cm, half_lx_cm, half_ly_cm, &
                           source_z_cm, x, y, z)

    !--- (B)  CosmoALEPH / Power-law / Guan / Frosin modes -------------------
    else
      call generate_muon(source_mode, radius_cm, half_lx_cm, half_ly_cm, &
                         source_z_cm, theta_max, angular_mode,            &
                         x, y, z, emu, cx, cy, cz, muon_charge)
      ! generate_muon returns flat local-frame (lu=x, lv=y, z=0) for disk/rect.
      ! Apply the same tilt transform as sample_position so that path B
      ! produces the correct tilted geometry (matches the visualiser).
      if ((source_mode == 1 .or. source_mode == 2) .and. src_tilt_rad > 1.0d-9) then
        tmp_pos = x   ! save lu
        tmp_dir = y   ! save lv
        x = tmp_pos * cos(src_tilt_rad)*cos(src_tilt_az_rad) &
          - tmp_dir  * sin(src_tilt_az_rad)
        y = tmp_pos * cos(src_tilt_rad)*sin(src_tilt_az_rad) &
          + tmp_dir  * cos(src_tilt_az_rad)
        z = -tmp_pos * sin(src_tilt_rad)
      end if
    end if

    !==========================================================================
    ! APPLY BOUNDING-BOX OFFSET + PLANE ROTATION
    ! Canonical sampling always produces a centred XY source with z_c = 0.
    ! Step 1: shift the canonical (x,y) to the bounding-box centre, set z to
    !         the fixed coordinate.
    ! Step 2: permute (x,y,z) and (cx,cy,cz) to map the canonical XY plane
    !         onto the chosen world plane.
    !
    ! Rotation table (verified with vertical muon (0,0,-1)):
    !   XY (plane=1): no change          → (0, 0,-1) ✓ downward
    !   XZ (plane=2): swap y↔z in pos,   → (0,-1, 0) ✓ travels -Y
    !                 swap cy↔cz in dir
    !   YZ (plane=3): (x,y,z)→(z,x,y),  → (-1,0, 0) ✓ travels -X
    !                 (cx,cy,cz)→(cz,cx,cy)
    !==========================================================================
    if (source_mode /= 3) then
      ! Flat source: apply centre offset; z tilt offset already set above
      x = x + centre_u_cm
      y = y + centre_v_cm
      if ((source_mode == 1 .or. source_mode == 2) .and. src_tilt_rad > 1.0d-9) then
        z = z + src_w_cm   ! add fixed W-offset to tilt offset
      else
        z = src_w_cm
      end if

      if (source_plane == 2) then
        ! XZ plane: world(x,y,z) = canonical(x, w, y)
        tmp_pos = y;  y = z;  z = tmp_pos   ! swap y ↔ z in position
        tmp_dir = cy; cy = cz; cz = tmp_dir  ! swap cy ↔ cz in direction
      else if (source_plane == 3) then
        ! YZ plane: world(x,y,z) = canonical(w, x, y)
        tmp_pos = x;  x = z;  z = y;  y = tmp_pos   ! (x,y,z)→(z,x,y) cyclic
        tmp_dir = cx; cx = cz; cz = cy; cy = tmp_dir  ! (cx,cy,cz)→(cz,cx,cy)
      end if
    else
      ! Hemisphere: source_z_cm inside generate_muon already added hemi_cz_cm
      ! (sample_position uses source_z_cm = hemi_cz_cm set above; no extra shift)
      continue
    end if

    ! Recompute theta/phi from world-frame direction after rotation
    theta = acos(max(-1.0d0, min(1.0d0, -cz)))
    if (abs(sin(theta)) < 1.0d-9) then
      phi = 0.0d0
    else
      phi = atan2(cy, cx)
      if (phi < 0.0d0) phi = phi + 2.0d0 * PI
    end if

    !==========================================================================
    ! DERIVED KINEMATICS
    !==========================================================================
    if (is_electron_mode()) then
      momentum = sqrt(max(emu**2 - ELECTRON_MASS**2, 0.0d0))
    else
      momentum = sqrt(emu**2 - MUON_MASS**2)
    end if
    px = momentum * cx
    py = momentum * cy
    pz = momentum * cz

    !==========================================================================
    ! DETECTOR INTERSECTION
    !==========================================================================
    hit      = .true.
    hit_flag = 1
    det_mask = 0

    if (use_detector == 1) then
      hit = .false.
      do i_det = 1, ndet
        hit_i = .false.
        if (det_shape_arr(i_det) == 1) then
          call ray_hits_cylinder(x, y, z, cx, cy, cz, &
                                 det_cyl_arr(i_det), hit_i, t_hit)
        else
          call ray_hits_aabb(x, y, z, cx, cy, cz, &
                             det_box_arr(i_det), hit_i, t_enter, t_exit)
        end if
        if (hit_i) then
          hit      = .true.
          det_mask = ior(det_mask, ishft(1, i_det-1))
        end if
      end do
      hit_flag = merge(1, 0, hit)
    end if

    !==========================================================================
    ! UPDATE GLOBAL TRY COUNTER — atomic capture gives unique event ID
    !==========================================================================
    !$OMP ATOMIC CAPTURE
    ntry = ntry + 1_8
    my_ntry = ntry
    !$OMP END ATOMIC

    !==========================================================================
    ! WRITE SURFACE FILE — only when detector filter is ON (logs all attempts,
    ! including misses).  When no detector filter is used, the write happens
    ! inside c_acc below so it is gated by  i < nmuons  and the file contains
    ! exactly nmuons rows (fixes the OMP over-write race).
    !==========================================================================
    if (write_surface == 1 .and. use_detector == 1) then
      !$OMP CRITICAL(c_surf)
      write(20,'(I10,1X,F13.4,1X,F13.4,1X,F13.4,1X,'// &
                'F13.6,1X,F13.6,1X,F13.6,1X,F13.6,1X,'// &
                'F13.9,1X,F13.9,1X,F13.6,1X,I4,1X,I2,1X,I8)') &
            my_ntry, x, y, z, momentum, px, py, pz, &
            theta, phi, emu, muon_charge, hit_flag, det_mask
      !$OMP END CRITICAL(c_surf)
    end if

    ! Skip the acceptance logic if this muon missed the detector
    if (use_detector == 1 .and. .not. hit) cycle

    !==========================================================================
    ! ACCEPTED MUON — increment counter, write selected/PHITS files
    !==========================================================================
    !$OMP CRITICAL(c_acc)
    if (i < int(nmuons, 8)) then
      i  = i + 1_8
      if (is_electron_mode()) then
        kf       = merge(-11, 11, muon_charge == 1)
        ekin_MeV = (emu - ELECTRON_MASS) * 1000.0d0
      else
        kf       = merge(-13, 13, muon_charge == 1)
        ekin_MeV = (emu - MUON_MASS) * 1000.0d0
      end if

      ! No detector filter: write surface file here (inside the i < nmuons gate)
      if (write_surface == 1 .and. use_detector == 0) then
        write(20,'(I10,1X,F13.4,1X,F13.4,1X,F13.4,1X,'// &
                  'F13.6,1X,F13.6,1X,F13.6,1X,F13.6,1X,'// &
                  'F13.9,1X,F13.9,1X,F13.6,1X,I4,1X,I2,1X,I8)') &
              my_ntry, x, y, z, momentum, px, py, pz, &
              theta, phi, emu, muon_charge, hit_flag, det_mask
      end if

      if (use_detector == 1) then
        write(21,'(I10,1X,F13.4,1X,F13.4,1X,F13.4,1X,'// &
                  'F13.6,1X,F13.6,1X,F13.6,1X,F13.6,1X,'// &
                  'F13.9,1X,F13.9,1X,F13.6,1X,I4,1X,I8)') &
              i, x, y, z, momentum, px, py, pz, &
              theta, phi, emu, muon_charge, det_mask
      end if

      if (save_phits == 1) then
        call write_phits_dump(22, kf, x, y, z, cx, cy, cz, &
                              ekin_MeV, 1.0d0, 0.0d0)
      end if

      if (mod(i, progress_interval) == 0_8) then
        write(*,'(A,I10,A,I16,A,F8.4,A)') &
          '  Saved', i, ' / tried', ntry, &
          '  (', 100d0*dble(i)/dble(ntry), '%)'
        flush(6)
      end if
    end if
    !$OMP END CRITICAL(c_acc)

  end do  ! infinite do — exits when i >= nmuons

  !$OMP END PARALLEL

  !===========================================================================
  ! CLOSE FILES & SUMMARY
  !===========================================================================
  if (write_surface == 1) close(20)
  if (use_detector  == 1) close(21)
  if (save_phits    == 1) close(22)

  ! OMP scheduling artefact: up to N_threads extra ntry increments occur when
  ! threads simultaneously pass the non-atomic early-exit check. With no
  ! detector filter every muon is accepted, so reset ntry = i for correct rate.
  if (use_detector == 0) ntry = i

  write(*,*)
  write(*,*) ' ============================================================'
  write(*,*) '                   COMPLETE'
  write(*,*) ' ============================================================'
  write(*,'(A,I12)') '  Saved:    ', i
  write(*,'(A,I16)') '  Tried:    ', ntry
  if (ntry > 0_8) &
    write(*,'(A,F8.4,A)') '  Rate:  ', 100d0*dble(i)/dble(ntry), ' %'
  write(*,'(A,A)')       '  File:  ', trim(output_all)
  write(*,*) ' ============================================================'

  stop

  !===========================================================================
  ! PARMA helper: sample position on generation surface
  ! Called from PARMA path inside the parallel region.
  ! Uses par_ranlux — thread-safe.
  !===========================================================================
contains

  subroutine sample_position(smode, r_cm, lx_cm, ly_cm, sz_cm, xo, yo, zo)
    integer, intent(in)  :: smode
    real(8), intent(in)  :: r_cm, lx_cm, ly_cm, sz_cm
    real(8), intent(out) :: xo, yo, zo
    real(8) :: r1, r2, r3, rr, ang, pol, lu, lv
    real(4) :: s1, s2, s3
    call par_ranlux(s1);  r1 = dble(s1)
    call par_ranlux(s2);  r2 = dble(s2)
    if (smode == 1) then
      rr  = r_cm * sqrt(r1)
      ang = r2 * 2.0d0 * PI
      lu  = rr * cos(ang)   ! local-frame U coordinate [cm]
      lv  = rr * sin(ang)   ! local-frame V coordinate [cm]
      if (src_tilt_rad < 1.0d-9) then
        ! No tilt: flat disk
        xo = lu
        yo = lv
        zo = sz_cm
      else
        ! Tilted disk using tangent vectors of the tilted plane:
        !   t1 = (cos(α)cos(φ), cos(α)sin(φ), -sin(α))  ← steepest-slope dir
        !   t2 = (-sin(φ),      cos(φ),         0      ) ← contour dir
        xo = lu * cos(src_tilt_rad)*cos(src_tilt_az_rad) &
           - lv * sin(src_tilt_az_rad)
        yo = lu * cos(src_tilt_rad)*sin(src_tilt_az_rad) &
           + lv * cos(src_tilt_az_rad)
        zo = sz_cm - lu * sin(src_tilt_rad)
      end if
    else if (smode == 2) then
      lu  = (r1 - 0.5d0) * 2.0d0 * lx_cm   ! local U in [−lx_cm, +lx_cm]
      lv  = (r2 - 0.5d0) * 2.0d0 * ly_cm   ! local V in [−ly_cm, +ly_cm]
      if (src_tilt_rad < 1.0d-9) then
        xo = lu;  yo = lv;  zo = sz_cm
      else
        xo = lu * cos(src_tilt_rad)*cos(src_tilt_az_rad) &
           - lv * sin(src_tilt_az_rad)
        yo = lu * cos(src_tilt_rad)*sin(src_tilt_az_rad) &
           + lv * cos(src_tilt_az_rad)
        zo = sz_cm - lu * sin(src_tilt_rad)
      end if
    else
      call par_ranlux(s3);  r3 = dble(s3)
      ang = r1 * 2.0d0 * PI
      pol = acos(r2)
      xo  = r_cm * sin(pol) * cos(ang)
      yo  = r_cm * sin(pol) * sin(ang)
      zo  = r_cm * cos(pol) + sz_cm
    end if
  end subroutine sample_position

end program ucmuon_gen_omp

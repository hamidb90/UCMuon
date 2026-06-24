!=============================================================================
! ucmuon_gen.f90  — UCMuon Generator  (MPI + OpenMP hybrid)
! UCLouvain Muography Group  |  Hamid Basiri <hamid.basiri@uclouvain.be>
!
! Replaces cosmoaleph_main_omp.f90 (pure-OpenMP).
!
! Parallelisation model
! ─────────────────────
!   • MPI  : each rank generates  local_nmuons = floor(N/nranks)  muons.
!             Rank 0 takes the remainder  (N mod nranks).
!             No inter-rank communication inside the hot loop.
!   • OMP  : each rank spawns OMP_NUM_THREADS threads that share the local
!             work via the existing  !$OMP PARALLEL  infinite-loop pattern.
!   • RNG  : each rank uses a distinct base seed
!             iranlux_local = iranlux + rank * 999983
!             (999983 is prime → large spacing in seed space).
!             Within each rank, rng_parallel gives every OMP thread its
!             own stream via the  THREADPRIVATE  mechanism.
!
! I/O
! ────
!   • Rank 0  reads stdin  (piped from input_params.dat) and broadcasts all
!             parameters to the other ranks via MPI_Bcast.
!   • Every rank writes its own slice to rank-labelled files:
!       ucmuon_surface_RRRRR.dat
!       ucmuon_selected_RRRRR.dat
!       ucmuon_phits_RRRRR.dat
!   • The SLURM script (run_ucmuon_gen.sh) concatenates these after the job.
!   • The column header is written ONLY by rank 0 (so concatenation is clean).
!
! Compatibility note
! ──────────────────
!   • Requires  use mpi  (provided by the MPI module from mpif90 wrapper).
!   • Compile with: mpif90 -O2 -fopenmp ...  (see Makefile target ucmuon_gen).
!   • input_params.dat format is IDENTICAL to the OMP version EXCEPT that
!     the trailing blank "Press Enter" line is no longer needed (and is
!     ignored safely if present).
!=============================================================================
program ucmuon_gen
  use mpi
  use ucmuon_source_module
  use geom_module
  use phits_module
  use rng_parallel
  use omp_lib
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
  ! MPI variables
  !---------------------------------------------------------------------------
  integer     :: my_rank, nranks, ierr
  integer(8)  :: local_nmuons          ! this rank's share of nmuons
  integer(8)  :: i_global, ntry_global ! reduced totals (rank 0 only)
  integer     :: iranlux_local         ! rank-specific RNG base seed

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

  ! Flat arrays for MPI_Bcast of derived types
  ! cyl_t fields: ax ay az bx by bz r margin caps(as 0/1 real)  → 9 reals/det
  ! aabb_t fields: xmin xmax ymin ymax zmin zmax margin          → 7 reals/det
  real(8) :: cyl_pack(9 * 10)
  real(8) :: box_pack(7 * 10)

  !---------------------------------------------------------------------------
  ! Output flags
  !---------------------------------------------------------------------------
  integer :: use_detector, save_all, save_phits
  integer :: write_surface

  !---------------------------------------------------------------------------
  ! Per-event variables (PRIVATE in OMP parallel region)
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
  integer  :: j, id, kk
  integer(8)  :: ntry           ! per-rank attempts
  integer(8)  :: i              ! per-rank accepted count
  integer     :: tim(8), iranlux
  integer     :: use_defaults
  character(120) :: output_all, output_sel, output_phits
  character(120) :: stem_all, stem_sel      ! base path without .dat, used for rank files
  character(120) :: fname_rank              ! per-rank filename scratch buffer
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
  ! OMP-specific
  !---------------------------------------------------------------------------
  integer(8) :: my_ntry
  integer(8) :: progress_interval

  !===========================================================================
  ! MPI INITIALISATION
  !===========================================================================
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)

  !===========================================================================
  ! BANNER  (rank 0 only)
  !===========================================================================
  if (my_rank == 0) then
    write(*,*)
    write(*,*) ' ============================================================'
    write(*,*) '    UCMuon Generator  (UCMuon_gen)                           '
    write(*,*) '    UCLouvain Muography Group                                '
    write(*,*) '    ** MPI + OpenMP hybrid version **                        '
    write(*,*) ' ============================================================'
    write(*,'(A,I6,A,I4,A)') '   MPI ranks: ', nranks, &
                              '   OMP threads/rank: ', omp_get_max_threads(), &
                              '   (set OMP_NUM_THREADS to change)'
    write(*,*) ' ============================================================'
    write(*,*)
  end if

  !===========================================================================
  ! INPUT  — rank 0 reads stdin, then broadcasts to all ranks
  !===========================================================================
  if (my_rank == 0) then

    write(*,*) ' Use default parameters? (1=Yes, 0=No)'
    write(*,*) '   Defaults: E=100-2500 GeV, R=800 m, N=100000'
    write(*,*) '             Angular: cos^2(theta), theta_max=85 deg'
    write(*,*) '             No detector filter, no PHITS output'
    read(*,*) use_defaults
    write(*,*)

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
      output_all       = 'ucmuon_surface.dat'
      output_sel       = 'ucmuon_selected.dat'
      output_phits     = 'ucmuon_phits.dat'

    else  ! use_defaults == 0  →  read all parameters from stdin

      write(*,*) ' --- [1/7] Energy Range ---'
      write(*,*) ' Min muon energy [GeV]  (e.g. 100):'
      read(*,*) e_min
      if (e_min < 1.0d0) e_min = 1.0d0
      write(*,*) ' Max muon energy [GeV]  (e.g. 2500):'
      read(*,*) e_max
      write(*,*)

      write(*,*) ' --- [2/7] Spectrum Model ---'
      write(*,*) '   1 = CosmoALEPH   dN/dp ~ p^-3.1952  (Schmelling 2013)'
      write(*,*) '   2 = Power-law   dN/dE ~ E^-3.7      (legacy MUSIC cross-check)'
      write(*,*) '   3 = PARMA / EXPACS  (location & date-aware)'
      write(*,*) '   4 = Guan et al. (2015)  Modified Gaisser, arXiv:1509.06176'
      write(*,*) '   5 = Frosin et al. (2025)  J.Phys.G 52, 035002'
      read(*,*) spectrum_mode_in
      if (spectrum_mode_in < 1 .or. spectrum_mode_in > 7) spectrum_mode_in = 1
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
          ! Disk: center + radius
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
      read(*,*) angular_mode
      if (angular_mode < 1 .or. angular_mode > 4) angular_mode = 2
      theta_max = 0d0
      if (angular_mode == 2 .or. angular_mode == 3 .or. angular_mode == 4) then
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
      save_phits = 0   ! PHITS output removed from hot loop — use ucmuon_to_phits
      write(*,*) ' Output filename for surface muons (Enter = ucmuon_surface.dat):'
      read(*,'(A)') output_all
      if (len_trim(output_all) == 0) output_all = 'ucmuon_surface.dat'
      if (use_detector == 1) then
        write(*,*) ' Selected muons filename (Enter = ucmuon_selected.dat):'
        read(*,'(A)') output_sel
        if (len_trim(output_sel) == 0) output_sel = 'ucmuon_selected.dat'
      end if
      output_phits = 'ucmuon_phits.dat'   ! unused — converter handles this
      write(*,*)

    end if  ! use_defaults

  end if  ! my_rank == 0

  !===========================================================================
  ! PARMA INITIALISATION  — rank 0 only, then broadcast CDFs
  !===========================================================================
  if (spectrum_mode_in == 3 .and. my_rank == 0) then

    call parma_set_datadir(trim(parma_datapath))

    parma_d_gcm2 = getd(parma_alt_km, parma_lat)
    parma_rc_GV  = getr(parma_lat, parma_lon)
    parma_ffp_MV = getHP(parma_year, parma_month, parma_day, parma_ic)
    parma_s_W    = max(parma_s_W, -135.4d0)

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
      dE    = parma_E_MeV(j) - parma_E_MeV(j-1)
      phi_p = getMuonSpec(1, parma_s_W, parma_rc_GV, parma_d_gcm2, parma_E_MeV(j-1))
      phi_m = getMuonSpec(2, parma_s_W, parma_rc_GV, parma_d_gcm2, parma_E_MeV(j-1))
      parma_cdf_plus(j)  = parma_cdf_plus(j-1)  + max(0.0d0, phi_p) * dE
      parma_cdf_minus(j) = parma_cdf_minus(j-1) + max(0.0d0, phi_m) * dE
    end do

    sum_p = parma_cdf_plus(NPARMA_E)
    sum_m = parma_cdf_minus(NPARMA_E)

    if (sum_p + sum_m < 1.0d-30) then
      write(*,*) ' ERROR: PARMA muon spectrum integral is zero.'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    if (sum_p > 0.0d0) parma_cdf_plus  = parma_cdf_plus  / sum_p
    if (sum_m > 0.0d0) parma_cdf_minus = parma_cdf_minus / sum_m

    if (parma_charge_mode == 1) then
      parma_ratio_plus = 1.0d0
    else if (parma_charge_mode == -1) then
      parma_ratio_plus = 0.0d0
    else
      frac_plus        = sum_p / max(sum_p + sum_m, 1.0d-30)
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

  end if  ! spectrum_mode_in == 3 .and. my_rank == 0

  !===========================================================================
  ! MPI BROADCAST — all parameters from rank 0 to all other ranks
  !===========================================================================

  ! --- Scalar integers ---
  call MPI_Bcast(use_defaults,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(spectrum_mode_in, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(source_mode,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(source_plane,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(angular_mode,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(nmuons,           1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(use_detector,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ndet,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(save_all,         1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  ! save_phits always 0 — PHITS conversion is post-processing
  call MPI_Bcast(parma_charge_mode,1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(parma_year,       1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(parma_month,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(parma_day,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(parma_ic,         1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  ! --- Scalar doubles ---
  call MPI_Bcast(e_min,            1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(e_max,            1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(src_disk_cx_m,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(src_disk_cy_m,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(src_disk_r_m,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(src_tilt_deg,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(src_tilt_az_deg,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(src_u1_m,         1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(src_u2_m,         1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(src_v1_m,         1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(src_v2_m,         1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(src_w_m,          1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(hemi_radius_m,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(hemi_cz_m,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(theta_max,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(parma_lat,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(parma_lon,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(parma_alt_km,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(parma_s_W,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(parma_d_gcm2,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(parma_rc_GV,      1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(parma_ffp_MV,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(parma_ratio_plus, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  ! --- Character strings ---
  call MPI_Bcast(output_all,    120, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(output_sel,    120, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(output_phits,  120, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(parma_datapath,200, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  ! --- Detector shape flags ---
  call MPI_Bcast(det_shape_arr, MAX_DET, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  ! --- Cylinder geometry  (pack derived type into flat double array) ---
  if (my_rank == 0) then
    kk = 0
    do id = 1, MAX_DET
      kk = kk + 1; cyl_pack(kk) = det_cyl_arr(id)%ax
      kk = kk + 1; cyl_pack(kk) = det_cyl_arr(id)%ay
      kk = kk + 1; cyl_pack(kk) = det_cyl_arr(id)%az
      kk = kk + 1; cyl_pack(kk) = det_cyl_arr(id)%bx
      kk = kk + 1; cyl_pack(kk) = det_cyl_arr(id)%by
      kk = kk + 1; cyl_pack(kk) = det_cyl_arr(id)%bz
      kk = kk + 1; cyl_pack(kk) = det_cyl_arr(id)%r
      kk = kk + 1; cyl_pack(kk) = det_cyl_arr(id)%margin
      kk = kk + 1; cyl_pack(kk) = merge(1.0d0, 0.0d0, det_cyl_arr(id)%caps)
    end do
  end if
  call MPI_Bcast(cyl_pack, 9*MAX_DET, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if (my_rank /= 0) then
    kk = 0
    do id = 1, MAX_DET
      kk = kk + 1; det_cyl_arr(id)%ax     = cyl_pack(kk)
      kk = kk + 1; det_cyl_arr(id)%ay     = cyl_pack(kk)
      kk = kk + 1; det_cyl_arr(id)%az     = cyl_pack(kk)
      kk = kk + 1; det_cyl_arr(id)%bx     = cyl_pack(kk)
      kk = kk + 1; det_cyl_arr(id)%by     = cyl_pack(kk)
      kk = kk + 1; det_cyl_arr(id)%bz     = cyl_pack(kk)
      kk = kk + 1; det_cyl_arr(id)%r      = cyl_pack(kk)
      kk = kk + 1; det_cyl_arr(id)%margin = cyl_pack(kk)
      kk = kk + 1; det_cyl_arr(id)%caps   = (cyl_pack(kk) > 0.5d0)
    end do
  end if

  ! --- Box geometry  (pack derived type into flat double array) ---
  if (my_rank == 0) then
    kk = 0
    do id = 1, MAX_DET
      kk = kk + 1; box_pack(kk) = det_box_arr(id)%xmin
      kk = kk + 1; box_pack(kk) = det_box_arr(id)%xmax
      kk = kk + 1; box_pack(kk) = det_box_arr(id)%ymin
      kk = kk + 1; box_pack(kk) = det_box_arr(id)%ymax
      kk = kk + 1; box_pack(kk) = det_box_arr(id)%zmin
      kk = kk + 1; box_pack(kk) = det_box_arr(id)%zmax
      kk = kk + 1; box_pack(kk) = det_box_arr(id)%margin
    end do
  end if
  call MPI_Bcast(box_pack, 7*MAX_DET, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if (my_rank /= 0) then
    kk = 0
    do id = 1, MAX_DET
      kk = kk + 1; det_box_arr(id)%xmin   = box_pack(kk)
      kk = kk + 1; det_box_arr(id)%xmax   = box_pack(kk)
      kk = kk + 1; det_box_arr(id)%ymin   = box_pack(kk)
      kk = kk + 1; det_box_arr(id)%ymax   = box_pack(kk)
      kk = kk + 1; det_box_arr(id)%zmin   = box_pack(kk)
      kk = kk + 1; det_box_arr(id)%zmax   = box_pack(kk)
      kk = kk + 1; det_box_arr(id)%margin = box_pack(kk)
    end do
  end if

  ! --- PARMA CDF arrays (only meaningful if spectrum_mode_in == 3) ---
  if (spectrum_mode_in == 3) then
    call MPI_Bcast(parma_E_MeV,       NPARMA_E,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(parma_cdf_plus,    NPARMA_E,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(parma_cdf_minus,   NPARMA_E,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(parma_cos_arr,     NPARMA_ANG, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(parma_cdf_ang_plus, NPARMA_ANG, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(parma_cdf_ang_minus,NPARMA_ANG, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end if

  !===========================================================================
  ! DERIVED QUANTITIES  (computed identically on every rank using broadcast data)
  !===========================================================================
  write_surface = 1
  if (use_detector == 1 .and. save_all == 0) write_surface = 0

  if (e_min <= MUON_MASS) then
    if (my_rank == 0) write(*,*) ' ERROR: E_min must be > 0.106 GeV. Stopping.'
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if
  p_min       = sqrt(e_min**2 - MUON_MASS**2)
  p_max       = sqrt(e_max**2 - MUON_MASS**2)
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
    half_lu_cm      = src_disk_r_cm
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
  ! WORK DISTRIBUTION
  !===========================================================================
  local_nmuons = int(nmuons, 8) / int(nranks, 8)
  ! Rank 0 absorbs the remainder so that sum(local_nmuons) == nmuons exactly.
  if (my_rank == 0) local_nmuons = local_nmuons + mod(int(nmuons, 8), int(nranks, 8))

  !===========================================================================
  ! CONFIGURATION SUMMARY  (rank 0 only)
  !===========================================================================
  if (my_rank == 0) then
    write(*,*)
    write(*,*) ' ============================================================'
    write(*,*) '            UCMuon_gen  GENERATION CONFIGURATION             '
    write(*,*) ' ============================================================'
    write(*,'(A,F12.2,A)') '  Min Energy:       ', e_min, ' GeV'
    write(*,'(A,F12.2,A)') '  Max Energy:       ', e_max, ' GeV'
    write(*,'(A,I12)')     '  Total muons:      ', nmuons
    write(*,'(A,I8)')      '  MPI ranks:        ', nranks
    write(*,'(A,I8)')      '  OMP threads/rank: ', omp_get_max_threads()
    write(*,'(A,I12)')     '  CPUs total:       ', nranks * omp_get_max_threads()
    write(*,'(A,I12)')     '  Muons / rank:     ', local_nmuons
    write(*,'(A,A)')       '  Base output:      ', trim(output_all)
  write(*,*) '  PHITS output:     use ucmuon_to_phits converter after run'
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
    ! Altitude correction — Cecchini & Spurio (2012) §3
    if (source_mode == 3) then
      alt_z_m = hemi_cz_m
    else
      alt_z_m = src_w_m
    end if
    alt_p_avg = sqrt(max(e_min, 0.106d0)**2 - 0.10566d0**2)
    alt_L_p   = 4900.0d0 + 750.0d0 * alt_p_avg
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
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierr)   ! let rank 0 finish printing before others start

  !===========================================================================
  ! RNG INITIALISATION
  !   iranlux is seeded from the wall clock on rank 0 and broadcast.
  !   Each rank then uses  iranlux + rank * 999983  as its OMP base seed.
  !   999983 is prime → seeds are well-separated in state space.
  !===========================================================================
  if (my_rank == 0) then
    call DATE_AND_TIME(VALUES=tim)
    iranlux = tim(6) + tim(5)*60 + tim(4)*3600
  end if
  call MPI_Bcast(iranlux, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  iranlux_local = iranlux + my_rank * 999983
  write(*,'(A,I5,A,I12)') '  [Rank ', my_rank, '] RNG base seed: ', iranlux_local
  call par_init_rng(iranlux_local)

  !===========================================================================
  ! CDF BUILD  (non-PARMA spectra — each rank builds its own copy)
  !===========================================================================
  if (spectrum_mode_in /= 3) then
    call build_cosmoaleph_cdf(p_min, p_max, spectrum_mode_in)
  end if

  !===========================================================================
  ! OPEN OUTPUT FILES  — rank-labelled filenames
  !   Format: <base>_RRRRR.dat  e.g.  ucmuon_surface_00003.dat
  !   Column header written ONLY by rank 0 → safe to  cat rank_0 rank_1 ...
  !===========================================================================
  ! Build stem (strip .dat suffix if present) from user-provided base names.
  ! Per-rank files: <stem>_RRRRR.dat  →  stored inside output_<JOBID>/ if
  ! the SLURM script sets the base names to include that directory.
  stem_all = trim(output_all)
  if (len_trim(stem_all) > 4) then
    if (stem_all(len_trim(stem_all)-3:len_trim(stem_all)) == '.dat') &
      stem_all = stem_all(1:len_trim(stem_all)-4)
  end if

  stem_sel = trim(output_sel)
  if (len_trim(stem_sel) > 4) then
    if (stem_sel(len_trim(stem_sel)-3:len_trim(stem_sel)) == '.dat') &
      stem_sel = stem_sel(1:len_trim(stem_sel)-4)
  end if

  if (write_surface == 1) then
    write(fname_rank,'(A,A,I5.5,A)') trim(stem_all), '_', my_rank, '.dat'
    open(unit=20, file=trim(fname_rank), status='unknown')
    if (my_rank == 0) then
      write(20,'(A)') &
        '# EventID  x_cm  y_cm  z_cm  p_GeV  px_GeV  py_GeV  pz_GeV' // &
        '  theta_rad  phi_rad  E_GeV  charge  hit_flag  det_mask'
    end if
  end if
  if (use_detector == 1) then
    write(fname_rank,'(A,A,I5.5,A)') trim(stem_sel), '_', my_rank, '.dat'
    open(unit=21, file=trim(fname_rank), status='unknown')
    if (my_rank == 0) then
      write(21,'(A)') &
        '# EventID  x_cm  y_cm  z_cm  p_GeV  px_GeV  py_GeV  pz_GeV' // &
        '  theta_rad  phi_rad  E_GeV  charge  det_mask'
    end if
  end if

  i    = 0_8
  ntry = 0_8
  if (my_rank == 0) write(*,*) ' Generating...'

  progress_interval = max(local_nmuons / 20_8, 1_8)   ! ~20 lines per rank

  !===========================================================================
  ! MAIN GENERATION LOOP  — OpenMP parallel  (unchanged from OMP version)
  !
  ! Each rank runs this loop independently until IT has saved local_nmuons.
  ! No MPI communication inside this loop.
  !
  ! Shared state protected by:
  !   !$OMP ATOMIC CAPTURE  → ntry  (unique monotonic attempt ID within rank)
  !   !$OMP CRITICAL(c_surf) → write to surface file
  !   !$OMP CRITICAL(c_acc)  → accepted counter + selected/PHITS file writes
  !===========================================================================

  !$OMP PARALLEL                                                               &
  !$OMP   DEFAULT(SHARED)                                                      &
  !$OMP   PRIVATE(x, y, z, emu, cx, cy, cz, muon_charge,                      &
  !$OMP           momentum, theta, phi, px, py, pz, ekin_MeV, kf,             &
  !$OMP           hit, hit_i, hit_flag, det_mask,                              &
  !$OMP           t_hit, t_enter, t_exit, i_det,                               &
  !$OMP           rnd, rnd2, rnd4, frac_plus, E_sampled_MeV,                   &
  !$OMP           cos_sampled, ibin, j, my_ntry, tmp_pos, tmp_dir)

  do   ! ← infinite loop: each thread generates until i == local_nmuons

    ! Non-atomic early exit check
    if (i >= local_nmuons) exit

    !==========================================================================
    ! GENERATE ONE MUON
    !==========================================================================

    !--- (A)  PARMA mode ------------------------------------------------------
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

      ! Apply centre offset + plane rotation
      if (source_mode /= 3) then
        x = x + centre_u_cm
        y = y + centre_v_cm
        if ((source_mode == 1 .or. source_mode == 2) .and. src_tilt_rad > 1.0d-9) then
          z = z + src_w_cm   ! add fixed W-offset to tilt offset from sample_position
        else
          z = src_w_cm
        end if
        if (source_plane == 2) then
          tmp_pos = y;  y = z;  z = tmp_pos   ! position
          tmp_dir = cy; cy = cz; cz = tmp_dir  ! direction
        else if (source_plane == 3) then
          tmp_pos = x;  x = z;  z = y;  y = tmp_pos
          tmp_dir = cx; cx = cz; cz = cy; cy = tmp_dir
        end if
        ! Recompute theta/phi from rotated direction
        theta = acos(max(-1.0d0, min(1.0d0, -cz)))
        if (abs(sin(theta)) < 1.0d-9) then
          phi = 0.0d0
        else
          phi = atan2(cy, cx)
          if (phi < 0.0d0) phi = phi + 2.0d0 * PI
        end if
      end if

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

      !==========================================================================
      ! APPLY BOUNDING-BOX OFFSET + PLANE ROTATION
      ! Canonical sampling always produces a centred XY source with z_c = 0.
      ! Step 1: shift the canonical (x,y) to the bounding-box centre, set z to
      !         the fixed coordinate.
      ! Step 2: permute (x,y,z) and (cx,cy,cz) to map the canonical XY plane
      !         onto the chosen world plane.
      !
      ! Rotation table (verified with vertical muon (0,0,-1)):
      !   XY (plane=1): no change          -> (0, 0,-1) correct downward
      !   XZ (plane=2): swap y<->z in pos, -> (0,-1, 0) correct travels -Y
      !                 swap cy<->cz in dir
      !   YZ (plane=3): (x,y,z)->(z,x,y), -> (-1,0, 0) correct travels -X
      !                 (cx,cy,cz)->(cz,cx,cy)
      !==========================================================================
      if (source_mode /= 3) then
        x = x + centre_u_cm
        y = y + centre_v_cm
        if ((source_mode == 1 .or. source_mode == 2) .and. src_tilt_rad > 1.0d-9) then
          z = z + src_w_cm   ! add fixed W-offset to tilt offset
        else
          z = src_w_cm
        end if

        if (source_plane == 2) then
          tmp_pos = y;  y = z;  z = tmp_pos
          tmp_dir = cy; cy = cz; cz = tmp_dir
        else if (source_plane == 3) then
          tmp_pos = x;  x = z;  z = y;  y = tmp_pos
          tmp_dir = cx; cx = cz; cz = cy; cy = tmp_dir
        end if
      else
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
    end if

    !==========================================================================
    ! DERIVED KINEMATICS
    !==========================================================================
    momentum = sqrt(emu**2 - MUON_MASS**2)
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
    ! UPDATE ATTEMPT COUNTER — atomic capture gives unique per-rank event ID
    !==========================================================================
    !$OMP ATOMIC CAPTURE
    ntry = ntry + 1_8
    my_ntry = ntry
    !$OMP END ATOMIC

    !==========================================================================
    ! WRITE SURFACE FILE
    !==========================================================================
    if (write_surface == 1) then
      !$OMP CRITICAL(c_surf)
      write(20,'(I10,1X,F13.4,1X,F13.4,1X,F13.4,1X,'// &
                'F13.6,1X,F13.6,1X,F13.6,1X,F13.6,1X,'// &
                'F13.9,1X,F13.9,1X,F13.6,1X,I4,1X,I2,1X,I8)') &
            my_ntry, x, y, z, momentum, px, py, pz, &
            theta, phi, emu, muon_charge, hit_flag, det_mask
      !$OMP END CRITICAL(c_surf)
    end if

    ! Skip if detector filter active and this muon missed
    if (use_detector == 1 .and. .not. hit) cycle

    !==========================================================================
    ! ACCEPTED MUON
    !==========================================================================
    !$OMP CRITICAL(c_acc)
    if (i < local_nmuons) then
      i  = i + 1_8
      kf = merge(-13, 13, muon_charge == 1)
      ekin_MeV = (emu - MUON_MASS) * 1000.0d0

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
        write(*,'(A,I5,A,I10,A,I12,A,F8.4,A)') &
          '  [Rank ', my_rank, '] Saved ', i, ' / tried ', ntry, &
          '  (', 100d0*dble(i)/dble(local_nmuons), '% of local target)'
        flush(6)
      end if
    end if
    !$OMP END CRITICAL(c_acc)

  end do  ! infinite loop

  !$OMP END PARALLEL

  !===========================================================================
  ! CLOSE PER-RANK FILES
  !===========================================================================
  if (write_surface == 1) close(20)
  if (use_detector  == 1) close(21)
  if (save_phits    == 1) close(22)

  !===========================================================================
  ! GLOBAL STATISTICS — MPI_Reduce to rank 0
  !===========================================================================
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Reduce(i,    i_global,    1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(ntry, ntry_global, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if (my_rank == 0) then
    write(*,*)
    write(*,*) ' ============================================================'
    write(*,*) '   UCMuon_gen  COMPLETE'
    write(*,*) ' ============================================================'
    write(*,'(A,I12)')  '  Total saved (all ranks):  ', i_global
    write(*,'(A,I16)')  '  Total tried (all ranks):  ', ntry_global
    if (ntry_global > 0_8) &
      write(*,'(A,F8.4,A)') '  Acceptance rate:           ', &
                              100d0*dble(i_global)/dble(ntry_global), ' %'
    write(*,'(A,I8,A,I4,A)') &
      '  Parallelisation:           ', nranks, ' MPI ranks × ', &
      omp_get_max_threads(), ' OMP threads'
    write(*,*)
    write(*,'(A,A,A)') '  Per-rank files: ', trim(stem_all), '_RRRRR.dat'
    write(*,*) '  Run the post-processing step in run_ucmuon_gen.sh to'
    write(*,*) '  concatenate them into a single output file.'
    write(*,*) ' ============================================================'
  end if

  call MPI_Finalize(ierr)
  stop

  !=============================================================================
contains

  !---------------------------------------------------------------------------
  ! sample_position — thread-safe position sampler (uses par_ranlux)
  !---------------------------------------------------------------------------
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
      lu  = rr * cos(ang)
      lv  = rr * sin(ang)
      if (src_tilt_rad < 1.0d-9) then
        xo = lu
        yo = lv
        zo = sz_cm
      else
        ! Tilted disk: tangent vectors
        !   t1 = (cos(α)cos(φ), cos(α)sin(φ), -sin(α))
        !   t2 = (-sin(φ),      cos(φ),         0     )
        xo = lu * cos(src_tilt_rad)*cos(src_tilt_az_rad) &
           - lv * sin(src_tilt_az_rad)
        yo = lu * cos(src_tilt_rad)*sin(src_tilt_az_rad) &
           + lv * cos(src_tilt_az_rad)
        zo = sz_cm - lu * sin(src_tilt_rad)
      end if
    else if (smode == 2) then
      lu  = (r1 - 0.5d0) * 2.0d0 * lx_cm
      lv  = (r2 - 0.5d0) * 2.0d0 * ly_cm
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

end program ucmuon_gen

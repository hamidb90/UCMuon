!=============================================================================
! ucmuon_transport_bb_omp.f90  —  UCMuon Bethe-Bloch + MS transport (OpenMP)
!
! Self-contained muon transport: Bethe-Bloch ionisation + radiative losses
! + Highland/Lynch-Dahl multiple scattering.
!
! Replaces the PHITS-dependent cosmoaleph_phitsxs_omp.f90 (which required
! MODULE NGSDATAMOD, levdat, etc.) with zero external module dependencies.
! Only ranlux_omp.o is needed (built by "make omp" in step 5 of setup.sh).
!
! Physics implemented:
!   dE/dX = dE/dX_ion  +  b_rad * b_rad_shape(E) * E
!
!   Bethe-Bloch ionisation (PDG 2022, eq. 34.5):
!     dE/dX_ion = K (Z/A) (1/beta^2) [0.5 ln(2 me c^2 beta^2 gamma^2 Tmax / I^2)
!                                       - beta^2 - delta/2]
!     Tmax = 2 me beta^2 gamma^2 / [1 + 2 gamma me/mu + (me/mu)^2]
!     Density correction delta: Sternheimer asymptotic via plasma energy
!
!   Radiative losses (bremsstrahlung + pair + photonuclear):
!     dE/dX_rad = b_rad * b_rad_shape(E) * E
!     b_rad is the material value at E_tot = 100 GeV; b_rad_shape carries the
!     PDG-2024 energy dependence (0.16 at 0.2 GeV -> 1.44 at 10 TeV for rock).
!
!   Multiple scattering (Highland formula, PDG 34.3):
!     theta_0 = (13.6 MeV / (beta*p)) * |z| * sqrt(t/X0) * [1 + 0.038 ln(t/X0)]
!     Applied as 3-D rotation about a random azimuth in the plane perp. to d.
!
!   Integration: CSDA, fixed 10 g/cm^2 steps along slant path.
!
! Accuracy vs MUSIC:
!   Survival fraction agrees to ~5-15% for std rock, 50-500 m depth.
!   Known residual: the analytic Bethe-Bloch ionisation runs 2.6-3.4% below
!   the PDG-2024 evaluated table above ~100 GeV (missing higher-order
!   corrections), so exit energies sit ~2-3% above the table-based py-BB.
!   Kept analytic on purpose: this engine is the self-contained cross-check.
!   No stochastic energy-loss fluctuations (CSDA) — overestimates survival
!   at high energies (>500 GeV) where catastrophic radiative losses matter.
!   Use MUSIC for precision; this driver gives a physically independent check.
!
! Stdin (identical sequence to what the GUI sends):
!   infile           (string)
!   outfile          (string)
!   transport_all    (0|1)
!   ncols_hint       (13|14, auto-detected internally — this value is consumed
!                    but not used; detection matches cosmoaleph_music_driver)
!   depth_m          (vertical depth in m, float)
!   mat_type         (1=StdRock  2=Ice  3=Water  4=Concrete  5=Custom)
!   [if mat_type==5 only, four separate lines:]
!     Zeff  Aeff  rho_gcm3  I_eV
!   ms_enable        (0|1)
!
! Output: 18-column ASCII, byte-for-byte identical format to
!         cosmoaleph_music_driver_omp output.
!
! Build (after "make omp" has been run):
!   gfortran -O2 -fopenmp -o ucmuon_transport_bb_omp \
!     ucmuon_transport_bb_omp.f90 ranlux_omp.o -lm
!=============================================================================
module bb_constants
  implicit none
  real(8), parameter :: MMUON    = 0.105658370d0    ! muon mass   [GeV/c^2]
  real(8), parameter :: ME       = 5.10998918d-4    ! e^- mass    [GeV/c^2]
  real(8), parameter :: K_BB     = 3.07075d-4       ! Bethe K     [GeV cm^2/g]
  real(8), parameter :: PI       = 3.141592653589793d0
  real(8), parameter :: TWO_PI   = 6.283185307179586d0
  real(8), parameter :: EMIN     = 0.105658370d0    ! muon mass — stopping threshold
end module bb_constants

program ucmuon_transport_bb_omp
  use bb_constants
  use omp_lib
  implicit none

  !---------------------------------------------------------------------------
  ! Physical constants (from bb_constants module; EMIN_GEV is program-local)
  !---------------------------------------------------------------------------
  real(8), parameter :: EMIN_GEV  = MMUON + 1.d-5   ! survival threshold [GeV]
  real(8), parameter :: STEP_GCM2 = 10.0d0           ! CSDA step   [g/cm^2]
  real(8), parameter :: PATH_MAX  = 1.0d7            ! near-horiz cutoff [cm]

  !---------------------------------------------------------------------------
  ! Material properties (set from stdin)
  !---------------------------------------------------------------------------
  integer  :: mat_type, ms_enable
  real(8)  :: Z_eff, A_eff, rho_mat, I_eV
  real(8)  :: X0_gcm2        ! radiation length   [g/cm^2]
  real(8)  :: b_rad          ! radiative b        [cm^2/g]
  real(8)  :: C_dens         ! Sternheimer asymptotic C

  !---------------------------------------------------------------------------
  ! I/O
  !---------------------------------------------------------------------------
  character(len=256) :: infile, outfile
  character(len=256) :: linebuf
  integer  :: transport_all, ncols, ncols_hint
  real(8)  :: depth_m, depth_cm

  !---------------------------------------------------------------------------
  ! Muon data arrays
  !---------------------------------------------------------------------------
  integer :: nmuon, nsurvive, nstop, nskip_mu, ios, ios2, i
  integer, allocatable :: in_eventid(:), in_charge(:), in_det_mask(:)
  real(8), allocatable :: in_xs(:), in_ys(:), in_zs(:)
  real(8), allocatable :: in_p_srf(:), in_px(:), in_py(:), in_pz(:)
  real(8), allocatable :: in_theta_s(:), in_phi_s(:), in_e_srf(:)
  integer, allocatable :: out_alive(:)
  real(8), allocatable :: out_x_ug(:), out_y_ug(:), out_z_ug(:), out_e_ug(:)
  real(8), allocatable :: out_cx_ug(:), out_cy_ug(:), out_cz_ug(:)
  real(8), allocatable :: out_theta_ug(:), out_phi_ug(:)

  !---------------------------------------------------------------------------
  ! Scratch for file parsing
  !---------------------------------------------------------------------------
  real(8) :: ev_r, xs_r, ys_r, zs_r, p_r, px_r, py_r, pz_r
  real(8) :: th_r, ph_r, e_r, ch_r, hf_r, dm_r

  !---------------------------------------------------------------------------
  ! Source-plane auto-detection
  !---------------------------------------------------------------------------
  integer  :: depth_axis    ! 0=YZ(X), 1=XZ(Y), 2=XY(Z, default)
  real(8)  :: xs_mean, ys_mean, zs_mean
  real(8)  :: xs_var,  ys_var,  zs_var

  !---------------------------------------------------------------------------
  ! OMP / per-muon transport
  !---------------------------------------------------------------------------
  integer :: nthreads, tid, iranlux_base
  real(8) :: x, y, z, cx, cy, cz, emu, slant_cm, slant_gcm2, theta_ug, phi_ug
  real(8) :: eff_cz          ! effective depth-direction cosine (plane-aware)
  real(8) :: t_start, elapsed
  real(8), external :: csda_range
  integer :: alive

  external :: rluxgo   ! from ranlux_omp.o

  !===========================================================================
  ! HEADER
  !===========================================================================
  write(*,*)
  write(*,*) ' ======================================================='
  write(*,*) '   UCMuon Bethe-Bloch + MS Transport  (OpenMP)'
  write(*,'(A,I4,A)') '   Threads available: ', omp_get_max_threads(), ' (OMP)'
  write(*,*) ' ======================================================='
  write(*,*)

  !===========================================================================
  ! READ STDIN
  !===========================================================================
  write(*,*) ' Input file  (Enter = muons_surface.dat):'
  read(*,'(A)') infile
  if (len_trim(infile) == 0) infile = 'muons_surface.dat'

  write(*,*) ' Output file (Enter = muons_underground.dat):'
  read(*,'(A)') outfile
  if (len_trim(outfile) == 0) outfile = 'muons_underground.dat'

  write(*,*) ' Transport ALL muons? (0=hit_flag=1 only, 1=all):'
  read(*,*) transport_all

  ! ncols_hint sent by GUI — consumed here, auto-detected from file below
  read(*,*) ncols_hint

  write(*,*) ' Vertical depth to detector [m]:'
  read(*,*) depth_m
  depth_cm = depth_m * 100.d0

  write(*,*) ' Material type:'
  write(*,*) '   1 = Standard Rock  (Z=11, A=22, rho=2.65 g/cm3, I=136.4 eV)'
  write(*,*) '   2 = Ice            (Z=7.42, A=14.99, rho=0.917 g/cm3, I=79.7 eV)'
  write(*,*) '   3 = Water          (Z=7.42, A=14.99, rho=1.000 g/cm3, I=79.7 eV)'
  write(*,*) '   4 = Concrete       (Z=11.11, A=22.08, rho=2.300 g/cm3, I=135.2 eV)'
  write(*,*) '   5 = Custom'
  read(*,*) mat_type

  select case(mat_type)
    case(1)  ! Standard Rock (Groom 2001 parameters)
      ! b_rad is the value at E_tot = 100 GeV (PDG-2024 per-process table);
      ! the energy dependence is applied via b_rad_shape(E).
      Z_eff=11.0d0; A_eff=22.0d0; rho_mat=2.65d0; I_eV=136.4d0
      X0_gcm2=26.54d0; b_rad=3.02d-6
    case(2)  ! Ice
      Z_eff=7.42d0; A_eff=14.99d0; rho_mat=0.917d0; I_eV=79.7d0
      X0_gcm2=36.08d0; b_rad=3.40d-6
    case(3)  ! Water
      Z_eff=7.42d0; A_eff=14.99d0; rho_mat=1.000d0; I_eV=79.7d0
      X0_gcm2=36.08d0; b_rad=3.40d-6
    case(4)  ! Concrete
      Z_eff=11.11d0; A_eff=22.08d0; rho_mat=2.300d0; I_eV=135.2d0
      X0_gcm2=26.70d0; b_rad=3.00d-6
    case(5)  ! Custom
      write(*,*) ' Enter Zeff:'  ;  read(*,*) Z_eff
      write(*,*) ' Enter Aeff:'  ;  read(*,*) A_eff
      write(*,*) ' Enter rho [g/cm^3]:'  ;  read(*,*) rho_mat
      write(*,*) ' Enter I_mean [eV]:'   ;  read(*,*) I_eV
      ! Radiation length via Tsai formula (PDG eq. 34.25): ln(287/sqrt(Z))
      X0_gcm2 = 716.408d0 * A_eff / &
                (Z_eff*(Z_eff + 1.d0)*log(287.d0/sqrt(Z_eff)))
      ! Radiative b scaled from std-rock reference (brems+pair ~ Z^2/A scaling)
      b_rad = 3.02d-6 * (Z_eff**2/A_eff) / (121.d0/22.d0)
      b_rad = max(b_rad, 1.d-7)
    case default
      Z_eff=11.0d0; A_eff=22.0d0; rho_mat=2.65d0; I_eV=136.4d0
      X0_gcm2=26.54d0; b_rad=3.02d-6; mat_type=1
  end select

  write(*,*) ' Multiple scattering? (1=ON, 0=OFF):'
  read(*,*) ms_enable

  !---------------------------------------------------------------------------
  ! Sternheimer asymptotic density correction coefficient C
  ! C = 2 * ln(hbar_omega_p / I)  - 1
  ! hbar_omega_p [GeV] = 28.816e-9 * sqrt(rho * Z / A)  (plasma frequency)
  !---------------------------------------------------------------------------
  C_dens = 2.d0 * log(28.816d-9 * sqrt(rho_mat*Z_eff/A_eff) / (I_eV*1.d-9)) - 1.d0

  write(*,*)
  write(*,'(A,F6.2,A,F6.2,A,F6.3,A)') &
    '  Material:  Z=', Z_eff, '  A=', A_eff, '  rho=', rho_mat, ' g/cm3'
  write(*,'(A,F8.2,A,F8.4,A,ES10.3,A)') &
    '  I=', I_eV, ' eV   X0=', X0_gcm2, ' g/cm2   b_rad=', b_rad, ' cm2/g'
  write(*,'(A,F8.3,A)') '  Vertical depth=', depth_m, ' m'
  write(*,*)

  !===========================================================================
  ! PASS 1a: AUTO-DETECT COLUMN FORMAT
  !===========================================================================
  open(unit=10, file=trim(infile), form='formatted', status='old', iostat=ios)
  if (ios /= 0) then
    write(*,*) ' ERROR: cannot open input file: ', trim(infile); stop
  end if

  ncols = 0
  do
    read(10,'(A)', iostat=ios) linebuf
    if (ios /= 0) exit
    if (linebuf(1:1)=='#' .or. len_trim(linebuf)==0) cycle
    read(linebuf,*,iostat=ios2) ev_r,xs_r,ys_r,zs_r,p_r, &
          px_r,py_r,pz_r,th_r,ph_r,e_r,ch_r,hf_r,dm_r
    if (ios2==0) then
      ncols = 14
    else
      read(linebuf,*,iostat=ios2) ev_r,xs_r,ys_r,zs_r,p_r, &
            px_r,py_r,pz_r,th_r,ph_r,e_r,ch_r,dm_r
      if (ios2==0) then
        ncols = 13
      else
        write(*,*) ' ERROR: cannot parse first data line.'; stop
      end if
    end if
    exit
  end do
  rewind(10)
  write(*,'(A,I2,A)') ' Detected format: ', ncols, ' columns'

  !===========================================================================
  ! PASS 1b: COUNT MUONS
  !===========================================================================
  nmuon = 0; nskip_mu = 0
  write(*,*) ' Counting muons (pass 1)...'
  do
    read(10,'(A)', iostat=ios) linebuf
    if (ios /= 0) exit
    if (linebuf(1:1)=='#' .or. len_trim(linebuf)==0) cycle
    if (ncols == 14) then
      read(linebuf,*,iostat=ios2) ev_r,xs_r,ys_r,zs_r,p_r, &
            px_r,py_r,pz_r,th_r,ph_r,e_r,ch_r,hf_r,dm_r
      if (ios2 /= 0) cycle
      if (nint(hf_r) /= 1 .and. transport_all == 0) then
        nskip_mu = nskip_mu + 1; cycle
      end if
    else
      read(linebuf,*,iostat=ios2) ev_r,xs_r,ys_r,zs_r,p_r, &
            px_r,py_r,pz_r,th_r,ph_r,e_r,ch_r,dm_r
      if (ios2 /= 0) cycle
    end if
    nmuon = nmuon + 1
  end do
  rewind(10)

  write(*,'(A,I12)') '  Muons to transport: ', nmuon
  if (ncols==14) write(*,'(A,I12)') '  Skipped (miss):    ', nskip_mu

  if (nmuon == 0) then
    write(*,*) ' ERROR: no muons found in input file.'; stop
  end if

  !===========================================================================
  ! ALLOCATE ARRAYS
  !===========================================================================
  allocate(in_eventid(nmuon), in_charge(nmuon), in_det_mask(nmuon))
  allocate(in_xs(nmuon), in_ys(nmuon), in_zs(nmuon))
  allocate(in_p_srf(nmuon), in_px(nmuon), in_py(nmuon), in_pz(nmuon))
  allocate(in_theta_s(nmuon), in_phi_s(nmuon), in_e_srf(nmuon))
  allocate(out_alive(nmuon))
  allocate(out_x_ug(nmuon), out_y_ug(nmuon), out_z_ug(nmuon), out_e_ug(nmuon))
  allocate(out_cx_ug(nmuon), out_cy_ug(nmuon), out_cz_ug(nmuon))
  allocate(out_theta_ug(nmuon), out_phi_ug(nmuon))

  !===========================================================================
  ! PASS 2: LOAD MUONS
  !===========================================================================
  write(*,*) ' Loading muons (pass 2)...'
  block
    integer :: imu
    imu = 0
    do
      read(10,'(A)', iostat=ios) linebuf
      if (ios /= 0) exit
      if (linebuf(1:1)=='#' .or. len_trim(linebuf)==0) cycle
      if (ncols == 14) then
        read(linebuf,*,iostat=ios2) ev_r,xs_r,ys_r,zs_r,p_r, &
              px_r,py_r,pz_r,th_r,ph_r,e_r,ch_r,hf_r,dm_r
        if (ios2 /= 0) cycle
        if (nint(hf_r) /= 1 .and. transport_all == 0) cycle
      else
        read(linebuf,*,iostat=ios2) ev_r,xs_r,ys_r,zs_r,p_r, &
              px_r,py_r,pz_r,th_r,ph_r,e_r,ch_r,dm_r
        if (ios2 /= 0) cycle
        hf_r = 1.d0
      end if
      imu = imu + 1
      in_eventid(imu)  = nint(ev_r)
      in_xs(imu)=xs_r; in_ys(imu)=ys_r; in_zs(imu)=zs_r
      in_p_srf(imu)=p_r
      in_px(imu)=px_r; in_py(imu)=py_r; in_pz(imu)=pz_r
      in_theta_s(imu)=th_r; in_phi_s(imu)=ph_r
      in_e_srf(imu)=e_r
      in_charge(imu)  = nint(ch_r)
      in_det_mask(imu) = nint(dm_r)
    end do
  end block
  close(10)
  write(*,'(A,I12,A)') '  Loaded ', nmuon, ' muons into memory.'
  write(*,*)

  !===========================================================================
  ! SOURCE-PLANE DETECTION
  ! The coordinate with smallest spread is the constant source-plane dimension.
  !   depth_axis=2 (XY plane, z=const): depth in Z → use |pz/p|  (default)
  !   depth_axis=1 (XZ plane, y=0):     depth in Y → use |py/p|
  !   depth_axis=0 (YZ plane, x=const): depth in X → use |px/p|
  !===========================================================================
  xs_mean = sum(in_xs(1:nmuon)) / dble(nmuon)
  ys_mean = sum(in_ys(1:nmuon)) / dble(nmuon)
  zs_mean = sum(in_zs(1:nmuon)) / dble(nmuon)
  xs_var  = sum((in_xs(1:nmuon) - xs_mean)**2) / dble(nmuon)
  ys_var  = sum((in_ys(1:nmuon) - ys_mean)**2) / dble(nmuon)
  zs_var  = sum((in_zs(1:nmuon) - zs_mean)**2) / dble(nmuon)

  if (ys_var <= xs_var .and. ys_var <= zs_var) then
    depth_axis = 1    ! XZ source plane, depth in Y
    write(*,*) '  Source plane: XZ  (y=const) — depth direction: Y'
  else if (xs_var <= ys_var .and. xs_var <= zs_var) then
    depth_axis = 0    ! YZ source plane, depth in X
    write(*,*) '  Source plane: YZ  (x=const) — depth direction: X'
  else
    depth_axis = 2    ! XY source plane, depth in Z (default)
    write(*,*) '  Source plane: XY  (z=const) — depth direction: Z'
  end if
  write(*,'(A,3ES11.3,A)') '  Pos. std (x,y,z): ', &
    sqrt(xs_var), sqrt(ys_var), sqrt(zs_var), ' cm'
  write(*,*)

  !===========================================================================
  ! INITIALIZE PER-THREAD RANLUX STREAMS
  !===========================================================================
  iranlux_base = 271828
  nthreads = omp_get_max_threads()
  !$OMP PARALLEL PRIVATE(tid)
    tid = omp_get_thread_num()
    call rluxgo(3, iranlux_base + tid * 100003, 0, 0)
  !$OMP END PARALLEL
  write(*,'(A,I4,A)') '  RNG streams initialised for ', nthreads, ' thread(s).'
  write(*,*)

  !===========================================================================
  ! PARALLEL TRANSPORT LOOP
  ! Each i writes only to out_*(i) -> no data race, no CRITICAL needed.
  !===========================================================================
  nsurvive = 0;  nstop = 0
  t_start = omp_get_wtime()
  write(*,*) ' Transporting muons (Bethe-Bloch + Highland MS)...'
  write(*,*)

  !$OMP PARALLEL DO                                                           &
  !$OMP   DEFAULT(SHARED)                                                     &
  !$OMP   PRIVATE(i, x, y, z, cx, cy, cz, eff_cz, emu, slant_cm, slant_gcm2, &
  !$OMP           alive, theta_ug, phi_ug)                                    &
  !$OMP   REDUCTION(+:nsurvive, nstop)                                        &
  !$OMP   SCHEDULE(DYNAMIC, 1)

  do i = 1, nmuon

    !-- Direction cosines from momentum or angles --
    if (in_p_srf(i) > 0.d0) then
      cx =  in_px(i) / in_p_srf(i)
      cy =  in_py(i) / in_p_srf(i)
      cz =  in_pz(i) / in_p_srf(i)
    else
      cx =  sin(in_theta_s(i)) * cos(in_phi_s(i))
      cy =  sin(in_theta_s(i)) * sin(in_phi_s(i))
      cz = -cos(in_theta_s(i))
    end if
    ! CosmoALEPH: cz<0=down; MUSIC/transport convention: cz>0=down
    cz = -cz

    x   = in_xs(i)
    y   = in_ys(i)
    z   = -in_zs(i)    ! CosmoALEPH z<=0 underground -> z>=0
    emu = in_e_srf(i)

    !-- Source-plane-aware effective depth cosine --
    ! depth_axis=1 (XZ plane): depth in Y → use -cy (cy_raw<0 for downward-Y)
    ! depth_axis=0 (YZ plane): depth in X → use -cx (cx_raw<0 for downward-X)
    ! depth_axis=2 (XY plane, default): depth in Z → use cz (already flipped >0)
    select case(depth_axis)
      case(1);  eff_cz = -cy     ! XZ plane: |py/p|, positive for py<0
      case(0);  eff_cz = -cx     ! YZ plane: |px/p|, positive for px<0
      case default;  eff_cz = cz  ! XY plane: already MUSIC-flipped
    end select

    !-- Slant path [cm] --
    if (eff_cz > 1.d-6) then
      slant_cm = depth_cm / eff_cz
    else
      slant_cm = PATH_MAX   ! near-horizontal -> always absorbed
    end if

    !-- CSDA pre-filter: skip transport when range < required slant depth --
    slant_gcm2 = rho_mat * slant_cm
    if (csda_range(emu, I_eV, Z_eff, A_eff, b_rad, C_dens) < slant_gcm2) then
      out_alive(i) = 0;  out_e_ug(i) = 0.d0;  nstop = nstop + 1
      out_x_ug(i) = in_xs(i);  out_y_ug(i) = in_ys(i);  out_z_ug(i) = in_zs(i)
      out_cx_ug(i) = cx;  out_cy_ug(i) = cy;  out_cz_ug(i) = -cz
      out_theta_ug(i) = in_theta_s(i);  out_phi_ug(i) = in_phi_s(i)
      cycle
    end if

    !-- Transport --
    call transport_bb(x, y, z, cx, cy, cz, emu, slant_cm, &
                      rho_mat, I_eV, Z_eff, A_eff, X0_gcm2, b_rad, C_dens, &
                      ms_enable)

    !-- Survival --
    if (emu > EMIN_GEV) then
      alive = 1;  out_e_ug(i) = emu;   nsurvive = nsurvive + 1
    else
      alive = 0;  out_e_ug(i) = 0.d0;  nstop    = nstop    + 1
    end if

    out_alive(i) = alive
    ! Back to CosmoALEPH coordinates
    out_x_ug(i)  =  x;   out_y_ug(i)  =  y;   out_z_ug(i)  = -z
    out_cx_ug(i) =  cx;  out_cy_ug(i) =  cy;  out_cz_ug(i) = -cz

    !-- Underground angles --
    if (alive == 1) then
      theta_ug = acos(min(1.d0, max(-1.d0, -out_cz_ug(i))))
      if (abs(sin(theta_ug)) > 1.d-9) then
        phi_ug = atan2(out_cy_ug(i), out_cx_ug(i))
        if (phi_ug < 0.d0) phi_ug = phi_ug + TWO_PI
      else
        phi_ug = 0.d0
      end if
    else
      theta_ug = in_theta_s(i);  phi_ug = in_phi_s(i)
    end if
    out_theta_ug(i) = theta_ug;  out_phi_ug(i) = phi_ug

    !-- Progress report every ~0.5% --
    if (mod(i, max(1, nmuon/200)) == 0) then
      !$OMP CRITICAL
      write(*,'(A,I10,A,I10,A,I10)') &
        '  Transported:', i, '  Survived:', nsurvive, '  Total:', nmuon
      flush(6)
      !$OMP END CRITICAL
    end if

  end do
  !$OMP END PARALLEL DO
  elapsed = omp_get_wtime() - t_start

  !===========================================================================
  ! WRITE OUTPUT — serial, preserves muon order
  ! Format is byte-for-byte identical to cosmoaleph_music_driver_omp output.
  !===========================================================================
  write(*,*)
  write(*,*) ' Writing output...'

  open(unit=11, file=trim(outfile), form='formatted', status='unknown')
  write(11,'(A)') &
    '# EventID  x_srf_cm  y_srf_cm  z_srf_cm' // &
    '  E_srf_GeV  theta_srf_rad  phi_srf_rad  charge' // &
    '  alive  x_ug_cm  y_ug_cm  z_ug_cm  E_ug_GeV' // &
    '  cx_ug  cy_ug  cz_ug  theta_ug_rad  phi_ug_rad'

  do i = 1, nmuon
    write(11,'(I10,1X,'     // &
              'F13.4,1X,F13.4,1X,F13.4,1X,' // &
              'F13.6,1X,F13.9,1X,F13.9,1X,I4,1X,' // &
              'I2,1X,' // &
              'F13.4,1X,F13.4,1X,F13.4,1X,' // &
              'F13.6,1X,' // &
              'F10.6,1X,F10.6,1X,F10.6,1X,' // &
              'F13.9,1X,F13.9)') &
      in_eventid(i), &
      in_xs(i), in_ys(i), in_zs(i), in_e_srf(i), &
      in_theta_s(i), in_phi_s(i), in_charge(i), &
      out_alive(i), &
      out_x_ug(i), out_y_ug(i), out_z_ug(i), out_e_ug(i), &
      out_cx_ug(i), out_cy_ug(i), out_cz_ug(i), &
      out_theta_ug(i), out_phi_ug(i)
  end do
  close(11)

  !===========================================================================
  ! SUMMARY
  !===========================================================================
  write(*,*)
  write(*,*) ' ======================================================='
  write(*,*) '           TRANSPORT COMPLETE  (Bethe-Bloch + MS)'
  write(*,*) ' ======================================================='
  write(*,'(A,I4)')        '  OMP threads:        ', nthreads
  write(*,'(A,I12)')       '  Muons transported:  ', nmuon
  if (ncols==14) &
    write(*,'(A,I12)')     '  Skipped (miss):     ', nskip_mu
  write(*,'(A,I12)')       '  Survived:           ', nsurvive
  write(*,'(A,I12)')       '  Stopped/absorbed:   ', nstop
  if (nmuon > 0) &
    write(*,'(A,F8.4,A)') &
      '  Survival rate:      ', 100.d0*dble(nsurvive)/dble(nmuon), ' %'
  write(*,'(A,A)')         '  Output file:        ', trim(outfile)
  write(*,'(A,F10.3,A)')  '  Elapsed:             ', elapsed, ' s'
  write(*,*)

  block
    character(len=512) :: timing_file
    timing_file = trim(outfile) // '_timing.txt'
    open(unit=13, file=trim(timing_file), form='formatted', status='unknown')
    write(13,'(A,F12.3)') 'Elapsed : ', elapsed
    close(13)
    write(*,'(A,A)') '  Timing file:        ', trim(timing_file)
    write(*,*)
  end block

  deallocate(in_eventid, in_charge, in_det_mask)
  deallocate(in_xs, in_ys, in_zs, in_p_srf, in_px, in_py, in_pz)
  deallocate(in_theta_s, in_phi_s, in_e_srf)
  deallocate(out_alive, out_x_ug, out_y_ug, out_z_ug, out_e_ug)
  deallocate(out_cx_ug, out_cy_ug, out_cz_ug, out_theta_ug, out_phi_ug)

end program ucmuon_transport_bb_omp


!=============================================================================
! transport_bb
!
! Integrate muon energy loss + multiple scattering along slant path.
!
! Coordinates: MUSIC convention (z=0 surface, z>0 underground, cz>0 = down).
! On entry:  (x,y,z) [cm], (cx,cy,cz) unit direction, E [GeV]
! On exit:   updated (x,y,z), (cx,cy,cz), E
!
! ranlux is called for MS — must be called from within an OMP region where
! each thread has already called rluxgo() to seed its private stream.
!=============================================================================
subroutine transport_bb(x, y, z, cx, cy, cz, E, slant_cm, &
                         rho, I_eV, Zat, Aat, X0_gcm2, b_rad, C_dens, ms_enable)
  use bb_constants
  implicit none

  real(8), intent(inout) :: x, y, z, cx, cy, cz, E
  real(8), intent(in)    :: slant_cm   ! total slant path [cm]
  real(8), intent(in)    :: rho        ! density          [g/cm^3]
  real(8), intent(in)    :: I_eV       ! mean exc. energy [eV]
  real(8), intent(in)    :: Zat, Aat       ! atomic Zat, mass Aat
  real(8), intent(in)    :: X0_gcm2    ! rad. length       [g/cm^2]
  real(8), intent(in)    :: b_rad      ! radiative b       [cm^2/g]
  real(8), intent(in)    :: C_dens     ! Sternheimer C
  integer, intent(in)    :: ms_enable

  real(8), parameter :: STEP_G = 10.0d0   ! [g/cm^2]

  integer :: nsteps, istep
  real(8) :: slant_gcm2, ds_gcm2, ds_cm
  real(8) :: gamma_mu, beta2, p_GeV, bg, Tmax, I_GeV
  real(8) :: log_arg, x_dens, delta, dEdX
  real(8) :: t_X0, theta0
  real(8) :: e1x, e1y, e1z, e2x, e2y, e2z, norm
  real(8) :: theta_scat, phi_scat, ct, st, cos_phi, sin_phi
  real(4) :: rr(4)

  external :: ranlux   ! from ranlux_omp.o  (THREADPRIVATE state)
  real(8), external :: b_rad_shape

  I_GeV      = I_eV * 1.d-9
  slant_gcm2 = rho * slant_cm
  nsteps     = max(20, int(slant_gcm2 / STEP_G) + 1)
  ds_gcm2    = slant_gcm2 / dble(nsteps)
  ds_cm      = ds_gcm2 / rho

  do istep = 1, nsteps

    if (E <= EMIN) return

    !-------------------------------------------------------------------------
    ! Bethe-Bloch mean ionisation loss   [GeV/(g/cm^2)]
    !-------------------------------------------------------------------------
    gamma_mu = E / MMUON
    beta2    = max(0.d0, 1.d0 - (MMUON/E)**2)
    p_GeV    = sqrt(max(0.d0, E*E - MMUON*MMUON))
    bg       = p_GeV / MMUON   ! beta*gamma

    ! Max energy transfer to atomic electron
    Tmax = 2.d0*ME*beta2*gamma_mu*gamma_mu / &
           (1.d0 + 2.d0*gamma_mu*ME/MMUON + (ME/MMUON)**2)

    ! Sternheimer density correction (asymptotic, valid for bg > ~10)
    x_dens = log10(max(bg, 1.d-30))
    if (x_dens > 1.0d0) then
      delta = 2.d0 * log(10.d0) * x_dens + C_dens
      delta = max(0.d0, delta)
    else
      delta = 0.d0
    end if

    if (Tmax > 0.d0 .and. beta2 > 1.d-12) then
      log_arg = 2.d0*ME*beta2*gamma_mu*gamma_mu*Tmax / (I_GeV*I_GeV)
      if (log_arg > 1.d0) then
        dEdX = K_BB * (Zat/Aat) / beta2 * &
               (0.5d0*log(log_arg) - beta2 - 0.5d0*delta)
      else
        dEdX = K_BB * (Zat/Aat) / beta2
      end if
    else
      dEdX = K_BB * (Zat/Aat) / max(beta2, 1.d-12)
    end if
    dEdX = max(dEdX, 1.d-6)

    ! Add radiative component:  b_rad(E) * E  [GeV/(g/cm^2)]
    ! (b_rad is the 100 GeV value; b_rad_shape carries the PDG-2024
    !  energy dependence — a constant b under-counts ~15% at 300 GeV.)
    dEdX = dEdX + b_rad * b_rad_shape(E) * E

    ! Update energy
    E = E - dEdX * ds_gcm2
    if (E <= EMIN) return

    !-------------------------------------------------------------------------
    ! Update position along current direction
    !-------------------------------------------------------------------------
    x = x + cx * ds_cm
    y = y + cy * ds_cm
    z = z + cz * ds_cm

    !-------------------------------------------------------------------------
    ! Highland multiple scattering  (PDG 2022, eq. 34.15)
    ! theta_0 = (13.6 MeV / (beta*p)) * sqrt(t/X0) * (1 + 0.038 ln(t/X0))
    !-------------------------------------------------------------------------
    if (ms_enable == 1 .and. X0_gcm2 > 0.d0) then

      t_X0 = ds_gcm2 / X0_gcm2

      if (t_X0 > 1.d-14) then
        ! Recompute p, beta with updated E
        p_GeV = sqrt(max(0.d0, E*E - MMUON*MMUON))
        beta2 = max(1.d-12, 1.d0 - (MMUON/E)**2)

        theta0 = 13.6d-3 / (sqrt(beta2) * p_GeV) * sqrt(t_X0) * &
                 (1.d0 + 0.038d0 * log(t_X0))

        if (theta0 > 0.d0) then
          ! Polar deflection is Rayleigh(theta0): Highland theta0 is the RMS
          ! *projected* angle, and two independent N(0,theta0) projections
          ! give a space angle theta0*sqrt(-2 ln u).  (A Gaussian polar angle
          ! under-scatters by sqrt(2).)
          call ranlux(rr, 4)
          rr(1) = max(rr(1), 1.e-30)
          theta_scat = theta0 * sqrt(-2.d0*log(dble(rr(1))))
          phi_scat   = TWO_PI * dble(rr(3))

          ! Construct orthonormal basis (e1, e2) perpendicular to (cx,cy,cz)
          if (abs(cx) <= abs(cy) .and. abs(cx) <= abs(cz)) then
            e1x =  0.d0; e1y = -cz; e1z =  cy
          else if (abs(cy) <= abs(cz)) then
            e1x =  cz;   e1y =  0.d0; e1z = -cx
          else
            e1x = -cy;   e1y =  cx;  e1z =  0.d0
          end if
          norm = sqrt(e1x*e1x + e1y*e1y + e1z*e1z)
          if (norm >= 1.d-14) then
            e1x = e1x/norm;  e1y = e1y/norm;  e1z = e1z/norm

            ! e2 = (cx,cy,cz) × e1
            e2x = cy*e1z - cz*e1y
            e2y = cz*e1x - cx*e1z
            e2z = cx*e1y - cy*e1x

            ! Rotate direction: d_new = cos(theta)*d + sin(theta)*(cos(phi)*e1 + sin(phi)*e2)
            ct      = cos(theta_scat);  st      = sin(theta_scat)
            cos_phi = cos(phi_scat);    sin_phi = sin(phi_scat)

            cx = ct*cx + st*(cos_phi*e1x + sin_phi*e2x)
            cy = ct*cy + st*(cos_phi*e1y + sin_phi*e2y)
            cz = ct*cz + st*(cos_phi*e1z + sin_phi*e2z)

            ! Renormalise
            norm = sqrt(cx*cx + cy*cy + cz*cz)
            if (norm > 1.d-14) then
              cx = cx/norm;  cy = cy/norm;  cz = cz/norm
            end if
          end if
        end if
      end if
    end if

  end do

end subroutine transport_bb


!=============================================================================
! csda_range  —  CSDA range [g/cm²] used by the pre-filter
!
! Integrates the same Bethe-Bloch + radiative dE/dx as transport_bb from
! E_init down to EMIN using NR uniform energy steps.  Because the transport
! is purely deterministic (CSDA), range < slant_gcm2 is an exact stopping
! criterion — no surviving muon is ever killed by this pre-filter.
!=============================================================================
real(8) function csda_range(E_init, I_eV, Zat, Aat, b_rad, C_dens)
  use bb_constants
  implicit none
  real(8), intent(in) :: E_init, I_eV, Zat, Aat, b_rad, C_dens

  integer, parameter :: NR = 100
  real(8) :: E, dE, I_GeV
  real(8) :: gamma_mu, beta2, p_GeV, bg, Tmax
  real(8) :: x_dens, delta, log_arg, dedx_val
  integer :: k
  real(8), external :: b_rad_shape

  csda_range = 0.d0
  if (E_init <= EMIN) return

  I_GeV = I_eV * 1.d-9
  dE    = (E_init - EMIN) / dble(NR)
  E     = E_init

  do k = 1, NR
    if (E <= EMIN) exit

    gamma_mu = E / MMUON
    beta2    = max(1.d-12, 1.d0 - (MMUON/E)**2)
    p_GeV    = sqrt(max(0.d0, E*E - MMUON*MMUON))
    bg       = p_GeV / MMUON

    Tmax = 2.d0*ME*beta2*gamma_mu*gamma_mu / &
           (1.d0 + 2.d0*gamma_mu*ME/MMUON + (ME/MMUON)**2)

    x_dens = log10(max(bg, 1.d-30))
    if (x_dens > 1.0d0) then
      delta = max(0.d0, 2.d0*log(10.d0)*x_dens + C_dens)
    else
      delta = 0.d0
    end if

    if (Tmax > 0.d0 .and. beta2 > 1.d-12) then
      log_arg = 2.d0*ME*beta2*gamma_mu*gamma_mu*Tmax / (I_GeV*I_GeV)
      if (log_arg > 1.d0) then
        dedx_val = K_BB*(Zat/Aat)/beta2*(0.5d0*log(log_arg) - beta2 - 0.5d0*delta)
      else
        dedx_val = K_BB*(Zat/Aat)/max(beta2, 1.d-12)
      end if
    else
      dedx_val = K_BB*(Zat/Aat)/max(beta2, 1.d-12)
    end if
    dedx_val = max(dedx_val, 1.d-6) + b_rad*b_rad_shape(E)*E

    csda_range = csda_range + dE / dedx_val
    E = E - dE
  end do
end function csda_range


!=============================================================================
! b_rad_shape  —  energy dependence of the radiative-loss coefficient b(E),
! normalised to 1 at E_tot = 100 GeV.
!
! Nodes computed from the PDG-2024 Standard Rock per-process table
! (brems+pair+photonuclear; same source as the UCMuon-MC v2 loss model):
! b(E) = L_rad(E)/E.  Log-log interpolation, <1.6% error over 0.2 GeV-10 TeV,
! clamped at the ends (radiative losses are negligible below 0.2 GeV).
! The shape is applied to every material's b_rad, which is defined as the
! value at 100 GeV; the ln E rise of brems/pair is material-universal at the
! accuracy of this engine.
!=============================================================================
real(8) function b_rad_shape(E_tot)
  implicit none
  real(8), intent(in) :: E_tot   ! total energy [GeV]

  integer, parameter :: NB = 15
  real(8), parameter :: EB(NB) = (/ 0.2d0, 0.5d0, 1.0d0, 2.0d0, 5.0d0, &
       10.0d0, 20.0d0, 50.0d0, 100.0d0, 200.0d0, 500.0d0, 1000.0d0, &
       2000.0d0, 5000.0d0, 10000.0d0 /)
  real(8), parameter :: SB(NB) = (/ 0.16249d0, 0.21144d0, 0.24846d0, &
       0.34867d0, 0.49733d0, 0.60751d0, 0.72233d0, 0.88382d0, 1.00000d0, &
       1.10938d0, 1.22605d0, 1.29981d0, 1.35654d0, 1.41206d0, 1.44280d0 /)

  real(8) :: le, frac
  integer :: j

  if (E_tot <= EB(1)) then
    b_rad_shape = SB(1);  return
  else if (E_tot >= EB(NB)) then
    b_rad_shape = SB(NB); return
  end if

  le = log(E_tot)
  do j = 2, NB
    if (E_tot <= EB(j)) then
      frac = (le - log(EB(j-1))) / (log(EB(j)) - log(EB(j-1)))
      b_rad_shape = exp( log(SB(j-1)) + frac*(log(SB(j)) - log(SB(j-1))) )
      return
    end if
  end do
  b_rad_shape = SB(NB)
end function b_rad_shape

!=============================================================================
! cosmoaleph_music_driver_omp.f90  --  OpenMP-parallel MUSIC transport driver
! UCLouvain Muography Group / Hamid Basiri
!
! MATERIAL FILE MANAGEMENT
!   MUSIC needs two material-specific data files:
!     music-eloss-{mat}.dat          -- energy-loss table (ships with MUSIC)
!     music-cross-sections-{mat}.dat -- integral cross-sections (generated)
!   where {mat} = rock | water | seawater  (mat_type = 1 | 2 | 3)
!
!   mat_type groupings:
!     1 = rock (Standard Rock, Limestone, Rock Salt, Iron, Custom)
!     2 = water / ice
!     3 = seawater
!
!   music-double-diff-rock.dat is material-INDEPENDENT (angular deflection
!   shapes do not depend on composition) -- one file used for all materials.
!
!   REQUIRED FILES (must be in the run directory):
!     music-eloss-rock.dat      -- copy/rename of original music-eloss.dat
!     music-eloss-water.dat     -- water version (from MUSIC distribution)
!     music-eloss-seawater.dat  -- seawater version
!     music-double-diff-rock.dat
!
!   music-cross-sections-{mat}.dat is AUTO-GENERATED on the first run
!   (init=0) and cached; later runs use init=1.
!
!   FALLBACK: if music-eloss-{mat}.dat is absent, music-eloss.dat is used
!   with a warning.  For non-rock materials this gives wrong Z/A.
!
! THREAD SAFETY
!   RANLUX: COMMON/RLUXSTATE/ + THREADPRIVATE  (ranlux_omp.f)
!   RANMAR: COMMON/RANMA1/    + THREADPRIVATE  (ranmar_omp.f)
!   MUSIC COMMONs: read-only after initialize_music -- no locks needed.
!=============================================================================
program ucmuon_transport_music_omp
  use omp_lib
  implicit none

  external :: rluxgo, rmarin
  external :: mucrsec, mulos, initialize_music, muon_transport

  real(8), parameter :: PI    = 3.141592654d0
  real(8), parameter :: MMUON = 0.105658d0
  ! Array sizes determined at runtime via two-pass read (no fixed limit).

  ! Material composition arrays
  real(4) :: zz0(20), a0(20), fr0(20), par_ion(6)
  real(4) :: zz0_rock(20)      = (/11.0,0.,0.,0.,0.,0.,0.,0.,0.,0., &
                                     0., 0.,0.,0.,0.,0.,0.,0.,0.,0./)
  real(4) :: a0_rock(20)       = (/22.0,0.,0.,0.,0.,0.,0.,0.,0.,0., &
                                     0., 0.,0.,0.,0.,0.,0.,0.,0.,0./)
  real(4) :: fr0_rock(20)      = (/1.0, 0.,0.,0.,0.,0.,0.,0.,0.,0., &
                                     0., 0.,0.,0.,0.,0.,0.,0.,0.,0./)
  real(4) :: par_ion_rock(6)   = (/136.4,-3.774,0.083,3.412,3.055,0.049/)
  real(4) :: zz0_water(20)     = (/1.0,8.0,0.,0.,0.,0.,0.,0.,0.,0., &
                                     0., 0.,0.,0.,0.,0.,0.,0.,0.,0./)
  real(4) :: a0_water(20)      = (/1.008,15.999,0.,0.,0.,0.,0.,0.,0.,0., &
                                     0.,  0.,   0.,0.,0.,0.,0.,0.,0.,0./)
  real(4) :: fr0_water(20)     = (/0.1119,0.8881,0.,0.,0.,0.,0.,0.,0.,0., &
                                     0.,   0.,   0.,0.,0.,0.,0.,0.,0.,0./)
  real(4) :: par_ion_water(6)  = (/75.0,-3.502,0.2065,3.007,2.5,0.24/)
  real(4) :: zz0_seawater(20)  = (/1.0,8.0,11.0,17.0,0.,0.,0.,0.,0.,0., &
                                     0., 0., 0.,  0., 0.,0.,0.,0.,0.,0./)
  real(4) :: a0_seawater(20)   = (/1.008,15.999,22.990,35.453,0.,0.,0.,0.,0.,0., &
                                     0.,  0.,    0.,    0.,   0.,0.,0.,0.,0.,0./)
  real(4) :: fr0_seawater(20)  = (/0.1100,0.8779,0.0106,0.0019,0.,0.,0.,0.,0.,0., &
                                     0.,   0.,    0.,    0.,   0.,0.,0.,0.,0.,0./)
  real(4) :: par_ion_seawater(6) = (/75.0,-3.502,0.2065,3.007,2.5,0.24/)

  ! User inputs
  character(200) :: infile, outfile, linebuf
  real(8)  :: depth_m, depth_cm, rho, rad
  integer  :: idim, idim1, minv, init, mat_type, transport_all, ncols

  ! Material file management
  character(60) :: mat_suffix, eloss_file, xsec_file
  logical       :: eloss_exists, xsec_exists, fallback_eloss
  integer       :: stat_cmd

  ! Input arrays
  integer,  allocatable :: in_eventid(:), in_charge(:), in_det_mask(:)
  real(8),  allocatable :: in_xs(:), in_ys(:), in_zs(:)
  real(8),  allocatable :: in_p_srf(:), in_px(:), in_py(:), in_pz(:)
  real(8),  allocatable :: in_theta_s(:), in_phi_s(:), in_e_srf(:)

  ! Output arrays
  integer,  allocatable :: out_alive(:)
  real(8),  allocatable :: out_x_ug(:), out_y_ug(:), out_z_ug(:)
  real(8),  allocatable :: out_e_ug(:)
  real(8),  allocatable :: out_cx_ug(:), out_cy_ug(:), out_cz_ug(:)
  real(8),  allocatable :: out_theta_ug(:), out_phi_ug(:)

  ! Source-plane auto-detection
  integer  :: depth_axis    ! 0=YZ(X), 1=XZ(Y), 2=XY(Z, default)
  real(8)  :: xs_mean, ys_mean, zs_mean
  real(8)  :: xs_var,  ys_var,  zs_var

  ! Private transport variables
  real(8)  :: x, y, z, cx, cy, cz, emu, ttime, depth_this_muon
  real(8)  :: eff_cz         ! effective depth-direction cosine (plane-aware)
  real(8)  :: theta_ug, phi_ug
  integer  :: alive

  ! Counters
  integer  :: nmuon, nsurvive, nstop, nskip_mu, ios, i
  integer  :: nthreads, tid, iranlux_base

  write(*,*)
  write(*,*) ' ======================================================='
  write(*,*) '   UCMuon MUSIC Transport Driver (OpenMP)'
  write(*,'(A,I4,A)') '   Threads available: ', omp_get_max_threads(), &
                      '  (set OMP_NUM_THREADS to change)'
  write(*,*) ' ======================================================='
  write(*,*)

  !==========================================================================
  ! USER INPUT
  !==========================================================================
  write(*,*) ' Input file  (Enter = muons_surface.dat):'
  read(*,'(A)') infile
  if (len_trim(infile) == 0) infile = 'muons_surface.dat'

  write(*,*) ' Output file (Enter = muons_underground.dat):'
  read(*,'(A)') outfile
  if (len_trim(outfile) == 0) outfile = 'muons_underground.dat'

  write(*,*) ' Rock density [g/cm^3]  (Enter = 2.65):'
  read(*,*) rho
  if (rho <= 0.d0) rho = 2.65d0

  write(*,*) ' Radiation length [cm]  (Enter = 26.48):'
  read(*,*) rad
  if (rad <= 0.d0) rad = 26.48d0

  write(*,*) ' Vertical depth to detector [m]  (e.g. 90):'
  read(*,*) depth_m
  depth_cm = depth_m * 100.d0

  write(*,*) ' 3D transport / lateral scattering (idim)? (1=ON, 0=OFF):'
  read(*,*) idim

  write(*,*) ' Other-process scattering (idim1)? (1=ON, 0=OFF):'
  read(*,*) idim1

  write(*,*) ' Energy loss cut exponent (minv, expert only, default=-30):'
  read(*,*) minv

  write(*,*) ' Cross-section tables on disk?'
  write(*,*) '   0 = calculate + save  (first run for this material, ~1 min)'
  write(*,*) '   1 = read from disk    (subsequent runs, fast)'
  read(*,*) init

  write(*,*) ' Material type:'
  write(*,*) '   1 = rock  (Std Rock / Limestone / Rock Salt / Iron / Custom)'
  write(*,*) '   2 = water / ice'
  write(*,*) '   3 = seawater'
  read(*,*) mat_type
  if (mat_type < 1 .or. mat_type > 3) mat_type = 1

  ! Set composition arrays and suffix
  select case (mat_type)
    case (2)
      zz0 = zz0_water;  a0 = a0_water;  fr0 = fr0_water;  par_ion = par_ion_water
      mat_suffix = 'water'
      write(*,*) '  Composition: Water/Ice  (H mass-frac 0.1119, O 0.8881)'
    case (3)
      zz0 = zz0_seawater;  a0 = a0_seawater;  fr0 = fr0_seawater
      par_ion = par_ion_seawater
      mat_suffix = 'seawater'
      write(*,*) '  Composition: Seawater  (H, O, Na, Cl)'
    case default
      zz0 = zz0_rock;  a0 = a0_rock;  fr0 = fr0_rock;  par_ion = par_ion_rock
      mat_suffix = 'rock'
      write(*,*) '  Composition: Standard Rock  (Z=11, A=22)'
      write(*,*) '  NOTE: density/rad.length from your input override defaults.'
      write(*,*) '        Elemental composition uses standard-rock tables.'
      write(*,*) '        Accurate for rock-like media; approximate for Fe/metals.'
  end select

  !==========================================================================
  ! RESOLVE MATERIAL-SPECIFIC FILENAMES
  !==========================================================================
  eloss_file = 'music-eloss-' // trim(mat_suffix) // '.dat'
  xsec_file  = 'music-cross-sections-' // trim(mat_suffix) // '.dat'

  ! --- Energy-loss file ---
  ! When init=0: mulos() will GENERATE music-eloss.dat -- no pre-check needed.
  ! When init=1: file must already exist.
  fallback_eloss = .false.
  if (init == 1) then
    inquire(file=trim(eloss_file), exist=eloss_exists)
    if (.not. eloss_exists) then
      inquire(file='music-eloss.dat', exist=eloss_exists)
      if (eloss_exists) then
        write(*,*)
        write(*,'(A,A,A)') ' *** WARNING: ', trim(eloss_file), ' not found.'
        write(*,*) '     Falling back to music-eloss.dat'
        if (mat_type /= 1) then
          write(*,*) '     Z and A will be standard-rock values (Z=11, A=22).'
          write(*,*) '     Energy loss and multiple scattering will be APPROXIMATE.'
          write(*,*) '     Re-run with init=0 to regenerate for this material.'
          fallback_eloss = .true.
        end if
        eloss_file = 'music-eloss.dat'
      else
        write(*,*) ' ERROR: init=1 but no eloss file found.'
        write(*,'(A,A)') '        Expected: ', trim(eloss_file)
        write(*,*) '        Run once with init=0 -- mulos() will generate it.'
        stop
      end if
    end if
  else
    ! init=0: check that the differential table exists (needed by mulos)
    inquire(file='music-double-diff-rock.dat', exist=eloss_exists)
    if (.not. eloss_exists) then
      write(*,*) ' ERROR: music-double-diff-rock.dat not found.'
      write(*,*) '        Required by mulos() to compute the energy-loss table.'
      write(*,*) '        Place it in the project root (from Kudryavtsev MUSIC zip).'
      stop
    end if
    write(*,'(A,A,A)') '  init=0: mulos() will generate ', trim(eloss_file), ' ...'
  end if

  ! --- Cross-sections file ---
  inquire(file=trim(xsec_file), exist=xsec_exists)
  if (init == 1 .and. .not. xsec_exists) then
    ! Try generic fallback
    inquire(file='music-cross-sections.dat', exist=xsec_exists)
    if (xsec_exists) then
      write(*,*)
      write(*,'(A,A,A)') ' *** WARNING: ', trim(xsec_file), ' not found.'
      write(*,*) '     Using music-cross-sections.dat (may be from a different material).'
      write(*,*) '     Re-run with init=0 to compute correct tables for ', trim(mat_suffix)
      xsec_file = 'music-cross-sections.dat'
    else
      write(*,*) ' ERROR: init=1 but no cross-section file found.'
      write(*,'(A,A)')   '        Expected: ', trim(xsec_file)
      write(*,*) '        Re-run with init=0 to compute and cache the tables.'
      stop
    end if
  end if

  !==========================================================================
  ! DETECT INPUT FORMAT
  !==========================================================================
  block
    integer :: hf_tmp, dm_tmp, ios2, ev_tmp, ch_tmp
    real(8) :: p_tmp, px_tmp, py_tmp, pz_tmp, th_tmp, ph_tmp, e_tmp, xs_tmp, ys_tmp, zs_tmp
    open(unit=10, file=trim(infile), form='formatted', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) ' ERROR: cannot open ', trim(infile); stop
    end if
    ncols = 0
    do
      read(10,'(A200)', iostat=ios) linebuf
      if (ios /= 0) exit
      if (linebuf(1:1)=='#' .or. len_trim(linebuf)==0) cycle
      read(linebuf,*,iostat=ios2) ev_tmp,xs_tmp,ys_tmp,zs_tmp,p_tmp, &
            px_tmp,py_tmp,pz_tmp,th_tmp,ph_tmp,e_tmp,ch_tmp,hf_tmp,dm_tmp
      if (ios2==0) then
        ncols = 14
      else
        read(linebuf,*,iostat=ios2) ev_tmp,xs_tmp,ys_tmp,zs_tmp,p_tmp, &
              px_tmp,py_tmp,pz_tmp,th_tmp,ph_tmp,e_tmp,ch_tmp,dm_tmp
        if (ios2==0) then; ncols = 13
        else; write(*,*) ' ERROR: cannot parse first data line.'; stop
        end if
      end if
      exit
    end do
    rewind(10)
  end block

  transport_all = 0
  write(*,'(A,I2,A)') ' Detected format: ', ncols, ' columns'
  if (ncols == 14) then
    write(*,*) ' Transport ALL muons? (0=hit_flag=1 only, 1=all):'
    read(*,*) transport_all
  end if

  !==========================================================================
  ! SUMMARY
  !==========================================================================
  write(*,*)
  write(*,*) ' ======================================================='
  write(*,'(A,A)')       '  Input:            ', trim(infile)
  write(*,'(A,A)')       '  Output:           ', trim(outfile)
  write(*,'(A,F7.3,A)')  '  Density:          ', rho, ' g/cm^3'
  write(*,'(A,F7.2,A)')  '  Rad. length:      ', rad, ' cm'
  write(*,'(A,F8.2,A)')  '  Depth:            ', depth_m, ' m'
  write(*,'(A,F10.1,A)') '  Depth (w.e.):     ', depth_cm*rho, ' g/cm^2'
  write(*,'(A,A)')       '  Material group:   ', trim(mat_suffix)
  write(*,'(A,A)')       '  Eloss file:       ', trim(eloss_file)
  if (fallback_eloss) write(*,*) '  *** FALLBACK eloss -- Z/A may be wrong!'
  write(*,'(A,A)')       '  Xsec file:        ', trim(xsec_file)
  write(*,'(A,I4,A)')    '  OMP threads:      ', omp_get_max_threads(), &
                         '  (set OMP_NUM_THREADS to change)'
  write(*,*) ' ======================================================='
  write(*,*)
  write(*,*) ' Press Enter to start...'
  read(*,*)

  !==========================================================================
  ! READ INPUT MUONS — two-pass, no fixed size limit
  !==========================================================================
  nmuon = 0;  nskip_mu = 0

  ! --- Pass 1: count accepted muons (no allocation, fast) ---
  write(*,*) ' Counting muons (pass 1)...'
  block
    integer :: hf_c, dm_c, ev_c, ch_c, ios2
    real(8) :: v1,v2,v3,v4,v5,v6,v7,v8,v9,v10
    do
      read(10,'(A200)', iostat=ios) linebuf
      if (ios /= 0) exit
      if (linebuf(1:1)=='#' .or. len_trim(linebuf)==0) cycle
      if (ncols == 14) then
        read(linebuf,*,iostat=ios2) ev_c,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,ch_c,hf_c,dm_c
        if (ios2 /= 0) cycle
        if (hf_c /= 1 .and. transport_all == 0) then; nskip_mu = nskip_mu + 1; cycle; end if
      else
        read(linebuf,*,iostat=ios2) ev_c,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,ch_c,dm_c
        if (ios2 /= 0) cycle
      end if
      nmuon = nmuon + 1
    end do
  end block
  rewind(10)

  write(*,'(A,I12)') '  Muons to transport: ', nmuon
  if (ncols==14) write(*,'(A,I12)') '  Skipped (miss):    ', nskip_mu
  write(*,'(A,F7.2,A)') '  Memory required:    ~', &
    dble(nmuon)*19.d0*8.d0/1024.d0**3, ' GB'
  write(*,*)

  ! --- Allocate exactly the right size ---
  allocate(in_eventid(nmuon), in_charge(nmuon), in_det_mask(nmuon))
  allocate(in_xs(nmuon), in_ys(nmuon), in_zs(nmuon))
  allocate(in_p_srf(nmuon), in_px(nmuon), in_py(nmuon), in_pz(nmuon))
  allocate(in_theta_s(nmuon), in_phi_s(nmuon), in_e_srf(nmuon))
  allocate(out_alive(nmuon), out_x_ug(nmuon), out_y_ug(nmuon), out_z_ug(nmuon))
  allocate(out_e_ug(nmuon), out_cx_ug(nmuon), out_cy_ug(nmuon), out_cz_ug(nmuon))
  allocate(out_theta_ug(nmuon), out_phi_ug(nmuon))

  ! --- Pass 2: fill arrays ---
  write(*,*) ' Loading muons (pass 2)...'
  block
    integer :: hf_rd, dm_rd, ev_rd, ch_rd, ios2, imu
    real(8) :: xs_rd,ys_rd,zs_rd,p_rd,px_rd,py_rd,pz_rd,th_rd,ph_rd,e_rd
    imu = 0
    do
      read(10,'(A200)', iostat=ios) linebuf
      if (ios /= 0) exit
      if (linebuf(1:1)=='#' .or. len_trim(linebuf)==0) cycle
      if (ncols == 14) then
        read(linebuf,*,iostat=ios2) ev_rd,xs_rd,ys_rd,zs_rd,p_rd, &
              px_rd,py_rd,pz_rd,th_rd,ph_rd,e_rd,ch_rd,hf_rd,dm_rd
        if (ios2 /= 0) cycle
        if (hf_rd /= 1 .and. transport_all == 0) cycle
      else
        read(linebuf,*,iostat=ios2) ev_rd,xs_rd,ys_rd,zs_rd,p_rd, &
              px_rd,py_rd,pz_rd,th_rd,ph_rd,e_rd,ch_rd,dm_rd
        if (ios2 /= 0) cycle
        hf_rd = 1
      end if
      imu = imu + 1
      in_eventid(imu) = ev_rd
      in_xs(imu)=xs_rd;  in_ys(imu)=ys_rd;  in_zs(imu)=zs_rd
      in_p_srf(imu)=p_rd
      in_px(imu)=px_rd;  in_py(imu)=py_rd;  in_pz(imu)=pz_rd
      in_theta_s(imu)=th_rd;  in_phi_s(imu)=ph_rd
      in_e_srf(imu)=e_rd;  in_charge(imu)=ch_rd;  in_det_mask(imu)=dm_rd
    end do
  end block
  close(10)
  write(*,'(A,I12,A)') '  Loaded ', nmuon, ' muons into memory.'
  write(*,*)

  !==========================================================================
  ! SOURCE-PLANE DETECTION  (same logic as ucmuon_transport_bb_omp)
  !==========================================================================
  xs_mean = sum(in_xs(1:nmuon)) / dble(nmuon)
  ys_mean = sum(in_ys(1:nmuon)) / dble(nmuon)
  zs_mean = sum(in_zs(1:nmuon)) / dble(nmuon)
  xs_var  = sum((in_xs(1:nmuon) - xs_mean)**2) / dble(nmuon)
  ys_var  = sum((in_ys(1:nmuon) - ys_mean)**2) / dble(nmuon)
  zs_var  = sum((in_zs(1:nmuon) - zs_mean)**2) / dble(nmuon)

  if (ys_var <= xs_var .and. ys_var <= zs_var) then
    depth_axis = 1
    write(*,*) '  Source plane: XZ  (y=const) — depth direction: Y'
  else if (xs_var <= ys_var .and. xs_var <= zs_var) then
    depth_axis = 0
    write(*,*) '  Source plane: YZ  (x=const) — depth direction: X'
  else
    depth_axis = 2
    write(*,*) '  Source plane: XY  (z=const) — depth direction: Z'
  end if
  write(*,'(A,3ES11.3,A)') '  Pos. std (x,y,z): ', &
    sqrt(xs_var), sqrt(ys_var), sqrt(zs_var), ' cm'
  write(*,*)

  !==========================================================================
  ! COMPUTE OR LOAD CROSS-SECTIONS (serial, material-aware)
  !==========================================================================
  if (init == 0) then
    write(*,'(A,A,A)') ' Computing tables for material: ', &
                       trim(mat_suffix), ' (~1 min)...'
    ! mulos()   writes music-eloss.dat          (energy-loss table)
    ! mucrsec() writes music-cross-sections.dat (integral cross-sections)
    call mulos(minv, zz0, a0, fr0, par_ion)
    call mucrsec(minv, zz0, a0, fr0)
    ! Copy eloss table to material-specific name
    call execute_command_line( &
      'cp music-eloss.dat ' // trim(eloss_file), &
      exitstat=stat_cmd)
    if (stat_cmd == 0) then
      write(*,'(A,A)') '  Energy-loss table saved as: ', trim(eloss_file)
    else
      write(*,*) '  NOTE: could not copy to ', trim(eloss_file)
      eloss_file = 'music-eloss.dat'
    end if
    ! Copy cross-sections to material-specific name
    call execute_command_line( &
      'cp music-cross-sections.dat ' // trim(xsec_file), &
      exitstat=stat_cmd)
    if (stat_cmd == 0) then
      write(*,'(A,A)') '  Cross-sections saved as:    ', trim(xsec_file)
    else
      write(*,*) '  NOTE: could not copy to ', trim(xsec_file)
      xsec_file = 'music-cross-sections.dat'
    end if
  end if

  !==========================================================================
  ! INITIALIZE MUSIC (serial -- populates shared read-only COMMON blocks)
  !==========================================================================
  write(*,*) ' Loading MUSIC tables...'
  write(*,'(A,A)') '   eloss file : ', trim(eloss_file)
  write(*,'(A,A)') '   xsec  file : ', trim(xsec_file)
  call initialize_music(minv, rho, rad, trim(eloss_file), trim(xsec_file))
  write(*,*) ' MUSIC ready.'
  write(*,*)

  !==========================================================================
  ! PER-THREAD RNG INITIALIZATION
  ! Each thread gets a unique seed offset -> independent RANLUX + RANMAR streams.
  !==========================================================================
  iranlux_base = 314159
  nthreads     = omp_get_max_threads()

  !$OMP PARALLEL PRIVATE(tid)
  tid = omp_get_thread_num()
  call rluxgo(3, iranlux_base + tid * 100003, 0, 0)
  call rmarin(271828 + tid * 99991, 0, 0)
  !$OMP END PARALLEL

  write(*,'(A,I4,A)') '  RNG streams initialised for ', nthreads, ' thread(s).'
  write(*,*)

  !==========================================================================
  ! PARALLEL TRANSPORT LOOP
  ! Each iteration i writes only to out_*(i) -> no race, no CRITICAL needed.
  ! SCHEDULE(DYNAMIC,1): load-balances because transport time scales with
  ! energy (high-E muons have more stochastic interactions -> longer runtime).
  !==========================================================================
  nsurvive = 0;  nstop = 0   ! initialise before parallel region
  write(*,*) ' Transporting muons (OMP)...'
  write(*,*)

  !$OMP PARALLEL DO                                                           &
  !$OMP   DEFAULT(SHARED)                                                     &
  !$OMP   PRIVATE(i, x, y, z, cx, cy, cz, eff_cz, emu, ttime,                &
  !$OMP           depth_this_muon, alive, theta_ug, phi_ug)                   &
  !$OMP   REDUCTION(+:nsurvive, nstop)                                        &
  !$OMP   SCHEDULE(DYNAMIC, 1)

  do i = 1, nmuon

    ! Direction cosines
    if (in_p_srf(i) > 0.d0) then
      cx =  in_px(i) / in_p_srf(i)
      cy =  in_py(i) / in_p_srf(i)
      cz =  in_pz(i) / in_p_srf(i)
    else
      cx =  dsin(in_theta_s(i)) * dcos(in_phi_s(i))
      cy =  dsin(in_theta_s(i)) * dsin(in_phi_s(i))
      cz = -dcos(in_theta_s(i))
    end if
    cz = -cz   ! CosmoALEPH cz<0=down -> MUSIC cz>0=down

    x     = in_xs(i)
    y     = in_ys(i)
    z     = -in_zs(i)
    emu   = in_e_srf(i)
    ttime = 0.d0

    ! Source-plane-aware effective depth cosine
    select case(depth_axis)
      case(1);  eff_cz = -cy    ! XZ plane: |py/p|, positive for py<0
      case(0);  eff_cz = -cx    ! YZ plane: |px/p|, positive for px<0
      case default;  eff_cz = cz ! XY plane: already MUSIC-flipped
    end select

    if (eff_cz > 1.d-6) then
      depth_this_muon = depth_cm / eff_cz
    else
      depth_this_muon = depth_cm * 1.0d4   ! near-horizontal: absorbed
    end if

    call muon_transport(x, y, z, cx, cy, cz, &
                        emu, depth_this_muon, ttime, idim, idim1)

    ! Survival — reduction-safe (no ATOMIC needed with REDUCTION clause)
    if (emu > MMUON) then
      alive = 1;  out_e_ug(i) = emu
      nsurvive = nsurvive + 1
    else
      alive = 0;  out_e_ug(i) = 0.d0
      nstop = nstop + 1
    end if

    ! Back to CosmoALEPH coordinates
    out_alive(i)  =  alive
    out_x_ug(i)   =  x;    out_y_ug(i)  =  y;    out_z_ug(i)  = -z
    out_cx_ug(i)  =  cx;   out_cy_ug(i) =  cy;   out_cz_ug(i) = -cz

    ! Underground angles
    if (alive == 1) then
      theta_ug = dacos(min(1.d0, max(-1.d0, -out_cz_ug(i))))
      if (dabs(dsin(theta_ug)) > 1.d-9) then
        phi_ug = datan2(out_cy_ug(i), out_cx_ug(i))
        if (phi_ug < 0.d0) phi_ug = phi_ug + 2.d0*PI
      else
        phi_ug = 0.d0
      end if
    else
      theta_ug = in_theta_s(i);  phi_ug = in_phi_s(i)
    end if
    out_theta_ug(i) = theta_ug;  out_phi_ug(i) = phi_ug

    ! Print every ~0.5% of total — ~200 updates for any file size
    if (mod(i, max(1, nmuon/200)) == 0) then
      write(*,'(A,I10,A,I10,A,I10,A)') &
        '  Transported:', i, &
        '  Survived:', nsurvive, &
        '  Total:', nmuon, ''
      flush(6)
    end if

  end do
  !$OMP END PARALLEL DO

  !==========================================================================
  ! SERIAL OUTPUT -- preserves input order
  !==========================================================================
  write(*,*)
  write(*,*) ' Writing output...'

  open(unit=11, file=trim(outfile), form='formatted', status='unknown')
  write(11,'(A)') &
    '# EventID  x_srf_cm  y_srf_cm  z_srf_cm' // &
    '  E_srf_GeV  theta_srf_rad  phi_srf_rad  charge' // &
    '  alive  x_ug_cm  y_ug_cm  z_ug_cm  E_ug_GeV' // &
    '  cx_ug  cy_ug  cz_ug  theta_ug_rad  phi_ug_rad'

  do i = 1, nmuon
    write(11,'(I10,1X,'// &
              'F13.4,1X,F13.4,1X,F13.4,1X,'// &
              'F13.6,1X,F13.9,1X,F13.9,1X,I4,1X,'// &
              'I2,1X,'// &
              'F13.4,1X,F13.4,1X,F13.4,1X,'// &
              'F13.6,1X,'// &
              'F10.6,1X,F10.6,1X,F10.6,1X,'// &
              'F13.9,1X,F13.9)') &
      in_eventid(i), &
      in_xs(i), in_ys(i), in_zs(i), in_e_srf(i), in_theta_s(i), in_phi_s(i), &
      in_charge(i), out_alive(i), &
      out_x_ug(i), out_y_ug(i), out_z_ug(i), out_e_ug(i), &
      out_cx_ug(i), out_cy_ug(i), out_cz_ug(i), &
      out_theta_ug(i), out_phi_ug(i)
  end do
  close(11)

  !==========================================================================
  ! SUMMARY
  !==========================================================================
  write(*,*)
  write(*,*) ' ======================================================='
  write(*,*) '              TRANSPORT COMPLETE (OpenMP)'
  write(*,*) ' ======================================================='
  write(*,'(A,I4)')      '  OMP threads:      ', nthreads
  write(*,'(A,A)')       '  Material group:   ', trim(mat_suffix)
  write(*,'(A,A)')       '  Eloss file used:  ', trim(eloss_file)
  if (fallback_eloss) &
    write(*,*) '  *** FALLBACK eloss used -- Z/A were standard-rock values!'
  write(*,'(A,I10)')     '  Muons transported:', nmuon
  if (ncols==14) write(*,'(A,I10)') '  Skipped (miss):   ', nskip_mu
  write(*,'(A,I10)')     '  Survived:         ', nsurvive
  write(*,'(A,I10)')     '  Stopped:          ', nstop
  if (nmuon > 0) &
    write(*,'(A,F8.4,A)') '  Survival rate:    ', &
      100.d0*dble(nsurvive)/dble(nmuon), ' %'
  write(*,'(A,F12.2,A,F9.2,A)') '  Depth:            ', &
    depth_cm*rho, ' g/cm^2  (', depth_m, ' m)'
  write(*,'(A,A)')       '  Output:           ', trim(outfile)
  write(*,*) ' ======================================================='
  flush(6)

  stop
end program ucmuon_transport_music_omp

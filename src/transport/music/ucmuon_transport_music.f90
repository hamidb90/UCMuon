!=============================================================================
! ucmuon_transport_music.f90  —  UCMuon MUSIC Transport Driver
! UCLouvain Muography Group  |  Hamid Basiri <hamid.basiri@uclouvain.be>
!
! Muon transport through rock/water using the MUSIC code
! (Kudryavtsev 2009, CPC 180, 339).
!
! Parallelisation model — MPI + OpenMP hybrid
! ────────────────────────────────────────────
!   MPI  : input muon list is split across ranks via MPI_Scatterv.
!          Rank 0 reads the file, distributes chunks; no inter-rank
!          communication during transport.
!   OMP  : each rank transports its chunk with an OpenMP parallel do loop
!          (SCHEDULE DYNAMIC to balance high/low-energy muons).
!   Work split:
!          base = floor(N / nranks)
!          rank 0 gets  base + (N mod nranks)  (absorbs remainder)
!          ranks 1..nranks-1 get  base  each
!
! MUSIC initialisation strategy
! ──────────────────────────────
!   init=0 (first run):  rank 0 calls mulos() + mucrsec() to generate
!          the energy-loss and cross-section table files on disk, then
!          MPI_Barrier ensures all ranks see the files before calling
!          initialize_music().
!   init=1 (tables exist): all ranks call initialize_music() directly.
!   MUSIC COMMON blocks are read-only after initialisation → thread-safe.
!
! Input file format (flexible, # comments stripped by SLURM script)
! ──────────────────────────────────────────────────────────────────
!   Compatible with ucmuon_selected.dat / ucmuon_surface.dat written
!   by ucmuon_gen.  Auto-detects 13-col (no hit_flag) or 14-col format.
!
! Output
! ──────
!   Per-rank files: ucmuon_underground_RRRRR.dat
!   SLURM script concatenates them into a single output file.
!   Format: 18 columns (identical to existing MUSIC driver output).
!
! Required files in run directory
! ────────────────────────────────
!   music-eloss-rock.dat          (or water / seawater variant)
!   music-cross-sections-rock.dat (auto-generated on first run if init=0)
!   music-double-diff-rock.dat    (needed only for init=0)
!   ucmuon_selected.dat           (or whatever input_file is set to)
!
! Build
! ─────
!   make ucmuon_transport_music     (see Makefile)
!   Requires: mpif90 wrapper + gfortran + OpenMP
!             music.o  music-crosssections.o  ranlux_omp.o  ranmar_omp.o
!             rnorml.o  corgen.o
!
! Usage
! ─────
!   sbatch run_ucmuon_transport.sh input_transport.dat
!=============================================================================
program ucmuon_transport_music
  use mpi
  use omp_lib
  implicit none

  !--- External MUSIC routines (from music.f / music-crosssections.f) ---------
  external :: rluxgo, rmarin
  external :: mucrsec, mulos, initialize_music, muon_transport

  !--- Physical constants ------------------------------------------------------
  real(8), parameter :: PI    = 3.141592653589793d0
  real(8), parameter :: MMUON = 0.105658370d0     ! muon mass [GeV/c²]

  !--- MPI variables -----------------------------------------------------------
  integer :: my_rank, nranks, ierr
  integer :: local_nmuon                           ! this rank's muon count
  integer, allocatable :: sendcounts(:), displs(:) ! for MPI_Scatterv

  !--- User parameters (read on rank 0, broadcast to all) ----------------------
  character(200) :: infile, outfile
  real(8)        :: rho, rad, depth_m, depth_cm
  integer        :: idim, idim1, minv, init_tables
  integer        :: mat_type, transport_all
  character(60)  :: mat_suffix, eloss_file, xsec_file
  integer        :: ncols                          ! 13 or 14 (auto-detected)

  !--- Material composition arrays (for MUSIC) ---------------------------------
  real(4) :: zz0(20), a0(20), fr0(20), par_ion(6)

  real(4), parameter :: zz0_rock(20)    = (/11.0,0.,0.,0.,0.,0.,0.,0.,0.,0., &
                                              0., 0.,0.,0.,0.,0.,0.,0.,0.,0./)
  real(4), parameter :: a0_rock(20)     = (/22.0,0.,0.,0.,0.,0.,0.,0.,0.,0., &
                                              0., 0.,0.,0.,0.,0.,0.,0.,0.,0./)
  real(4), parameter :: fr0_rock(20)    = (/1.0, 0.,0.,0.,0.,0.,0.,0.,0.,0., &
                                              0., 0.,0.,0.,0.,0.,0.,0.,0.,0./)
  real(4), parameter :: par_ion_rock(6) = (/136.4,-3.774,0.083,3.412,3.055,0.049/)

  real(4), parameter :: zz0_water(20)   = (/1.0,8.0,0.,0.,0.,0.,0.,0.,0.,0., &
                                              0., 0.,0.,0.,0.,0.,0.,0.,0.,0./)
  real(4), parameter :: a0_water(20)    = (/1.008,15.999,0.,0.,0.,0.,0.,0.,0.,0., &
                                              0.,  0.,   0.,0.,0.,0.,0.,0.,0.,0./)
  real(4), parameter :: fr0_water(20)   = (/0.1119,0.8881,0.,0.,0.,0.,0.,0.,0.,0., &
                                              0.,   0.,   0.,0.,0.,0.,0.,0.,0.,0./)
  real(4), parameter :: par_ion_water(6)= (/75.0,-3.502,0.2065,3.007,2.5,0.24/)

  real(4), parameter :: zz0_sw(20)      = (/1.0,8.0,11.0,17.0,0.,0.,0.,0.,0.,0., &
                                              0., 0., 0.,  0., 0.,0.,0.,0.,0.,0./)
  real(4), parameter :: a0_sw(20)       = (/1.008,15.999,22.990,35.453,0.,0.,0.,0.,0.,0., &
                                              0.,  0.,    0.,    0.,   0.,0.,0.,0.,0.,0./)
  real(4), parameter :: fr0_sw(20)      = (/0.1100,0.8779,0.0106,0.0019,0.,0.,0.,0.,0.,0., &
                                              0.,   0.,    0.,    0.,   0.,0.,0.,0.,0.,0./)
  real(4), parameter :: par_ion_sw(6)   = (/75.0,-3.502,0.2065,3.007,2.5,0.24/)

  !--- Global arrays (rank 0 only during read, then scattered) -----------------
  integer,  allocatable :: g_eventid(:), g_charge(:), g_det_mask(:)
  real(8),  allocatable :: g_xs(:), g_ys(:), g_zs(:)
  real(8),  allocatable :: g_p(:), g_px(:), g_py(:), g_pz(:)
  real(8),  allocatable :: g_theta(:), g_phi(:), g_e(:)

  !--- Local arrays (each rank's chunk after scatter) --------------------------
  integer,  allocatable :: l_eventid(:), l_charge(:), l_det_mask(:)
  real(8),  allocatable :: l_xs(:), l_ys(:), l_zs(:)
  real(8),  allocatable :: l_p(:), l_px(:), l_py(:), l_pz(:)
  real(8),  allocatable :: l_theta(:), l_phi(:), l_e(:)

  !--- Output arrays (per rank, written to per-rank file) ----------------------
  integer,  allocatable :: out_alive(:)
  real(8),  allocatable :: out_x_ug(:), out_y_ug(:), out_z_ug(:)
  real(8),  allocatable :: out_e_ug(:)
  real(8),  allocatable :: out_cx_ug(:), out_cy_ug(:), out_cz_ug(:)
  real(8),  allocatable :: out_theta_ug(:), out_phi_ug(:)

  !--- Per-event transport variables (private in OMP loop) ---------------------
  real(8)  :: x, y, z, cx, cy, cz, emu, ttime, slant_depth
  real(8)  :: theta_ug, phi_ug
  integer  :: alive

  !--- Counters ----------------------------------------------------------------
  integer  :: nmuon_total           ! total muons across all ranks
  integer  :: nsurvive, nstop       ! per-rank counters (OMP REDUCTION)
  integer  :: nsurvive_all, nstop_all ! reduced totals (rank 0)
  integer  :: ios, ios2, i
  integer  :: nthreads, tid, iranlux_base
  integer  :: nskip_mu              ! muons skipped (hit_flag /= 1)

  !--- Source-plane detection ---------------------------------------------------
  integer  :: depth_axis
  real(8)  :: eff_cz, xs_mean, ys_mean, zs_mean, xs_var, ys_var, zs_var

  !--- Misc --------------------------------------------------------------------
  character(200) :: linebuf
  character(200) :: outfile_rank    ! per-rank output filename
  integer        :: stat_cmd
  logical        :: file_exists
  integer        :: imu             ! loop index for file read (rank 0)
  real(8)        :: t_start, t_end  ! timing

  !===========================================================================
  ! MPI INITIALISATION
  !===========================================================================
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nranks,  ierr)

  allocate(sendcounts(0:nranks-1), displs(0:nranks-1))

  !===========================================================================
  ! BANNER
  !===========================================================================
  if (my_rank == 0) then
    write(*,*)
    write(*,*) ' ============================================================'
    write(*,*) '   UCMuon Transport  —  MUSIC engine'
    write(*,*) '   UCLouvain Muography Group'
    write(*,*) '   ** MPI + OpenMP hybrid version **'
    write(*,*) ' ============================================================'
    write(*,'(A,I6,A,I4,A)') '   MPI ranks: ', nranks, &
                              '   OMP threads/rank: ', omp_get_max_threads(), &
                              '   (set OMP_NUM_THREADS to change)'
    write(*,*) ' ============================================================'
    write(*,*)
  end if

  !===========================================================================
  ! INPUT — rank 0 reads stdin, broadcasts all parameters
  !===========================================================================
  if (my_rank == 0) then

    write(*,*) ' --- [1/6] Input / Output files ---'
    write(*,*) ' Input file (Enter = ucmuon_selected.dat):'
    read(*,'(A)') infile
    if (len_trim(infile) == 0) infile = 'ucmuon_selected.dat'

    write(*,*) ' Output file prefix (Enter = ucmuon_underground):'
    read(*,'(A)') outfile
    if (len_trim(outfile) == 0) outfile = 'ucmuon_underground'

    write(*,*) ' --- [2/6] Material Properties ---'
    write(*,*) ' Material type:'
    write(*,*) '   1 = Rock     (Standard Rock  Z=11, A=22)'
    write(*,*) '   2 = Water / Ice'
    write(*,*) '   3 = Seawater  (H, O, Na, Cl)'
    read(*,*) mat_type
    if (mat_type < 1 .or. mat_type > 3) mat_type = 1

    write(*,*) ' Rock/material density [g/cm³]  (e.g. 2.65 std rock, 1.7 volcanic):'
    read(*,*) rho
    if (rho <= 0.d0) rho = 2.65d0

    write(*,*) ' Radiation length [cm]  (e.g. 26.48 std rock, 36.08 water):'
    read(*,*) rad
    if (rad <= 0.d0) rad = 26.48d0

    write(*,*) ' --- [3/6] Geometry ---'
    write(*,*) ' Vertical depth to detector [m]:'
    read(*,*) depth_m
    depth_cm = depth_m * 100.d0

    write(*,*) ' --- [4/6] Transport Physics ---'
    write(*,*) ' 3D lateral transport / multiple scattering?  (1=ON  0=OFF):'
    write(*,*) '   1 = recommended for inclined muons (slower)'
    write(*,*) '   0 = 1D depth-only (faster, good for near-vertical)'
    read(*,*) idim

    write(*,*) ' Other-process scattering (idim1)?  (1=ON  0=OFF):'
    write(*,*) '   1 = include nuclear interactions (standard)'
    write(*,*) '   0 = suppress (testing only)'
    read(*,*) idim1

    write(*,*) ' Energy-loss cut exponent (minv)  [expert, default -30]:'
    write(*,*) '   -30 is the standard MUSIC value — do not change unless'
    write(*,*) '   you know what you are doing.'
    read(*,*) minv

    write(*,*) ' --- [5/6] Cross-section Tables ---'
    write(*,*) ' Cross-section table status:'
    write(*,*) '   0 = compute + save to disk  (first run, ~1 min)'
    write(*,*) '   1 = read from disk          (subsequent runs, fast)'
    write(*,*) '   NOTE: files are named  music-eloss-{mat}.dat'
    write(*,*) '         and  music-cross-sections-{mat}.dat'
    read(*,*) init_tables

    write(*,*) ' --- [6/6] Transport scope ---'
    write(*,*) ' Transport ALL muons?'
    write(*,*) '   0 = only muons with hit_flag=1 (detector-aimed, recommended)'
    write(*,*) '   1 = all muons in file (use if file is already pre-filtered)'
    write(*,*) '   Note: auto-detected 13-col format always transports all.'
    read(*,*) transport_all

  end if  ! my_rank == 0

  !===========================================================================
  ! BROADCAST SCALAR PARAMETERS
  !===========================================================================
  call MPI_Bcast(mat_type,      1,   MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(idim,          1,   MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(idim1,         1,   MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(minv,          1,   MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(init_tables,   1,   MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(transport_all, 1,   MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(rho,           1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(rad,           1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(depth_m,       1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(depth_cm,      1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(infile,        200, MPI_CHARACTER,        0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(outfile,       200, MPI_CHARACTER,        0, MPI_COMM_WORLD, ierr)

  !===========================================================================
  ! SET MATERIAL COMPOSITION ARRAYS (all ranks)
  !===========================================================================
  select case (mat_type)
    case (2)
      zz0=zz0_water; a0=a0_water; fr0=fr0_water; par_ion=par_ion_water
      mat_suffix='water'
    case (3)
      zz0=zz0_sw; a0=a0_sw; fr0=fr0_sw; par_ion=par_ion_sw
      mat_suffix='seawater'
    case default
      zz0=zz0_rock; a0=a0_rock; fr0=fr0_rock; par_ion=par_ion_rock
      mat_suffix='rock'
  end select

  eloss_file = 'music-eloss-'          // trim(mat_suffix) // '.dat'
  xsec_file  = 'music-cross-sections-' // trim(mat_suffix) // '.dat'

  if (my_rank == 0) then
    write(*,*)
    write(*,*) ' ============================================================'
    write(*,*) '        UCMuon Transport — CONFIGURATION'
    write(*,*) ' ============================================================'
    write(*,'(A,A)')       '  Input file:        ', trim(infile)
    write(*,'(A,A)')       '  Output prefix:     ', trim(outfile)
    write(*,'(A,A)')       '  Material:          ', trim(mat_suffix)
    write(*,'(A,F7.3,A)')  '  Density:           ', rho, ' g/cm³'
    write(*,'(A,F7.2,A)')  '  Radiation length:  ', rad, ' cm'
    write(*,'(A,F8.2,A)')  '  Vertical depth:    ', depth_m, ' m'
    write(*,'(A,F10.1,A)') '  Depth (w.e.):      ', depth_cm*rho, ' g/cm²'
    write(*,'(A,I2,A,I2)') '  idim=', idim, '  idim1=', idim1
    write(*,'(A,I4)')      '  minv:              ', minv
    write(*,'(A,A)')       '  Eloss file:        ', trim(eloss_file)
    write(*,'(A,A)')       '  Xsec file:         ', trim(xsec_file)
    write(*,'(A,I6,A,I4)') '  MPI ranks:         ', nranks, &
                            '   OMP threads/rank: ', omp_get_max_threads()
    write(*,*) ' ============================================================'
    write(*,*)
  end if

  !===========================================================================
  ! MUSIC TABLE INITIALISATION (rank 0 generates if init=0, all ranks load)
  !===========================================================================
  if (init_tables == 0) then
    !--- Rank 0 generates tables; other ranks wait at the Barrier below -------
    if (my_rank == 0) then
      write(*,*) ' init=0: computing energy-loss and cross-section tables...'
      write(*,*) '  (this takes ~1 minute for Standard Rock)'

      ! music-double-diff-rock.dat: look in root first, then data/
      inquire(file='music-double-diff-rock.dat', exist=file_exists)
      if (.not. file_exists) then
        inquire(file='data/music-double-diff-rock.dat', exist=file_exists)
        if (file_exists) then
          call execute_command_line( &
            'cp data/music-double-diff-rock.dat music-double-diff-rock.dat', &
            exitstat=stat_cmd)
          file_exists = (stat_cmd == 0)
          if (file_exists) &
            write(*,*) '  Copied data/music-double-diff-rock.dat to run directory.'
        end if
      end if
      if (.not. file_exists) then
        write(*,*) ' ERROR: music-double-diff-rock.dat not found.'
        write(*,*) '        Checked ./ and ./data/'
        write(*,*) '        This file ships with the MUSIC distribution.'
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end if

      call mulos(minv, zz0, a0, fr0, par_ion)
      call mucrsec(minv, zz0, a0, fr0)

      ! Rename generic output files to material-specific names in root.
      ! Tables MUST stay in the project root — music.f opens them by name
      ! without path.  Never move them to a subdirectory.
      call execute_command_line('cp music-eloss.dat ' // trim(eloss_file), &
                                exitstat=stat_cmd)
      if (stat_cmd /= 0) eloss_file = 'music-eloss.dat'

      call execute_command_line('cp music-cross-sections.dat ' // trim(xsec_file), &
                                exitstat=stat_cmd)
      if (stat_cmd /= 0) xsec_file = 'music-cross-sections.dat'

      write(*,'(A,A)') '  Tables written to root: ', trim(eloss_file)
      write(*,'(A,A)') '                          ', trim(xsec_file)
    end if  ! my_rank == 0

    ! Wait for rank 0 to finish writing, then broadcast the final filenames
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Bcast(eloss_file, 60, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(xsec_file,  60, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  else  ! init_tables == 1 — read existing tables from disk
    if (my_rank == 0) then
      ! Tables should be in project root (where binary runs).
      ! Also check data/ as a fallback (user may have placed them there).
      inquire(file=trim(eloss_file), exist=file_exists)
      if (.not. file_exists) then
        inquire(file='data/' // trim(eloss_file), exist=file_exists)
        if (file_exists) then
          ! Copy from data/ to root so music.f can open it by name
          call execute_command_line('cp data/' // trim(eloss_file) // &
                                    ' ' // trim(eloss_file), exitstat=stat_cmd)
          if (stat_cmd == 0) then
            write(*,'(A,A)') '  Copied from data/ to root: ', trim(eloss_file)
          else
            write(*,'(A,A)') ' ERROR: could not copy from data/: ', trim(eloss_file)
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
          end if
        else
          write(*,'(A,A)') ' ERROR: eloss file not found: ', trim(eloss_file)
          write(*,*) '  Checked ./ and ./data/'
          write(*,*) '  Run once with init_tables=0 to generate it.'
          call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
      end if
      inquire(file=trim(xsec_file), exist=file_exists)
      if (.not. file_exists) then
        inquire(file='data/' // trim(xsec_file), exist=file_exists)
        if (file_exists) then
          call execute_command_line('cp data/' // trim(xsec_file) // &
                                    ' ' // trim(xsec_file), exitstat=stat_cmd)
          if (stat_cmd == 0) then
            write(*,'(A,A)') '  Copied from data/ to root: ', trim(xsec_file)
          else
            write(*,'(A,A)') ' ERROR: could not copy from data/: ', trim(xsec_file)
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
          end if
        else
          write(*,'(A,A)') ' ERROR: xsec file not found: ', trim(xsec_file)
          write(*,*) '  Checked ./ and ./data/'
          write(*,*) '  Run once with init_tables=0 to generate it.'
          call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
      end if
    end if  ! my_rank == 0

    ! Broadcast filenames (unchanged but ensures all ranks agree)
    call MPI_Bcast(eloss_file, 60, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(xsec_file,  60, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  end if

  ! All ranks load MUSIC tables from disk (read-only COMMON blocks after this)
  if (my_rank == 0) write(*,*) ' Loading MUSIC tables (all ranks)...'
  call initialize_music(minv, rho, rad, trim(eloss_file), trim(xsec_file))
  if (my_rank == 0) write(*,*) ' MUSIC ready.'
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  !===========================================================================
  ! READ INPUT FILE — rank 0 only, two-pass (count then load)
  !===========================================================================
  nmuon_total = 0
  ncols       = 0
  nskip_mu    = 0

  if (my_rank == 0) then
    write(*,*)
    write(*,*) ' Reading input file...'

    open(unit=10, file=trim(infile), form='formatted', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') ' ERROR: cannot open input file: ', trim(infile)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    !--- Auto-detect column format (13 or 14) ---
    block
      integer  :: hf_tmp, dm_tmp, ev_tmp, ch_tmp, ios3
      real(8)  :: p_t,px_t,py_t,pz_t,th_t,ph_t,e_t,xs_t,ys_t,zs_t
      do
        read(10,'(A200)', iostat=ios) linebuf
        if (ios /= 0) exit
        if (linebuf(1:1)=='#' .or. len_trim(linebuf)==0) cycle
        read(linebuf,*,iostat=ios3) ev_tmp,xs_t,ys_t,zs_t,p_t, &
              px_t,py_t,pz_t,th_t,ph_t,e_t,ch_tmp,hf_tmp,dm_tmp
        if (ios3==0) then
          ncols = 14
        else
          read(linebuf,*,iostat=ios3) ev_tmp,xs_t,ys_t,zs_t,p_t, &
                px_t,py_t,pz_t,th_t,ph_t,e_t,ch_tmp,dm_tmp
          if (ios3==0) then
            ncols = 13
          else
            write(*,*) ' ERROR: cannot parse first data line of ', trim(infile)
            write(*,*) '        Expected 13 or 14 columns.'
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
          end if
        end if
        exit
      end do
      rewind(10)
    end block

    write(*,'(A,I2,A)') '  Detected format: ', ncols, ' columns'
    if (ncols == 13) transport_all = 1  ! 13-col has no hit_flag → transport all

    !--- Pass 1: count muons to transport ---
    write(*,*) '  Pass 1: counting muons...'
    block
      integer :: hf_c, dm_c, ev_c, ch_c, ios3
      real(8) :: v1,v2,v3,v4,v5,v6,v7,v8,v9,v10
      do
        read(10,'(A200)', iostat=ios) linebuf
        if (ios /= 0) exit
        if (linebuf(1:1)=='#' .or. len_trim(linebuf)==0) cycle
        if (ncols == 14) then
          read(linebuf,*,iostat=ios3) ev_c,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10, &
                ch_c,hf_c,dm_c
          if (ios3 /= 0) cycle
          if (hf_c /= 1 .and. transport_all == 0) then
            nskip_mu = nskip_mu + 1;  cycle
          end if
        else
          read(linebuf,*,iostat=ios3) ev_c,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10, &
                ch_c,dm_c
          if (ios3 /= 0) cycle
        end if
        nmuon_total = nmuon_total + 1
      end do
    end block
    rewind(10)

    write(*,'(A,I12)') '  Muons to transport: ', nmuon_total
    if (ncols==14 .and. transport_all==0) &
      write(*,'(A,I12)') '  Skipped (missed):   ', nskip_mu
    write(*,'(A,F7.3,A)') '  Memory on rank 0:  ~', &
      dble(nmuon_total)*14.d0*8.d0/1024.d0**3, ' GB'
    write(*,*)

    if (nmuon_total == 0) then
      write(*,*) ' ERROR: no muons found in input file.'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    !--- Allocate global arrays (rank 0 only) ---
    allocate(g_eventid(nmuon_total), g_charge(nmuon_total), g_det_mask(nmuon_total))
    allocate(g_xs(nmuon_total), g_ys(nmuon_total), g_zs(nmuon_total))
    allocate(g_p(nmuon_total), g_px(nmuon_total), g_py(nmuon_total), g_pz(nmuon_total))
    allocate(g_theta(nmuon_total), g_phi(nmuon_total), g_e(nmuon_total))

    !--- Pass 2: load muons into arrays ---
    write(*,*) '  Pass 2: loading muons...'
    block
      integer :: hf_rd, dm_rd, ev_rd, ch_rd, ios3
      real(8) :: xs_rd,ys_rd,zs_rd,p_rd,px_rd,py_rd,pz_rd,th_rd,ph_rd,e_rd
      imu = 0
      do
        read(10,'(A200)', iostat=ios) linebuf
        if (ios /= 0) exit
        if (linebuf(1:1)=='#' .or. len_trim(linebuf)==0) cycle
        if (ncols == 14) then
          read(linebuf,*,iostat=ios3) ev_rd,xs_rd,ys_rd,zs_rd,p_rd, &
                px_rd,py_rd,pz_rd,th_rd,ph_rd,e_rd,ch_rd,hf_rd,dm_rd
          if (ios3 /= 0) cycle
          if (hf_rd /= 1 .and. transport_all == 0) cycle
        else
          read(linebuf,*,iostat=ios3) ev_rd,xs_rd,ys_rd,zs_rd,p_rd, &
                px_rd,py_rd,pz_rd,th_rd,ph_rd,e_rd,ch_rd,dm_rd
          if (ios3 /= 0) cycle
          hf_rd = 1
        end if
        imu = imu + 1
        g_eventid(imu) = ev_rd
        g_xs(imu)=xs_rd;  g_ys(imu)=ys_rd;  g_zs(imu)=zs_rd
        g_p(imu)=p_rd
        g_px(imu)=px_rd;  g_py(imu)=py_rd;  g_pz(imu)=pz_rd
        g_theta(imu)=th_rd;  g_phi(imu)=ph_rd
        g_e(imu)=e_rd;  g_charge(imu)=ch_rd;  g_det_mask(imu)=dm_rd
      end do
    end block
    close(10)
    write(*,'(A,I12,A)') '  Loaded ', nmuon_total, ' muons.'

    ! Auto-detect source plane by finding the coordinate with near-zero spread
    xs_mean = sum(g_xs(1:nmuon_total)) / dble(nmuon_total)
    ys_mean = sum(g_ys(1:nmuon_total)) / dble(nmuon_total)
    zs_mean = sum(g_zs(1:nmuon_total)) / dble(nmuon_total)
    xs_var  = sum((g_xs(1:nmuon_total) - xs_mean)**2) / dble(nmuon_total)
    ys_var  = sum((g_ys(1:nmuon_total) - ys_mean)**2) / dble(nmuon_total)
    zs_var  = sum((g_zs(1:nmuon_total) - zs_mean)**2) / dble(nmuon_total)
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

  end if  ! my_rank == 0

  !===========================================================================
  ! BROADCAST total muon count + ncols to all ranks
  !===========================================================================
  call MPI_Bcast(nmuon_total,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ncols,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(transport_all,1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(depth_axis,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  !===========================================================================
  ! COMPUTE SEND COUNTS AND DISPLACEMENTS FOR MPI_Scatterv
  ! Rank 0 gets floor(N/nranks) + (N mod nranks); others get floor(N/nranks)
  !===========================================================================
  block
    integer :: base, rem, r
    base = nmuon_total / nranks
    rem  = mod(nmuon_total, nranks)
    sendcounts(0) = base + rem
    do r = 1, nranks-1
      sendcounts(r) = base
    end do
    displs(0) = 0
    do r = 1, nranks-1
      displs(r) = displs(r-1) + sendcounts(r-1)
    end do
    local_nmuon = sendcounts(my_rank)
  end block

  if (my_rank == 0) then
    write(*,*)
    write(*,'(A,I6,A,I8,A,I8)') '  Work split:  base=', nmuon_total/nranks, &
      '  rank0 gets ', sendcounts(0), '  others get ', sendcounts(1)
    write(*,*)
  end if

  !===========================================================================
  ! ALLOCATE LOCAL ARRAYS (each rank allocates only its chunk)
  !===========================================================================
  allocate(l_eventid(local_nmuon), l_charge(local_nmuon), l_det_mask(local_nmuon))
  allocate(l_xs(local_nmuon), l_ys(local_nmuon), l_zs(local_nmuon))
  allocate(l_p(local_nmuon), l_px(local_nmuon), l_py(local_nmuon), l_pz(local_nmuon))
  allocate(l_theta(local_nmuon), l_phi(local_nmuon), l_e(local_nmuon))

  allocate(out_alive(local_nmuon))
  allocate(out_x_ug(local_nmuon), out_y_ug(local_nmuon), out_z_ug(local_nmuon))
  allocate(out_e_ug(local_nmuon))
  allocate(out_cx_ug(local_nmuon), out_cy_ug(local_nmuon), out_cz_ug(local_nmuon))
  allocate(out_theta_ug(local_nmuon), out_phi_ug(local_nmuon))

  ! Non-rank-0 need dummy global arrays of size 1 for MPI_Scatterv call
  if (my_rank /= 0) then
    allocate(g_eventid(1), g_charge(1), g_det_mask(1))
    allocate(g_xs(1), g_ys(1), g_zs(1))
    allocate(g_p(1), g_px(1), g_py(1), g_pz(1))
    allocate(g_theta(1), g_phi(1), g_e(1))
  end if

  !===========================================================================
  ! SCATTER INPUT DATA — 11 real(8) arrays + 3 integer arrays
  !===========================================================================
  call MPI_Scatterv(g_xs,      sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_xs,      local_nmuon,         MPI_DOUBLE_PRECISION, &
                    0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_ys,      sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_ys,      local_nmuon,         MPI_DOUBLE_PRECISION, &
                    0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_zs,      sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_zs,      local_nmuon,         MPI_DOUBLE_PRECISION, &
                    0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_p,       sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_p,       local_nmuon,         MPI_DOUBLE_PRECISION, &
                    0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_px,      sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_px,      local_nmuon,         MPI_DOUBLE_PRECISION, &
                    0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_py,      sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_py,      local_nmuon,         MPI_DOUBLE_PRECISION, &
                    0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_pz,      sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_pz,      local_nmuon,         MPI_DOUBLE_PRECISION, &
                    0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_theta,   sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_theta,   local_nmuon,         MPI_DOUBLE_PRECISION, &
                    0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_phi,     sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_phi,     local_nmuon,         MPI_DOUBLE_PRECISION, &
                    0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_e,       sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_e,       local_nmuon,         MPI_DOUBLE_PRECISION, &
                    0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_eventid, sendcounts, displs, MPI_INTEGER, &
                    l_eventid, local_nmuon,         MPI_INTEGER, &
                    0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_charge,  sendcounts, displs, MPI_INTEGER, &
                    l_charge,  local_nmuon,         MPI_INTEGER, &
                    0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_det_mask,sendcounts, displs, MPI_INTEGER, &
                    l_det_mask,local_nmuon,         MPI_INTEGER, &
                    0, MPI_COMM_WORLD, ierr)

  ! Free global arrays — no longer needed on any rank
  deallocate(g_eventid, g_charge, g_det_mask)
  deallocate(g_xs, g_ys, g_zs, g_p, g_px, g_py, g_pz, g_theta, g_phi, g_e)

  !===========================================================================
  ! PER-THREAD RNG INITIALISATION
  ! Each rank gets a different base seed; each thread within a rank is further
  ! offset.  Seeds are derived from the prime 99991 for good separation.
  !===========================================================================
  iranlux_base = 314159 + my_rank * 99991
  nthreads     = omp_get_max_threads()

  !$OMP PARALLEL PRIVATE(tid)
  tid = omp_get_thread_num()
  call rluxgo(3, iranlux_base + tid * 100003, 0, 0)
  call rmarin(271828  + my_rank*99991 + tid*77777, 0, 0)
  !$OMP END PARALLEL

  write(*,'(A,I5,A,I4,A)') '  [Rank ', my_rank, '] RNG ready for ', &
                             nthreads, ' thread(s).'
  flush(6)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  !===========================================================================
  ! TRANSPORT LOOP — MPI × OMP
  ! Each rank transports its local_nmuon muons with OMP parallel do.
  ! SCHEDULE(DYNAMIC,1): balances work because high-E muons are slower.
  ! No CRITICAL needed — each thread writes only to out_*(i) (private index).
  !===========================================================================
  if (my_rank == 0) then
    write(*,*)
    write(*,*) ' Transporting muons (MPI × OMP)...'
    write(*,*)
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  t_start = MPI_Wtime()

  nsurvive = 0;  nstop = 0

  !$OMP PARALLEL DO                                              &
  !$OMP   DEFAULT(SHARED)                                        &
  !$OMP   PRIVATE(i, x, y, z, cx, cy, cz, emu, ttime,           &
  !$OMP           slant_depth, eff_cz, alive, theta_ug, phi_ug)  &
  !$OMP   REDUCTION(+:nsurvive, nstop)                           &
  !$OMP   SCHEDULE(DYNAMIC, 1)

  do i = 1, local_nmuon

    !--- Direction cosines from momentum components or theta/phi ---
    if (l_p(i) > 0.d0) then
      cx =  l_px(i) / l_p(i)
      cy =  l_py(i) / l_p(i)
      cz =  l_pz(i) / l_p(i)
    else
      cx =  sin(l_theta(i)) * cos(l_phi(i))
      cy =  sin(l_theta(i)) * sin(l_phi(i))
      cz = -cos(l_theta(i))
    end if
    ! UCMuon convention: cz<0 = downward
    ! MUSIC convention:  cz>0 = downward  → flip sign
    cz = -cz

    x     =  l_xs(i)
    y     =  l_ys(i)
    z     = -l_zs(i)   ! UCMuon z<0 underground → MUSIC z>0
    emu   =  l_e(i)
    ttime =  0.d0

    !--- Slant depth using the correct depth-direction cosine ---
    select case(depth_axis)
      case(1);  eff_cz = -cy    ! XZ plane: depth in Y, cy<0 → eff_cz>0
      case(0);  eff_cz = -cx    ! YZ plane: depth in X
      case default; eff_cz = cz  ! XY plane: cz already flipped positive
    end select
    if (eff_cz > 1.d-6) then
      slant_depth = depth_cm / eff_cz
    else
      slant_depth = depth_cm * 1.0d4  ! near-horizontal → effectively absorbed
    end if

    !--- MUSIC transport (modifies x,y,z,cx,cy,cz,emu in place) ---
    call muon_transport(x, y, z, cx, cy, cz, emu, slant_depth, ttime, idim, idim1)

    !--- Survival check ---
    if (emu > MMUON) then
      alive = 1;  out_e_ug(i) = emu;  nsurvive = nsurvive + 1
    else
      alive = 0;  out_e_ug(i) = 0.d0; nstop    = nstop    + 1
    end if

    !--- Back to UCMuon coordinate convention ---
    out_alive(i) = alive
    out_x_ug(i)  =  x;    out_y_ug(i)  =  y;    out_z_ug(i)  = -z
    out_cx_ug(i) =  cx;   out_cy_ug(i) =  cy;   out_cz_ug(i) = -cz

    !--- Underground angles ---
    if (alive == 1) then
      theta_ug = acos(min(1.d0, max(-1.d0, -out_cz_ug(i))))
      if (abs(sin(theta_ug)) > 1.d-9) then
        phi_ug = atan2(out_cy_ug(i), out_cx_ug(i))
        if (phi_ug < 0.d0) phi_ug = phi_ug + 2.d0*PI
      else
        phi_ug = 0.d0
      end if
    else
      theta_ug = l_theta(i);  phi_ug = l_phi(i)
    end if
    out_theta_ug(i) = theta_ug;  out_phi_ug(i) = phi_ug

    !--- Progress report every ~5% ---
    if (mod(i, max(1, local_nmuon/20)) == 0) then
      write(*,'(A,I5,A,I8,A,I8,A,F6.1,A)') &
        '  [Rank ', my_rank, '] Transported: ', i, &
        '  Survived: ', nsurvive, &
        '  (', 100.d0*dble(i)/dble(local_nmuon), '%)'
      flush(6)
    end if

  end do
  !$OMP END PARALLEL DO

  t_end = MPI_Wtime()

  !===========================================================================
  ! WRITE PER-RANK OUTPUT FILE
  ! Format: 18 columns, identical to existing MUSIC driver output.
  ! Header line written only by rank 0 so merged file has exactly one header.
  !===========================================================================
  write(outfile_rank,'(A,A,I5.5,A)') trim(outfile), '_', my_rank, '.dat'
  open(unit=20, file=trim(outfile_rank), form='formatted', status='unknown')

  if (my_rank == 0) then
    write(20,'(A)') &
      '# EventID  x_srf_cm  y_srf_cm  z_srf_cm' // &
      '  E_srf_GeV  theta_srf_rad  phi_srf_rad  charge' // &
      '  alive  x_ug_cm  y_ug_cm  z_ug_cm  E_ug_GeV' // &
      '  cx_ug  cy_ug  cz_ug  theta_ug_rad  phi_ug_rad'
  end if

  do i = 1, local_nmuon
    write(20,'(I10,1X,'// &
              'F13.4,1X,F13.4,1X,F13.4,1X,'// &
              'F13.6,1X,F13.9,1X,F13.9,1X,I4,1X,'// &
              'I2,1X,'// &
              'F13.4,1X,F13.4,1X,F13.4,1X,'// &
              'F13.6,1X,'// &
              'F10.6,1X,F10.6,1X,F10.6,1X,'// &
              'F13.9,1X,F13.9)') &
      l_eventid(i), &
      l_xs(i), l_ys(i), l_zs(i), l_e(i), l_theta(i), l_phi(i), l_charge(i), &
      out_alive(i), &
      out_x_ug(i), out_y_ug(i), out_z_ug(i), out_e_ug(i), &
      out_cx_ug(i), out_cy_ug(i), out_cz_ug(i), &
      out_theta_ug(i), out_phi_ug(i)
  end do
  close(20)

  write(*,'(A,I5,A,A)') '  [Rank ', my_rank, '] Output: ', trim(outfile_rank)
  flush(6)

  !===========================================================================
  ! GLOBAL STATISTICS — MPI_Reduce to rank 0
  !===========================================================================
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Reduce(nsurvive, nsurvive_all, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(nstop,    nstop_all,    1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if (my_rank == 0) then
    write(*,*)
    write(*,*) ' ============================================================'
    write(*,*) '   UCMuon Transport — MUSIC — COMPLETE'
    write(*,*) ' ============================================================'
    write(*,'(A,I12)')    '  Total transported:  ', nmuon_total
    write(*,'(A,I12)')    '  Survived:           ', nsurvive_all
    write(*,'(A,I12)')    '  Stopped:            ', nstop_all
    if (nmuon_total > 0) &
      write(*,'(A,F8.4,A)') '  Survival rate:      ', &
        100.d0*dble(nsurvive_all)/dble(nmuon_total), ' %'
    write(*,'(A,F8.2,A,F8.2,A)') '  Depth (vertical):   ', depth_m, ' m  (', &
        depth_cm*rho, ' g/cm²)'
    write(*,'(A,A)')      '  Material:           ', trim(mat_suffix)
    write(*,'(A,I6,A,I4)')  '  Parallelisation:    ', nranks, &
        ' MPI ranks × ', omp_get_max_threads(), ' OMP threads'
    write(*,'(A,F10.2,A)') '  Wall time:          ', t_end - t_start, ' s'
    write(*,*)
    write(*,*) '  Per-rank output files: ', trim(outfile), '_RRRRR.dat'
    write(*,*) '  Run post-processing in SLURM script to merge.'
    write(*,*) ' ============================================================'
  end if

  call MPI_Finalize(ierr)
  stop
end program ucmuon_transport_music

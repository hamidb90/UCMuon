!=============================================================================
! ucmuon_transport_bb.f90  —  UCMuon Bethe-Bloch Transport Driver
! UCLouvain Muography Group  |  Hamid Basiri <hamid.basiri@uclouvain.be>
!
! Muon transport using deterministic CSDA Bethe-Bloch energy loss +
! Highland multiple scattering.  Fully self-contained — no external
! table files required (unlike MUSIC).
!
! Physics model
! dE/dX = dE/dX_ion + b_rad * b_rad_shape(E) * E
! Bethe-Bloch ionisation (PDG 2022, eq. 34.5) + Sternheimer density correction.
! Radiative: dE/dX_rad = b_rad * b_rad_shape(E) * E; b_rad is the material
! value at E_tot = 100 GeV, b_rad_shape carries the PDG-2024 energy
! dependence (see ucmuon_transport_bb_omp.f90).
! Multiple scattering: Highland formula (PDG 34.3), 3D rotation per step.
! Integration: CSDA, fixed 10 g/cm2 steps along slant path.
!
! Parallelisation: same MPI+OMP pattern as ucmuon_transport_music.
! Rank 0 reads, MPI_Scatterv distributes, each rank transports with OMP DO.
!
! Build:  make ucmuon_transport_bb
! Usage:  sbatch run_ucmuon_transport.sh input_transport_bb.dat
!=============================================================================
program ucmuon_transport_bb
  use mpi
  use omp_lib
  implicit none

  external :: rluxgo, ranlux

  real(8), parameter :: PI      = 3.141592653589793d0
  real(8), parameter :: TWO_PI  = 6.283185307179586d0
  real(8), parameter :: MMUON   = 0.105658370d0
  real(8), parameter :: ME      = 5.10998918d-4
  real(8), parameter :: K_BB    = 3.07075d-4
  real(8), parameter :: EMIN    = 0.105658370d0 + 1.d-5

  integer :: my_rank, nranks, ierr, local_nmuon
  integer, allocatable :: sendcounts(:), displs(:)

  character(200) :: infile, outfile
  real(8)        :: depth_m, depth_cm
  integer        :: mat_type, ms_enable, transport_all, ncols
  real(8)        :: Z_eff, A_eff, rho_mat, I_eV, X0_gcm2, b_rad, C_dens

  integer,  allocatable :: g_eventid(:), g_charge(:), g_det_mask(:)
  real(8),  allocatable :: g_xs(:), g_ys(:), g_zs(:)
  real(8),  allocatable :: g_p(:), g_px(:), g_py(:), g_pz(:)
  real(8),  allocatable :: g_theta(:), g_phi(:), g_e(:)

  integer,  allocatable :: l_eventid(:), l_charge(:), l_det_mask(:)
  real(8),  allocatable :: l_xs(:), l_ys(:), l_zs(:)
  real(8),  allocatable :: l_p(:), l_px(:), l_py(:), l_pz(:)
  real(8),  allocatable :: l_theta(:), l_phi(:), l_e(:)

  integer,  allocatable :: out_alive(:)
  real(8),  allocatable :: out_x_ug(:), out_y_ug(:), out_z_ug(:), out_e_ug(:)
  real(8),  allocatable :: out_cx_ug(:), out_cy_ug(:), out_cz_ug(:)
  real(8),  allocatable :: out_theta_ug(:), out_phi_ug(:)

  real(8)  :: x, y, z, cx, cy, cz, emu, slant_cm, theta_ug, phi_ug
  integer  :: alive, nmuon_total, nsurvive, nstop, nsurvive_all, nstop_all
  integer  :: nskip_mu, ios, ios2, imu, i, nthreads, tid, iranlux_base
  character(200) :: linebuf, outfile_rank
  real(8)        :: t_start, t_end, hbar_omega_p
  integer  :: depth_axis
  real(8)  :: eff_cz, xs_mean, ys_mean, zs_mean, xs_var, ys_var, zs_var

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nranks,  ierr)
  allocate(sendcounts(0:nranks-1), displs(0:nranks-1))

  !===========================================================================
  if (my_rank == 0) then
    write(*,*)
    write(*,*) ' ============================================================'
    write(*,*) '   UCMuon Transport  —  Bethe-Bloch + MS engine'
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
  ! INPUT
  !===========================================================================
  if (my_rank == 0) then
    write(*,*) ' --- [1/5] Input / Output files ---'
    write(*,*) ' Input file (Enter = ucmuon_selected.dat):'
    read(*,'(A)') infile
    if (len_trim(infile) == 0) infile = 'ucmuon_selected.dat'

    write(*,*) ' Output prefix (Enter = ucmuon_underground_bb):'
    read(*,'(A)') outfile
    if (len_trim(outfile) == 0) outfile = 'ucmuon_underground_bb'

    write(*,*) ' --- [2/5] Material ---'
    write(*,*) ' Material type:'
    write(*,*) '   1 = Standard Rock  (Z=11,   A=22,    rho=2.65, I=136.4 eV)'
    write(*,*) '   2 = Ice            (Z=7.42,  A=14.99, rho=0.917, I=79.7 eV)'
    write(*,*) '   3 = Water          (Z=7.42,  A=14.99, rho=1.000, I=79.7 eV)'
    write(*,*) '   4 = Concrete       (Z=11.11, A=22.08, rho=2.300, I=135.2 eV)'
    write(*,*) '   5 = Custom         (you supply Z, A, rho, I)'
    read(*,*) mat_type

    select case (mat_type)
      case (1)
        Z_eff=11.0d0; A_eff=22.0d0; rho_mat=2.65d0; I_eV=136.4d0
        ! b_rad = value at E_tot = 100 GeV; energy dependence via b_rad_shape
        X0_gcm2=26.54d0; b_rad=3.02d-6
      case (2)
        Z_eff=7.42d0; A_eff=14.99d0; rho_mat=0.917d0; I_eV=79.7d0
        X0_gcm2=36.08d0; b_rad=3.40d-6
      case (3)
        Z_eff=7.42d0; A_eff=14.99d0; rho_mat=1.000d0; I_eV=79.7d0
        X0_gcm2=36.08d0; b_rad=3.40d-6
      case (4)
        Z_eff=11.11d0; A_eff=22.08d0; rho_mat=2.300d0; I_eV=135.2d0
        X0_gcm2=26.70d0; b_rad=3.00d-6
      case (5)
        write(*,*) ' Zeff:';       read(*,*) Z_eff
        write(*,*) ' Aeff:';       read(*,*) A_eff
        write(*,*) ' Density [g/cm3]:';   read(*,*) rho_mat
        write(*,*) ' Mean excitation energy I [eV]:'; read(*,*) I_eV
        X0_gcm2 = 716.408d0*A_eff / (Z_eff*(Z_eff+1.d0)*(11.319d0-log(Z_eff)))
        b_rad   = max(3.02d-6*(Z_eff**2/A_eff)/(121.d0/22.d0), 1.d-7)
      case default
        write(*,*) ' Unknown mat_type — defaulting to Standard Rock.'
        Z_eff=11.0d0; A_eff=22.0d0; rho_mat=2.65d0; I_eV=136.4d0
        X0_gcm2=26.54d0; b_rad=3.02d-6; mat_type=1
    end select

    write(*,*) ' --- [3/5] Geometry ---'
    write(*,*) ' Vertical depth to detector [m]:'
    read(*,*) depth_m
    depth_cm = depth_m * 100.d0

    write(*,*) ' --- [4/5] Physics options ---'
    write(*,*) ' Multiple scattering (Highland)?  (1=ON  0=OFF):'
    read(*,*) ms_enable

    write(*,*) ' --- [5/5] Transport scope ---'
    write(*,*) ' Transport all muons? (0=hit_flag=1 only  1=all):'
    read(*,*) transport_all
  end if

  !===========================================================================
  ! BROADCAST
  !===========================================================================
  call MPI_Bcast(mat_type,      1,   MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ms_enable,     1,   MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(transport_all, 1,   MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(depth_m,       1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(depth_cm,      1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(Z_eff,         1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(A_eff,         1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(rho_mat,       1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(I_eV,          1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(X0_gcm2,       1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(b_rad,         1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(infile,        200, MPI_CHARACTER,        0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(outfile,       200, MPI_CHARACTER,        0, MPI_COMM_WORLD, ierr)

  ! Sternheimer C — same formula on all ranks (inputs already broadcast)
  hbar_omega_p = 28.816d-9 * sqrt(rho_mat * Z_eff / A_eff)
  if (hbar_omega_p > 0.d0 .and. I_eV > 0.d0) then
    C_dens = 2.d0*log(hbar_omega_p / (I_eV*1.d-9)) - 1.d0
  else
    C_dens = 0.d0
  end if

  if (my_rank == 0) then
    write(*,*)
    write(*,*) ' ============================================================'
    write(*,*) '        UCMuon Transport BB — CONFIGURATION'
    write(*,*) ' ============================================================'
    write(*,'(A,A)')       '  Input:             ', trim(infile)
    write(*,'(A,A)')       '  Output prefix:     ', trim(outfile)
    write(*,'(A,F7.3,A,F6.3)') '  Z / A:             ', Z_eff, ' / ', A_eff
    write(*,'(A,F7.3,A)')  '  Density:           ', rho_mat, ' g/cm3'
    write(*,'(A,F7.2,A)')  '  Rad. length X0:    ', X0_gcm2, ' g/cm2'
    write(*,'(A,E10.3)')   '  b_rad:             ', b_rad
    write(*,'(A,F7.2,A)')  '  I_mean:            ', I_eV, ' eV'
    write(*,'(A,F8.2,A)')  '  Vertical depth:    ', depth_m, ' m'
    write(*,'(A,F10.1,A)') '  Depth (w.e.):      ', depth_cm*rho_mat, ' g/cm2'
    write(*,'(A,I2,A,I6,A,I4)') '  MS: ', ms_enable, &
      '   MPI ranks: ', nranks, '  OMP threads/rank: ', omp_get_max_threads()
    write(*,*) ' ============================================================'
    write(*,*)
  end if

  !===========================================================================
  ! READ INPUT — rank 0, two-pass
  !===========================================================================
  nmuon_total = 0; ncols = 0; nskip_mu = 0

  if (my_rank == 0) then
    write(*,*) ' Reading input file...'
    open(unit=10, file=trim(infile), form='formatted', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') ' ERROR: cannot open: ', trim(infile)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    ! Format detection
    block
      integer  :: hf_t, dm_t, ev_t, ch_t, ios3
      real(8)  :: p_t,px_t,py_t,pz_t,th_t,ph_t,e_t,xs_t,ys_t,zs_t
      do
        read(10,'(A200)', iostat=ios) linebuf
        if (ios /= 0) exit
        if (linebuf(1:1)=='#' .or. len_trim(linebuf)==0) cycle
        read(linebuf,*,iostat=ios3) ev_t,xs_t,ys_t,zs_t,p_t, &
              px_t,py_t,pz_t,th_t,ph_t,e_t,ch_t,hf_t,dm_t
        if (ios3==0) then
          ncols = 14
        else
          read(linebuf,*,iostat=ios3) ev_t,xs_t,ys_t,zs_t,p_t, &
                px_t,py_t,pz_t,th_t,ph_t,e_t,ch_t,dm_t
          if (ios3==0) then; ncols = 13
          else
            write(*,*) ' ERROR: cannot parse first data line.'
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
          end if
        end if
        exit
      end do
      rewind(10)
    end block

    write(*,'(A,I2,A)') '  Format: ', ncols, ' columns'
    if (ncols == 13) transport_all = 1

    ! Pass 1: count
    block
      integer :: hf_c, dm_c, ev_c, ch_c, ios3
      real(8) :: v1,v2,v3,v4,v5,v6,v7,v8,v9,v10
      do
        read(10,'(A200)', iostat=ios) linebuf
        if (ios /= 0) exit
        if (linebuf(1:1)=='#' .or. len_trim(linebuf)==0) cycle
        if (ncols == 14) then
          read(linebuf,*,iostat=ios3) ev_c,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,ch_c,hf_c,dm_c
          if (ios3 /= 0) cycle
          if (hf_c /= 1 .and. transport_all == 0) then; nskip_mu=nskip_mu+1; cycle; end if
        else
          read(linebuf,*,iostat=ios3) ev_c,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,ch_c,dm_c
          if (ios3 /= 0) cycle
        end if
        nmuon_total = nmuon_total + 1
      end do
    end block
    rewind(10)

    write(*,'(A,I12)') '  Muons to transport: ', nmuon_total
    if (nmuon_total == 0) then
      write(*,*) ' ERROR: no muons in input file.'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    allocate(g_eventid(nmuon_total), g_charge(nmuon_total), g_det_mask(nmuon_total))
    allocate(g_xs(nmuon_total), g_ys(nmuon_total), g_zs(nmuon_total))
    allocate(g_p(nmuon_total), g_px(nmuon_total), g_py(nmuon_total), g_pz(nmuon_total))
    allocate(g_theta(nmuon_total), g_phi(nmuon_total), g_e(nmuon_total))

    ! Pass 2: load
    block
      integer :: hf_r, dm_r, ev_r, ch_r, ios3
      real(8) :: xs_r,ys_r,zs_r,p_r,px_r,py_r,pz_r,th_r,ph_r,e_r
      imu = 0
      do
        read(10,'(A200)', iostat=ios) linebuf
        if (ios /= 0) exit
        if (linebuf(1:1)=='#' .or. len_trim(linebuf)==0) cycle
        if (ncols == 14) then
          read(linebuf,*,iostat=ios3) ev_r,xs_r,ys_r,zs_r,p_r, &
                px_r,py_r,pz_r,th_r,ph_r,e_r,ch_r,hf_r,dm_r
          if (ios3 /= 0) cycle
          if (hf_r /= 1 .and. transport_all == 0) cycle
        else
          read(linebuf,*,iostat=ios3) ev_r,xs_r,ys_r,zs_r,p_r, &
                px_r,py_r,pz_r,th_r,ph_r,e_r,ch_r,dm_r
          if (ios3 /= 0) cycle; hf_r = 1
        end if
        imu = imu + 1
        g_eventid(imu)=ev_r; g_xs(imu)=xs_r; g_ys(imu)=ys_r; g_zs(imu)=zs_r
        g_p(imu)=p_r; g_px(imu)=px_r; g_py(imu)=py_r; g_pz(imu)=pz_r
        g_theta(imu)=th_r; g_phi(imu)=ph_r; g_e(imu)=e_r
        g_charge(imu)=ch_r; g_det_mask(imu)=dm_r
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
    write(*,*)
  end if  ! my_rank == 0

  call MPI_Bcast(nmuon_total,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ncols,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(transport_all,1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(depth_axis,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  !===========================================================================
  ! WORK SPLIT
  !===========================================================================
  block
    integer :: base, rem, r
    base = nmuon_total / nranks
    rem  = mod(nmuon_total, nranks)
    sendcounts(0) = base + rem
    do r = 1, nranks-1; sendcounts(r) = base; end do
    displs(0) = 0
    do r = 1, nranks-1; displs(r) = displs(r-1) + sendcounts(r-1); end do
    local_nmuon = sendcounts(my_rank)
  end block

  !===========================================================================
  ! ALLOCATE LOCAL + OUTPUT ARRAYS
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

  if (my_rank /= 0) then
    allocate(g_eventid(1),g_charge(1),g_det_mask(1))
    allocate(g_xs(1),g_ys(1),g_zs(1))
    allocate(g_p(1),g_px(1),g_py(1),g_pz(1))
    allocate(g_theta(1),g_phi(1),g_e(1))
  end if

  !===========================================================================
  ! SCATTER
  !===========================================================================
  call MPI_Scatterv(g_xs,       sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_xs,       local_nmuon, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_ys,       sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_ys,       local_nmuon, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_zs,       sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_zs,       local_nmuon, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_p,        sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_p,        local_nmuon, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_px,       sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_px,       local_nmuon, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_py,       sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_py,       local_nmuon, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_pz,       sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_pz,       local_nmuon, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_theta,    sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_theta,    local_nmuon, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_phi,      sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_phi,      local_nmuon, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_e,        sendcounts, displs, MPI_DOUBLE_PRECISION, &
                    l_e,        local_nmuon, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_eventid,  sendcounts, displs, MPI_INTEGER, &
                    l_eventid,  local_nmuon, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_charge,   sendcounts, displs, MPI_INTEGER, &
                    l_charge,   local_nmuon, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Scatterv(g_det_mask, sendcounts, displs, MPI_INTEGER, &
                    l_det_mask, local_nmuon, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  deallocate(g_eventid,g_charge,g_det_mask,g_xs,g_ys,g_zs)
  deallocate(g_p,g_px,g_py,g_pz,g_theta,g_phi,g_e)

  !===========================================================================
  ! RNG INIT
  !===========================================================================
  iranlux_base = 271828 + my_rank * 99991
  nthreads     = omp_get_max_threads()
  !$OMP PARALLEL PRIVATE(tid)
  tid = omp_get_thread_num()
  call rluxgo(3, iranlux_base + tid * 100003, 0, 0)
  !$OMP END PARALLEL

  write(*,'(A,I5,A,I4,A)') '  [Rank ',my_rank,'] RNG ready for ',nthreads,' thread(s).'
  flush(6)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  !===========================================================================
  ! TRANSPORT LOOP
  !===========================================================================
  if (my_rank == 0) then
    write(*,*)
    write(*,*) ' Transporting muons (BB — MPI x OMP)...'
    write(*,*)
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  t_start = MPI_Wtime()

  nsurvive = 0; nstop = 0

  !$OMP PARALLEL DO                                        &
  !$OMP   DEFAULT(SHARED)                                  &
  !$OMP   PRIVATE(i,x,y,z,cx,cy,cz,emu,slant_cm,eff_cz,  &
  !$OMP           alive,theta_ug,phi_ug)                   &
  !$OMP   REDUCTION(+:nsurvive,nstop)                      &
  !$OMP   SCHEDULE(DYNAMIC,1)
  do i = 1, local_nmuon
    if (l_p(i) > 0.d0) then
      cx = l_px(i)/l_p(i); cy = l_py(i)/l_p(i); cz = l_pz(i)/l_p(i)
    else
      cx =  sin(l_theta(i))*cos(l_phi(i))
      cy =  sin(l_theta(i))*sin(l_phi(i))
      cz = -cos(l_theta(i))
    end if
    cz = -cz                    ! UCMuon cz<0=down → BB cz>0=down
    x =  l_xs(i); y =  l_ys(i); z = -l_zs(i)
    emu = l_e(i)

    select case(depth_axis)
      case(1);  eff_cz = -cy   ! XZ plane: depth in Y, cy<0 → eff_cz>0
      case(0);  eff_cz = -cx   ! YZ plane: depth in X
      case default; eff_cz = cz ! XY plane: cz already flipped positive
    end select
    if (eff_cz > 1.d-6) then
      slant_cm = min(depth_cm/eff_cz, 1.0d7)
    else
      slant_cm = 1.0d7
    end if

    call transport_bb(x, y, z, cx, cy, cz, emu, slant_cm, &
                      rho_mat, I_eV, Z_eff, A_eff, X0_gcm2, b_rad, C_dens, ms_enable)

    if (emu > EMIN) then
      alive=1; out_e_ug(i)=emu; nsurvive=nsurvive+1
    else
      alive=0; out_e_ug(i)=0.d0; nstop=nstop+1
    end if

    out_alive(i)=alive
    out_x_ug(i)=x; out_y_ug(i)=y; out_z_ug(i)=-z
    out_cx_ug(i)=cx; out_cy_ug(i)=cy; out_cz_ug(i)=-cz

    if (alive == 1) then
      theta_ug = acos(min(1.d0,max(-1.d0,-out_cz_ug(i))))
      if (abs(sin(theta_ug)) > 1.d-9) then
        phi_ug = atan2(out_cy_ug(i),out_cx_ug(i))
        if (phi_ug < 0.d0) phi_ug = phi_ug + TWO_PI
      else
        phi_ug = 0.d0
      end if
    else
      theta_ug = l_theta(i); phi_ug = l_phi(i)
    end if
    out_theta_ug(i) = theta_ug; out_phi_ug(i) = phi_ug

    if (mod(i, max(1,local_nmuon/20)) == 0) then
      write(*,'(A,I5,A,I8,A,I8,A,F6.1,A)') &
        '  [Rank ',my_rank,'] Transported: ',i, &
        '  Survived: ',nsurvive,'  (',100.d0*dble(i)/dble(local_nmuon),'%)'
      flush(6)
    end if
  end do
  !$OMP END PARALLEL DO

  t_end = MPI_Wtime()

  !===========================================================================
  ! OUTPUT — 18-column format identical to MUSIC driver
  !===========================================================================
  write(outfile_rank,'(A,A,I5.5,A)') trim(outfile),'_',my_rank,'.dat'
  open(unit=20, file=trim(outfile_rank), form='formatted', status='unknown')
  if (my_rank == 0) &
    write(20,'(A)') '# EventID  x_srf_cm  y_srf_cm  z_srf_cm' // &
      '  E_srf_GeV  theta_srf_rad  phi_srf_rad  charge' // &
      '  alive  x_ug_cm  y_ug_cm  z_ug_cm  E_ug_GeV' // &
      '  cx_ug  cy_ug  cz_ug  theta_ug_rad  phi_ug_rad'

  do i = 1, local_nmuon
    write(20,'(I10,1X,F13.4,1X,F13.4,1X,F13.4,1X,' // &
              'F13.6,1X,F13.9,1X,F13.9,1X,I4,1X,' // &
              'I2,1X,F13.4,1X,F13.4,1X,F13.4,1X,' // &
              'F13.6,1X,F10.6,1X,F10.6,1X,F10.6,1X,' // &
              'F13.9,1X,F13.9)') &
      l_eventid(i), l_xs(i),l_ys(i),l_zs(i),l_e(i),l_theta(i),l_phi(i),l_charge(i), &
      out_alive(i), out_x_ug(i),out_y_ug(i),out_z_ug(i),out_e_ug(i), &
      out_cx_ug(i),out_cy_ug(i),out_cz_ug(i), &
      out_theta_ug(i),out_phi_ug(i)
  end do
  close(20)
  write(*,'(A,I5,A,A)') '  [Rank ',my_rank,'] Output: ',trim(outfile_rank)
  flush(6)

  !===========================================================================
  ! STATISTICS
  !===========================================================================
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Reduce(nsurvive,nsurvive_all,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_Reduce(nstop,   nstop_all,   1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  if (my_rank == 0) then
    write(*,*)
    write(*,*) ' ============================================================'
    write(*,*) '   UCMuon Transport — Bethe-Bloch — COMPLETE'
    write(*,*) ' ============================================================'
    write(*,'(A,I12)')    '  Total transported:  ', nmuon_total
    write(*,'(A,I12)')    '  Survived:           ', nsurvive_all
    write(*,'(A,I12)')    '  Stopped:            ', nstop_all
    if (nmuon_total > 0) &
      write(*,'(A,F8.4,A)') '  Survival rate:      ', &
        100.d0*dble(nsurvive_all)/dble(nmuon_total),' %'
    write(*,'(A,F8.2,A,F8.1,A)') '  Depth (vertical):   ',depth_m,' m  (', &
        depth_cm*rho_mat,' g/cm2)'
    write(*,'(A,I6,A,I4)') '  Parallelisation:    ',nranks, &
        ' MPI ranks x ',omp_get_max_threads(),' OMP threads'
    write(*,'(A,F10.2,A)') '  Wall time:          ',t_end-t_start,' s'
    write(*,*) ' ============================================================'
  end if

  call MPI_Finalize(ierr)
  stop

contains

  !=============================================================================
  ! transport_bb — CSDA Bethe-Bloch + Highland MS
  ! Coordinates: z>0 = underground (MUSIC convention on entry and exit).
  ! Thread-safe: all state on stack + THREADPRIVATE ranlux from ranlux_omp.f.
  !=============================================================================
  subroutine transport_bb(x, y, z, cx, cy, cz, E, slant_cm, &
                           rho, I_eV_in, Zat, Aat, X0, b, C_d, ms)
    implicit none
    real(8), intent(inout) :: x, y, z, cx, cy, cz, E
    real(8), intent(in)    :: slant_cm, rho, I_eV_in, Zat, Aat, X0, b, C_d
    integer, intent(in)    :: ms

    real(8), parameter :: MMUON_ = 0.105658370d0
    real(8), parameter :: ME_    = 5.10998918d-4
    real(8), parameter :: K_     = 3.07075d-4
    real(8), parameter :: TWO_P_ = 6.283185307179586d0
    real(8), parameter :: EMIN_  = 0.105658370d0
    real(8), parameter :: STEP_  = 10.0d0

    integer :: nsteps, istep
    real(8) :: slant_gcm2, ds_gcm2, ds_cm
    real(8) :: gamma_mu, beta2, p_GeV, bg, Tmax, I_GeV
    real(8) :: log_arg, x_dens, delta, dEdX
    real(8) :: t_X0, theta0
    real(8) :: e1x, e1y, e1z, e2x, e2y, e2z, norm
    real(8) :: theta_s, phi_s, ct, st, cp, sp
    real(4) :: rr(4)
    external :: ranlux

    I_GeV      = I_eV_in * 1.d-9
    slant_gcm2 = rho * slant_cm
    nsteps     = max(20, int(slant_gcm2 / STEP_) + 1)
    ds_gcm2    = slant_gcm2 / dble(nsteps)
    ds_cm      = ds_gcm2 / rho

    do istep = 1, nsteps
      if (E <= EMIN_) return

      gamma_mu = E / MMUON_
      beta2    = max(0.d0, 1.d0 - (MMUON_/E)**2)
      p_GeV    = sqrt(max(0.d0, E*E - MMUON_*MMUON_))
      bg       = p_GeV / MMUON_

      Tmax = 2.d0*ME_*beta2*gamma_mu*gamma_mu / &
             (1.d0 + 2.d0*gamma_mu*ME_/MMUON_ + (ME_/MMUON_)**2)

      x_dens = log10(max(bg, 1.d-30))
      delta  = merge(max(0.d0, 2.d0*log(10.d0)*x_dens + C_d), 0.d0, x_dens > 1.d0)

      if (Tmax > 0.d0 .and. beta2 > 1.d-12) then
        log_arg = 2.d0*ME_*beta2*gamma_mu*gamma_mu*Tmax / (I_GeV*I_GeV)
        if (log_arg > 1.d0) then
          dEdX = K_*(Zat/Aat)/beta2*(0.5d0*log(log_arg) - beta2 - 0.5d0*delta)
        else
          dEdX = K_*(Zat/Aat)/beta2
        end if
      else
        dEdX = K_*(Zat/Aat)/max(beta2,1.d-12)
      end if
      ! b is the 100 GeV value; b_rad_shape carries the energy dependence
      dEdX = max(dEdX, 1.d-6) + b*b_rad_shape(E)*E

      E = E - dEdX*ds_gcm2
      if (E <= EMIN_) return

      x = x + cx*ds_cm; y = y + cy*ds_cm; z = z + cz*ds_cm

      if (ms == 1 .and. X0 > 0.d0) then
        t_X0 = ds_gcm2 / X0
        if (t_X0 > 1.d-14) then
          p_GeV  = sqrt(max(0.d0, E*E - MMUON_*MMUON_))
          beta2  = max(1.d-12, 1.d0 - (MMUON_/E)**2)
          theta0 = 13.6d-3/(sqrt(beta2)*p_GeV) * sqrt(t_X0) * &
                   (1.d0 + 0.038d0*log(t_X0))
          if (theta0 <= 0.d0) cycle

          ! Rayleigh(theta0) polar deflection: Highland theta0 is the RMS
          ! projected angle; a Gaussian polar angle under-scatters by sqrt(2).
          call ranlux(rr, 4)
          rr(1) = max(rr(1), 1.e-30)
          theta_s = theta0*sqrt(-2.d0*log(dble(rr(1))))
          phi_s   = TWO_P_*dble(rr(3))

          if (abs(cx) <= abs(cy) .and. abs(cx) <= abs(cz)) then
            e1x=0.d0; e1y=-cz; e1z= cy
          else if (abs(cy) <= abs(cz)) then
            e1x= cz; e1y=0.d0; e1z=-cx
          else
            e1x=-cy; e1y= cx;  e1z=0.d0
          end if
          norm = sqrt(e1x*e1x+e1y*e1y+e1z*e1z)
          if (norm < 1.d-14) cycle
          e1x=e1x/norm; e1y=e1y/norm; e1z=e1z/norm
          e2x=cy*e1z-cz*e1y; e2y=cz*e1x-cx*e1z; e2z=cx*e1y-cy*e1x

          ct=cos(theta_s); st=sin(theta_s); cp=cos(phi_s); sp=sin(phi_s)
          cx=ct*cx+st*(cp*e1x+sp*e2x)
          cy=ct*cy+st*(cp*e1y+sp*e2y)
          cz=ct*cz+st*(cp*e1z+sp*e2z)
          norm=sqrt(cx*cx+cy*cy+cz*cz)
          if (norm > 1.d-14) then; cx=cx/norm; cy=cy/norm; cz=cz/norm; end if
        end if
      end if
    end do
  end subroutine transport_bb

  !===========================================================================
  ! b_rad_shape — PDG-2024 rock radiative-b energy dependence, normalised
  ! to 1 at E_tot = 100 GeV (mirror of ucmuon_transport_bb_omp.f90).
  !===========================================================================
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

end program ucmuon_transport_bb

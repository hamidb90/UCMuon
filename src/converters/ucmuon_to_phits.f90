!=============================================================================
! ucmuon_to_phits.f90  —  UCMuon → PHITS s-type=17 dump converter
! UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>
!
! Fast compiled converter — handles GB-scale files efficiently.
! Reads from stdin, writes to stdout (pipe-friendly).
!
! Usage
! -----
!   ./ucmuon_to_phits gen       < ucmuon_selected.dat   > output_phits.dat
!   ./ucmuon_to_phits transport < ucmuon_underground.dat > output_phits.dat
!
! Or via the SLURM scripts (called automatically after each job).
!
! Input formats
! -------------
!   gen mode       : 13-col or 14-col generator output (auto-detected)
!   transport mode : 18-col transport output (alive==1 rows only)
!
! Output format  (PHITS s-type=17 ASCII dump)
! --------------------------------------------
!   10 columns per muon, D-exponent notation (Fortran 1pd24.15):
!     kf  x[cm]  y[cm]  z[cm]  u  v  w  Ekin[MeV]  wt  time[ns]
!
!   kf codes (PDG): mu- = 13,  mu+ = -13
!
!   Use in PHITS input:
!     [Source]
!     s-type = 17
!     file   = ucmuon_selected_phits.dat
!     dump   = -10
!     1 2 3 4 5 6 7 8 9 10
!
! Performance
! -----------
!   Fortran formatted I/O on a 1GB file: ~15-30 s (vs ~100 s in Python).
!   No MPI or OMP needed — purely serial post-processing.
!=============================================================================
program ucmuon_to_phits
  implicit none

  real(8), parameter :: MMUON_GEV = 0.105658370d0
  real(8), parameter :: GEV_TO_MEV = 1000.0d0

  character(16)  :: mode_arg
  character(512) :: linebuf
  integer :: ios, ncols, nwritten, nskipped

  ! Generator columns
  integer  :: id_g, charge_g, hit_flag_g, det_mask_g
  real(8)  :: x_g, y_g, z_g, p_g, px_g, py_g, pz_g
  real(8)  :: theta_g, phi_g, e_g

  ! Transport columns
  integer  :: id_t, charge_t, alive_t
  real(8)  :: x_srf, y_srf, z_srf, e_srf, theta_srf, phi_srf
  real(8)  :: x_ug, y_ug, z_ug, e_ug, cx_ug, cy_ug, cz_ug
  real(8)  :: theta_ug, phi_ug

  ! PHITS output fields
  integer  :: kf
  real(8)  :: x, y, z, u, v, w, ekin_mev, norm

  !--- Read mode argument -------------------------------------------------------
  if (command_argument_count() < 1) then
    write(0,*) 'Usage: ucmuon_to_phits gen|transport < input.dat > output.dat'
    stop 1
  end if
  call get_command_argument(1, mode_arg)
  mode_arg = adjustl(mode_arg)

  if (trim(mode_arg) /= 'gen' .and. trim(mode_arg) /= 'transport') then
    write(0,*) 'ERROR: mode must be "gen" or "transport", got: ', trim(mode_arg)
    stop 1
  end if

  nwritten = 0
  nskipped = 0
  ncols    = 0

  !===========================================================================
  ! GENERATOR MODE  (13 or 14 columns)
  !===========================================================================
  if (trim(mode_arg) == 'gen') then

    do
      read(*, '(A512)', iostat=ios) linebuf
      if (ios /= 0) exit
      linebuf = adjustl(linebuf)
      if (len_trim(linebuf) == 0) cycle
      if (linebuf(1:1) == '#')   cycle

      ! Auto-detect column count from first data line
      if (ncols == 0) then
        read(linebuf, *, iostat=ios) id_g, x_g, y_g, z_g, p_g, &
              px_g, py_g, pz_g, theta_g, phi_g, e_g, charge_g, hit_flag_g, det_mask_g
        if (ios == 0) then
          ncols = 14
        else
          read(linebuf, *, iostat=ios) id_g, x_g, y_g, z_g, p_g, &
                px_g, py_g, pz_g, theta_g, phi_g, e_g, charge_g, det_mask_g
          if (ios == 0) then
            ncols = 13
          else
            write(0,*) 'ERROR: cannot parse first data line as 13 or 14 columns.'
            stop 1
          end if
        end if
      else
        if (ncols == 14) then
          read(linebuf, *, iostat=ios) id_g, x_g, y_g, z_g, p_g, &
                px_g, py_g, pz_g, theta_g, phi_g, e_g, charge_g, hit_flag_g, det_mask_g
        else
          read(linebuf, *, iostat=ios) id_g, x_g, y_g, z_g, p_g, &
                px_g, py_g, pz_g, theta_g, phi_g, e_g, charge_g, det_mask_g
        end if
        if (ios /= 0) then
          nskipped = nskipped + 1
          cycle
        end if
      end if

      ! Direction cosines from momentum components
      if (p_g > 0.0d0) then
        u = px_g / p_g
        v = py_g / p_g
        w = pz_g / p_g
      else
        u = 0.0d0;  v = 0.0d0;  w = -1.0d0
      end if

      ! Normalise (safety)
      norm = sqrt(u*u + v*v + w*w)
      if (norm > 1.0d-14) then
        u = u/norm;  v = v/norm;  w = w/norm
      end if

      ! PDG kf: mu+ = -13, mu- = +13
      kf = merge(-13, 13, charge_g == 1)

      ! Kinetic energy in MeV
      ekin_mev = max(0.0d0, (e_g - MMUON_GEV) * GEV_TO_MEV)

      ! Write PHITS record
      write(*, '(10(1pd24.15))') dble(kf), x_g, y_g, z_g, &
                                  u, v, w, ekin_mev, 1.0d0, 0.0d0
      nwritten = nwritten + 1
    end do

  !===========================================================================
  ! TRANSPORT MODE  (18 columns, alive==1 only)
  !===========================================================================
  else

    do
      read(*, '(A512)', iostat=ios) linebuf
      if (ios /= 0) exit
      linebuf = adjustl(linebuf)
      if (len_trim(linebuf) == 0) cycle
      if (linebuf(1:1) == '#')   cycle

      read(linebuf, *, iostat=ios) &
        id_t, x_srf, y_srf, z_srf, e_srf, theta_srf, phi_srf, charge_t, &
        alive_t, x_ug, y_ug, z_ug, e_ug, cx_ug, cy_ug, cz_ug, theta_ug, phi_ug

      if (ios /= 0) then
        nskipped = nskipped + 1
        cycle
      end if

      ! Skip stopped muons
      if (alive_t /= 1) then
        nskipped = nskipped + 1
        cycle
      end if

      u = cx_ug;  v = cy_ug;  w = cz_ug

      ! Normalise direction
      norm = sqrt(u*u + v*v + w*w)
      if (norm > 1.0d-14) then
        u = u/norm;  v = v/norm;  w = w/norm
      end if

      kf = merge(-13, 13, charge_t == 1)
      ekin_mev = max(0.0d0, (e_ug - MMUON_GEV) * GEV_TO_MEV)

      write(*, '(10(1pd24.15))') dble(kf), x_ug, y_ug, z_ug, &
                                  u, v, w, ekin_mev, 1.0d0, 0.0d0
      nwritten = nwritten + 1
    end do

  end if

  !===========================================================================
  ! SUMMARY  (to stderr so it doesn't pollute the PHITS file)
  !===========================================================================
  write(0, '(A,I12)') '  ucmuon_to_phits: written  ', nwritten
  if (nskipped > 0) &
    write(0, '(A,I12)') '  ucmuon_to_phits: skipped  ', nskipped
  write(0, '(A)')      '  Use in PHITS: s-type=17, dump=-10, 1 2 3 4 5 6 7 8 9 10'

  stop
end program ucmuon_to_phits

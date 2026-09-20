!=============================================================================
! parma_vertical_spec.f90 — helper for make_fig02_flux_spectra.py
!
! Prints the PARMA/EXPACS vertical differential muon intensity
!   I(T, 0 deg) = [Phi_mu+(E) + Phi_mu-(E)] * F_ang(E, cos=1)
! in cm^-2 s^-1 sr^-1 GeV^-1 on a log-spaced kinetic-energy grid.
! Mirrors the call sequence of src/generator/ucmuon_gen_omp.f90 (mode 3).
!
! Build (from repo root):
!   gfortran -O2 -o manuscript/scripts/_parma_spec \
!     src/parma/parma_path_module.f90 src/parma/parma_subroutines.f90 \
!     manuscript/scripts/parma_vertical_spec.f90
!
! stdin (one value per line):
!   datapath  lat[deg]  lon[deg]  alt[km]  year  month  day  W  Tmin[GeV]  Tmax[GeV]  npts
! stdout: two columns  T[GeV]  I[cm^-2 s^-1 sr^-1 GeV^-1]
!=============================================================================
program parma_vertical_spec
  use parma_path
  implicit none
  real(8), external :: getd, getr, getHP, getMuonSpec, getSpecAngFinal
  character(200) :: datapath
  real(8) :: lat, lon, alt_km, s_W, t_min, t_max
  integer :: year, month, day, npts, j
  integer :: ic   ! getHP error-code output — must be a variable
  real(8) :: d_gcm2, rc_GV, ffp_MV, T_GeV, E_MeV, phi_p, phi_m, ang, intens

  read(*,'(A)') datapath
  read(*,*) lat
  read(*,*) lon
  read(*,*) alt_km
  read(*,*) year
  read(*,*) month
  read(*,*) day
  read(*,*) s_W
  read(*,*) t_min
  read(*,*) t_max
  read(*,*) npts

  call parma_set_datadir(trim(datapath))
  ic = 0
  d_gcm2 = getd(alt_km, lat)
  rc_GV  = getr(lat, lon)
  ffp_MV = getHP(year, month, day, ic)
  s_W    = max(s_W, -135.4d0)

  write(0,'(A,F8.2,A,F8.3,A,F8.1)') '# depth=', d_gcm2, ' g/cm2  rc=', &
        rc_GV, ' GV  FFP=', ffp_MV

  do j = 1, npts
    T_GeV = t_min * exp( dble(j-1)/dble(npts-1) * log(t_max/t_min) )
    E_MeV = T_GeV * 1.0d3
    phi_p = max(0.0d0, getMuonSpec(1, s_W, rc_GV, d_gcm2, E_MeV))
    phi_m = max(0.0d0, getMuonSpec(2, s_W, rc_GV, d_gcm2, E_MeV))
    ang   = max(0.0d0, getSpecAngFinal(4, s_W, rc_GV, d_gcm2, E_MeV, &
                                       0.0d0, 1.0d0))
    ! (phi_p+phi_m) is /cm2/s/MeV; ang is /sr; *1e3 -> per GeV
    intens = (phi_p + phi_m) * ang * 1.0d3
    write(*,'(ES14.6, 2X, ES14.6)') T_GeV, intens
  end do
end program parma_vertical_spec

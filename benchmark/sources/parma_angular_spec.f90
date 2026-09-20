!=============================================================================
! parma_angular_spec.f90 — helper for make_parma_spectrum.py
!
! Prints the PARMA/EXPACS muon angular distribution factor
!   F_ang(E, cos theta)   [1/sr]
! at a fixed kinetic energy, on a zenith-angle grid theta = 0..theta_max.
! This is the angular factor used by the surface generator (mode 3); the
! ratio F_ang(E,theta)/F_ang(E,0) is the zenith dependence I(p,theta)/I(p,0)
! plotted in fig_zenith_dependence.  Muon = particle id 4 in PARMA.
!
! stdin (one value per line):
!   datapath  lat[deg]  lon[deg]  alt[km]  year  month  day  W  T[GeV]
!   theta_max[deg]  ntheta
! stdout: two columns  theta[deg]  F_ang[1/sr]
!=============================================================================
program parma_angular_spec
  use parma_path
  implicit none
  real(8), external :: getd, getr, getHP, getSpecAngFinal
  character(200) :: datapath
  real(8) :: lat, lon, alt_km, s_W, t_gev, theta_max
  integer :: year, month, day, nth, j
  integer :: ic
  real(8) :: d_gcm2, rc_GV, ffp_MV, e_mev, theta, cth, ang
  real(8), parameter :: PI = 3.141592653589793d0

  read(*,'(A)') datapath
  read(*,*) lat
  read(*,*) lon
  read(*,*) alt_km
  read(*,*) year
  read(*,*) month
  read(*,*) day
  read(*,*) s_W
  read(*,*) t_gev
  read(*,*) theta_max
  read(*,*) nth

  call parma_set_datadir(trim(datapath))
  ic = 0
  d_gcm2 = getd(alt_km, lat)
  rc_GV  = getr(lat, lon)
  ffp_MV = getHP(year, month, day, ic)
  s_W    = max(s_W, -135.4d0)
  e_mev  = t_gev * 1.0d3

  do j = 1, nth
    theta = (dble(j-1)/dble(nth-1)) * theta_max
    cth   = cos(theta * PI / 180.0d0)
    ang   = max(0.0d0, getSpecAngFinal(4, s_W, rc_GV, d_gcm2, e_mev, &
                                       0.0d0, cth))
    write(*,'(F10.5, 2X, ES14.6)') theta, ang
  end do
end program parma_angular_spec

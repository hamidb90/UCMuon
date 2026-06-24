!=============================================================================
! ucmuon_source_module.f90
! CosmoALEPH muon flux parametrization  — OpenMP version
!
! Changes from v10:
!   - uses rng_parallel (par_ranlux) instead of legacy RANLUX
!   - all  call RANLUX(yfl, 1)  →  call par_ranlux(yfl)
!   - removed:  external :: RANLUX
!   - added:    use rng_parallel
!
! Spectrum options:
!   1 = CosmoALEPH:  dN/dp = 10^3.8467 * p^(-3.1952)  [default]
!   2 = Power-law:   dN/dE ∝ E^(-3.7), E = E_min * u^(-1/2.7)
!                    (Kudryavtsev/MUSIC convention, unbounded above E_min)
!   3 = PARMA/EXPACS (handled externally in cosmoaleph_main_omp.f90)
!   4 = Guan et al. (2015), arXiv:1509.06176
!   5 = Frosin et al. (2025), J. Phys. G 52, 035002
!   6 = Bugaev/Gaisser (1990): Gaisser pion+kaon formula, no atm. correction
!   7 = Reyna–Bugaev (2006): log-polynomial p^3*F_vert, arXiv:hep-ph/0604145
!   8 = Cosmic electrons: dN/dE ∝ E^−3.0, sea level, 10 MeV–1 GeV
!=============================================================================
module ucmuon_source_module
  use rng_parallel          ! <-- OMP change: replaces  external :: RANLUX
  implicit none
  private


  integer,  parameter :: NGRID     = 300
  real(8),  parameter :: PI        = 3.141592654d0
  real(8),  parameter :: MUON_MASS = 0.10566d0
  real(8),  parameter :: A_COSMO   = 3.8467d0
  real(8),  parameter :: B_COSMO   = -3.1952d0

  real(8), parameter :: GUAN_P1 =  0.102573d0
  real(8), parameter :: GUAN_P2 = -0.068287d0
  real(8), parameter :: GUAN_P3 =  0.958633d0
  real(8), parameter :: GUAN_P4 =  0.0407253d0
  real(8), parameter :: GUAN_P5 =  0.817285d0
  real(8), parameter :: GUAN_DENOM = 0.99144315d0
  real(8), parameter :: GUAN_EPI = 115.0d0
  real(8), parameter :: GUAN_EK  = 850.0d0
  real(8), parameter :: GUAN_KF  = 0.054d0
  real(8), parameter :: GUAN_PRE = 0.14d0
  real(8), parameter :: GUAN_IDX = -2.7d0
  real(8), parameter :: GUAN_A   = 3.64d0
  real(8), parameter :: GUAN_B   = 1.29d0
  real(8), parameter :: FROSIN_A = 3.512d0
  real(8), parameter :: FROSIN_B = 1.388d0
  ! Mode 6: Bugaev/Gaisser — standard Gaisser formula, no atmospheric correction
  real(8), parameter :: BUGAEV_A = 0.0d0
  real(8), parameter :: BUGAEV_B = 1.0d0
  ! Mode 7: Reyna–Bugaev (arXiv:hep-ph/0604145, Eq. 6-7)
  !   I_V(p) = C0 * p^-(C1 + C2*z + C3*z^2 + C4*z^3),  z = log10(p)
  real(8), parameter :: REYNA_C0 =  0.00253d0
  real(8), parameter :: REYNA_C1 =  0.2455d0
  real(8), parameter :: REYNA_C2 =  1.288d0
  real(8), parameter :: REYNA_C3 = -0.2555d0
  real(8), parameter :: REYNA_C4 =  0.0209d0
  real(8), parameter :: ELECTRON_MASS = 0.000511d0  ! GeV/c^2
  real(8), parameter :: EL_IDX        = -3.0d0      ! electron spectral index

  real(8), dimension(NGRID) :: p_cdf, cdf
  real(8) :: p_min_stored   = 0.0d0
  real(8) :: p_max_stored   = 0.0d0
  logical :: use_analytical = .false.

  integer :: spectrum_mode = 1

  public :: PI, MUON_MASS, ELECTRON_MASS
  public :: build_cosmoaleph_cdf
  public :: generate_muon
  public :: is_electron_mode


contains


  !==========================================================================
  subroutine build_cosmoaleph_cdf(p_min, p_max, mode)
    real(8), intent(in) :: p_min, p_max
    integer, intent(in) :: mode
    real(8) :: flux(NGRID), partial(NGRID), dp, I_total
    real(8) :: E_j
    integer :: j

    p_min_stored  = p_min
    p_max_stored  = p_max
    spectrum_mode = mode

    write(*,*) 'Building spectrum CDF...'

    ! Mono-energetic: all particles at exactly p_min — skip CDF entirely
    if (p_min >= p_max) then
      use_analytical = .false.
      write(*,*) '  Spectrum mode:     Mono-energetic (delta function)'
      write(*,*) '  CDF ready.'
      write(*,*)
      return
    end if

    if (spectrum_mode == 2) then
      use_analytical = .true.
      I_total = p_min**(-2.7d0) / 2.7d0
      write(*,'(A,ES12.4)') '  Integrated flux:   ', I_total
      write(*,*) '  Spectrum mode:     Power-law E^-3.7 (legacy MUSIC cross-check)'
      write(*,*) '  CDF ready.'
      write(*,*)
      return
    end if

    if (spectrum_mode == 4) then
      use_analytical = .false.
      do j = 1, NGRID
        p_cdf(j) = p_min * exp( dble(j-1)/dble(NGRID-1) * log(p_max/p_min) )
      end do
      do j = 1, NGRID
        E_j      = sqrt(p_cdf(j)**2 + MUON_MASS**2)
        flux(j)  = guan_flux(E_j, 1.0d0, GUAN_A, GUAN_B) * p_cdf(j) / E_j
      end do
      partial(1) = 0.0d0
      do j = 2, NGRID
        dp         = p_cdf(j) - p_cdf(j-1)
        partial(j) = partial(j-1) + 0.5d0*(flux(j)+flux(j-1))*dp
      end do
      I_total = partial(NGRID)
      do j = 1, NGRID
        cdf(j) = partial(j) / I_total
      end do
      write(*,'(A,ES12.4)') '  Integrated flux:   ', I_total
      write(*,*) '  Spectrum mode:     Guan et al. (2015)'
      write(*,*) '  CDF ready.'
      write(*,*)
      return
    end if

    if (spectrum_mode == 5) then
      use_analytical = .false.
      do j = 1, NGRID
        p_cdf(j) = p_min * exp( dble(j-1)/dble(NGRID-1) * log(p_max/p_min) )
      end do
      do j = 1, NGRID
        E_j      = sqrt(p_cdf(j)**2 + MUON_MASS**2)
        flux(j)  = guan_flux(E_j, 1.0d0, FROSIN_A, FROSIN_B) * p_cdf(j) / E_j
      end do
      partial(1) = 0.0d0
      do j = 2, NGRID
        dp         = p_cdf(j) - p_cdf(j-1)
        partial(j) = partial(j-1) + 0.5d0*(flux(j)+flux(j-1))*dp
      end do
      I_total = partial(NGRID)
      do j = 1, NGRID
        cdf(j) = partial(j) / I_total
      end do
      write(*,'(A,ES12.4)') '  Integrated flux:   ', I_total
      write(*,*) '  Spectrum mode:     Frosin et al. (2025)'
      write(*,*) '  CDF ready.'
      write(*,*)
      return
    end if

    if (spectrum_mode == 6) then
      use_analytical = .false.
      do j = 1, NGRID
        p_cdf(j) = p_min * exp( dble(j-1)/dble(NGRID-1) * log(p_max/p_min) )
      end do
      do j = 1, NGRID
        E_j     = sqrt(p_cdf(j)**2 + MUON_MASS**2)
        flux(j) = guan_flux(E_j, 1.0d0, BUGAEV_A, BUGAEV_B) * p_cdf(j) / E_j
      end do
      partial(1) = 0.0d0
      do j = 2, NGRID
        dp         = p_cdf(j) - p_cdf(j-1)
        partial(j) = partial(j-1) + 0.5d0*(flux(j)+flux(j-1))*dp
      end do
      I_total = partial(NGRID)
      do j = 1, NGRID
        cdf(j) = partial(j) / I_total
      end do
      write(*,'(A,ES12.4)') '  Integrated flux:   ', I_total
      write(*,*) '  Spectrum mode:     Bugaev/Gaisser (1990)  dN/dE ~ E^-2.7*(pion+kaon)'
      write(*,*) '  CDF ready.'
      write(*,*)
      return
    end if

    if (spectrum_mode == 7) then
      use_analytical = .false.
      do j = 1, NGRID
        p_cdf(j) = p_min * exp( dble(j-1)/dble(NGRID-1) * log(p_max/p_min) )
      end do
      do j = 1, NGRID
        flux(j) = reyna_flux(p_cdf(j), 1.0d0)
      end do
      partial(1) = 0.0d0
      do j = 2, NGRID
        dp         = p_cdf(j) - p_cdf(j-1)
        partial(j) = partial(j-1) + 0.5d0*(flux(j)+flux(j-1))*dp
      end do
      I_total = partial(NGRID)
      do j = 1, NGRID
        cdf(j) = partial(j) / I_total
      end do
      write(*,'(A,ES12.4)') '  Integrated flux:   ', I_total
      write(*,*) '  Spectrum mode:     Reyna-Bugaev (2006)  log-poly p^3*F_vert'
      write(*,*) '  CDF ready.'
      write(*,*)
      return
    end if

    if (spectrum_mode == 8) then
      use_analytical = .true.
      write(*,*) '  Spectrum mode:     Cosmic electrons  dN/dE ~ E^-3.0'
      write(*,*) '  CDF ready.'
      write(*,*)
      return
    end if

    ! Mode 1: CosmoALEPH
    if (p_min < 10.0d0) then
      use_analytical = .true.
      I_total = (10d0**A_COSMO) &
              * (p_max**(B_COSMO+1d0) - p_min**(B_COSMO+1d0)) &
              / (B_COSMO + 1d0)
      write(*,'(A,ES12.4)') '  Integrated flux:   ', I_total
      write(*,*) '  Spectrum mode:     CosmoALEPH (analytical)'
      write(*,*) '  CDF ready.'
      write(*,*)
      return
    end if

    use_analytical = .false.
    do j = 1, NGRID
      p_cdf(j) = p_min * exp( dble(j-1)/dble(NGRID-1) * log(p_max/p_min) )
    end do
    do j = 1, NGRID
      flux(j) = 10d0**A_COSMO * p_cdf(j)**B_COSMO
    end do
    partial(1) = 0d0
    do j = 2, NGRID
      dp = p_cdf(j) - p_cdf(j-1)
      partial(j) = partial(j-1) + 0.5d0*(flux(j)+flux(j-1))*dp
    end do
    I_total = partial(NGRID)
    do j = 1, NGRID
      cdf(j) = partial(j) / I_total
    end do
    write(*,'(A,ES12.4)') '  Integrated flux:   ', I_total
    write(*,*) '  Spectrum mode:     CosmoALEPH (table)'
    write(*,*) '  CDF ready.'
    write(*,*)
  end subroutine build_cosmoaleph_cdf



  !==========================================================================
  ! generate_muon — unchanged interface, OMP-safe via par_ranlux
  !==========================================================================
  subroutine generate_muon(source_mode, radius_cm, half_lx_cm, half_ly_cm, &
                            source_z_cm,                                      &
                            theta_max, angular_mode,                         &
                            x, y, z, emu, cx, cy, cz, muon_charge)
    integer, intent(in)  :: source_mode
    real(8), intent(in)  :: radius_cm
    real(8), intent(in)  :: half_lx_cm
    real(8), intent(in)  :: half_ly_cm
    real(8), intent(in)  :: source_z_cm
    real(8), intent(in)  :: theta_max
    integer, intent(in)  :: angular_mode

    real(8), intent(out) :: x, y, z, emu, cx, cy, cz
    integer, intent(out) :: muon_charge

    real(8) :: p, u, r, tc, theta, phi, s_th
    real(8) :: charge_ratio, pos_frac
    real(4) :: yfl      ! real(4) buffer — par_ranlux writes real(4)

    !--- Momentum sampling ---
    call par_ranlux(yfl)           ! <-- OMP change
    u = dble(yfl)
    p = sample_momentum(u)

    !--- Charge and total energy (particle-type dependent) ---
    if (spectrum_mode == 8) then
      call par_ranlux(yfl)
      muon_charge = merge(1, -1, dble(yfl) < 0.5d0)
      emu = sqrt(p*p + ELECTRON_MASS*ELECTRON_MASS)
    else
      charge_ratio = charge_ratio_from_p(p)
      pos_frac = charge_ratio / (1d0 + charge_ratio)
      call par_ranlux(yfl)         ! <-- OMP change
      muon_charge = merge(1, -1, dble(yfl) < pos_frac)
      emu = sqrt(p*p + MUON_MASS*MUON_MASS)
    end if

    !--- Position ---
    if (source_mode == 2) then
      call par_ranlux(yfl)         ! <-- OMP change
      x = half_lx_cm * (2d0*dble(yfl) - 1d0)
      call par_ranlux(yfl)         ! <-- OMP change
      y = half_ly_cm * (2d0*dble(yfl) - 1d0)
      z = source_z_cm

    else if (source_mode == 3) then
      call par_ranlux(yfl)         ! <-- OMP change
      phi = 2d0*PI*dble(yfl)
      call par_ranlux(yfl)         ! <-- OMP change
      u   = dble(yfl)
      theta = acos(1d0 - u)
      s_th  = sin(theta)
      x     = radius_cm * s_th * cos(phi)
      y     = radius_cm * s_th * sin(phi)
      z     = radius_cm * cos(theta) + source_z_cm

    else
      call par_ranlux(yfl)         ! <-- OMP change
      tc = 2d0*PI*dble(yfl)
      call par_ranlux(yfl)         ! <-- OMP change
      r  = radius_cm * sqrt(dble(yfl))
      x  = r*cos(tc);  y = r*sin(tc);  z = source_z_cm
    end if

    !--- Direction ---
    select case (angular_mode)
      case (1)
        theta = 0d0
      case (2)
        theta = sample_cos2(theta_max)
      case (3)
        call par_ranlux(yfl)       ! <-- OMP change
        theta = dble(yfl) * theta_max
      case (4)
        theta = sample_guan_angle(emu, theta_max)
      case (5)
        theta = sample_cos3(theta_max)
      case default
        theta = 0d0
    end select

    call par_ranlux(yfl)           ! <-- OMP change
    phi = 2d0*PI*dble(yfl)

    cx =  sin(theta)*cos(phi)
    cy =  sin(theta)*sin(phi)
    cz = -cos(theta)
  end subroutine generate_muon



  !==========================================================================
  ! Internal: momentum sampler
  !==========================================================================
  function sample_momentum(u) result(p)
    real(8), intent(in) :: u
    real(8) :: p
    integer :: j
    real(8) :: frac

    ! Mono-energetic: return fixed momentum regardless of u
    if (p_min_stored >= p_max_stored) then
      p = p_min_stored
      return
    end if

    if (spectrum_mode == 8) then
      ! Exact inverse CDF for dN/dE ∝ E^EL_IDX = E^-3: α = EL_IDX+1 = -2
      block
        real(8) :: alpha, A_lo, A_hi
        alpha = EL_IDX + 1.0d0
        A_lo  = p_min_stored ** alpha
        A_hi  = p_max_stored ** alpha
        p     = (A_lo + u * (A_hi - A_lo)) ** (1.0d0 / alpha)
      end block
      return
    end if

    if (spectrum_mode == 2) then
      p = p_min_stored * (u**(-1.0d0/2.7d0))
      if (p > p_max_stored) p = p_max_stored
      return
    end if

    if (use_analytical) then
      block
        real(8) :: alpha, A_lo, A_hi
        alpha = B_COSMO + 1.0d0
        A_lo  = p_min_stored ** alpha
        A_hi  = p_max_stored ** alpha
        p     = (A_lo + u * (A_hi - A_lo)) ** (1.0d0 / alpha)
      end block
      return
    end if

    if (u <= cdf(1)) then
      p = p_cdf(1); return
    end if
    do j = 2, NGRID
      if (u <= cdf(j)) then
        frac = (u - cdf(j-1)) / (cdf(j) - cdf(j-1))
        p    = p_cdf(j-1) + frac*(p_cdf(j)-p_cdf(j-1))
        return
      end if
    end do
    p = p_cdf(NGRID)
  end function sample_momentum



  !==========================================================================
  ! CosmoALEPH charge ratio
  !==========================================================================
  function charge_ratio_from_p(p) result(ratio)
    real(8), intent(in) :: p
    real(8) :: ratio
    if      (p <=  112d0) then; ratio = 1.252d0
    else if (p <=  141d0) then; ratio = 1.293d0
    else if (p <=  178d0) then; ratio = 1.259d0
    else if (p <=  224d0) then; ratio = 1.271d0
    else if (p <=  282d0) then; ratio = 1.239d0
    else if (p <=  355d0) then; ratio = 1.348d0
    else if (p <=  447d0) then; ratio = 1.541d0
    else if (p <=  562d0) then; ratio = 1.373d0
    else if (p <=  708d0) then; ratio = 1.243d0
    else if (p <=  891d0) then; ratio = 1.547d0
    else if (p <= 1122d0) then; ratio = 1.785d0
    else if (p <= 1413d0) then; ratio = 1.361d0
    else if (p <= 1778d0) then; ratio = 0.648d0
    else;                       ratio = 1.495d0
    end if
  end function charge_ratio_from_p



  !==========================================================================
  ! is_electron_mode — returns .true. when spectrum_mode == 8
  !==========================================================================
  function is_electron_mode() result(flag)
    logical :: flag
    flag = (spectrum_mode == 8)
  end function is_electron_mode


  !==========================================================================
  ! cos^2(theta) angular sampler
  !==========================================================================
  function sample_cos2(theta_max) result(theta)
    real(8), intent(in) :: theta_max
    real(8) :: theta
    real(4) :: yfl
    real(8) :: pdf, u
    real(8), parameter :: MAX_PDF = 0.385d0
    logical :: accepted

    accepted = .false.
    do while (.not. accepted)
      call par_ranlux(yfl)         ! <-- OMP change
      theta = dble(yfl) * theta_max
      pdf   = cos(theta)**2 * sin(theta)
      call par_ranlux(yfl)         ! <-- OMP change
      u     = dble(yfl)
      if (u*MAX_PDF < pdf) accepted = .true.
    end do
  end function sample_cos2


  !==========================================================================
  ! sample_cos3 — cos³θ sampler (Reyna-Bugaev angular distribution)
  ! Exact inverse CDF: cosθ = (1 − u·(1 − cos⁴θ_max))^(1/4)
  !==========================================================================
  function sample_cos3(theta_max) result(theta)
    real(8), intent(in) :: theta_max
    real(8) :: theta
    real(4) :: yfl
    real(8) :: u, cos4_max
    call par_ranlux(yfl)
    u        = dble(yfl)
    cos4_max = cos(theta_max)**4
    theta    = acos((1.0d0 - u*(1.0d0 - cos4_max))**0.25d0)
  end function sample_cos3


  !==========================================================================
  ! sample_guan_angle — self-consistent angular sampler for Guan/Frosin
  !==========================================================================
  function sample_guan_angle(E_GeV, theta_max) result(theta)
    real(8), intent(in) :: E_GeV, theta_max
    real(8)             :: theta

    integer, parameter :: NANG = 100
    real(8) :: cos_arr(NANG), pdf_arr(NANG), cdf_arr(NANG)
    real(8) :: cos_min, dc, c, frac, total
    real(8) :: a_p, b_p
    real(4) :: yfl
    real(8) :: u
    integer :: k

    if (spectrum_mode == 5) then
      a_p = FROSIN_A;  b_p = FROSIN_B
    else
      a_p = GUAN_A;    b_p = GUAN_B
    end if

    cos_min = cos(theta_max)
    do k = 1, NANG
      cos_arr(k) = cos_min + dble(k-1)/dble(NANG-1) * (1.0d0 - cos_min)
    end do

    do k = 1, NANG
      pdf_arr(k) = guan_flux(E_GeV, cos_arr(k), a_p, b_p) * cos_arr(k)
    end do

    cdf_arr(1) = 0.0d0
    do k = 2, NANG
      dc         = cos_arr(k) - cos_arr(k-1)
      cdf_arr(k) = cdf_arr(k-1) + 0.5d0*(pdf_arr(k)+pdf_arr(k-1))*dc
    end do
    total = cdf_arr(NANG)
    if (total > 0.0d0) then
      cdf_arr = cdf_arr / total
    else
      theta = 0.0d0;  return
    end if

    call par_ranlux(yfl);  u = dble(yfl)   ! <-- OMP change
    c = cos_arr(NANG)
    do k = 2, NANG
      if (u <= cdf_arr(k)) then
        if (cdf_arr(k) > cdf_arr(k-1)) then
          frac = (u - cdf_arr(k-1)) / (cdf_arr(k) - cdf_arr(k-1))
        else
          frac = 0.0d0
        end if
        c = cos_arr(k-1) + frac*(cos_arr(k) - cos_arr(k-1))
        exit
      end if
    end do
    c     = max(cos_min, min(1.0d0, c))
    theta = acos(c)
  end function sample_guan_angle


  !==========================================================================
  ! reyna_flux — Reyna–Bugaev (2006) dΦ/dp, arXiv:hep-ph/0604145, Eq. 6-7
  ! Returns differential flux [cm⁻²s⁻¹sr⁻¹(GeV/c)⁻¹] at momentum p [GeV/c]
  ! and zenith angle cosine cos_th.
  ! Formula: Φ = C0·x^-(C1+C2·log₁₀x+C3·log₁₀²x+C4·log₁₀³x),  x = p·cos*
  ! Integrates to 7.0e-3 cm⁻²s⁻¹sr⁻¹ above 1 GeV at cosθ=1 (PDG value).
  !==========================================================================
  function reyna_flux(p_GeV, cos_th) result(phi)
    real(8), intent(in) :: p_GeV, cos_th
    real(8) :: phi, cs, p_eff, lp, n_exp
    cs    = guan_cos_star(cos_th)
    p_eff = p_GeV * cs
    if (p_eff <= 0.0d0) then; phi = 0.0d0; return; end if
    lp    = log10(p_eff)
    n_exp = REYNA_C1 + REYNA_C2*lp + REYNA_C3*lp**2 + REYNA_C4*lp**3
    phi   = REYNA_C0 * p_eff**(-n_exp)
    if (phi < 0.0d0) phi = 0.0d0
  end function reyna_flux


  !==========================================================================
  ! guan_cos_star
  !==========================================================================
  function guan_cos_star(cos_th) result(cs)
    real(8), intent(in) :: cos_th
    real(8)             :: cs, numer
    numer = cos_th**2 + GUAN_P1**2 &
          + GUAN_P2 * cos_th**GUAN_P3 &
          + GUAN_P4 * cos_th**GUAN_P5
    cs = sqrt(max(0.0d0, numer)) / GUAN_DENOM
    cs = max(0.0d0, min(1.0d0, cs))
  end function guan_cos_star


  !==========================================================================
  ! guan_flux
  !==========================================================================
  function guan_flux(E_GeV, cos_th, a_par, b_par) result(phi)
    real(8), intent(in) :: E_GeV, cos_th, a_par, b_par
    real(8)             :: phi, cs, E_eff, pion_t, kaon_t
    cs     = guan_cos_star(cos_th)
    E_eff  = E_GeV * (1.0d0 + a_par / (E_GeV * cs**b_par))
    pion_t = 1.0d0     / (1.0d0 + 1.1d0 * E_GeV * cs / GUAN_EPI)
    kaon_t = GUAN_KF   / (1.0d0 + 1.1d0 * E_GeV * cs / GUAN_EK)
    phi    = GUAN_PRE  * E_eff**GUAN_IDX * (pion_t + kaon_t)
  end function guan_flux


end module ucmuon_source_module

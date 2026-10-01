!=============================================================================
! intensity_dump.f90 -- Fortran side of the cross-implementation flux test.
!
! Prints the same grid, vertical integrated flux and horizontal-surface
! quadrature as intensity_dump.cc, from src/generator/ucmuon_source_module.f90.
! The vertical integrated flux is the "Integrated flux" line that
! build_cosmoaleph_cdf itself prints, i.e. the number the GUI rate uses.
!
!   gfortran -O2 -fopenmp src/generator/rng_parallel.f90 \
!     src/generator/ucmuon_source_module.f90 tests/flux/intensity_dump.f90
!=============================================================================
program intensity_dump
  use ucmuon_source_module
  implicit none
  integer,  parameter :: modes(5) = [1, 4, 5, 6, 7]
  real(8),  parameter :: ps(8) = [1d0, 3d0, 10d0, 30d0, 100d0, 300d0, 1000d0, 2000d0]
  real(8),  parameter :: cs(5) = [1d0, 0.8660254037844387d0, 0.5d0, &
                                  0.25881904510252074d0, 0.1d0]
  real(8),  parameter :: emins(3) = [1d0, 10d0, 100d0]
  real(8),  parameter :: EMAX = 1500d0, THETA_MAX_DEG = 85d0
  integer :: i, j, k
  real(8) :: p_min, p_max

  do i = 1, size(modes)
    do j = 1, size(ps)
      do k = 1, size(cs)
        write(*,'(A,I2,3ES26.17)') 'I ', modes(i), ps(j), cs(k), &
             spectrum_intensity(modes(i), ps(j), cs(k))
      end do
    end do
  end do

  do i = 1, size(modes)
    do j = 1, size(emins)
      p_min = sqrt(emins(j)**2 - MUON_MASS**2)
      p_max = sqrt(EMAX**2 - MUON_MASS**2)
      write(*,'(A,I2,ES26.17)') 'DUMPCDF ', modes(i), emins(j)
      call build_cosmoaleph_cdf(p_min, p_max, modes(i))
      write(*,'(A,I2,2ES26.17)') 'Q ', modes(i), emins(j), &
           horizontal_rate(modes(i), p_min, p_max)
    end do
  end do

contains

  ! Same algorithm as horizontalRate() in intensity_dump.cc.
  function horizontal_rate(mode, p_min, p_max) result(total)
    integer, intent(in) :: mode
    real(8), intent(in) :: p_min, p_max
    integer, parameter :: NP = 2000, NC = 400
    real(8) :: total, c_lo, lr, c, inner, prev_p, prev_f, p, f
    integer :: jj, kk
    c_lo = cos(THETA_MAX_DEG * PI / 180d0)
    lr = log(p_max / p_min)
    total = 0d0
    do kk = 0, NC - 1
      c = c_lo + (kk + 0.5d0) * (1d0 - c_lo) / NC
      inner = 0d0
      prev_p = p_min
      prev_f = spectrum_intensity(mode, p_min, c)
      do jj = 1, NP - 1
        p = p_min * exp(lr * jj / (NP - 1))
        f = spectrum_intensity(mode, p, c)
        inner = inner + 0.5d0 * (f + prev_f) * (p - prev_p)
        prev_p = p
        prev_f = f
      end do
      total = total + c * inner * (1d0 - c_lo) / NC
    end do
    total = 2d0 * PI * total
  end function horizontal_rate

end program intensity_dump

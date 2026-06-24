!=======================================================================
! corgen.f90
! CORSET / CORGEN - correlated Gaussian random vector generator
! Called by MUSIC's MULTIPLE subroutine (N=2 always in MUSIC).
! Calls RNORML (rnorml.f) which calls RANMAR (ranmar.f).
!=======================================================================

subroutine corset(v, c, n)
! Cholesky decomposition of covariance matrix V.
! v(n,n) real(4) input  - symmetric positive-definite covariance
! c(n,n) real(4) output - lower-triangular Cholesky factor
  implicit none
  integer, intent(in)  :: n
  real(4), intent(in)  :: v(n,n)
  real(4), intent(out) :: c(n,n)
  integer :: i, j, k
  real(4) :: s

  do i = 1, n
    do j = 1, i
      s = v(i,j)
      do k = 1, j-1
        s = s - c(i,k) * c(j,k)
      end do
      if (i == j) then
        if (s <= 0.0) s = 1.0e-30
        c(i,j) = sqrt(s)
      else
        if (c(j,j) > 0.0) then
          c(i,j) = s / c(j,j)
        else
          c(i,j) = 0.0
        end if
      end if
    end do
    do j = i+1, n
      c(i,j) = 0.0
    end do
  end do

end subroutine corset


subroutine corgen(c, x, n)
! Draw one correlated Gaussian N-vector: x = C * z, z ~ N(0,I)
! c(n,n) real(4) input  - Cholesky factor from corset
! x(n)   real(4) output - correlated Gaussian random numbers
  implicit none
  integer, intent(in)  :: n
  real(4), intent(in)  :: c(n,n)
  real(4), intent(out) :: x(n)
  real(4) :: z(10)
  integer :: i, j

  call rnorml(z, n)

  do i = 1, n
    x(i) = 0.0
    do j = 1, i
      x(i) = x(i) + c(i,j) * z(j)
    end do
  end do

end subroutine corgen

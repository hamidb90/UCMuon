!=============================================================================
! ucmuon_rng.f90  -  uniform and Gaussian pseudorandom number generators
!
! Independent implementation of three published generators, written from
! the algorithm descriptions in the papers:
!   RANLUX  Luscher's subtract-with-borrow generator with luxury decimation:
!           M. Luscher, Comput. Phys. Commun. 79 (1994) 100, with the calling
!           conventions of F. James, Comput. Phys. Commun. 79 (1994) 111
!   RANMAR  the universal generator of G. Marsaglia, A. Zaman and
!           W. W. Tsang, Stat. Probab. Lett. 9 (1990) 35, seeded from one
!           integer as in F. James, Comput. Phys. Commun. 60 (1990) 329
!   RNORML  standard normal deviates by the ratio-of-uniforms method with
!           the quadratic bounds of J. L. Leva, ACM Trans. Math. Softw. 18
!           (1992) 449, using uniforms from RANMAR
! The output is identical, bit for bit, to that of the CERNLIB routines
! (RANLUX V115, RANMAR V113, RNORML) that UCMuon used up to 1.4.1; this is
! checked by tests/rng/test_rng.py against tests/rng/reference_vectors.txt.
!
! Entry points (external subroutines, implicit interface, default REAL):
!   call RANLUX(RVEC, LENV)         LENV uniform deviates in (0,1)
!   call RLUXGO(LUX, INS, K1, K2)   luxury level LUX = 0..4 (p = 24, 48, 97,
!                                   223, 389) or p = LUX for 24..2000; seed
!                                   INS (0: default); skip K1 + K2*10**9
!   call RANMAR(RVEC, LENV)         LENV uniform deviates in (0,1)
!   call RMARIN(IJKL, NTOT, NTOT2)  seed IJKL; skip NTOT + NTOT2*10**9
!   call RNORML(DEVIAS, NDEV)       NDEV standard normal deviates
! Without a seeding call, RANLUX starts from seed 314159265 at level 3 and
! RANMAR from IJKL = 54217137.  Every OpenMP thread owns its own state
! (threadprivate module variables; one per generator, so that a call costs
! a single thread-local lookup).
!
! MIT License.  Copyright (c) 2026 Hamid Basiri, UCLouvain Muography Group
!=============================================================================
module ucmuon_rng
  implicit none

  integer, parameter :: mask24 = 16777215          ! 2**24 - 1
  real,    parameter :: two_m24 = 2.0**(-24)
  real,    parameter :: two_m48 = 2.0**(-48)

  !---------------------------------------------------------------------------
  ! RANLUX.  Integers x in [0, 2**24) with x(n) = x(n-10) - x(n-24) - c(n-1),
  ! plus 2**24 and borrow c(n) = 1 when negative.  Of every p numbers the
  ! first 24 are returned and the other p-24 are thrown away.
  !---------------------------------------------------------------------------
  integer, parameter :: lux_pmax = 2000
  integer, parameter :: lux_seed_default = 314159265
  integer, parameter :: lux_ptab(0:4) = [24, 48, 97, 223, 389]

  type :: lux_state
    integer :: x(-23:lux_pmax) = 0   ! (1:p) current block, (-23:0) the 24 before
    integer :: c = 0                 ! borrow bit
    integer :: p = 223               ! block length p
    integer :: k = 25                ! next element of r to hand out
    real    :: r(24) = 0.0           ! the 24 kept numbers of the block, scaled
    logical :: unset = .true.        ! true until this thread's state is seeded
  end type lux_state

  !---------------------------------------------------------------------------
  ! RANMAR.  Lagged Fibonacci U(n) = U(n-97) - U(n-33) mod 1 on 24-bit
  ! fractions, minus the arithmetic sequence C(n) = C(n-1) - 7654321/2**24
  ! mod 16777213/2**24.  Both are held as integers in units of 2**-24 and
  ! advanced 97 steps at a time.
  !---------------------------------------------------------------------------
  integer, parameter :: mar_seed_default = 54217137
  integer, parameter :: mar_cstep = 7654321, mar_cmod = 16777213

  type :: mar_state
    integer :: u(97) = 0             ! u(k) = U(n-98+k), n = next index
    integer :: c = 0                 ! C(n-1)
    integer :: k = 98                ! next element of r to hand out
    real    :: r(97) = 0.0           ! outputs of the current 97 steps
    logical :: unset = .true.
  end type mar_state

  type(lux_state) :: luxst
  type(mar_state) :: marst
  !$omp threadprivate(luxst, marst)

contains

  !---------------------------------------------------------------------------
  ! RANLUX internals
  !---------------------------------------------------------------------------

  ! Fill x(1:m) with the next m raw numbers (m <= lux_pmax), then move the
  ! last 24 of them to x(-23:0) unless keep is set.
  subroutine lux_raw(m, keep)
    integer, intent(in) :: m
    logical, intent(in) :: keep
    integer :: n, t, c
    c = luxst%c
    do n = 1, m
      t = luxst%x(n-10) - luxst%x(n-24) - c
      c = ishft(t, -31)              ! 1 if t < 0
      luxst%x(n) = iand(t, mask24)   ! t + 2**24 if t < 0
    end do
    luxst%c = c
    if (.not. keep) luxst%x(-23:0) = luxst%x(m-23:m)
  end subroutine lux_raw

  ! Generate one block of p numbers and prepare the 24 that are returned.
  ! A number below 2**-12 is completed with x(n-9) * 2**-48, so that it
  ! still carries about 24 significant bits (and is never zero).
  subroutine lux_block()
    integer :: n, p
    p = luxst%p
    call lux_raw(p, .true.)
    do n = 1, 24
      if (luxst%x(n) >= 4096) then
        luxst%r(n) = real(luxst%x(n)) * two_m24
      else if (luxst%x(n) > 0 .or. luxst%x(n-9) > 0) then
        luxst%r(n) = real(luxst%x(n)) * two_m24 + real(luxst%x(n-9)) * two_m48
      else
        luxst%r(n) = two_m48
      end if
    end do
    luxst%x(-23:0) = luxst%x(p-23:p)
    luxst%k = 1
  end subroutine lux_block

  ! Throw away nskip raw numbers.
  subroutine lux_skip(nskip)
    integer(8), intent(in) :: nskip
    integer(8) :: left
    integer :: m
    left = nskip
    do while (left > 0)
      m = int(min(left, int(lux_pmax, 8)))
      call lux_raw(m, .false.)
      left = left - m
    end do
  end subroutine lux_skip

  ! Initial 24 numbers from the multiplicative congruential generator
  ! s <- 40014 s mod 2147483563 (L'Ecuyer), low 24 bits of each value.  The
  ! first value becomes x(0), the newest; the borrow starts at 1 only if
  ! the oldest, x(-23), is zero.
  subroutine lux_seed(iseed)
    integer, intent(in) :: iseed
    integer(8) :: s
    integer :: n
    s = int(iseed, 8)
    do n = 0, -23, -1
      s = mod(40014_8 * s, 2147483563_8)
      luxst%x(n) = int(iand(s, int(mask24, 8)))
    end do
    luxst%c = merge(1, 0, luxst%x(-23) == 0)
    luxst%k = 25
    luxst%unset = .false.
  end subroutine lux_seed

  !---------------------------------------------------------------------------
  ! RANMAR internals
  !---------------------------------------------------------------------------

  ! Advance the lagged Fibonacci state by 97 steps, in place.
  subroutine mar_lag97()
    integer :: k
    do k = 1, 33
      marst%u(k) = iand(marst%u(k) - marst%u(k+64), mask24)
    end do
    do k = 34, 97
      marst%u(k) = iand(marst%u(k) - marst%u(k-33), mask24)
    end do
  end subroutine mar_lag97

  ! Next 97 outputs.  An exact zero is replaced by U(m) * 2**-48, where
  ! U(m) is the latest value with m = 95 (mod 97), counting from the seed
  ! (u(96) of this block or, earlier in the block, of the previous one).
  subroutine mar_block()
    integer :: k, c, w, last96, fill
    last96 = marst%u(96)
    call mar_lag97()
    c = marst%c
    do k = 1, 97
      c = c - mar_cstep
      if (c < 0) c = c + mar_cmod
      w = iand(marst%u(k) - c, mask24)
      if (w > 0) then
        marst%r(k) = real(w) * two_m24
      else
        fill = merge(marst%u(96), last96, k >= 96)
        marst%r(k) = real(max(fill, 1)) * two_m48
      end if
    end do
    marst%c = c
    marst%k = 1
  end subroutine mar_block

  ! Seed from one integer: split it into the four seeds of the Marsaglia-
  ! Zaman-Tsang initialisation (a 3-lag Fibonacci product mod 179 and a
  ! congruential sequence mod 169), which give the 24 bits of each U.
  subroutine mar_seed(iseed)
    integer, intent(in) :: iseed
    integer :: hi, lo, f1, f2, f3, f4, g, n, b, bit, word
    hi = iseed / 30082
    lo = iseed - 30082 * hi
    f1 = mod(hi / 177, 177) + 2
    f2 = mod(hi, 177) + 2
    f3 = mod(lo / 169, 178) + 1
    g  = mod(lo, 169)
    do n = 97, 1, -1
      word = 0
      bit = 8388608                   ! 2**23
      do b = 1, 24
        f4 = mod(mod(f1 * f2, 179) * f3, 179)
        f1 = f2
        f2 = f3
        f3 = f4
        g = mod(53 * g + 1, 169)
        if (mod(g * f4, 64) >= 32) word = word + bit
        bit = bit / 2
      end do
      marst%u(n) = word
    end do
    marst%c = 362436
    marst%k = 98
    marst%unset = .false.
  end subroutine mar_seed

  ! Skip nskip outputs.
  subroutine mar_skip(nskip)
    integer(8), intent(in) :: nskip
    integer(8) :: nblk, i
    integer :: rest
    if (nskip <= 0) return
    nblk = nskip / 97
    rest = int(nskip - 97 * nblk)
    do i = 1, nblk
      call mar_lag97()
    end do
    marst%c = int(modulo(int(marst%c, 8) - mod(97 * nblk, int(mar_cmod, 8)) &
                         * mar_cstep, int(mar_cmod, 8)))
    if (rest > 0) then
      call mar_block()
      marst%k = rest + 1
    end if
  end subroutine mar_skip

  ! One uniform from the RANMAR stream (for RNORML).
  real function mar_next()
    if (marst%unset) call mar_seed(mar_seed_default)
    if (marst%k > 97) call mar_block()
    mar_next = marst%r(marst%k)
    marst%k = marst%k + 1
  end function mar_next

end module ucmuon_rng

!=============================================================================
! RANLUX(RVEC, LENV): LENV uniform deviates
!=============================================================================
subroutine ranlux(rvec, lenv)
  use ucmuon_rng
  implicit none
  integer, intent(in) :: lenv
  real, intent(out) :: rvec(*)
  integer :: i, n
  if (luxst%unset) then
    luxst%p = lux_ptab(3)
    call lux_seed(lux_seed_default)
  end if
  i = 0
  do while (i < lenv)
    if (luxst%k > 24) call lux_block()
    n = min(lenv - i, 25 - luxst%k)
    rvec(i+1:i+n) = luxst%r(luxst%k:luxst%k+n-1)
    luxst%k = luxst%k + n
    i = i + n
  end do
end subroutine ranlux

!=============================================================================
! RLUXGO(LUX, INS, K1, K2): luxury, seed, then skip K1 + K2*10**9 numbers
!   LUX < 0: level 3;  0..4: level;  24..2000: p = LUX;  otherwise level 4.
!   The skip counts every number generated, returned or thrown away (the
!   count James' RLUXAT reports); if it does not end within the 24 returned
!   numbers of a block, output resumes at the start of a new block.
!=============================================================================
subroutine rluxgo(lux, ins, k1, k2)
  use ucmuon_rng
  implicit none
  integer, intent(in) :: lux, ins, k1, k2
  integer(8) :: nskip
  integer :: phase

  if (lux < 0) then
    luxst%p = lux_ptab(3)
  else if (lux <= 4) then
    luxst%p = lux_ptab(lux)
  else if (lux >= 24 .and. lux <= lux_pmax) then
    luxst%p = lux
  else
    write(*, '(A,I0,A)') ' RLUXGO: illegal luxury level ', lux, ', using level 4'
    luxst%p = lux_ptab(4)
  end if

  if (ins < 0) then
    write(*, '(A,I0,A)') ' RLUXGO: negative seed ', ins, ', using the default seed'
    call lux_seed(lux_seed_default)
  else if (ins == 0) then
    call lux_seed(lux_seed_default)
  else
    call lux_seed(ins)
  end if

  nskip = int(k1, 8) + 1000000000_8 * int(k2, 8)
  if (nskip > 0) then
    phase = int(mod(nskip, int(luxst%p, 8)))
    if (phase < 24) then
      call lux_skip(nskip - phase)
      if (phase > 0) then
        call lux_block()
        luxst%k = phase + 1
      end if
    else
      call lux_skip(nskip)
    end if
  end if
end subroutine rluxgo

!=============================================================================
! RANMAR(RVEC, LENV): LENV uniform deviates
!=============================================================================
subroutine ranmar(rvec, lenv)
  use ucmuon_rng
  implicit none
  integer, intent(in) :: lenv
  real, intent(out) :: rvec(*)
  integer :: i, n
  if (marst%unset) call mar_seed(mar_seed_default)
  i = 0
  do while (i < lenv)
    if (marst%k > 97) call mar_block()
    n = min(lenv - i, 98 - marst%k)
    rvec(i+1:i+n) = marst%r(marst%k:marst%k+n-1)
    marst%k = marst%k + n
    i = i + n
  end do
end subroutine ranmar

!=============================================================================
! RMARIN(IJKL, NTOT, NTOT2): seed RANMAR, then skip NTOT + NTOT2*10**9
!=============================================================================
subroutine rmarin(ijkl, ntot, ntot2)
  use ucmuon_rng
  implicit none
  integer, intent(in) :: ijkl, ntot, ntot2
  call mar_seed(ijkl)
  call mar_skip(int(ntot, 8) + 1000000000_8 * int(ntot2, 8))
end subroutine rmarin

!=============================================================================
! RNORML(DEVIAS, NDEV): NDEV standard normal deviates (Leva 1992).
! Each trial takes two RANMAR uniforms u and w, sets v = 1.7156 (w - 1/2)
! and accepts v/u if (u, v) lies inside the ratio-of-uniforms region
! v**2 <= -4 u**2 ln(u).  Two quadratic curves around the boundary decide
! almost every trial without the logarithm.  The final test is evaluated
! in single precision in this order (u**2 first), as in the reference.
!=============================================================================
subroutine rnorml(devias, ndev)
  use ucmuon_rng
  implicit none
  integer, intent(in) :: ndev
  real, intent(out) :: devias(*)
  real, parameter :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472
  real, parameter :: r_in = 0.27597, r_out = 0.27846
  real :: u, v, x, y, q
  integer :: i
  do i = 1, ndev
    do
      u = mar_next()
      v = 1.7156 * (mar_next() - 0.5)
      x = u - s
      y = abs(v) - t
      q = x * x + y * (a * y - b * x)
      if (q < r_in) exit
      if (q > r_out) cycle
      if (v * v <= -4.0 * u**2 * log(u)) exit
    end do
    devias(i) = v / u
  end do
end subroutine rnorml

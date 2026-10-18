!=============================================================================
! rng_dump.f90  -  dump the streams of the Fortran random number routines
!
! Prints "case position bits" lines, where bits is the IEEE single-precision
! bit pattern of the deviate (as a signed 32-bit integer), so two builds can
! be compared exactly with diff.  tests/rng/reference_vectors.txt holds the
! output of this program linked against the CERNLIB RANLUX (Mathlib V115),
! RANMAR (Kernlib V113) and RNORML routines that UCMuon shipped up to 1.4.1;
! tests/rng/test_rng.py checks src/common/ucmuon_rng.f90 against it.
!
! Usage:  rng_dump MODE     MODE = ldef | mdef | ndef | main | omp
!   ldef, mdef, ndef  first calls of RANLUX / RANMAR / RNORML without any
!                     seeding call (each needs a fresh process)
!   main              seeded streams, luxury levels, skips, vector lengths
!   omp               four OpenMP threads seeded as the transport engines do
!=============================================================================
program rng_dump
  use omp_lib
  implicit none
  external :: ranlux, rluxgo, ranmar, rmarin, rnorml
  character(16) :: mode
  integer :: i, k, s, l
  integer, parameter :: lseeds(6) = [1, 314159, 371831, 12345678, 2147483646, 0]
  integer, parameter :: mseeds(6) = [1, 54217137, 271828, 371819, 900000000, 0]
  integer, parameter :: luxes(10) = [0, 1, 2, 3, 4, -1, 5, 30, 223, 2001]
  real(4) :: v(1000), z(7)
  real(4) :: tv(200, 0:3)
  character(32) :: tag

  call get_command_argument(1, mode)
  select case (trim(mode))

  case ('ldef')
    do i = 1, 300
      call ranlux(v, 1)
      call put('ranlux_default', i, v(1))
    end do

  case ('mdef')
    do i = 1, 300
      call ranmar(v, 1)
      call put('ranmar_default', i, v(1))
    end do

  case ('ndef')
    do i = 1, 150
      call rnorml(z, 2)
      call put('rnorml_default', 2*i-1, z(1))
      call put('rnorml_default', 2*i,   z(2))
    end do

  case ('main')
    ! RANLUX: every luxury level and seed, singles then long vector draws
    do l = 1, size(luxes)
      do s = 1, size(lseeds)
        if (l > 5 .and. s > 2) cycle
        write(tag, '(A,I0,A,I0)') 'ranlux_l', luxes(l), '_s', lseeds(s)
        call rluxgo(luxes(l), lseeds(s), 0, 0)
        call stream_lux(trim(tag))
      end do
    end do
    ! RLUXGO skip: K1 numbers are skipped after seeding
    call rluxgo(3, 314159, 1000, 0)
    do i = 1, 50
      call ranlux(v, 1)
      call put('ranlux_l3_s314159_skip1000', i, v(1))
    end do
    ! RANMAR: seeds, then the NTOT skip of RMARIN
    do s = 1, size(mseeds)
      write(tag, '(A,I0)') 'ranmar_s', mseeds(s)
      call rmarin(mseeds(s), 0, 0)
      call stream_mar(trim(tag))
    end do
    call rmarin(271828, 1000, 0)
    do i = 1, 50
      call ranmar(v, 1)
      call put('ranmar_s271828_skip1000', i, v(1))
    end do
    ! RNORML in pairs (as CORGEN calls it for MUSIC), then the RANMAR draws
    ! that follow, which show how many uniforms RNORML consumed
    call rmarin(271828, 0, 0)
    do i = 1, 500
      call rnorml(z, 2)
      call put('rnorml_s271828_n2', 2*i-1, z(1))
      call put('rnorml_s271828_n2', 2*i,   z(2))
    end do
    do i = 1, 10
      call ranmar(v, 1)
      call put('rnorml_s271828_then_ranmar', i, v(1))
    end do
    call rmarin(371819, 0, 0)
    k = 0
    do i = 1, 300
      call rnorml(z, 1)
      k = k + 1
      call put('rnorml_s371819_n1_n7', k, z(1))
    end do
    do i = 1, 100
      call rnorml(z, 7)
      do s = 1, 7
        k = k + 1
        call put('rnorml_s371819_n1_n7', k, z(s))
      end do
    end do

  case ('omp')
    ! Seeded inside a parallel region exactly as the BB and MUSIC engines do
    call omp_set_num_threads(4)
    !$omp parallel private(k, i, v) num_threads(4)
    k = omp_get_thread_num()
    call rluxgo(3, 314159 + k * 100003, 0, 0)
    do i = 1, 100
      call ranlux(v, 1)
      tv(i, k) = v(1)
    end do
    call rmarin(271828 + k * 99991, 0, 0)
    do i = 1, 50
      call rnorml(v, 2)
      tv(100 + 2*i - 1, k) = v(1)
      tv(100 + 2*i,     k) = v(2)
    end do
    !$omp end parallel
    do k = 0, 3
      write(tag, '(A,I0)') 'omp_thread', k
      do i = 1, 200
        call put(trim(tag), i, tv(i, k))
      end do
    end do

  case default
    write(*, '(A)') 'usage: rng_dump ldef|mdef|ndef|main|omp'
    stop 2
  end select

contains

  subroutine put(name, pos, x)
    character(*), intent(in) :: name
    integer, intent(in) :: pos
    real(4), intent(in) :: x
    write(*, '(A,1X,I0,1X,I0)') name, pos, transfer(x, 0)
  end subroutine put

  ! 200 single draws, then vector draws of 4, 7 and 997 (the stream must not
  ! depend on how it is cut), recording a few positions up to 100 004
  subroutine stream_lux(name)
    character(*), intent(in) :: name
    call stream(name, 1)
  end subroutine stream_lux

  subroutine stream_mar(name)
    character(*), intent(in) :: name
    call stream(name, 2)
  end subroutine stream_mar

  subroutine stream(name, which)
    character(*), intent(in) :: name
    integer, intent(in) :: which
    integer :: pos, j, n, lens(3), il
    real(4) :: w(997)
    lens = [4, 7, 997]
    pos = 0
    do j = 1, 200
      call draw(which, w, 1)
      pos = pos + 1
      call put(name, pos, w(1))
    end do
    il = 0
    do while (pos < 100004)
      il = mod(il, 3) + 1
      n = lens(il)
      call draw(which, w, n)
      do j = 1, n
        pos = pos + 1
        if (pos <= 260 .or. pos == 1000 .or. (pos >= 10000 .and. pos <= 10004) &
            .or. (pos >= 100000 .and. pos <= 100004)) call put(name, pos, w(j))
      end do
    end do
  end subroutine stream

  subroutine draw(which, w, n)
    integer, intent(in) :: which, n
    real(4), intent(out) :: w(n)
    if (which == 1) then
      call ranlux(w, n)
    else
      call ranmar(w, n)
    end if
  end subroutine draw

end program rng_dump

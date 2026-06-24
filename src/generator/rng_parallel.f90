!=============================================================================
! rng_parallel.f90  —  Thread-local RNG for OpenMP parallelisation
! UCLouvain CosmoALEPH / Muography Group
!=============================================================================
!
! INTERFACE NOTE
!   The dummy argument is a plain SCALAR (real(4) :: r), NOT an array.
!   This avoids the rank-mismatch error that gfortran raises when an explicit
!   Fortran-90 module interface is present and a scalar is passed to an
!   assumed-size array dummy (real(4)(*)).
!
!   Every call site uses:   call par_ranlux(yfl)      (n argument dropped)
!   Replaces:               call RANLUX(yfl, 1)
!
! THREAD-SAFETY
!   tl_seed is !$OMP THREADPRIVATE — each thread owns its own private copy.
!
!=============================================================================
module rng_parallel
  use omp_lib
  implicit none
  private

  integer(8), save :: tl_seed = 1234567890123456789_8
  !$OMP THREADPRIVATE(tl_seed)

  public :: par_ranlux
  public :: par_init_rng

contains

  subroutine par_init_rng(base_seed)
    integer, intent(in) :: base_seed
    integer   :: tid
    integer(8) :: s
    !$OMP PARALLEL PRIVATE(tid, s)
    tid = omp_get_thread_num()
    s = int(base_seed, 8) * 6364136223846793005_8 &
      + int(tid + 1,   8) * 2654435761_8          &
      + 1442695040888963407_8
    s = s * 6364136223846793005_8 + 1442695040888963407_8
    s = s * 6364136223846793005_8 + 1442695040888963407_8
    s = s * 6364136223846793005_8 + 1442695040888963407_8 + int(tid, 8)
    s = s * 6364136223846793005_8 + 1442695040888963407_8
    tl_seed = s
    !$OMP END PARALLEL
    write(*,'(A,I4,A)') &
      '  RNG streams initialised for ', omp_get_max_threads(), ' thread(s).'
  end subroutine par_init_rng

  ! Generate ONE real(4) in (0,1).  Scalar dummy — no rank mismatch.
  subroutine par_ranlux(r)
    real(4), intent(out) :: r
    integer(8) :: bits
    real(4), parameter :: INV24 = 5.96046448e-8_4
    tl_seed = tl_seed * 6364136223846793005_8 + 1442695040888963407_8
    bits    = ibits(tl_seed, 40, 24)
    if (bits == 0_8) bits = 1_8
    r = real(bits, 4) * INV24
  end subroutine par_ranlux

end module rng_parallel

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
  public :: ucmuon_base_seed

contains

  !---------------------------------------------------------------------------
  ! Base seed for one run.
  !
  ! If the environment variable UCMUON_SEED holds an integer, it is used
  ! verbatim, which makes a run exactly reproducible (needed by the validation
  ! suite, and by any paper that wants to quote its seed).  Otherwise a seed is
  ! derived from the wall clock at millisecond resolution and mixed with the
  ! process id, so that jobs launched at the same instant (SLURM arrays,
  ! repeated GUI runs) still get distinct streams.
  !
  ! History: this used to be  tim(6) + tim(5)*60 + tim(4)*3600, which reads
  ! DATE_AND_TIME's VALUES as if it were (.., hour, minute, second).  It is
  ! actually (year, month, day, UTC-offset-minutes, hour, minute, second, ms),
  ! so the expression evaluated to minute + hour*60 + utc_offset*3600: it never
  ! touched seconds or milliseconds and stayed constant for a whole wall-clock
  ! minute.  Two runs started in the same minute produced byte-identical
  ! output, silently duplicating events across concurrent jobs.
  !---------------------------------------------------------------------------
  function ucmuon_base_seed() result(seed)
    integer :: seed
    integer :: tim(8), ios, ln, st
    integer(8) :: s8
    character(32) :: env

    call get_environment_variable('UCMUON_SEED', env, ln, st)
    if (st == 0 .and. ln > 0) then
      read(env, *, iostat=ios) seed
      if (ios == 0) then
        write(*,'(A,I0,A)') '  RNG base seed: ', seed, '   (UCMUON_SEED)'
        return
      end if
      write(*,'(A)') '  WARNING: UCMUON_SEED is not an integer — using clock seed.'
    end if

    call DATE_AND_TIME(VALUES=tim)
    s8 = int(tim(8), 8)                &   ! millisecond
       + int(tim(7), 8) * 1000_8       &   ! second
       + int(tim(6), 8) * 60000_8      &   ! minute
       + int(tim(5), 8) * 3600000_8        ! hour
    s8   = s8 * 1000003_8 + int(getpid(), 8)
    seed = int(mod(abs(s8), 2147483647_8), 4)
    write(*,'(A,I0,A)') '  RNG base seed: ', seed, '   (clock+pid)'
  end function ucmuon_base_seed

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

!=============================================================================
! parma_path_module.f90 — configurable data directory for standalone PARMA
! UCLouvain Muography Group
!
! Companion module to parma_subroutines.f90 (EXPACS/PARMA v4.19, Sato JAEA).
! See parma_subroutines.f90 for full license terms and citation requirements.
!
! Usage:
!   call parma_set_datadir('/path/to/data/EXPACS/parma/')
!   (default: '' → opens files relative to CWD, i.e. 'input/...' as shipped)
!=============================================================================
module parma_path
  implicit none
  character(200), save :: parma_datadir = ''
  ! .true. silences getHP's "no W index for this date" warning while the
  ! generator walks back to the latest measured day (it reports that itself).
  logical, save :: parma_quiet = .false.
contains
  subroutine parma_set_datadir(path)
    character(*), intent(in) :: path
    integer :: n
    parma_datadir = trim(path)
    n = len_trim(parma_datadir)
    if (n > 0 .and. parma_datadir(n:n) /= '/') then
      parma_datadir = trim(parma_datadir)//'/'
    end if
  end subroutine parma_set_datadir

  ! .true. if the date exists. PARMA 4.19's getHP stops on one that does not
  ! (30 February), so callers check first.
  logical function parma_valid_date(y, m, d)
    integer, intent(in) :: y, m, d
    integer :: mdays
    parma_valid_date = .false.
    if (m < 1 .or. m > 12) return
    select case (m)
    case (4, 6, 9, 11)
      mdays = 30
    case (2)
      mdays = 28
      if ((mod(y, 4) == 0 .and. mod(y, 100) /= 0) .or. mod(y, 400) == 0) mdays = 29
    case default
      mdays = 31
    end select
    parma_valid_date = (d >= 1 .and. d <= mdays)
  end function parma_valid_date
end module parma_path

!=============================================================================
! parma_path_module.f90 — configurable data directory for standalone PARMA
! UCLouvain Muography Group
!
! Companion module to parma_subroutines.f90 (EXPACS/PARMA v4.10, Sato JAEA).
! See parma_subroutines.f90 for full license terms and citation requirements.
!
! Usage:
!   call parma_set_datadir('/path/to/data/EXPACS/parma/')
!   (default: '' → opens files relative to CWD, i.e. 'input/...' as shipped)
!=============================================================================
module parma_path
  implicit none
  character(200), save :: parma_datadir = ''
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
end module parma_path

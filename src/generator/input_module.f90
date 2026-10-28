!=============================================================================
! input_module.f90 - read and check the generator's input parameters
!
! The generators read one value per line from standard input. Every value is
! read here with its name, so that a missing or unreadable line stops the run
! with a message naming the parameter, and every range check stops it with a
! message and exit status 2. Up to UCMuon 1.4.1 bad values were replaced
! silently (an angular mode of 9 became 2, e_min > e_max gave every muon the
! energy e_min, a negative muon count gave an empty file, all with status 0).
!
! The MPI generator reads on rank 0 only; it points input_abort at a routine
! that calls MPI_Abort, so that the other ranks stop as well.
!=============================================================================
module input_module
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none
  private
  public :: get_int, get_real, get_real2, get_real3, get_str, need, &
            input_abort, i2s, r2s

  abstract interface
    subroutine abort_iface()
    end subroutine abort_iface
  end interface
  procedure(abort_iface), pointer :: input_abort => null()

  integer :: n_read = 0      ! values read so far, to locate a missing one

contains

  subroutine fail(msg)
    character(*), intent(in) :: msg
    write(error_unit, '(A)') ' INPUT ERROR: ' // trim(msg)
    flush(error_unit)
    if (associated(input_abort)) call input_abort()
    stop 2
  end subroutine fail

  subroutine check_read(ios, name, what)
    integer,      intent(in) :: ios
    character(*), intent(in) :: name, what
    n_read = n_read + 1
    if (ios < 0) call fail('the input ends before ' // name // ' (value ' // &
                           i2s(n_read) // ' of the input)')
    if (ios > 0) call fail(name // ' is not ' // what // ' (value ' // &
                           i2s(n_read) // ' of the input)')
  end subroutine check_read

  subroutine get_int(name, v)
    character(*), intent(in)  :: name
    integer,      intent(out) :: v
    integer :: ios
    read(*, *, iostat=ios) v
    call check_read(ios, name, 'an integer')
  end subroutine get_int

  subroutine get_real(name, v)
    character(*), intent(in)  :: name
    real(8),      intent(out) :: v
    integer :: ios
    read(*, *, iostat=ios) v
    call check_read(ios, name, 'a number')
  end subroutine get_real

  subroutine get_real2(name, a, b)
    character(*), intent(in)  :: name
    real(8),      intent(out) :: a, b
    integer :: ios
    read(*, *, iostat=ios) a, b
    call check_read(ios, name, 'two numbers')
  end subroutine get_real2

  subroutine get_real3(name, a, b, c)
    character(*), intent(in)  :: name
    real(8),      intent(out) :: a, b, c
    integer :: ios
    read(*, *, iostat=ios) a, b, c
    call check_read(ios, name, 'three numbers')
  end subroutine get_real3

  ! A text line; an empty line is allowed (the caller supplies the default).
  subroutine get_str(name, s)
    character(*), intent(in)  :: name
    character(*), intent(out) :: s
    integer :: ios
    read(*, '(A)', iostat=ios) s
    call check_read(ios, name, 'a text line')
  end subroutine get_str

  subroutine need(ok, msg)
    logical,      intent(in) :: ok
    character(*), intent(in) :: msg
    if (.not. ok) call fail(msg)
  end subroutine need

  function i2s(i) result(s)
    integer, intent(in) :: i
    character(:), allocatable :: s
    character(24) :: buf
    write(buf, '(I0)') i
    s = trim(buf)
  end function i2s

  function r2s(x) result(s)
    real(8), intent(in) :: x
    character(:), allocatable :: s
    character(32) :: buf
    write(buf, '(G0.6)') x
    s = trim(adjustl(buf))
  end function r2s

end module input_module

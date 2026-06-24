!=============================================================================
! phits_module.f90
! Writes one muon per line in PHITS s-type=17 ASCII dump format.
!
! Record fields (IDs 1..10):
!   kf  x  y  z  u  v  w  e  wt  time
!
! In PHITS input use:
!   [Source]
!   s-type = 17
!   file   = muons_for_phits.dat
!   dump   = -10
!   1 2 3 4 5 6 7 8 9 10
!
! kf codes (PDG convention):  mu- = 13,  mu+ = -13
! Energy unit: kinetic energy in MeV  (e = (E_total - m_mu) * 1000)
!=============================================================================
module phits_module
  implicit none
  private
  public :: write_phits_dump

contains

  subroutine write_phits_dump(iu, kf, x,y,z, u,v,w, ekin_MeV, wt, time_ns)
    integer, intent(in) :: iu, kf
    real(8), intent(in) :: x,y,z, u,v,w, ekin_MeV, wt, time_ns
    ! Format matches PHITS dump-a: (30(1p1d24.15))
    ! No leading spaces, D-exponent, 24 wide, 15 decimal places
    write(iu,'(10(1pd24.15))') dble(kf), x,y,z, u,v,w, ekin_MeV, wt, time_ns
  end subroutine write_phits_dump

end module phits_module

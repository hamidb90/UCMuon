!=============================================================================
! parma_stub.f90: stand-ins for the PARMA routines the generator calls, linked
! when PARMA is not installed (src/parma/parma_subroutines.f90 absent).
! UCLouvain Muography Group, MIT License 2026
!
! PARMA (T. Sato, JAEA) is not shipped with UCMuon; tools/install_parma.py
! downloads it from JAEA (setup.sh and install.ps1 offer to run it). Without
! it every spectrum except mode 3 works, and mode 3 stops here, before any
! muon is generated, with a message saying how to install it.
!=============================================================================
subroutine parma_missing()
  implicit none
  write(*,'(A)') ''
  write(*,'(A)') ' ERROR: spectrum mode 3 needs PARMA, which is not installed.'
  write(*,'(A)') '        PARMA (JAEA, non-commercial use) is downloaded from JAEA by'
  write(*,'(A)') '          python3 tools/install_parma.py'
  write(*,'(A)') '        then rebuild the generator (make local), or simply re-run'
  write(*,'(A)') '        setup.sh (Windows: install_windows.bat), which does both.'
  error stop 3
end subroutine parma_missing

real(8) function getMuonSpec(ip, s, r1, d1, e)
  implicit none
  integer :: ip
  real(8) :: s, r1, d1, e
  getMuonSpec = 0.0d0
  call parma_missing()
end function getMuonSpec

real(8) function getSpecAngFinal(ip, s, r, d, e, g, ang)
  implicit none
  integer :: ip
  real(8) :: s, r, d, e, g, ang
  getSpecAngFinal = 0.0d0
  call parma_missing()
end function getSpecAngFinal

real(8) function getd(alti, cido)
  implicit none
  real(8) :: alti, cido
  getd = 0.0d0
  call parma_missing()
end function getd

real(8) function getr(cido, ckei)
  implicit none
  real(8) :: cido, ckei
  getr = 0.0d0
  call parma_missing()
end function getr

real(8) function getFFPfromW(s)
  implicit none
  real(8) :: s
  getFFPfromW = 0.0d0
  call parma_missing()
end function getFFPfromW

real(8) function getHP(iy0, im0, id0, ic)
  implicit none
  integer :: iy0, im0, id0, ic
  getHP = 0.0d0
  call parma_missing()
end function getHP

!-------------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                     !
!-------------------------------------------------------------------------
!bop
!
! !Module:  m_kind
!
! !Description: This module declares a real kind with decimal
!precision of at least 8 digits. The portability is allowed.
!\\
!\\
! !Interface:
!
module m_kind
!
  implicit none
!
! !Public data members:
!
  integer, public, parameter :: rkind = selected_real_kind(14)
  integer, public :: m_seed ! Initialize ranlux code
  integer, public, parameter :: kacum = 500 ! Accumulator statistical error 
!
  !Revision history:
! 05Aug 2015 Asdrubal Lozada
!
!eop
!------------------------------------------------------------------------
end module m_kind

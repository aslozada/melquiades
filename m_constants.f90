!---------------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                       !
!---------------------------------------------------------------------------
!bop
!
! !Module:  m_constants
!
! !Description: This module defines some constant values.
!\\
!\\
! !Interface:
!
module m_constants
!
! !Uses:
!
  use m_kind
  implicit none
!
! !Revision history:
!  06Aug 2015 Asdrubal Lozada
!
!eop
!--------------------------------------------------------------------------
    real(rkind), parameter :: kB = 0.0019872041_rkind ! kcal/mol/K
    real(rkind), parameter :: pi = 3.141592653589793_rkind
    real(rkind), parameter :: c_kcalmol = 332.06574_rkind ! Coulumb to kcal/mol
    real(rkind), parameter :: c_k = 273.15_rkind ! Celsius to Kelvin degree
    real(rkind), parameter :: pv_conv = 0.000014576_rkind ! relation P x V
    real(rkind), parameter :: d_con = 1.6605393_rkind
    integer,     parameter :: ndim = 3 ! Space dimension

end module m_constants

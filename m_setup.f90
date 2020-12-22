!---------------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                       !
!---------------------------------------------------------------------------
!bop
!
! !Module:  m_init
!
! !Description: This module contains routines for initializing the variables
!and activating the Linked-Cell List 
!\\
!\\
! !Interface:
!
module m_setup
  use m_kind
  use m_simtype
  use m_boxtype
  use m_read
  use m_init
  use m_zeros
  use m_configuration
  use m_metropolis
  use m_dealloc
  use m_error
  use m_unit
!
 implicit none

 public :: r_setup

 contains

subroutine r_setup( edge, t, y , x)
  implicit none
  type(simulation), intent(inout) :: y  
  type(box), pointer :: x
  type(temporary), intent(inout) :: t
  real(rkind) :: engconf, virconf
  real(rkind), dimension(:), intent(inout) :: edge

  
  call r_configuration(edge, engconf, virconf, t, y, x, 0, .true.)  
  call r_metropolis(edge, engconf, virconf, t, y, x)
  call r_configuration(edge, engconf, virconf, t, y, x, 1, .true.)

end subroutine r_setup  

end module m_setup  

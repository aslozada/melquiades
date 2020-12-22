!---------------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                       !
!---------------------------------------------------------------------------
!bop
!
! !Module:  m_cells
!
! !Description: This module contains a routine for defining the 
!maximum number of neighbors in the Linked-Cells List.
!\\
!\\
! !Interface:
!
module m_cells
!
! !Uses:
!
  use m_kind
  use m_simtype
  implicit none
!
! !Public member functions:
!
  public :: r_neighs
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
!eop
!-------------------------------------------------------------------------
  contains
!    
!bop
!
! !Iroutine:  r_neighs
!
! !Description: This routine builds a list of neighbors cells as function 
!of total number particles.
!
  subroutine r_neighs( y )
!
! !Input parameters:
!
  type(simulation), intent(inout) :: y
!
! !REVISION HISTORY:
!  06 August 2015 by Asdrubal Lozada
!
! !Remarks: Decreases in performance can be produced by reducing
!of total number particles. Is advised in this case use a asintotic
!$\Theta(n^{2})$ algorithm.
! 
!eop
!-----------------------------------------------------------------------
!boc
  if(y%m_mxmol < 27) then
    y%m_viz = y%m_mxmol
  else
    y%m_viz = 27
  end if
!eoc

 end subroutine r_neighs 
end module m_cells

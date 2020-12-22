!----------------------------------------------------------------------
!          MELQUIADES: Metropolis Monte Carlo Program
!----------------------------------------------------------------------
!bop
!
! !Module:  m_pbc
!
! !Description: This module contains a routine to
!define the periodical boundary conditions.
!\\
!\\
! !Interface:
!
module m_pbc
!
! !Uses:
  use m_kind
  use m_simtype
  use m_boxtype
  implicit none
!
! !Public member functions:
!
  private
  public :: r_pbc
!
! !Revision history:
!  06Aug 2015 Asdrubal Lozada
!
!eop
!---------------------------------------------------------------------
  contains
!bop
!
! !Iroutine: r_pbc
!
! !Description: This routine implements the periodic boundary
!conditions criterium and minimum imagen convention.
!\\
!\\
! !Interface:
  subroutine r_pbc( com, edge )
  implicit none
!
! !Input parameters:   
  real(rkind), dimension(:), intent(inout) :: com ! Center of mass position
  real(rkind), dimension(:), intent(in)    :: edge ! Lenght of edge box
!
! !Revision history:
!  06Aug 2015 Asdrubal Lozada
!
!eop
!------------------------------------------------------------------
!
!boc
  com(1) = com(1) - edge(1) * dnint(com(1) / edge(1))
  com(2) = com(2) - edge(2) * dnint(com(2) / edge(2))
  com(3) = com(3) - edge(3) * dnint(com(3) / edge(3))
!eoc

  end subroutine r_pbc
end module m_pbc

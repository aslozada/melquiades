!----------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program
!----------------------------------------------------------------------
!bop
!
! !Module: m_shift
!
! !Description: This module contains a routine for volume shifted!\\
!\\
! !Interface:
!
  module m_shift
!  
! !Uses
  use m_kind
  use m_boxtype
  use m_simtype
  use m_unit
!
! !Public member functions:
  private
  public :: r_shift
!
! ! Revision history:
! ! 12Nov 2015 Asdrubal Lozada
!
!eop
!-----------------------------------------------
  contains
!
!bop
!
! !Iroutine: r_shift
!
! !Description: Shifted volume
!\\
!\\
! !Interface:
  subroutine r_shift( edge, iter, y, x )
  implicit none
! Dummy arguments
  type(box), pointer :: x
  type(simulation), intent(inout) :: y
  real(rkind) :: dv, vol, edgen
  real(rkind), dimension(:), intent(inout) :: edge
  integer, intent(inout) :: iter

 

    if( y%m_drho > 0.0_rkind .and. y%m_drho <= 1.0_rkind ) then

     iter = iter - 0.05

     dv = y%m_vol * (1.0_rkind - y%m_drho ) / y%m_drho
     vol = y%m_vol + iter * dv
     edgen = dlog(vol) / 3.0_rkind
     edgen = dexp(edgen)

     edge(1) = edgen
     edge(2) = edgen
     edge(3) = edgen

     x%m_edge(1) = edge(1)
     x%m_edge(2) = edge(1)
     x%m_edge(3) = edge(1)

     y%m_dens = real(y%m_mxmol,rkind) * (1.0_rkind/vol)
     y%m_voli = 1.0_rkind/vol

      write(*,*)'edge var: ', y%m_drho, dv, edge(1) 


   else 
    
     continue

   end if  


  end subroutine r_shift


  end module m_shift

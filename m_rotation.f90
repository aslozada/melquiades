!----------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program
!----------------------------------------------------------------------
!bop
!
! !Module: m_rotation
!
! !Description: This module contains a routine that applies a transformation
!matrix over cartesian coordinates in the molecules.
!\\
!\\
! !Interface:
!
  module m_rotation
!
! !Uses:
    use m_kind
    use m_simtype
    use m_boxtype
    use m_constants
    use m_random
!
! !Public member functions:
!
  private
  public :: r_euler
!
! !Revision history
! 08Aug 2015 Asdrubal Lozada
!
!eop
!----------------------------------------------------------------------
  contains
!
!bop
!
! !Iroutine: r_euler
!
! !Description: This routine performs random rotation about space-fixed 
!axes. The method is based in the Euler's rotation theorem.
!\\
!\\
! !Interface:
  subroutine r_euler( j, y, x )
!    
  implicit none
!
! !Input parameters:
    type(simulation), intent(inout) :: y
    type(box), pointer :: x    
    integer, intent(in) :: j
!
! !Revision history:
! 08Aug 2015 Asdrubal Lozada
!
!eop
!----------------------------------------------------------------------
! Local variables
  integer :: iaxis
  real(rkind) :: angle, dgamax, dgamma
  real(rkind) :: cosdg, sindg
  real(rkind) :: exiold, eyiold, eziold
  real(rkind) :: exinew, eyinew, ezinew
  integer :: k
!  integer,parameter :: iseed = 86456

  angle = y%m_rotats
  iaxis = int( 3.0_rkind * m_rand()) + 1
  dgamax = pi * angle / 180.0_rkind
  dgamma = (2.0_rkind * m_rand() - 1.0) * dgamax

  cosdg = dcos(dgamma)
  sindg = dsin(dgamma)

  do k = 1, x%m_ns(j)

   exiold = x%m_rot(1,k)
   eyiold = x%m_rot(2,k)
   eziold = x%m_rot(3,k)

   if( iaxis == 1 ) then
     exinew = exiold
     eyinew = cosdg * eyiold + sindg * eziold
     ezinew = cosdg * eziold  -  sindg * eyiold
   else if( iaxis == 2 ) then
     exinew = cosdg * exiold - sindg * eziold
     eyinew = eyiold
     ezinew = cosdg * eziold + sindg * exiold
   else if(iaxis == 3) then
     exinew = cosdg * exiold + sindg * eyiold
     eyinew   = cosdg * eyiold -  sindg * exiold
     ezinew = eziold
   end if

   x%m_rot(1,k) = exinew
   x%m_rot(2,k) = eyinew
   x%m_rot(3,k) = ezinew

  end do 

  end subroutine r_euler
end module m_rotation

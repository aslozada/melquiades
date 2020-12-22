!---------------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                       !
!---------------------------------------------------------------------------
!bop
!
! !Module:  m_error
!
! !Description: This module contains routines for evaluating failures
!while reading files and the management memory.
!\\
!\\
! !Interface:
!
module m_error
!
  implicit none
!
! !Public member functions:
!
  public :: r_ealloc
  public :: r_eread
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
!eop
!--------------------------------------------------------------------------
!
  contains
!
!bop
!
! !Iroutine: r_ealloc
!
! !Description: Returns a error message if management memory fails.
!\\
!\\
! !Interface:
!
  subroutine r_ealloc( var_name, error )
!
  implicit none
!
! !Input parameters:
!
  character(len=*), intent(in) :: var_name
  integer, intent(in) :: error
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
!eop
!---------------------------------------------------------------------------
! Local variables
!boc
  character(len=33) :: msg
  msg = "Space requested not possible for "

  if( error /= 0 ) then
    write(*,'(a)') msg//" "//trim(var_name)
    stop
  end if
!eoc

  end subroutine r_ealloc
  
!bop
!
! !Iroutine: r_eread
!
! !Description: Returns a error message if the reading fails.
!\\
!\\
! !Interface:
!
  subroutine r_eread( ios, arg )
!
  implicit none
!
! !Input parameters:
!
  integer, intent(in) :: ios
  character(len=*), intent(in) :: arg
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
!eop
!---------------------------------------------------------------------------
  if(ios>0) then
    write(*,'("*** Check the input file ***")')
    write(*,*) 'Error at line-argument ', trim(arg)
    write(*,'("***                      ***")')
    stop
  end if

  end subroutine r_eread

end module m_error

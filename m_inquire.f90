!---------------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                       !
!---------------------------------------------------------------------------
!bop
!
! !Module:  m_inquire
!
! !Description: This module contains a routine for program initialization
! from command line and inquire by input files.
!\\
!\\
! !Interface:
!
module m_inquire
!
! !Uses:
!
  use m_unit, only : infile, m_inpfile
  implicit none
!
! !Public member functions:
!
  public :: r_comline
  public :: r_request
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
! !Iroutine: r_comline
!
! !Description: Read arguments from command line (the input file)
!\\
!\\
! !Interface:
!
  subroutine r_comline
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
!eop
!--------------------------------------------------------------------------
  implicit none
! Local variables:

  integer, parameter :: iarg = 1
  character(len=40) :: arg 
!boc
  if( iargc() == 1 ) then
    call getarg( iarg, arg)

    m_inpfile  =  arg   
    call r_request( infile, m_inpfile )
!eoc  
  else
    write(*,'("MELQUIADES requires a input file : melquiades.x input.dat")')
    stop
  end if

  end subroutine r_comline
!
!--------------------------------------------------------------------------------
!bop
!
! !Iroutine:  r_request
!
! !Description: This routine request a name of input file and returns a error
! message if this not there.
!\\
!\\
! !Interface:
!
  subroutine r_request( iunit, file_name )
  implicit none
!
! !Input parameters:
!
  integer, intent(in) :: iunit ! Unit number
  character(len=40), intent(in) :: file_name ! name input file
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
!eop
!-----------------------------------------------------------------------------------
!  Local variables
  logical :: lopen, lexists
  character(len=40) :: line
!boc
  inquire(unit=iunit,opened=lopen)
  if( .not. lopen ) then
    line = file_name(:len_trim(file_name))
    inquire(file=line,exist=lexists)
!eoc    
    if(lexists) then
      open(unit=iunit, file=line, status='old', form='formatted')
      rewind(unit=iunit)
    else
      write(*,*)'File ', line(:len_trim(line)), ' does not exist and&
                    & MELQUIADES is unable to continue.'
      stop
      rewind(iunit)
      return
    endif ! lexist
  end if  ! lopen

  end subroutine r_request

end module m_inquire

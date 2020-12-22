!----------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                  !
!----------------------------------------------------------------------
!bop
!
! !Module: m_store
!
! !Description: This modules contains a routine for writting:
!1) A like-formatted simulation box file; 
!2) A file in conventional xyz file format.
!\\
!\\
! !Interface: 
!
  module m_store
!
! !Uses:
  use m_kind
  use m_simtype
  use m_boxtype
  use m_unit
  
!  
! !Public member functions:
!
  private

  public :: r_saves
  public :: r_newbox
!
! !Revision history:
! 08Aug 2015 Asdrubal Lozada
!
!eop
!----------------------------------------------------------------------
  contains
!
!bop
!
! !Iroutine:  r_saves
!
! !Description: This routine made a format xyz file
!\\
!\\
! !Interface:
  subroutine r_saves(energy, virial, y, x )
!
  implicit none
!
! !Input parameters:
  type(simulation), intent(inout) :: y
  type(box), pointer :: x
  real(rkind), intent(inout) :: energy
  real(rkind), intent(inout) :: virial
!
! !Revision history
! 08Aug 2015 Asdrubal Lozada
!
!eop
!----------------------------------------------------------------------
! Local variables
  integer :: j, i
  character(len=2) :: label
  real(rkind) :: xx,yy,zz
  character(len=10) :: ival, rkval
  character(len=10) :: caption
   
  caption = 'Energy = '

  write(ival,'(i10)') y%m_mxatms  
  write(rkval,'(f10.1)') energy
  write(ixyz,'(a10)') adjustl(ival)
  write(ixyz,'(2a10)') caption, adjustl(rkval)

  do j = 1, y%m_mxmol
   do i = 1, x%m_ns(j)
    label =  x%m_symbol(i,j)
     xx = x%m_site(1,i,j) + x%m_cmass(1,j)
     yy = x%m_site(2,i,j) + x%m_cmass(2,j)
     zz = x%m_site(3,i,j) + x%m_cmass(3,j)

     write(ixyz,'(a2,2x,f10.5,2x,f10.5,2x,f10.5)')label,xx,yy,zz

    end do
  end do
  end subroutine r_saves
!
!bop
!
! !Iroutine: r_newbox
!
! !Description: This routine update boxfile at the end of Markov chain.
!\\
!\\
!
! !Interface:
  subroutine r_newbox( y, x )
!
! !Input parameters:
   type(simulation), intent(inout) :: y
   type(box), pointer :: x

!
! !Revision history:
! 08Aug 2015 Asdrubal Lozada 
!
!eop
!----------------------------------------------------------------------
! Local Variable    
  integer :: i, j
  character(len=40) :: line
  
  line = trim(y%m_boxfile)  
  open(inbox, file=line, status='old')

  rewind(inbox)
  write(inbox,*) y%m_ntf
  write(inbox,*) x%m_nmol(:), x%m_nsite(:), x%m_mass(:), x%m_edge(:)
  write(inbox,*) m_seed
  
  do j = 1, y%m_mxmol
   write(inbox,*) trim(x%m_molname(j)), x%m_idtype(j), x%m_cmass(:,j)
     do i = 1, x%m_ns(j)
      write(inbox,*) x%m_symbol(i,j), x%m_idpar(i,j), x%m_site(:,i,j)
     end do
  end do
    
  end subroutine r_newbox
end module m_store

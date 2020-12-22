!---------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlos Program
!---------------------------------------------------------------------
!bop
!
! !Module: m_uplist
!
! !Description: This module contains the routines needed for updating
!list of neighbors in box simulation.
!\\
!\\
! !Interface:
!
  module m_uplist
! 
! !Uses:
!
  use m_kind
  use m_simtype
  use m_boxtype 
  use m_precells
!
! !Public member functions: 
!
  private
  public :: r_update
  public :: r_up2date
  public :: r_up3date
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
! !Iroutine: r_update
!
! !Description: This routine updates the list of neighbors if cell index 
!was changed.
!\\
!\\
! !Interface:
  subroutine r_update( com_new, i_d, a, y, x )
!
  implicit none
!
! !Input parameters:
  type(simulation), intent(inout) :: y
  type(box), pointer :: x
  integer,intent(in) :: a
  integer,intent(in) :: i_d
  real(rkind), dimension(:),intent(in) :: com_new
!
! !Revision history:
! 08Aug 2015 Asdrubal Lozada
!
!eop
!----------------------------------------------------------------------
! Local variables
  real(rkind), dimension(3) :: com_old
  integer :: jold, jnew, i

  if( a == 1 ) then
    com_old(1) = x%m_cmass(1,i_d)
    com_old(2) = x%m_cmass(2,i_d)
    com_old(3) = x%m_cmass(3,i_d)

    x%m_cmass(1,i_d) = com_new(1)
    x%m_cmass(2,i_d) = com_new(2)
    x%m_cmass(3,i_d) = com_new(3)

    do i=1,x%m_ns(i_d)
     x%m_site(1,i,i_d) = x%m_rot(1,i)
     x%m_site(2,i,i_d) = x%m_rot(2,i)
     x%m_site(3,i,i_d) = x%m_rot(3,i)
    end do

    jold = f_idcell(com_old,y,x)
    jnew = f_idcell(com_new,y,x)
    
    if( jnew /= jold ) then   
      call r_headinit(y,x)
      call r_linkscell(y,x)
    end if ! index cell

   else
     x%m_cmass(1,i_d) = x%m_rcom(1,i_d)
     x%m_cmass(2,i_d) = x%m_rcom(2,i_d)
     x%m_cmass(3,i_d) = x%m_rcom(3,i_d)
   
    do i=1,x%m_ns(i_d)
     x%m_site(1,i,i_d) = x%m_rsite(1,i)
     x%m_site(2,i,i_d) = x%m_rsite(2,i)
     x%m_site(3,i,i_d) = x%m_rsite(3,i)
    end do
  end if ! a

  end subroutine r_update
  
  subroutine r_up2date( com_new, i_d, a, y, x )
  implicit none

  type(simulation), intent(inout) :: y
  type(box), pointer :: x
  integer,intent(in) :: a
  integer,intent(in) :: i_d
  real(rkind), dimension(:),intent(in) :: com_new

! Local variables  
  real(rkind), dimension(3) :: com_old
  integer :: jold, jnew, i

  if( a == 1 ) then
    com_old(1) = x%m_cmass(1,i_d)
    com_old(2) = x%m_cmass(2,i_d)
    com_old(3) = x%m_cmass(3,i_d)

    x%m_cmass(1,i_d) = com_new(1)
    x%m_cmass(2,i_d) = com_new(2)
    x%m_cmass(3,i_d) = com_new(3)

    do i=1,x%m_ns(i_d)
     x%m_site(1,i,i_d) = x%m_rot(1,i)
     x%m_site(2,i,i_d) = x%m_rot(2,i)
     x%m_site(3,i,i_d) = x%m_rot(3,i)
    end do

    jold = f_idcell(com_old,y,x)
    jnew = f_idcell(com_new,y,x)
 
    if( jnew /= jold ) then
      call r_headinit2(y,x)
      call r_linkscell2(y,x)
   end if

 else
   x%m_cmass(1,i_d) = x%m_rcom(1,i_d)
   x%m_cmass(2,i_d) = x%m_rcom(2,i_d)
   x%m_cmass(3,i_d) = x%m_rcom(3,i_d)
  
   do i=1,x%m_ns(i_d)
    x%m_site(1,i,i_d) = x%m_rsite(1,i)
    x%m_site(2,i,i_d) = x%m_rsite(2,i)
    x%m_site(3,i,i_d) = x%m_rsite(3,i)
   end do

 end if ! a

  y%m_solute = .true.

 end subroutine r_up2date

 subroutine r_up3date( com_new, i_d, a, y, x )
  implicit none

  type(simulation), intent(inout) :: y
  type(box), pointer :: x
  integer,intent(in) :: a
  integer,intent(in) :: i_d
  real(rkind), dimension(:),intent(in) :: com_new
! Local variable  
  real(rkind), dimension(3) :: com_old
  integer :: jold, jnew, i

  if( a == 1 ) then
    com_old(1) = x%m_cmass(1,i_d)
    com_old(2) = x%m_cmass(2,i_d)
    com_old(3) = x%m_cmass(3,i_d)

    x%m_cmass(1,i_d) = com_new(1)
    x%m_cmass(2,i_d) = com_new(2)
    x%m_cmass(3,i_d) = com_new(3)

    do i=1,x%m_ns(i_d)
     x%m_site(1,i,i_d) = x%m_rot(1,i)
     x%m_site(2,i,i_d) = x%m_rot(2,i)
     x%m_site(3,i,i_d) = x%m_rot(3,i)
    end do

    jold = f_idcell(com_old,y,x)
    jnew = f_idcell(com_new,y,x)
 
    if( jnew /= jold ) then
      call r_headinit3(y,x)
      call r_linkscell3(y,x)
    end if

  else
    x%m_cmass(1,i_d) = x%m_rcom(1,i_d)
    x%m_cmass(2,i_d) = x%m_rcom(2,i_d)
    x%m_cmass(3,i_d) = x%m_rcom(3,i_d)

   do i=1,x%m_ns(i_d)
    x%m_site(1,i,i_d) = x%m_rsite(1,i)
    x%m_site(2,i,i_d) = x%m_rsite(2,i)
    x%m_site(3,i,i_d) = x%m_rsite(3,i)
   end do

 end if 

  y%m_solute = .true.

 end subroutine r_up3date

end module m_uplist

 

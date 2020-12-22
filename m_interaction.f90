!----------------------------------------------------------------------
!         MELQUIADES : Metropolis Monte Carlo Program
!----------------------------------------------------------------------
!bop
!
! !Module:  m_interaction
!
! !Description: This module contains a routine for
! interaction energy calculation.
!\\
!\\
! !Interface:
!
  module m_interaction
!
! !Uses:
  use m_kind
  use m_simtype
  use m_boxtype
  use m_precells
  use m_constants, only :  c_kcalmol, pi
  use m_zeros
  use m_pbc
  use m_ljones
!
! !Public member functions:
!
  private
!  
  public :: r_interactions
!
! !Revision history:
! 06Aug 2015  Asdrubal Lozada
!
!eop
!----------------------------------------------------------------------
  contains
!
!bop
!
! !Iroutine: r_interactions
!
! !Description: This routine uses a neighbors list for intermolecular
!energy calculation. Through a scheme of case selection different types
!of potential can be used. In the current version this potential functions
!are: The Lennard-Jones, Buckingham, Yukawa  and Coulombic potentials.
!
!\\
!\\
! !Interface:
 subroutine r_interactions( com, edge, i, jb, ids, energy, virial, t, y, x )
!
  implicit none
!--------------------
!
! !Input arguments:

  type(simulation), intent(inout) :: y
  type(box), pointer :: x
  type(temporary), intent(inout) :: t
  real(rkind), dimension(:), intent(in)  :: com
  real(rkind), dimension(:), intent(in)  :: edge
  real(rkind), intent(inout) :: energy, virial
  integer, intent(in) :: i, jb, ids 
!-------------------------------------
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
!eop
!----------------------------------------------------------------
!  Local variables
  real(rkind) :: rix, riy, riz
  real(rkind) :: rjx, rjy, rjz
  real(rkind), dimension(3) :: rij

  real(rkind), dimension(6) :: v
!---------------------------------  
  integer :: icell
  integer :: k, j, iter
  real(rkind) :: rijsq, rm, rp

   energy = 0.0_rkind
   virial = 0.0_rkind
   v = 0.0_rkind

   rix = com(1)
   riy = com(2) 
   riz = com(3)

   icell = f_idcell( com, y, x )  

   call r_idneighs( icell, y, x )  

   do k = 1, y%m_viz
  
    if( .not. y%m_solute ) then
     iter = x%m_ncell(k)
     j = x%m_head(iter)
   else
     iter = x%m_ncelsa(k)
     j = x%m_hedsb(iter)
   end if 

   do while(j /= 0)

    if( j /= i .and. j >= jb ) then


      rjx = x%m_cmass(1,j); rjy = x%m_cmass(2,j); rjz = x%m_cmass(3,j)

      rij(1) = rix - rjx
      rij(2) = riy - rjy
      rij(3) = riz - rjz
    
   if( y%m_pbcs ) then 
    call r_pbc(rij,edge)
   end if

   rijsq = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
   rm = rijsq 

   if( .not.y%m_solute ) then
     rp = y%m_rpair
   else
     rp = y%m_rsol
   end if

   if( rm<=rp ) then
!     
!boc       
   select case(y%m_potens)
   case('m_ljc')  

     call p_ljones(rij, v, i, j, ids, t, x) 
!!   case('m_bc')    
!!     call p_ljones(rij, v, i, j, ids, t, x) 
!!     call p_buckhm(rij, v, i, j, ids, t, x) 
!!   case('m_ljy')
!!     call p_ljones(rij, v, i, j, ids, t, x) 
!!     call p_yukawa(rij, v, i, j, ids, t, x) 
   case default
    write(*,*) 'Test without potential'
  end select
!eoc
 end if
end if 

   if( .not.y%m_solute ) then
    j = x%m_list(j)
   else
    j = x%m_listb(j)
   end if 

   end do  ! While
 end do  ! k

  energy  = v(1) + v(2) + v(3)
  virial  = v(4) + v(5) + v(6)

 end subroutine r_interactions 

end module m_interaction

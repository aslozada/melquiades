!----------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program
!----------------------------------------------------------------------
!bop
!
! !Module: m_isobaric
!
! !Description: This module contains a routine for calculating configurational
!energy in the isobaric-isothermal ensemble.
!\\
!\\
! !Interface:
!
module m_isobaric
  use m_kind
  use m_simtype
  use m_boxtype
  use m_random
  use m_constants
  use m_zeros
  use m_configuration
! 
  implicit none
!
! !Public member functions:
!
  private
!  
  public :: r_npt
!
! !Revision history:
! 10Aug 2015 Asdrubal Lozada
!
!eop
!----------------------------------------------------------------------
  contains
!
!bop
!
! !Iroutine: r_npt
!
! !Description: In this routine random volume changes are performed.
!\\
!\\
! !Interface:
!
  subroutine r_npt( enconf, viconf, apv, t, y, x )
!
  implicit none
!
! !Input parameters

   type(simulation), intent(inout) :: y
   type(box), pointer :: x
   type(temporary), intent(inout) :: t
!   
! !Ouptup parameters:   

   real(rkind), intent(inout) :: enconf  
   real(rkind), intent(inout) :: viconf
   real(rkind), intent(inout) :: apv
!
! !Revision history:
! 10Aug 2015 Asdrubal Lozada
!
!eop
!----------------------------------------------------------------------
!  Local variables
   real(rkind) :: vol_new
   real(rkind) :: ivol_new
   real(rkind) :: del_vol
!--------------------------------
!  Save first values
   real(rkind) :: edge_old
   real(rkind) :: edge_new
   real(rkind) :: redge
!----------------------------------
   real(rkind) :: ener_old
   real(rkind) :: ener_new
   real(rkind) :: virl_old
   real(rkind) :: virl_new
   real(rkind) :: en_dvol
   real(rkind) :: vi_dvol
   real(rkind) :: dpv
   real(rkind) :: dvl 
   real(rkind) :: delthb   
!-----------------------------
   real(rkind), dimension(3) :: edgen
!---------------------------------------------
   integer     :: accept   
   integer     :: i
   integer     :: istate  
   integer     :: vmove
   integer     :: nk
!-----------------------------------------------
   accept = 0
   apv = 0.0_rkind
   vmove = 0

! Save coordiantes of center of mass

   do i = 1, y%m_mxmol
    x%m_svcom(1,i) = x%m_cmass(1,i)
    x%m_svcom(2,i) = x%m_cmass(2,i)
    x%m_svcom(3,i) = x%m_cmass(3,i)
   end do   

!boc
   vol_new = y%m_vol + (2.0_rkind * m_rand() - 1.0_rkind) * y%m_maxvol
!eoc

   edge_old = dlog(y%m_vol) / 3.0_rkind  
   edge_old = dexp(edge_old) 
  
   edge_new = dlog(vol_new) / 3.0_rkind
   edge_new = dexp(edge_new)
 
   redge = edge_new / edge_old  

   edgen(1) = x%m_edge(1) * redge
   edgen(2) = x%m_edge(2) * redge
   edgen(3) = x%m_edge(3) * redge

   vol_new  = edgen(1) * edgen(2) * edgen(3)
   y%m_voli = 1.0_rkind / vol_new
   del_vol  = vol_new - y%m_vol

   ener_old = enconf
   virl_old = viconf
!--------------------------------

  do i = 1, y%m_mxmol
   x%m_cmass(1,i) = x%m_cmass(1,i) * redge
   x%m_cmass(2,i) = x%m_cmass(2,i) * redge
   x%m_cmass(3,i) = x%m_cmass(3,i) * redge
  end do


 call r_configuration(edgen,ener_new, virl_new, t, y, x, istate, .false.)

  en_dvol = ener_new - ener_old

   
   dpv = pv_conv * del_vol
   dvl = real(y%m_mxmol, rkind) * dlog(vol_new / y%m_vol)
   delthb = y%m_beta * (en_dvol + dpv) - dvl


 if( delthb < 10.0_rkind ) then  
   if( delthb <= 0.0_rkind ) then      
     accept = 1
     vmove = vmove + 1
   else if( exp(-delthb) > m_rand() ) then
     accept = 1
     vmove = vmove + 1
   end if
 end if

 if( accept == 1 ) then
   enconf  = ener_new

 end if

!------------------------------------------------
! Update/Restart coordinates in box simulation
 if( accept == 1 ) then
   x%m_edge(1) = edgen(1)
   x%m_edge(2) = edgen(2)
   x%m_edge(3) = edgen(3)

   apv = dpv
   y%m_vol = vol_new
   y%m_dens = real(y%m_mxmol,rkind) * y%m_voli
 else if( accept == 0 ) then
 
   y%m_voli = 1.0_rkind / y%m_vol 

   do i = 1, y%m_mxmol  
    x%m_cmass(1,i) = x%m_svcom(1,i)
    x%m_cmass(2,i) = x%m_svcom(2,i)
    x%m_cmass(3,i) = x%m_svcom(3,i)
   end do
 end if

   t%m_amv = t%m_amv + 1

  end subroutine r_npt
 end module m_isobaric

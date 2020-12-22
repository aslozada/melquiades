!---------------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                       !
!---------------------------------------------------------------------------
!bop
!
! !Module:  m_init
!
! !Description: This module contains routines for initializing the variables
!and activating the Linked-Cell List 
!\\
!\\
! !Interface:
!
module m_init
!
! !Uses:
!
  use m_kind
  use m_simtype
  use m_boxtype
  use m_constants
  use m_cutoffs
  use m_cells
  use m_precells
  use m_zeros
  use m_longs
  use m_unit ! Unit to average file
  implicit none
!
! !Public member functions:
!
  public :: r_thermo
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
!eop
!----------------------------------------------------------------------
  contains
!
!bop
!
! !Iroutine: r_thermo
!
! !Description: This routine modifies the thermodynamic
!parameters to an appropiate scale and initializes the Linked-Cells List.
!\\
!\\
! !Interface:
!
  subroutine r_thermo( edge, t, y, x)
  implicit none
!
! !Input parameters:
!
  type(simulation), intent(inout) :: y
  type(box), pointer :: x
  type(temporary), intent(inout) :: t
  real(rkind), dimension(:), intent(out) :: edge
  real(rkind) :: one
  character(len=30) :: rkval
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
! 02Nov 2015 Asdrubal Lozada
!
!eop
!----------------------------------------------------------------------
!   Local variables
  integer :: i

  y%m_vol   = x%m_edge(1) * x%m_edge(2) * x%m_edge(3) 
  y%m_voli = 1.0_rkind / y%m_vol
  y%m_dens  = real(y%m_mxmol,rkind) * y%m_voli 
  
  
  y%m_weight = 0.0_rkind

  do i = 1, y%m_ntf   
  y%m_weight = y%m_weight + real(x%m_mass(i)*x%m_nmol(i),rkind)
  end do

  y%m_densi = (y%m_weight * y%m_voli) * d_con  

!
 write(rkval,'(e20.14)') y%m_densi
 write(*,'(" Initial density                : ",a20," g/cm**3 ")') adjustl(rkval)
!
  
!
!boc
  y%m_temper = y%m_temper + c_k   
  y%m_beta = 1.0_rkind / (kB * y%m_temper) 
  y%m_ppcon = y%m_press * pv_conv
!eoc

  edge(1) = x%m_edge(1)
  edge(2) = x%m_edge(2)
  edge(3) = x%m_edge(3)

  call r_cuts( y, x )

  call r_long ( y, x )

  call r_neighs( y )
  call r_intcell( y, x )


  if(.not.y%m_solute) then
    call r_headinit( y, x )  
    call r_linkscell( y, x )
 end if

  if( y%m_solute ) then
    call r_headinit2( y, x )
    call r_linkscell2( y, x )
    call r_headinit3( y, x )
    call r_linkscell3( y, x )
  end if
  
  call r_cero( t )

!-------------  
    ! Open file to average
  if( y%m_averag ) then 
   open(inaver, file='average', status='unknown', form='formatted')
    read(inaver,'(a)')
    read(inaver,*) y%m_eold
   close(inaver)
  end if

! Open file to correlation
  if( y%m_corr ) then
    open(incorr,file='correlation',status='unknown', form='formatted')
     read(incorr,'(a)')
     read(incorr,*) y%m_varold
    close(incorr)
  end if  
!-------------

 end subroutine r_thermo
 
end module m_init

!----------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program
!----------------------------------------------------------------------
!bop
!
! !Module: m_sampler
!
! !Description: This module contains a routine for writting output values.
!\\
!\\
! !Interface:
!
  module m_sampler
!
! !Uses:
  use m_kind
  use m_simtype
  use m_boxtype
  use m_zeros
  use m_constants
  use m_configuration
!
! !Public member functions:
!
  private 
  public :: r_samples
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
! !Iroutine: r_samples
!
! !Description: This routine saves some values temporarily.
!\\
!\\
! !Interface:
  subroutine r_samples( iunit, i, iter, engconf, virconf, enp, rhos, press, apv, t, y )
!
  implicit none
!
! !Input parameters:
  type(simulation), intent(inout) :: y
  type(temporary), intent(inout) :: t
  integer, intent(in) :: iunit
  integer, intent(in) :: i
!
! !Output parameters:
  real(rkind), intent(inout) :: engconf
  real(rkind), intent(inout) :: virconf
  real(rkind), intent(out) :: enp, press, rhos
  real(rkind), intent(in) :: apv
  

!
! !Local variables
   integer :: kfreq  ! Frequency to sampling correlation 
   integer, intent(inout) :: iter
   integer :: store
!
! !Revision history:
! 08Aug 2015 Asdrubal Lozada
!
!eop
!----------------------------------------------------------------------


  if( y%m_mxmol /=0 ) then

  enp   =  (engconf) / real(y%m_mxmol, rkind) ! Energy per molecule
  press =  y%m_dens * y%m_temper + virconf/y%m_vol  ! Pressure
  rhos  =  y%m_weight * y%m_voli * d_con


!-----------------------

! Accumulative volume

  t%m_avol = t%m_avol + y%m_vol
  t%m_avos = t%m_avos + y%m_vol * y%m_vol
 
! Accumulative energy per molecule  

  t%m_aeng = t%m_aeng + enp
  t%m_aens = t%m_aens + enp * enp

!  Accumulative total energy

  t%m_aeto = t%m_aeto + engconf
  t%m_aets = t%m_aets + engconf * engconf

! Accumulative enthalpy

  t%m_aent = t%m_aent + engconf + apv
  t%m_aehs = t%m_aehs + (engconf + apv) * (engconf + apv)
  t%m_aalp = t%m_aalp + y%m_vol * ( engconf + apv )

  
! Accumulative density
  
  t%m_arho = t%m_arho + y%m_dens
  t%m_arhs = t%m_arhs + y%m_dens * y%m_dens


  else
  
  enp = 0.0_rkind
  press = 0.0_rkind

  end if    

  if( y%m_averag ) then
      write(iunit,*) i, enp
  end if    
 

!!!  if(.not.y%m_averag ) then    
!-----------
! Define frequency to sampling
!!!  if(y%m_nsteps > 5 .and. y%m_nsteps < 100 ) then
!!!    kfreq = int(y%m_nsteps / 10 )
!!!  elseif(y%m_nsteps > 100 ) then
!!!    kfreq = int(y%m_nsteps / 2000)
!!!  else
!!!    kfreq = 5
!!!  end if  
!-----------

!!!   if( mod(i,kfreq) == 0) then
!!   if( i == y%m_nsteps) then

 !!!     iter = iter + 1
 !!!     store = (kfreq*int(i/kfreq))

 !!!     write(iunit,*) enp, store
 !!!  end if   
  
 !!! end if

  end subroutine r_samples
 
 end module m_sampler

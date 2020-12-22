  !----------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program
!----------------------------------------------------------------------
!bop
!
! !Module:  m_averages
!
! !Description: This module contains a routine for averages calculation.
!\\
!\\
! !Interface:
!
module m_averages
!
! !Uses:
!
  use m_kind
  use m_simtype
  use m_boxtype
  use m_constants
  use m_unit
  use m_zeros
  use m_fluctuations
  use m_configuration
  implicit none
!
! !Public member functions:
!
  private 
  public :: r_averages
  public :: r_rnormal
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
! !Iroutine:  r_averages
!
! !Description: This routine calculates the average of quantities if
!the input flag "average" have a true value.
!\\
!\\
! !Interface:
!
  subroutine r_averages( t, y, x )
!    
  implicit none
!
! !Input parameters:
  type(simulation), intent(inout) :: y
  type(temporary), intent(inout) :: t
  type(box), pointer :: x
!
! !Revision history:
! 10Aug 2015 Asdrubal Lozada
!
!eop
!----------------------------------------------------------------------
!  Local variables

   character(len=40) :: title
   real(rkind) :: tconf, acmold
   real(rkind) :: etold, etsqo
   real(rkind) :: esvold, esvqol
   real(rkind) :: denold, dsqold
   real(rkind) :: achold, achsqo
   real(rkind) :: volold, vsqold, alfold
   real(rkind) :: avd, avdsq, avd2
   real(rkind) :: avh, avhsq
   real(rkind) :: c1, c2, c3
   real(rkind) :: avol, avolsq, volrad

   real(rkind) :: aeng, aengsq
   real(rkind) :: aetot, aetotsq
   real(rkind) :: cp, kappa, alpha
   real(rkind) :: cpt, kappat, alphat
   real(rkind) :: cvap, avmol

   real(rkind) :: tmp1, tmp2

   integer :: nk, kk, ll
   integer :: n, m, iter, error
   
    if( y%m_averag) then

     open(inaver, file='average', status='unknown', form='formatted')
      
      rewind(inaver)
      
      read(inaver,'(a)') title
      read(inaver, *) acmold
      read(inaver, *) etold, etsqo
      read(inaver, *) esvold, esvqol
      read(inaver, *) denold, dsqold
      read(inaver, *) achold, achsqo
      read(inaver, *) volold, vsqold
      read(inaver, *) alfold

   else 

   acmold = 0.0_rkind    
   etold  = 0.0_rkind
   etsqo  = 0.0_rkind
   esvold = 0.0_rkind
   esvqol = 0.0_rkind
   denold = 0.0_rkind
   dsqold = 0.0_rkind
   achold = 0.0_rkind
   achsqo = 0.0_rkind
   volold = 0.0_rkind
   vsqold = 0.0_rkind
   alfold = 0.0_rkind
   end if    

   
    tconf = real(t%m_amv, rkind) + acmold
!-----------------
! Zero fluctuations
  if( y%m_averag ) then

   x%m_sehpy = 0.0_rkind
   x%m_serho = 0.0_rkind
   x%m_sevol = 0.0_rkind
   x%m_sete  = 0.0_rkind

  end if    

  call r_fluctuations( t, y, x )
!-----------------
    aeng    = t%m_aeng / real(t%m_amv, rkind)
    aengsq  = t%m_aens / real(t%m_amv, rkind)
    
    aetot   = t%m_aeto / real(t%m_amv, rkind)
    aetotsq = t%m_aets / real(t%m_amv, rkind)

!---------------------------------------------- 
   write(*,*)"--------------------------------------------------------------------------------"
   write(*,*)"                                Outputs                                         "
   write(*,*)"--------------------------------------------------------------------------------"
   write(*,'(" Total energy in this run      :         ",e20.10," kcal/mol ")') aetot 
   write(*,'(" Total energy per molecule     :         ",e20.10," kcal/mol ")') aeng
!---------------
   write(*,'(" Attempted to move molecules   :   ", i20)') y%m_nsteps 
   write(*,'(" Accepted moves molecules      :   ", i20)') t%m_imove
   write(*,*)"--------------------------------------------------------------------------------"

   cp = (t%m_aens / real(t%m_amv, rkind)) - ((t%m_aent / real(t%m_amv, rkind))**2)
   cp = cp / (kB * y%m_temper * y%m_temper * real(y%m_mxmol))
   cp = cp * 1000.0_rkind

   kappa = (t%m_avos / real(t%m_amv, rkind)) - ((t%m_avol / real(t%m_amv, rkind))**2)
   kappa = kappa / (kB * y%m_temper * (t%m_avol/real(t%m_amv,rkind)))
   kappa = kappa * 1.45836_rkind * 0.00001_rkind
 
   tmp1 = t%m_avol / real(t%m_amv, rkind)
   tmp2 = t%m_aent / real(t%m_amv, rkind)

   alpha = (t%m_aalp / real(t%m_amv, rkind)) - (tmp1 * tmp2)
   alpha = alpha / (kB * y%m_temper * y%m_temper * tmp1)

   write(*,'(" Specific heat capacity [Cp]   :         ",e20.10," cal/mol K")') cp
   write(*,'(" Isothermal compressibility[κ] :         ",e20.10," atm⁻¹")') kappa
   write(*,'(" Thermal expansion coeff.[α]   :         ",e20.10," K⁻¹")') alpha
   write(*,*)"--------------------------------------------------------------------------------"

   t%m_aeto  = t%m_aeto + etold
   t%m_aets  = t%m_aets + etsqo

   t%m_aeng  = t%m_aeng + esvold
   t%m_aens  = t%m_aens + esvqol

   aeng  = t%m_aeng / tconf
   aengsq = t%m_aens / tconf

   t%m_arho = t%m_arho + denold
   t%m_arhs = t%m_arhs + dsqold 

   avd    = t%m_arho / tconf
   avdsq  = t%m_arhs / tconf

   t%m_aent = t%m_aent + achold
   t%m_aehs = t%m_aehs + achsqo 
 
   avh   = t%m_aent / tconf
   avhsq = t%m_aehs / tconf 

   c1 = 1.0_rkind / real(t%m_amv)
   c2 = 1.0_rkind / real(t%m_amv * t%m_amv) 
   
   c3 = t%m_avol * c1

   t%m_avol = t%m_avol + volold
   t%m_avos = t%m_avos + vsqold

   avol    = t%m_avol / tconf
   volrad  = avol 
   avolsq  = t%m_avos / tconf
  
   cpt =  (avhsq - avh * avh) / (real(y%m_mxmol, rkind) * kB * y%m_temper * y%m_temper)
   cpt = cpt * 1000.0_rkind

   kappat = (avolsq - avol * avol ) /(kB * y%m_temper * avol)
   kappat = kappat * 1.45836_rkind * 0.00001_rkind
   
   tmp1 = (t%m_aent + alfold) / tconf

   t%m_aalp = t%m_aalp + alfold

   alphat = (tmp1 - avol * avh) / (kB * y%m_temper * y%m_temper * avol)

   cvap = - aeng + 0.0019891_rkind * y%m_temper

   write(*,*)"                               <Averages>                                        "

   avd = y%m_weight * avd * d_con / real(y%m_mxmol)
   avd2 = (y%m_weight/avol) * d_con

   avmol = volrad / (real(y%m_mxmol, rkind) * d_con) 

   open(inaver, file='average', status='unknown', form ='formatted')
   
   rewind(inaver)

   write(inaver,'(a)') y%m_titles
   write(inaver,*) tconf
   write(inaver,*) t%m_aeto, t%m_aets
   write(inaver,*) t%m_aeng, t%m_aens
   write(inaver,*) t%m_arho, t%m_arhs
   write(inaver,*) t%m_aent, t%m_aehs
   write(inaver,*) t%m_avol, t%m_avos
   write(inaver,*) t%m_aalp

   close(inaver)

   write(*,'(" Total number configurations   :   ",i20)') int(tconf)
   write(*,'(" Average energy per molecule   :   ",f20.15,2x," ±",f20.15," kcal/mol")') aeng,sget
   write(*,'(" Enthalpy per molecule         :   ",f20.15,2x," ±",f20.15" kcal/mol")') avh/real(y%m_mxmol, rkind), sgh
   write(*,'(" Density                       :   ",f20.15,2x," ±",f20.15," g/cm³")') avd, sgrho
   write(*,'(" Density [~vol⁻¹]              :   ",f20.15,2x," ±",f20.15," g/cm³")') avd2, sgrho
   write(*,'(" Molar volume (cm ** 3)        :   ",f20.15,2x," ±",f20.15," cm³")') avmol, sgvol
   write(*,*)"Box volume                    :   ", volrad
   write(*,*)"-----------------------------------------------------------------------------------------------"
   write(*,*)"                             <Fluctuations>                                        "
   write(*,'(" Specific heat capacity [Cp]   :   ",e20.10," cal/mol K")') cpt
   write(*,'(" Isothermal compressibility[κ] :   ",e20.10," atm⁻¹")') kappat
   write(*,'(" Thermal expansion coeff.[α]   :   ",e20.10," K⁻¹")') alphat
   write(*,*)"--------------------------------------------------------------------------------"
   write(*,'(" Heat vaporization [ΔHv]       :   ",f20.15," kcal/mol")') cvap
   write(*,*)"-----------------------------------------------------------------------------------------------"

    call r_rnormal( y ) 
   
   end subroutine r_averages

    subroutine r_rnormal( y )
    implicit none
    type(simulation), intent(inout) :: y
    integer :: i, id, nblock, nb20, ii
    real(rkind), dimension(y%m_nsteps) :: acum
    real(rkind) :: average, saverage
    integer, parameter :: nn = 10000
    real(rkind) :: sum, sums
 
 
    average = 0.0_rkind
    saverage = 0.0_rkind
    nblock = 0
    rewind(intemp)
 
    open(inplot,status='scratch')
 
    do i = 1, y%m_nsteps
     read(intemp,*,end=100) id, acum(i)
      average = average + acum(i)
      saverage = saverage + acum(i)*acum(i)
      nblock = nblock + 1
    end do
 
100 id = 0
    nb20 = nblock / nn
 
   do ii = 1, nn
    sum = 0.0_rkind
    sums = 0.0_rkind
 
    do i = 1, nb20
     id = id + 1
       sum = sum + acum(id)
       sums = sums + acum(id) * acum(id)
    end do
 
     write(inplot,*) ii, sum/(nblock/real(nn)),&
       & dsqrt((sum-sums)**2/(nblock/real(nn)))
   end do
 
  end subroutine r_rnormal

  end module m_averages

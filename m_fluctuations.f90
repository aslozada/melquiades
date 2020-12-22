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
module m_fluctuations
!
! !Uses:
!
  use m_kind
  use m_simtype
  use m_boxtype
  use m_constants
  use m_unit
  use m_zeros
  implicit none
  real(rkind), public :: sgh, sgrho, sgvol, sget
!
! !Public members functions:
!
  private
  public :: r_fluctuations
!
! !Revisionn history:
! 16Nov 2015 Asdrubal Lozada
!
!eop
!---------------------------------------------
  contains
!
!bop
!
! !Iroutine: r_flutuations
!
! !Description: This routine calcluates the statistical error
!
!\\
!\\
! !Interface:
!
  subroutine r_fluctuations( t, y, x)
!
  implicit none
!
! !Input parameters:
  type(simulation), intent(inout) :: y
  type(temporary), intent(inout) :: t
  type(box) :: x
!
! !Revision history:
! 16Nov 2015 Asdrubal Lozada
!
!eop
!------------------------------------------------------
! Local variables
  integer :: idone, k
  real(rkind) :: am, cv
  real(rkind) :: cnorm
  
  open(influt,file='error',status='unknown',form='formatted')
  
  rewind(influt)
  
  if( .not.y%m_averag ) then
    idone = 0
  else

   read(influt,*) idone

   do k = 1, idone
     read(influt,*) x%m_sehpy(k), x%m_serho(k), x%m_sevol(k), x%m_sete(k)
   end do
  end if  

   idone = idone + 1

   am = 1.0_rkind / t%m_amv
!----
   x%m_sehpy(idone) = (t%m_aehs*am-(t%m_aent*am)**2)/(y%m_mxmol**2)
   sgh = 0.0_rkind
   
   cv = y%m_weight * d_con / y%m_mxmol
   cv = cv * cv
  
   x%m_serho(idone) = (t%m_arhs*am -(t%m_arho*am)**2) * cv
   sgrho = 0.0_rkind

   cv = y%m_mxmol * d_con
   cv = cv * cv

   x%m_sevol(idone) = (t%m_avos*am -(t%m_avol*am)**2) / cv
   sgvol = 0.0_rkind

   x%m_sete(idone) = (t%m_aens*am -(t%m_aeng*am)**2)
   sget = 0.0_rkind 

   if( y%m_averag ) then

   do k = 1, idone

    sgrho = sgrho + x%m_serho(k)
    sget  = sget  + x%m_sete(k)
    sgvol = sgvol + x%m_sevol(k)
    sgh   = sgh   + x%m_sehpy(k)

  end do

   cnorm = real( (idone * (idone-1)), rkind )
   cnorm = 1.0_rkind / cnorm

   sgrho = sgrho * cnorm
   sgrho = 1.0_rkind * dsqrt(dabs(sgrho))

   sget  = sget * cnorm
   sget  = 1.0_rkind * dsqrt(dabs(sget))

   sgvol = sgvol * cnorm
   sgvol = 1.0_rkind * dsqrt(dabs(sgvol))

   sgh   = sgh * cnorm
   sgh   = 1.0_rkind * dsqrt(dabs(sgh))

   end if


   rewind(influt)

   write(influt,*) idone

   do k = 1, idone
    write(influt,*) x%m_sehpy(k), x%m_serho(k), x%m_sevol(k), x%m_sete(k)
   end do

   close(influt)

  end subroutine r_fluctuations

  end module m_fluctuations

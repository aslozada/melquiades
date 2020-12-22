!----------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                  !
!---------------------------------------------------------------------
!bop
!
! !Module:  m_configuration
!
! !Description: This module contains a routine for
! the configurational energy calculation.
!\\
!\\
! !Interface:
!
  module m_configuration
!
! !Uses:
  use m_kind
  use m_simtype
  use m_boxtype
  use m_precells
  use m_zeros
  use m_interaction
!
  implicit none
!
! !Public member functions:
!
  private
!  
  public  :: r_configuration
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
! !Iroutine: r_configuration
!
! !Description: In this routine the configurational energy values
!are calculated by input the center of mass coordinates sequencially.
!\\
!\\
! !Interface:
!
  subroutine r_configuration( edge, engconf, virconf, t, y, x, istate, iprint )
! return from this routine.
!    
  implicit none
!
! !Input parameters:
  type(simulation), intent(inout) :: y
  type(box), pointer :: x
  type(temporary), intent(inout) :: t
  real(rkind), dimension(:), intent(in) :: edge
  real(rkind), intent(out) :: engconf, virconf
  integer, intent(in) :: istate ! It's a flag for writting
  logical :: iprint
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
! !Remarks:
! The internal function cpu_time is called in order to measure
!the consumed time.
!
!eop
!--------------------------------------------------------
! Local variables
  real(rkind), dimension(3) :: com
  real(rkind) :: energy, virial, pressure
  integer :: ids, jb, jb2
  integer :: i, l, k, nk, m, n, iter, error
  real(rkind) :: apv
!--------------------------------------------------------------
! Check time execution
  real(rkind) :: start_time, stop_time
  real(rkind) :: init_time, end_time

!---------------------------------------------------
   ids = 0
   jb = 0
   jb2 = 0
   iter = 0

   engconf  = 0.0_rkind  ! Cofigurational energy
   virconf  = 0.0_rkind  ! Virial
   pressure = 0.0_rkind  

    if( iprint ) then
     write(*,'("                Configurational Energy                 ")') 
    end if


   if ( .not. y%m_solute ) then

   call cpu_time(start_time)

   do i = 1,  y%m_mxmol-1
   
   com(1) = x%m_cmass(1,i)
   com(2) = x%m_cmass(2,i) 
   com(3) = x%m_cmass(3,i)

   jb = i + 1

   call r_interactions( com, edge, i, jb, ids, energy, virial, t, y, x )

   engconf  = engconf + energy
   virconf  = virconf  + virial
   pressure = y%m_dens * y%m_temper + virial/y%m_vol


  end do ! i

  call cpu_time(stop_time)

   engconf = engconf !+ y%m_vlrc * y%m_voli

   if ( iprint ) then

   if(istate == 0 )  then
     write(*,7000) engconf 
!--------------------------------------------------------------------     
     write(*,*)"-------------------------------------------------------------------------------"
     write(*,7100) virconf
     write(*,7200) pressure
   else
     write(*,7001) engconf
!--------------------------------------------------------------------     
     write(*,*)"-------------------------------------------------------------------------------"
     write(*,7101) virconf
   end if 
    write(*,7002) stop_time - start_time
    write(*,*)'-------------------------------------------------------------------'
   end if ! iprint 
!!*  else

!!*  call cpu_time(init_time)
!!*  call cpu_time(start_time)

!!*  do i = 1, y%m_mxmol-1

!!*   logic = .false.

!!*    do l = 1, y%m_iexcl
!!*     if(x%m_idtype(i) == x%m_extype(l)) then
!!*      logic = .true.
!!*      exit
!!*     end if ! id_type
!!*    end do ! l
   
!!*    if(.not.logic) then
!!*      y%m_solute = .false.

!!*      com(1) = x%m_cmass(1,i)
!!*      com(2) = x%m_cmass(2,i)
!!*      com(3) = x%m_cmass(3,i)

!!*      jb = i + 1

!!*      call r_interactions(com, edge, i, jb, ids, e1, t, y, x)
   
!!*      y%m_enerne = y%m_enerne + e1(1)
!!*      y%m_virial = y%m_virial + e1(2)
     
!!*    end if ! logic
!!*  end do ! i

!!*   call cpu_time(stop_time)

!!*   engconf1 = y%m_enerne

!!*   if (istate == 0) then   
!!*     write(*,8000) engconf1
!!*   else
!!*     write(*,8001) engconf1
!!*   end if
!!*   write(*,8002) stop_time - start_time
    
!!*   y%m_solute = .true.

!!*  call cpu_time(start_time)

!!*  do i = 1, y%m_mxmol-1
  
!!*   logic2 = .false.

!!*   do l = 1, y%m_iexcl
!!*    if(x%m_idtype(i) == x%m_extype(l)) then  
!!*      logic2 = .true.
!!*      exit
!!*    end if ! id_type

!!*  end do ! l

!!*  if( logic2 )  then

!!*     y%m_solute = .true.
!!*     com(1) = x%m_cmass(1,i)
!!*     com(2) = x%m_cmass(2,i)
!!*    com(3) = x%m_cmass(3,i)

!!*    jb2 = i + 1

!!*    call r_interactions(com, edge, i, jb2, ids, e2, t, y, x)

!!*   y%m_eneret = y%m_eneret + e2(1)
!!*    y%m_virial = y%m_virial + e2(2)

!!*   end if ! logic2
!!*  end do !i  

!!*  call cpu_time(stop_time)
!!*  call cpu_time(end_time)

!!*  engconf2 = y%m_eneret 
 
!!*  if (istate == 0) then 
!!*    write(*,9000) engconf2
!!*  else
!!*    write(*,9001) engconf2
!!*  end if
!!*  write(*,9002) stop_time - start_time
!!*  write(*,*)'------------------------------------------------------------------------'

!!*  engconf = y%m_enerto  ! Total configurational eniergy
!!*  engconf = engconf1 + engconf2
 
!!*  if (istate == 0) then 
!!*    write(*,1000) engconf
!!*  else
!!*    write(*,2000) engconf
!!*  end if   
!!*   write(*,3000) end_time - init_time
!!*  write(*,*)'------------------------------------------------------------------------'

 end if ! Solute

7000 format(' Initial configurational energy   : ',' Total',e20.10,' kcal/mol')
7001 format(' Final configurational energy     : ',' Total',e20.10,' kcal/mol')
7002 format(' Time for configuration energy    : ',e20.10,' seconds')
7100 format(' Initial configurational virial   : ',e20.10)
7101 format(' Final configurational virial     : ',' Total',e20.10,' kcal/mol')
7200 format(' Initial pressure                 : ',e20.10)
!!*8000 format(' Initial configuration energy. Non-excluded types : ',e20.10,' kcal/mol')
!!*8001 format(' Final configuration energy. Non-excluded types   : ',e20.10,' kcal/mol')
!!*8002 format(' Time for configuration energy. Non-Excluded types: ',e20.10,' seconds')
!!*9000 format(' Initial configuration energy. Excluded types     : ',e20.10,' kcal/mol')
!!*9001 format(' Final configuration energy. Excluded types       : ',e20.10,' kcal/mol')
!!*9002 format(' Time for configuration energy. Excluded types    : ',e20.10,' seconds')
!!*1000 format(' Total initial configurational energy    : ',e20.10,' kcal/mol')
!!*2000 format(' Total final configurational energy      : ',e20.10,' kcal/mol')
!!*3000 format(' Total time configurational energy       : ',e20.10,' seconds')
   
 end subroutine r_configuration

end module m_configuration

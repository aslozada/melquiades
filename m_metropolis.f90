!----------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program
!----------------------------------------------------------------------
!bop
! !Module:  m_metropolis
!
! !Description: This module is the core of program. Here the Metropolis Monte Carlo algorithm is
!implemented.
!
!A general scheme in the calculation routines is:
!
!\begin{center}
!\begin{verbatim}
!Inital position in 'Markov chain' {Loop in total number of step}
!    call translation components (atoms - molecules - groups)
!     Evaluate energy
!    call rotation components
!     Evaluate energy       
!    call volume modification
!     Evaluate energy       
!    call exchage components
!     Evaluate energy       
!    call internal degree evaluation
!     Evaluate energy       
!      ....
!    call other process 
!     Evaluate energy      
! End position in 'Markov chain' {Loop in number step}  
!\end{verbatim}
!\end{center}
!\\
!\\
! !Interface:
!
  module m_metropolis
!
! !Uses:
!
    use m_kind
    use m_unit
    use m_simtype
    use m_boxtype
    use m_random
    use m_rotation
    use m_pbc
    use m_zeros
    use m_interaction
    use m_uplist
    use m_averages
    use m_sampler
    use m_store
    use m_isobaric
    use m_plots
    use m_configuration
    use m_correlation
    use m_shift
    use m_error
    implicit none
!
! !Public member functions:
!
  private
! 
  public :: r_markov, r_metropolis
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
! !Remarks:
! It's advised to avoid the writting of formatted files during run time.
!
!eop
!----------------------------------------------------------
  contains
!
!bop
!
! !Iroutine: r_metropolis
!
! !Description: In this routine a Markov chain is builded.
!\\
!\\
! !Interface:
  subroutine r_metropolis( edge, engconf, virconf, t, y, x )
!    
  implicit none

! !Input parameters:

  type(simulation), intent(inout) :: y
  type(box), pointer :: x
  type(temporary), intent(inout) :: t
  real(rkind), dimension(:), intent(inout) :: edge
  real(rkind), intent(inout) :: engconf
  real(rkind), intent(inout) :: virconf
!
! !Revision history:
! 10Aug 2015 Asdrubal Lozada
!
!eop
!----------------------------------------------------------------------
  integer :: naccn, attempt, step
  integer :: ncycl, nmoves, nsamp
  integer, parameter :: iseed = 86456

  integer :: acm ! accumulate
  real(rkind) :: enp, press
  real(rkind) :: rho
  real(rkind) :: apv

  real(rkind) :: start_time, stop_time  
  integer :: ifreq, jfreq, iter, icfreq, iter0
!---------
! Verify variables
  integer :: istate
  integer :: var
  real(rkind) :: tol_loc
  
!---------
  istate = 0
  iter = 0
  iter0 = 11

  naccn = 0
  attempt = 0
  nmoves = 1
  nsamp = y%m_mxmol
  acm = 0
  var = 2
  tol_loc = 0.001_rkind
!--------------------

! Frequeny for screening 
   if( y%m_nsteps > 10 ) then
     ifreq = int(y%m_nsteps / 10)
   else
     ifreq = 5
   end if 

! Frequency for writting data
   if( y%m_nsteps > 5 .and. y%m_nsteps < 100 ) then
     jfreq = int(y%m_nsteps / 10)
   elseif(y%m_nsteps > 100) then
     jfreq = int(y%m_nsteps / 1000)
   else
     jfreq = 5
   end if 

! Frequency to finite difference
  
   if( y%m_nsteps > 5 .and. y%m_nsteps < 100 ) then
     icfreq = int(y%m_nsteps / 11)
   elseif(y%m_nsteps > 100) then
     icfreq = int(y%m_nsteps / 11)
   else
     icfreq = 5
   end if 


! Temporary files
  open(intemp,status='scratch')
  open(ixyz,file='file.xyz',status='unknown')

   write(*,*)'--------------------------------------------------------------------------------'
   write(*,*)'                Start of Markov Chain                                           '
   if(y%m_averag) then
    write(*,*)' Sampling                                                                      '
    write(*,*)' Cycle                Energy per molecule[kcal/mol]           Density[g/cm³]   '
   else
    write(*,*)' Thermalization                                                                '
    write(*,*)' Cycle                Energy per molecule[kcal/mol]           Density[g/cm³]   '
    end if
   write(*,*)'--------------------------------------------------------------------------------'
 
   ncycl = y%m_nsteps


 call cpu_time(start_time)


  do step = 1, ncycl ! Start Markov chain

   call r_markov(edge, engconf, virconf, apv, step, t, y, x)

!*! [without save averages] 
  if( .not.y%m_averag ) then   
    call r_samples(intemp, step, iter, engconf, virconf, enp, rho, press, apv, t ,y )
!*    

  if( mod(step,ifreq) == 0 ) then
     write(*,'(1x,i10,a8,i10,a7,2e20.9)')step, ' out of ', ncycl, '   :   ', enp, rho
  end if
 end if

! Sampling for averages      
  if( y%m_averag ) then  
    call r_samples(intemp, step, iter, engconf, virconf, enp, rho, press, apv, t, y )
     if( mod(step,ifreq) == 0 ) then
     write(*,'(1x,i10,a8,i10,a7,2e20.9)')step, ' out of ', ncycl, '   :   ', enp, rho
     end if

     if( mod(step,jfreq) == 0 ) then
      call r_saves(engconf, virconf,y,x) 
     end if
   end if

!!!      call r_saves(engconf,y,x) 

!-----------------------------------------
! Shift volume

  if( y%m_ensbls=='nvt' .and. y%m_corr) then
    if( y%m_varold < tol_loc) then
      if( mod(step,icfreq) == 0) then
        call r_shift( edge, iter0, y, x)
      end if  
    end if   
  end if  
!------------------------------------------

  end do ! End Markov chain
  
  if( mod(step,jfreq) == 0 ) then
      call r_saves(engconf, virconf,y,x) 
  end if

  call cpu_time(stop_time)

   write(*,*)'--------------------------------------------------------------------------'
   write(*,*)'                End of Markov Chain                                       '
   write(*,*)'--------------------------------------------------------------------------'
   write(*,5000) (stop_time - start_time)
   write(*,*)'--------------------------------------------------------------------------'
   write(*,*)'--------------------------------------------------------------------------'

 !  if( y%m_corr .and. .not.y%m_averag) then
    call r_corr( y, var, iter,  1 ) 
!   end if

 !! if( y%m_averag ) then 
    call r_averages(t, y, x)
 !! end if

  call r_saves( engconf, virconf, y, x )   
  call r_newbox( y, x )
 
  if(y%m_plot) then
   call r_plot( y )
  end if

5000 format(' Total time in Markov chain : ',e20.10,' seconds')

end subroutine r_metropolis
!
!bop
!
! !Iroutine: r_markov
!
! !Description: This routine.
!\\
!\\
! !Interface:
  subroutine r_markov( edge, engconf, virconf, apv, steps, t, y, x )
!    
  implicit none
!
! !Input parameters:
  type(simulation), intent(inout) :: y
  type(box), pointer :: x
  type(temporary), intent(inout) :: t
  real(rkind), dimension(:), intent(in) :: edge
  real(rkind), intent(inout) :: engconf
  real(rkind), intent(inout) :: virconf
  real(rkind), intent(out) :: apv
  integer, intent(in) :: steps
!
! !Revision history:
! 10Aug 2015 Asdrubal Lozada
! 
!eop
!---------------------------------------------------------------------
! Local variables
  integer        :: i_d, i
  real(rkind), dimension(3) :: com_old
  real(rkind), dimension(3) :: com_new

  real(rkind) :: e_old, e_new
  real(rkind) :: v_old, v_new

  integer        :: ido, idn
  real(rkind)    :: deltv, deltvb
  real(rkind)    :: deltw
  integer        :: accept
  integer        :: jb
  logical        :: logic, logic2
  integer        :: l

  real(rkind)    :: eng_temp
  real(rkind)    :: xi
!--------------------------------------------------
    deltv    = 0.0_rkind
    deltvb   = 0.0_rkind
    deltw    = 0.0_rkind
    apv      = 0.0_rkind
    ido = 0
    idn = 1
    jb = 1

 if(.not.y%m_solute) then

     xi = m_rand()

!100 i_d =  int(y%m_mxmol * m_rand()) + 1
100 i_d =  int(y%m_mxmol * xi) + 1

    if(i_d > y%m_mxmol) goto 100

    accept = 0
!--------------------------------------------------
    com_old(1) = x%m_cmass(1,i_d)
    com_old(2) = x%m_cmass(2,i_d)
    com_old(3) = x%m_cmass(3,i_d)
!---------------------------------------------------
! Backup com position
    x%m_rcom(1,i_d) = x%m_cmass(1,i_d)  
    x%m_rcom(2,i_d) = x%m_cmass(2,i_d) 
    x%m_rcom(3,i_d) = x%m_cmass(3,i_d)
!----------------------------------------------------
! Start target site
    do i = 1, x%m_ns(i_d)
     x%m_rot(1,i) = x%m_site(1,i,i_d)
     x%m_rot(2,i) = x%m_site(2,i,i_d)
     x%m_rot(3,i) = x%m_site(3,i,i_d)
   end do
!----------------------
! Backup site position
   do i = 1, x%m_ns(i_d)
    x%m_rsite(1,i) = x%m_site(1,i,i_d)
    x%m_rsite(2,i) = x%m_site(2,i,i_d)
    x%m_rsite(3,i) = x%m_site(3,i,i_d)
   end do
!-------------------------
  call r_interactions( com_old, edge, i_d, jb, ido, e_old, v_old, t, y, x )
!---------------------------------------------------
! Transformation matrix in R3 apply
!-------------------------------------------------------------------------------
    com_new(1) = (com_old(1) + (2.0d0 * m_rand() - 1.0d0) * y%m_transl)
    com_new(2) = (com_old(2) + (2.0d0 * m_rand() - 1.0d0) * y%m_transl)
    com_new(3) = (com_old(3) + (2.0d0 * m_rand() - 1.0d0) * y%m_transl)
!------------------------------------------------------------------------------
!   Pick up central image

   if( y%m_pbcs ) then  

  call r_pbc(com_new, edge)

   end if
!---------------------------------------------------------------------
  call r_euler(i_d,y,x)
!--------------------------------------------------------------------
  call r_interactions( com_new, edge, i_d, jb, idn, e_new, v_new, t, y, x )
!-------------------------------------------------------------
!boc
    deltv = e_new - e_old ! Energy variation
    deltw = v_new - v_old ! Virial variation
    deltvb = deltv * y%m_beta
!-----------------------------

     if(deltv < 10.0d0) then
       if(deltv <= 0.0d0) then    
        accept = 1
        t%m_move = t%m_move + 1
       else  if(exp(-deltvb) > m_rand()) then
        accept = 1
        t%m_move = t%m_move + 1
       end if
      end if
!eoc
      if(accept == 1) then
       engconf = engconf + deltv
       virconf = virconf + deltw
!-----------------
      end if

      t%m_imove = t%m_imove + accept
      t%m_amv = t%m_amv + 1
!-------------------------------------------------------------
   call r_update(com_new, i_d, accept, y, x)
!-------------------------------------------------------------
  if(y%m_ensbls == 'npt') then
   if( mod (steps, y%m_ifvol) == 0 ) then
    call r_npt(engconf, virconf, apv, t, y, x)
   end if
  end if

! End conventional calculation
!--------------------------------------------------------------
!*  else 
!---
!*200 i_d =  int(y%m_mxmol * m_rand()) + 1
!*    if(i_d > y%m_mxmol) goto 200
!*
!*    logic2 = .false.
!--------------------------------------------------
!*  do l = 1, y%m_iexcl
!*   if(x%m_idtype(i_d) == x%m_extype(l)) then
!*    logic2 = .true.
!*    exit
!*   end if
!*  end do ! l 
!---------------------------------------------------
!*   if(logic2) then

!*   com_old(1) = x%m_cmass(1,i_d)
!*   com_old(2) = x%m_cmass(2,i_d)
!*   com_old(3) = x%m_cmass(3,i_d)
!---------------------------------------------------
! Backup com position
!*    x%m_rcom(1,i_d) = x%m_cmass(1,i_d)  
!*    x%m_rcom(2,i_d) = x%m_cmass(2,i_d) 
!*    x%m_rcom(3,i_d) = x%m_cmass(3,i_d)
!----------------------------------------------------
!*   do i = 1, x%m_ns(i_d)
!*    x%m_rot(1,i) = x%m_site(1,i,i_d)
!*    x%m_rot(2,i) = x%m_site(2,i,i_d)
!*    x%m_rot(3,i) = x%m_site(3,i,i_d)  
!*   end do ! i
!----------------------
! Backup site position
!*   do i = 1, x%m_ns(i_d)
!*    x%m_rsite(1,i) = x%m_site(1,i,i_d)
!*    x%m_rsite(2,i) = x%m_site(2,i,i_d)
!*    x%m_rsite(3,i) = x%m_site(3,i,i_d)
!*   end do
!-------------------------
!*   call interactions(com_old, edge, i_d, jb, ido, e_old, t, y, x)
!----
! Transformation matrix in R3 apply
!*   com_new(1) = (com_old(1) + (2.0d0 * m_rand() - 1.0d0) * y%m_transl)
!*   com_new(2) = (com_old(2) + (2.0d0 * m_rand() - 1.0d0) * y%m_transl)
!*   com_new(3) = (com_old(3) + (2.0d0 * m_rand() - 1.0d0) * y%m_transl)
!-----------
! Pick up central image
!*   call pbc(com_new, edge)
!------
!*   call rotation(i_d, y, x)
!------
!*   call interactions(com_new, edge, i_d, jb, idn,e_new, t, y, x)
!---------
!*   deltv = e_new(1) - e_old(1) ! Energy variation
!*   deltvb = deltv * y%m_beta
!-------  
!*   if(deltv < 75.0d0) then
!*    if(deltv <= 0.0d0) then
!*     accept = 1
!*   else if(exp(-deltvb) > m_rand()) then
!*     accept = 1
!*    end if
!*   end if
   
!*   if(accept == 1) then
!*    engconf  = engconf + deltv
!*    enthalpy = enthalpy + (e_new(1) + y%m_ppcon * y%m_vol)-(e_old(1) + y%m_ppcon * y%m_vol) 
!     y%m_virial = y%m_virial + (vir_new-vir_old)
!*   end if
     
!*    call update3(com_new, i_d, accept, y, x)

!*   end if ! logic
!----------------------
! -----------------------
! Non-exclude types
!--------------------------------------
!*   logic = .false.
!--------------------------------------------------
!*  do l = 1, y%m_iexcl
!*   if(x%m_idtype(i_d) == x%m_extype(l)) then
!*    logic = .true.
!*    exit
!*   end if
!*  end do ! l 
!---------------------------------------------------
!*   if(.not.logic) then
!*   y%m_solute = .false.

!*   com_old(1) = x%m_cmass(1,i_d)
!*   com_old(2) = x%m_cmass(2,i_d)
!*   com_old(3) = x%m_cmass(3,i_d)
!---------------------------------------------------
! Backup com position
!*    x%m_rcom(1,i_d) = x%m_cmass(1,i_d)  
!*    x%m_rcom(2,i_d) = x%m_cmass(2,i_d) 
!*    x%m_rcom(3,i_d) = x%m_cmass(3,i_d)
!----------------------------------------------------
!*   do i = 1, x%m_ns(i_d)
!*    x%m_rot(1,i) = x%m_site(1,i,i_d)
!*    x%m_rot(2,i) = x%m_site(2,i,i_d)
!*    x%m_rot(3,i) = x%m_site(3,i,i_d)  
!*   end do ! i
!----------------------
! Backup site position
!*   do i = 1, x%m_ns(i_d)
!*    x%m_rsite(1,i) = x%m_site(1,i,i_d)
!*    x%m_rsite(2,i) = x%m_site(2,i,i_d)
!*    x%m_rsite(3,i) = x%m_site(3,i,i_d)
!*   end do
!-------------------------
!*   call interactions(com_old, edge, i_d, jb, ido, e_old, t, y, x)
!----
! Transformation matrix in R3 apply
!*   com_new(1) = (com_old(1) + (2.0d0 * m_rand() - 1.0d0) * y%m_transl)
!*   com_new(2) = (com_old(2) + (2.0d0 * m_rand() - 1.0d0) * y%m_transl)
!*   com_new(3) = (com_old(3) + (2.0d0 * m_rand() - 1.0d0) * y%m_transl)
!-----------
! Pick up central image
!*   call pbc(com_new, edge)
!------
!*   call rotation(i_d, y, x)
!------
!*   call interactions(com_new, edge, i_d, jb, idn, e_new, t, y, x)

!---------
!*   deltv = e_new(1) - e_old(1) ! Energy variation
!*   deltvb = deltv * y%m_beta
!-------  
!*   if(deltv < 75.0d0) then
!*    if(deltv <= 0.0d0) then
!*     accept = 1
!*    else if(exp(-deltvb) > m_rand()) then
!*     accept = 1
!*    end if
!*   end if
   
!*   if(accept == 1) then
!*    engconf = engconf + deltv
!*    enthalpy = enthalpy + (e_new(1) + y%m_ppcon * y%m_vol)-(e_old(1) + y%m_ppcon * y%m_vol) 
!     y%m_virial = y%m_virial + (vir_new-vir_old)
!*   end if
     
!*    call update2(com_new, i_d, accept, y, x)

!*   end if ! logic
!-------------------------------------------------------------
!*  if(y%m_ensbls == 'npt') then
!*   if( mod (steps, y%m_ifvol) == 0 ) then
!*    call  npt(engconf, enthalpy, t, y, x)
!*   end if
!*  end if
!----------------------------------------------------------
!*   y%m_solute = .true.

 end if  ! Active

!------------------------------------------------------------
    end  subroutine r_markov

    end module m_metropolis

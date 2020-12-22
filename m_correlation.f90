!----------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program
!----------------------------------------------------------------------
!bop
!
! !Module: m_correlation
!
! !Description: This module contains a routine for correlation calculation!\\
!\\
! !Interface:
!
  module m_correlation
!
! !Uses
  use m_kind
  use m_simtype
  use m_unit
!
! !Public member functions:
!
  private
  public :: r_corr
!
! ! Revision history:
! 11Nov 2015 Asdrubal Lozada
!
!eop
!------------------------------------------------------------------------
  contains
!
!bop
!
! !Iroutine: r_corr
!
! !Description: Correlation eneregy calculation.
!\\
!\\
! !Interface:
   subroutine r_corr( y, var, iter, stats )
   implicit none
!  Dummy arguments
   type(simulation), intent(inout) :: y
   integer, intent(in) :: var ! Total number of variables to be used
   integer, intent(in) :: stats ! Define variable to correlation 
   real(rkind), allocatable, dimension(:) :: acc
   real(rkind), allocatable, dimension(:) :: corr
   integer, allocatable, dimension(:) :: istore
   integer :: id, i, j, ntotal, error, tau
   real(rkind) aver, dev
   real(rkind) :: threshold, corr0, tole
   integer, intent(in) :: iter

   aver = 0.0_rkind
   dev  = 0.0_rkind

   ntotal = 0
   threshold = 0.001_rkind
   tole = 0.005_rkind
   tau = int(iter/5.0_rkind) ! Value sugest by V.S. Pugachev "I.Teória Probabilidades" 

   allocate(acc(0:iter),istore(0:iter),stat=error)
    if( error /= 0) then
     write(*,*) "Allocation Error in m_correlation"
     stop
    end if    
   allocate(corr(0:iter),stat=error)
    if( error /= 0) then
     write(*,*) "Allocation Error in m_correlation"
     stop
    end if    

   rewind(intemp)


   do i = 0, iter-1  
     read(intemp,*,end=100) acc(i), istore(i)
     corr(i) = 0.0_rkind
     ntotal = ntotal + 1
   end do
  
100 continue

     do i = 0, iter-1
       aver = aver + acc(i)
       dev  = dev + acc(i) * acc(i)
    end do     
 
     aver = aver/real(ntotal,rkind)
     dev  = (dev/real(ntotal,rkind)) - aver * aver
!----------------------------------------------------
     open(incorr, file='correlation',status='unknown', form='formatted')

     rewind(incorr)

     write(incorr,'(a)')
     write(incorr,*) dev
     close(incorr) 
!----------------------------------------------------

!!*   if( dev > 0.0_rkind ) dev = dsqrt(dev)

    if(dev <= threshold) then 

    do i = 0, tau
      do j = 1, iter-tau
       corr(i) = corr(i) + (acc(j) - aver) *(acc(j+i)-aver)

       if( abs(corr(i)) <= tole ) then         

!        write(*,*) istore(j+i)
       end if
      end do  
    end do 

   corr0 = corr(0)

   do i =0, tau
    corr(i) = corr(i)/corr0
    write(39,*) i, corr(i)
   end do

   end if

  deallocate(acc,corr,istore)

  close(intemp)


   end subroutine r_corr

 end module m_correlation

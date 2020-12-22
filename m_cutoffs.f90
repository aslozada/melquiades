!---------------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                       !
!---------------------------------------------------------------------------
!bop
!
! !Module:  m_cutoffs
!
! !Description: This module contains routines for calculating cut-off values.
!\\
!\\
! !Interface:
!
module m_cutoffs
!
! !Uses:
!
  use m_kind
  use m_simtype
  use m_boxtype
  implicit none
!
! !Public member functions:
!
  public :: r_cuts
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
!eop
!------------------------------------------------------------------------
  contains
!
!bop
!
! !Iroutine: r_cuts

! !Description: This routine calculates the total cut-off value as
!function of the molecular topology.
!\begin{equation*}
!r_{c} = cut-off + max(d_{1},d_{2},d_{3},\cdots)
!\end{equation*}
!\vspace{-\topsep}
!\\
!\\
! !Interface:
!
  subroutine r_cuts( y, x )
!
  implicit none
!  
! !Input parameters:  
  type(simulation), intent(inout) :: y
  type(box), pointer :: x
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
!eop
!---------------------------------------------------------------------
!   Local variables
  real(rkind) :: dist, mx_dma
  real(rkind) :: dist2, mx_dma2
  integer :: i, k, l
  integer :: num ! temporary iterator
  logical :: logic ! temporary logic variable

  y%m_rcutsq = y%m_cutoff * y%m_cutoff
  x%m_fn  =  0.0_rkind
  x%m_fns =  0.0_rkind
  num = 0

  if(.not.y%m_solute) then
    do k = 1, y%m_ntf
      num = num + x%m_nmol(k)
      do i = 1, x%m_nsite(k)
        dist = x%m_site(1,i,num)*x%m_site(1,i,num)+&
          & x%m_site(2,i,num)*x%m_site(2,i,num)+&
          & x%m_site(3,i,num)*x%m_site(3,i,num)

        if(dist >= x%m_fn(k)) x%m_fn(k) = dist
      end do ! i
      x%m_fn(k) = dsqrt(x%m_fn(k))
    end do !  k
!boc
  mx_dma = x%m_fn(1)
    if( y%m_ntf > 1) then
      do k = 2, y%m_ntf
       if( x%m_fn(k) > mx_dma) mx_dma = x%m_fn(k)
      end do
    end if

 y%m_rpair = (y%m_cutoff + 2 * mx_dma) * (y%m_cutoff + 2 * mx_dma)
 y%m_cutoff = y%m_cutoff + 2 * mx_dma
!eoc 
 else
   do k = 1, y%m_ntf 
    logic = .false.
    num = num + x%m_nmol(k)
     do l = 1, y%m_iexcl
      if( x%m_idtype(num) == x%m_extype(l) ) then
        logic = .true.
        exit
      end if ! id_type
     end do ! l

  if( .not.logic ) then
    do i = 1, x%m_nsite(k)
    dist = x%m_site(1,i,num)*x%m_site(1,i,num)+&
      & x%m_site(2,i,num)*x%m_site(2,i,num)+&
      & x%m_site(3,i,num)*x%m_site(3,i,num)

    if( dist >= x%m_fn(k) ) x%m_fn(k) = dist
    end do ! i

   x%m_fn(k) = sqrt(x%m_fn(k))
  end if 
  
  do i = 1, x%m_nsite(k)
  dist2 = x%m_site(1,i,num)*x%m_site(1,i,num)+&
    & x%m_site(2,i,num)*x%m_site(2,i,num)+&
    & x%m_site(3,i,num)*x%m_site(3,i,num)

    if(dist2 >= x%m_fns(k)) x%m_fns(k) = dist2
  end do ! i   

   x%m_fns(k) = sqrt(x%m_fns(k))

  end do ! k

  mx_dma = x%m_fn(1)
  mx_dma2 = x%m_fns(1)

 if(  y%m_ntf > 1 ) then
  do k = 2, y%m_ntf
   if( x%m_fn(k) > mx_dma ) mx_dma = x%m_fn(k)
   if( x%m_fns(k) > mx_dma2 ) mx_dma2 = x%m_fns(k)
  end do
 end if

 y%m_rpair = (y%m_cutoff + 2 * mx_dma) * (y%m_cutoff + 2 * mx_dma)
 y%m_rsol  = (y%m_cutsol + 2 * mx_dma2) * &
               & (y%m_cutsol + 2 * mx_dma2)

 end if ! main condition

 end subroutine r_cuts

 end module m_cutoffs

!----------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program
!----------------------------------------------------------------------
!bop
!
! !Module: m_ljones
!
! !Description: This module contains a routine for intermolecular
!energy calculation using the Lennard-Jones + Coulombic potential model:
!  \begin{equation*}
!   V(r_{ij}) = \displaystyle
!   \sum^{i<j}4\epsilon_{ij}\left[\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12}-\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{6}\right]+\
!    \displaystyle \sum^{i<j}\frac{q_{i}q_{j}}{r_{ij}}
!      \end{equation*}
!        \vspace{-\topsep} 
!\\
!\\
! !Interface:
!
module m_longs
!
! !Uses:
  use m_kind
  use m_simtype
  use m_boxtype   
  use m_constants, only : c_kcalmol, pi
  use m_zeros
  implicit none
!
! !Public member functions:
!
  public ::r_long
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada 
!
!eop
!-------------------------------------------------------------------
  contains 
!
!bop
!
! !Iroutine:  p_ljones
!
! !Description: This routine defines the Lennard-Jones and Coulombic potentials.
!The Berthelot-Lorentz's rules are used in order to obtain the combined parameters, e.g.,
!\begin{equation*}
!\epsilon_{ij} = \sqrt{\epsilon_{ii}\epsilon_{jj}}
!\end{equation*}
! \vspace{-\topsep} 
!\\
!\\
! !Interface:
  subroutine r_long(y, x)
   implicit none
   type(simulation), intent(inout) :: y
   type(box), pointer :: x
   integer :: i, j, m, n
   real(rkind) :: ai_lj, bi_lj
   real(rkind) :: aj_lj, bj_lj
   real(rkind) :: aij_lj, bij_lj
   integer :: itera, iterb
   integer :: iterc, iterd
   real(rkind) :: rc3, rc9
   real(rkind) :: s3, s9
   real(rkind) :: sum_ii, sum_ij

   ai_lj = 0.0_rkind; bi_lj = 0.0_rkind
   aj_lj = 0.0_rkind; bj_lj = 0.0_rkind
   aij_lj = 0.0_rkind; bij_lj = 0.0_rkind
   itera = 0; iterb = 0
   iterc = 0; iterd = 0
   sum_ii = 0.0_rkind; sum_ij = 0.0_rkind

   x%s3_ii = 0.0_rkind; x%s9_ii = 0.0_rkind
   x%s3_ij = 0.0_rkind; x%s9_ij = 0.0_rkind
   x%lrc_ii = 0.0_rkind; x%lrc_ij = 0.0_rkind

   do i = 1, y%mi_ntf
     do j = i, y%mi_ntf

      if( j == i) then

          itera = itera + 1
        
          do m = 1, x%mi_nsite(i)
             ai_lj = 4.0_rkind * x%mi_epsilon(m,i) * (x%mi_sigma(m,i) ** 12)
             bi_lj = 4.0_rkind * x%mi_epsilon(m,i) * (x%mi_sigma(m,i) ** 6)
            do n = 1, x%mi_nsite(j)
             aj_lj = 4.0_rkind * x%mi_epsilon(n,j) * (x%mi_sigma(n,j) ** 12)
             bj_lj = 4.0_rkind * x%mi_epsilon(n,j) * (x%mi_sigma(n,j) ** 6)

             if(ai_lj * aj_lj /= 0.0_rkind) then
                 aij_lj = dsqrt(ai_lj * aj_lj)
                 bij_lj = -dsqrt(bi_lj * bj_lj)

                 x%s9_ii(itera) = x%s9_ii(itera) + aij_lj
                 x%s3_ii(itera) = x%s3_ii(itera) - bij_lj

             else
                 aij_lj = 0.0_rkind
                 bij_lj = 0.0_rkind
             end if    
             
            end do
          end do  


      elseif( j /= i) then

          iterb = iterb + 1

          do m = 1, x%mi_nsite(i)
             ai_lj = 4.0_rkind * x%mi_epsilon(m,i) * (x%mi_sigma(m,i) ** 12)
             bi_lj = 4.0_rkind * x%mi_epsilon(m,i) * (x%mi_sigma(m,i) ** 6)
            do n = 1, x%mi_nsite(j)
             aj_lj = 4.0_rkind * x%mi_epsilon(n,j) * (x%mi_sigma(n,j) ** 12)
             bj_lj = 4.0_rkind * x%mi_epsilon(n,j) * (x%mi_sigma(n,j) ** 6)

             if(ai_lj * aj_lj /= 0.0_rkind) then
                 aij_lj = dsqrt(ai_lj * aj_lj)
                 bij_lj = -dsqrt(bi_lj * bj_lj)

                 x%s9_ij(iterb) = x%s9_ij(iterb) + aij_lj
                 x%s3_ij(iterb) = x%s3_ij(iterb) - bij_lj

             else
                 aij_lj = 0.0_rkind
                 bij_lj = 0.0_rkind
             end if    
            end do
          end do  

      end if    

     end do
  end do  

  rc3 = y%m_cutoff * y%m_cutoff * y%m_cutoff
  rc9 = rc3 * rc3 * rc3

  s9 = 2.0_rkind * pi / ( 9.0_rkind * rc9 )
  s3 = 2.0_rkind * pi / ( 3.0_rkind * rc3 )


  do i = 1, y%mi_ntf

   iterc = iterc + 1 

   x%lrc_ii(iterc) = (s9 * x%s9_ii(iterc) - s3 * x%s3_ii(iterc)) * x%mi_nmol(i) *x%mi_nmol(i)
   sum_ii = sum_ii + x%lrc_ii(iterc)

  end do
  

  if( y%mi_ntf > 1 ) then

  do i = 1, y%mi_ntf-1
    do j = i+1, y%mi_ntf

     iterd = iterd + 1

   x%lrc_ij(iterd) = 2.0_rkind * (s9 * x%s9_ij(iterd) - s3 * x%s3_ij(iterd)) * x%mi_nmol(i) *x%mi_nmol(j)
     
   sum_ij = sum_ij + x%lrc_ij(iterd)
    end do
  end do  

  end if 

! Long-range corrections value

  y%m_vlrc = sum_ii + sum_ij
  


 end subroutine r_long

 end module m_longs

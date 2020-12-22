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
module m_ljones
!
! !Uses:
  use m_kind
  use m_simtype
  use m_boxtype   
  use m_constants, only : c_kcalmol, pi
  use m_zeros
  use m_longs
  implicit none
!
! !Public member functions:
!
  public :: p_ljones
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
  subroutine p_ljones( rij, v, i, j, ids, t, x )
! 
  implicit none
  !
! !Input parameters:
  type(box), pointer :: x
  type(temporary), intent(inout) :: t
  real(rkind), dimension(:), intent(in) :: rij
  real(rkind), dimension(:), intent(inout) :: v
  integer, intent(in) :: i, j, ids

!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
!eop
!------------------------------------------------------------
! Local variables
  real(rkind) :: alj_i, blj_i, qc_i, alj_j, blj_j, qc_j
  real(rkind) :: av_i, bv_i, qv_i, av_j, bv_j, qv_j
  real(rkind) :: alj_ij, blj_ij, qc_ij, av_ij, bv_ij, qv_ij
  real(rkind) :: sig6, sig12, vir6, vir12, s6, s12
  real(rkind) :: r6, r_i, t6, t12, rmx, rmy, rmz
  real(rkind) :: rnx, rny, rnz, rmnx, rmny, rmnz, rmnsq, rs
  integer :: m, n
!--------------------------------------------------------------

  alj_i  = 0.0_rkind; alj_j = 0.0_rkind; sig6  = 0.0_rkind; t12 = 0.0_rkind
  blj_i  = 0.0_rkind; blj_j = 0.0_rkind; sig12 = 0.0_rkind; rmx = 0.0_rkind
  qc_i  = 0.0_rkind;  qc_j = 0.0_rkind; vir6  = 0.0_rkind; rmy = 0.0_rkind
  av_i   = 0.0_rkind; av_j  = 0.0_rkind; vir12 = 0.0_rkind; rmz = 0.0_rkind
  bv_i   = 0.0_rkind; bv_j  = 0.0_rkind; s6    = 0.0_rkind; rnx = 0.0_rkind
  qv_i  = 0.0_rkind;  qv_j = 0.0_rkind; s12   = 0.0_rkind; rny = 0.0_rkind
  alj_ij = 0.0_rkind; av_ij = 0.0_rkind; r6    = 0.0_rkind; rnz = 0.0_rkind 
  blj_ij = 0.0_rkind; bv_ij = 0.0_rkind; r_i   = 0.0_rkind; rs  = 0.0_rkind
  qc_ij  = 0.0_rkind; qv_ij = 0.0_rkind; t6    = 0.0_rkind; rmnsq  = 0.0_rkind 


    do m = 1, x%m_ns(i)
        t%m_ski = x%m_sigma(m,i)
        t%m_eki = x%m_epsilon(m,i)
        t%m_qki = x%m_charge(m,i)

        alj_i = 4.0_rkind * t%m_eki * (t%m_ski ** 12)
        blj_i = 4.0_rkind * t%m_eki * (t%m_ski ** 6)
        qc_i   = t%m_qki

        av_i = 48.0_rkind * t%m_eki * (t%m_ski ** 12) / 3.0_rkind
        bv_i = 48.0_rkind * t%m_eki * 0.5_rkind * (t%m_ski ** 6) / 3.0_rkind
        qv_i = t%m_qki

          if( ids == 1 ) then
             rmx = x%m_rot(1,m)
             rmy = x%m_rot(2,m)
             rmz = x%m_rot(3,m)   
          else   
             rmx = x%m_site(1,m,i)
             rmy = x%m_site(2,m,i)
             rmz = x%m_site(3,m,i)     
          end if 

       do n = 1, x%m_ns(j)
            t%m_skj = x%m_sigma(n,j)
            t%m_ekj = x%m_epsilon(n,j)
            t%m_qkj = x%m_charge(n,j)

!boc
            alj_j = 4.0_rkind * t%m_ekj * (t%m_skj ** 12)
            blj_j = 4.0_rkind * t%m_ekj * (t%m_skj ** 6)
            qc_j  = t%m_qkj
!eoc   

            av_j = 48.0_rkind * t%m_ekj * (t%m_skj ** 12) / 3.0_rkind
            bv_j = 48.0_rkind * t%m_ekj * 0.5_rkind * (t%m_skj ** 6) / 3.0_rkind
            qv_j = t%m_qkj

            rnx = x%m_site(1,n,j)
            rny = x%m_site(2,n,j)
            rnz = x%m_site(3,n,j)

            rmnx = rij(1) + (rmx - rnx)
            rmny = rij(2) + (rmy - rny)
            rmnz = rij(3) + (rmz - rnz) 

            rmnsq = rmnx * rmnx + rmny * rmny + rmnz * rmnz
            rs = rmnsq

!boc

! Berthelot-Lorentz rule

            if(alj_i * alj_j /= 0.0_rkind)  then
              alj_ij = dsqrt(alj_i * alj_j)
              blj_ij = -dsqrt(blj_i * blj_j)
              qc_ij = (qc_i * qc_j) * c_kcalmol
            else
              alj_ij = 0.0_rkind
              blj_ij = 0.0_rkind
              qc_ij = (qc_i * qc_j) * c_kcalmol
            end if

            if(av_i * av_j /= 0.0_rkind)  then
              av_ij = dsqrt(av_i * av_j)
              bv_ij = -dsqrt(bv_i *bv_j)
              qv_ij = (qv_i * qv_j) * c_kcalmol
            else
              av_ij = 0.0_rkind
              bv_ij = 0.0_rkind
              qv_ij = (qv_i * qv_j) * c_kcalmol
            end if    
!eoc
              sig6  = blj_ij
              sig12 = alj_ij 

              vir6 = av_ij
              vir12 = bv_ij 

            if(sig6 /= 0.0_rkind .or. sig12 /= 0.0_rkind) then
              r6 = rs * rs * rs
              r_i = 1.0_rkind / r6

              t6 = sig6 * r_i
              v(1) = v(1) + t6
              t12 = sig12 * r_i * r_i
              v(2) = v(2) + t12
              s6 = vir6 * r_i
              v(4) = v(4) + s6
              s12 = vir12 * r_i * r_i
              v(5) = v(5) + s12
            end if 
   
            if(qc_ij /= 0.0_rkind) then
             qc_ij = qc_ij / dsqrt(rs)
             v(3) = v(3) + qc_ij
             qv_ij = qv_ij / dsqrt(rs)
             v(6) = v(6) + qv_ij
            end if 

       end do  ! n
    end do ! m

!---------------------------------------------------

 end subroutine p_ljones

 end module m_ljones

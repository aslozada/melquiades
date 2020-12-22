!--------------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                      ! 
!--------------------------------------------------------------------------
!bop
!
! !Module: m_zeros
!
! !Description: In this module the global variables are 
!initialized.
!\\
!\\
! !Interface:
!
module m_zeros
!
! !Uses:
!
  use m_kind
  use m_simtype
  use m_boxtype
  implicit none
!
! !Public data members:
!
  public :: temporary
!  
! !Public member functions:
!
  public :: r_cero
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
!eop
!------------------------------------------------------------------------

 type temporary
  real(rkind) :: m_ski    ! ~ van der Waals sigma i
  real(rkind) :: m_skj    ! ~ van der Waals sigma j
  real(rkind) :: m_eki    ! ~ van der Waals epsilon i
  real(rkind) :: m_ekj    ! ~ van der Waals epsilon j
  real(rkind) :: m_qki    ! ~ coulombic charge i
  real(rkind) :: m_qkj    ! ~ coulombic charge j
  real(rkind) :: m_aki    ! ~ buckingham a i
  real(rkind) :: m_akj    ! ~ buckingham a j
  real(rkind) :: m_bki    ! ~ buckingham b i
  real(rkind) :: m_bkj    ! ~ buckingham b j
  real(rkind) :: m_cki    ! ~ buckingham c i
  real(rkind) :: m_ckj    ! ~ buckingham c j 
  real(rkind) :: m_gki    ! ~ coupling yukawa g i
  real(rkind) :: m_gkj    ! ~ coupling yukawa g j
  real(rkind) :: m_mki    ! ~ inverse scope yukawa km i
  real(rkind) :: m_mkj    ! ~ inverse scope yukawa km j
  real(rkind) :: m_bo     ! ~ Harmonic equilibrium position
  real(rkind) :: m_kb     ! ~ Harmonic restoring constant
  real(rkind) :: m_ao     ! ~ Harmonic equilibrium angle
  real(rkind) :: m_ka     ! ~ Harmonic restoring angle
  real(rkind) :: m_c0     ! ~ expansion fourier serie constant
  real(rkind) :: m_c1     ! ~ expansion fourier serie constant
  real(rkind) :: m_c2     ! ~ expansion fourier serie constant
  real(rkind) :: m_c3     ! ~ expansion fourier serie constant
  real(rkind) :: m_c4     ! ~ expansion fourier serie constant
  real(rkind) :: m_c5     ! ~ expansion fourier serie constant
  integer      :: m_move   ! ~ moves counts in Markov chain
  integer      :: m_imove  ! ~ accepted moves 
  real(rkind) :: m_avol   ! ~ accumulative volume
  real(rkind) :: m_avos   ! ~ accumulative square volume
  real(rkind) :: m_aent   ! ~ accumulative enthalpy value
  real(rkind) :: m_aehs   ! ~ accumulative square enthalpy value
  real(rkind) :: m_aalp   ! ~ accumulative alpha value
  real(rkind) :: m_arho   ! ~ accumulative density value
  real(rkind) :: m_arhs   ! ~ accumulative square density value
  real(rkind) :: m_aeng   ! ~ accumulative energy per molecule
  real(rkind) :: m_aens   ! ~ accumulative square energy per molecule
  real(rkind) :: m_aeto   ! ~ accumulative total energy
  real(rkind) :: m_aets   ! ~ accumulative square total energy
  integer     :: m_amv    ! ~ accumulative moves  
  
 end type temporary

  contains
!bop
!
! !Iroutine:  r_cero
!
! !Description: Add zero value to variables.
!\\
!\\
! !Interface:
  subroutine r_cero( t )
!
! !Input/Output parameters:
  type(temporary), intent(inout) :: t
!
! !Revision history:
! 08Aug 2015 Asdrubal Lozada
!
!eop
!----------------------------------------------------------------------

! Local variables

   t%m_ski   = 0.0_rkind  ! ~ van der Waals sigma i
   t%m_skj   = 0.0_rkind  ! ~ van der Waals sigma j
   t%m_eki   = 0.0_rkind  ! ~ van der Waals epsilon i
   t%m_ekj   = 0.0_rkind  ! ~ van der Waals epsilon j
   t%m_qki   = 0.0_rkind  ! ~ coulombic charge i
   t%m_qkj   = 0.0_rkind  ! ~ coulombic charge j
   t%m_aki   = 0.0_rkind  ! ~ buckingham a i
   t%m_akj   = 0.0_rkind  ! ~ buckingham a j
   t%m_bki   = 0.0_rkind  ! ~ buckingham b i
   t%m_bkj   = 0.0_rkind  ! ~ buckingham b j
   t%m_cki   = 0.0_rkind  ! ~ buckingham c i
   t%m_ckj   = 0.0_rkind  ! ~ buckingham c j 
   t%m_gki   = 0.0_rkind  ! ~ coupling yukawa g i
   t%m_gkj   = 0.0_rkind  ! ~ coupling yukawa g j
   t%m_mki   = 0.0_rkind  ! ~ inverse scope yukawa km i
   t%m_mkj   = 0.0_rkind  ! ~ inverse scope yukawa km j
   t%m_bo    = 0.0_rkind  ! ~ Harmonic equilibrium position
   t%m_kb    = 0.0_rkind  ! ~ Harmonic restoring constant
   t%m_ao    = 0.0_rkind  ! ~ Harmonic equilibrium angle
   t%m_ka    = 0.0_rkind  ! ~ Harmonic restoring angle
   t%m_c0    = 0.0_rkind  ! ~ expansion fourier serie constant
   t%m_c1    = 0.0_rkind  ! ~ expansion fourier serie constant
   t%m_c2    = 0.0_rkind  ! ~ expansion fourier serie constant
   t%m_c3    = 0.0_rkind  ! ~ expansion fourier serie constant
   t%m_c4    = 0.0_rkind  ! ~ expansion fourier serie constant
   t%m_c5    = 0.0_rkind  ! ~ expansion fourier serie constant
   t%m_move  = 0          ! ~ moves in Markov chain
   t%m_avol  = 0.0_rkind
   t%m_avos  = 0.0_rkind
   t%m_aent  = 0.0_rkind
   t%m_aalp  = 0.0_rkind
   t%m_aehs  = 0.0_rkind
   t%m_arho  = 0.0_rkind
   t%m_arhs  = 0.0_rkind
   t%m_aeng  = 0.0_rkind
   t%m_aens  = 0.0_rkind
   t%m_aeto  = 0.0_rkind
   t%m_aets  = 0.0_rkind
   t%m_amv   = 0
   t%m_imove = 0
  

  end subroutine r_cero
 end module m_zeros

!---------------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                       !
!---------------------------------------------------------------------------
!bop
!
! !Module: m_boxtype
!
! !Description: The derived types definition. This section refers to 
! array or pointer variables.
!\\
!\\
! !Interface:
!
module m_boxtype
!
! !Uses:
!
  use m_kind
  implicit none
!
! !Public types: 
!
!-------------------------------------------

  type box
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
! !Remarks:
! Natural alinement in derived type is advised: 
! rkind -- integer -- logical -- character  --> padding
!----------------------------------------------------------------------------
!eop
!----------------------------------------------------------------------------
  ! Global variables for simulation

   real(rkind), dimension(:,:,:), pointer :: m_site    ! Atomic coordinates vector
   real(rkind), dimension(:,:),   pointer :: m_cmass   ! Center of mass coordinate vector
   real(rkind), dimension(:,:),   pointer :: m_sigma   ! Sigma parameter
   real(rkind), dimension(:,:),   pointer :: m_epsilon ! Epsilon parameter
   real(rkind), dimension(:,:),   pointer :: m_charge  ! Coulombic charge
   real(rkind), dimension(:,:),   pointer :: m_bucka   ! Parameter A
   real(rkind), dimension(:,:),   pointer :: m_buckr   ! Parameter Rho
   real(rkind), dimension(:,:),   pointer :: m_buckc   ! Parameter C
   real(rkind), dimension(:,:),   pointer :: m_yuka    ! Parameter A
   real(rkind), dimension(:,:),   pointer :: m_rot     ! Rotation molecule
   real(rkind), dimension(:,:),   pointer :: m_rcom    ! Restart position of com
   real(rkind), dimension(:,:),   pointer :: m_rsite   ! Restatt position of sitc_rsite,errore
   real(rkind), dimension(:,:),   pointer :: m_svcom   ! Stores centre of mass
   real(rkind), dimension(:),     pointer :: m_mass    ! Mass of molecules
   real(rkind), dimension(:),     pointer :: m_edge    ! Lenght box vector
   real(rkind), dimension(:),     pointer :: m_hedge   ! Half lenght box vector
   real(rkind), dimension(:),     pointer :: m_fn      ! Topology: furthest neighbor
   real(rkind), dimension(:),     pointer :: m_fns     ! Topology: furthest neighbor in solute
   real(rkind), dimension(:),     pointer :: m_nedge
   real(rkind), dimension(:),     pointer :: m_celli   ! Inverse cell in xyz direction 
   real(rkind), dimension(:),     pointer :: m_cellsi  ! Inverse cell in xyz direction [solute]

!--------------------------------------------------
! Statistical Error

   real(rkind), dimension(:),     pointer :: m_sehpy   ! Stastistical error enthalpy
   real(rkind), dimension(:),     pointer :: m_serho   ! Statistical error density
   real(rkind), dimension(:),     pointer :: m_sevol   ! Statistical error volume
   real(rkind), dimension(:),     pointer :: m_sete    ! Statistical error total energy
!-----------------------------------------------------------------------------------------------------------
   integer,     dimension(:,:), pointer :: m_sconstr   ! Constrained coordinates sites
   integer,     dimension(:,:), pointer :: m_idpar     ! Id parameters
   integer,     dimension(:),   pointer :: m_nmol      ! Total number of molecules
   integer,     dimension(:),   pointer :: m_nsite     ! Number of sites in molecules
   integer,     dimension(:),   pointer :: m_extype    ! Exclusion array
   integer,     dimension(:),   pointer :: m_ns        ! Dummy variable for transfer sites
   integer,     dimension(:),   pointer :: m_idtype    ! Groups id identificator 
   integer,     dimension(:),   pointer :: m_tconstr   ! Constrained coordinates
   integer,     dimension(:),   pointer :: m_rconstr   ! Array sites
   integer,     dimension(:),   pointer :: m_ndiv      ! Number of partitions in simulation box
   integer,     dimension(:),   pointer :: m_cell      ! Number of cells in R3 coordinates
   integer,     dimension(:),   pointer :: m_cells     ! Number of cells in R3 coordinates [solute]
   integer,     dimension(:),   pointer :: m_ncell     ! Store id-neighbors in R3
   integer,     dimension(:),   pointer :: m_ncelsa    ! Store id neighbors in R3 [solute: case 1]
   integer,     dimension(:),   pointer :: m_ncelsb    ! Store id neighbors in R3 [solute: case 2]
   integer,     dimension(:),   pointer :: m_head      ! Head in linked cell list
   integer,     dimension(:),   pointer :: m_list      ! List in linked cell list
   integer,     dimension(:),   pointer :: m_hedsa     ! Head in linked cell list [solute: case 1]
   integer,     dimension(:),   pointer :: m_lista     ! List in linked cell list [solute: case 1]
   integer,     dimension(:),   pointer :: m_hedsb     ! Head in linked cell list [solute: case 2]
   integer,     dimension(:),   pointer :: m_listb     ! List in linked cell list [solute: case 2]
!------------------------------------------------------------------------------------------------------
   character(len=1),  dimension(:),   pointer :: m_idlabel    ! Groups id identificator 
   character(len=40), dimension(:),   pointer :: m_molname ! Groups name identificator
   character(len=2),  dimension(:,:), pointer :: m_symbol  ! Atomic symbol
!------------------------------------------------------------------------------------
! Info variables
  character(len=2), dimension(:,:), pointer :: mi_symbol
  integer, dimension(:,:), pointer :: mi_param
  character(len=5), dimension(:), pointer :: mi_typename
  integer, dimension(:), pointer :: mi_idtype
  integer, dimension(:), pointer :: mi_nsite
  integer, dimension(:), pointer :: mi_nmol
  real(rkind), dimension(:,:),   pointer :: mi_sigma   ! Sigma parameter info
  real(rkind), dimension(:,:),   pointer :: mi_epsilon ! Epsilon parameter info
  real(rkind), dimension(:,:),   pointer :: mi_charge  ! Coulombic charge info
  real(rkind), dimension(:), pointer :: s3_ii          ! Integral in lrc
  real(rkind), dimension(:), pointer :: s9_ii          ! Integral in lrc
  real(rkind), dimension(:), pointer :: s3_ij
  real(rkind), dimension(:), pointer :: s9_ij
  real(rkind), dimension(:), pointer :: lrc_ii         ! Long-range ii
  real(rkind), dimension(:), pointer :: lrc_ij         ! Long-range ij
  
!------------------------------------------------------------------------------------
  
 end type box

  end module m_boxtype

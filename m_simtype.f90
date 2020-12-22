!---------------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                       !
!---------------------------------------------------------------------------
!bop
!
! !Module: m_simtype
!
! !Description: The derived type definition. This section refers to 
!non-array or pointer variables.
!\\
!\\
! !Interface:
!
module m_simtype
!
! !Uses:
!
  use m_kind
  implicit none
!
! !Public data members:
!
  public :: simulation
!
! !Revision history:
! 06Aug Asdrubal Lozada
!
! !Remarks:
! Natural alinement in derived type is advised: 
! rkind -- integer -- logical -- character  --> padding
!---------------------------------------------------
!eop
!----------------------------------------------------------------------------
  type simulation
    ! Variables in input file 
    real(rkind)  :: m_temper  ! Temperature
    real(rkind)  :: m_press  ! Pressure
    real(rkind)  :: m_cutoff  ! Cut-off radius in conventional calculations  
    real(rkind)  :: m_transl  ! Maximum value for translation of com 
    real(rkind)  :: m_rotats  ! Maximun angle for rotation of sites
    real(rkind)  :: m_maxvol  ! Maximun value for volume displacement at the npt ensemble
    real(rkind)  :: m_cutsol  ! Cut-off radius for 2th neighbor list type
    real(rkind)  :: m_vol     ! Volume of simulation box
    real(rkind)  :: m_voli    ! Inverse volume of simulation box [1/vol]
    real(rkind)  :: m_dens    ! Total density value in simulation box
    real(rkind)  :: m_densi   ! Total inverse density value in simulation box
    real(rkind)  :: m_weight  ! Molecular weight
    real(rkind)  :: m_beta    ! Beta value [1/kB*T]
    real(rkind)  :: m_ppcon   ! Constants of transformation pressure values
    real(rkind)  :: m_rcutsq  ! Square cut-off radius value
    real(rkind)  :: m_drho    ! Density variation
    real(rkind)  :: m_rpair   ! Pair cut-off radius value [minimal bound]
    real(rkind)  :: m_rsol    ! Pair cut-off radius value size-depended [minimal bound]
    real(rkind)  :: m_enercc  ! Energy value [Conventional calculation]
    real(rkind)  :: m_enerne  ! Energy value [Convential calulation solute depended]
    real(rkind)  :: m_eneret  ! Energy value [Excluded types]
    real(rkind)  :: m_enerto  ! Eneryg value [enerne + enerto]
    real(rkind)  :: m_virial  ! Virial value [Conventional calculation]
    real(rkind)  :: m_yukawa  ! Parameter lambda in yukawa potential  
    real(rkind)  :: m_vlrc    ! Long-range corrections
    real(rkind)  :: m_eold    ! Old energy value
    real(rkind)  :: m_varold  ! Old Variance value
!------------------------------------------------------------------------------------------
!---------------------
    integer      :: mi_ntf
    integer      :: m_ipair
    integer      :: m_siter  
!---------------------
    integer      :: m_nsteps  ! Total number steps for markov chain building 
    integer      :: m_ifvol   ! Frequecy for moviment of volumen
    integer      :: m_ntf     ! Total mumber of components
    integer      :: m_mxmol   ! Total mumber of molecules
    integer      :: m_mxns    ! Maximum number of sites in molecule
    integer      :: m_mxatms  ! Total number of particles in box
    integer      :: m_iexcl   ! Total number of excluded types
    integer      :: m_ncellc  ! Total number of cells :: nc1*nc2*nc3 [with conventional cut-off]
    integer      :: m_ncelln  ! Total number of cells [non conventional cut-off ]
    integer      :: m_viz     ! Number of neighbors
    integer      :: m_icons   ! Total types constraint
    integer      :: m_carga   ! Total charge input orca 
    integer      :: m_multi   ! Spin multiplicity    
!-------------------------------------------------------------------------------------------------    
    logical      :: m_averag  ! Logical key for sampling
    logical      :: m_solute  ! Logical key for activation 
    logical      :: m_plot    ! Scritp plot   
    logical      :: m_orca    ! Junk ORCA (pipe)
    logical      :: m_pbcs    ! Active periodical bounds conditions
    logical      :: m_corr    ! Active correlation use
!--------------------------------------------------------------------------------------------------
    character(len=40) :: m_boxfile ! Linked to unit 20
    character(len=40) :: m_parfile ! Linked to unit 30
    character(len=40) :: m_inffile ! Linked to unit 35
    character(len=40) :: m_intfile ! Linken to unit 40
    character(len=40) :: m_excfile ! Linked to unit 50
    character(len=40) :: m_iorfile ! Linked to unit 60
    character(len=40) :: m_titles  ! Title of work
    character(len=6)  :: m_torca   ! Type ORCA calculation
    character(len=5)  :: m_potens  ! Choice intermolecular potential type 
    character(len=3)  :: m_ensbls  ! Thermodinamic ensemble type
    character(len=40) :: m_pairs   ! Linked to unit 38
!---------------------------------------------------------------------------------------------------    
  end type simulation

end module m_simtype

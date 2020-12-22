!---------------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                       !
!---------------------------------------------------------------------------
!bop
!
! !Module:  m_read
!
!Description: This module contains routines for reading input files and 
!memory managament.
!\\
!\\
! !Interface:
!
module m_read
!
! !Uses:
!
  use m_kind
  use m_unit
  use m_inquire
  use m_simtype
  use m_boxtype
  use m_constants
  use m_error
  use m_head
!!  use m_aleph
  implicit none
!
! !Public member functions:
!
  public :: r_sim
  public :: r_box
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
!eop
!--------------------------------------------------------------------------
  contains
!    
!bop
!
! !Iroutine: r_sim
!
! !Description: This routine uses a quasi-parse function for reading
!the input file.
!
!\\
! !Interface: 
!
  subroutine r_sim( y )
!
  implicit none

! !Input parameters:
!
  type(simulation), intent(inout) :: y
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
! !Remarks:
! Quasi-parse function. It's case sensitive: lowercase
!letters should be used.
! 
!eop
!--------------------------------------------------------------------
!   Local variables
  integer :: loc, io, aloc
  character(len=80) :: arg1, arg2
  character(len=8)  :: ival, rkval

  y%m_solute = .false.  ! Default condition

  loc  = 0
  io   = 0
  aloc = 0

  do while(io == 0)
  read(infile,'(a)',iostat=io) arg1
  if(io == 0) then
    aloc = aloc + 1
    loc = scan(arg1, '   ')
    arg2 = arg1(1:loc)
    arg1 = arg1(loc+1:)

  select case(arg2)
  case('title')
    read(arg1,*,iostat=io) y%m_titles
  case('coordinates')
    read(arg1,*,iostat=io) y%m_boxfile
    call r_request( inbox, y%m_boxfile )
  case('parameters')
    read(arg1,*,iostat=io) y%m_parfile
    call r_request( inparm,y%m_parfile )
  case('info')
    read(arg1,*,iostat=io) y%m_inffile
    call r_request( ininfo,y%m_inffile)
  case('ensemble')
    read(arg1,*,iostat=io) y%m_ensbls
    
    if(y%m_ensbls =='nvt'.or.y%m_ensbls =='npt'.or.y%m_ensbls =='mvt') then
      continue
    else
      write(*,'("Undefined ensemble type. Check the input file.")')
      stop
    end if

  case('potential')
    read(arg1,*,iostat=io) y%m_potens
  case('temperature')
    read(arg1,*,iostat=io) y%m_temper
    call r_eread( io,arg1 )
  case('pressure')
    read(arg1,*,iostat=io) y%m_press
    call r_eread( io,arg1 )
  case('cutoff')
    read(arg1,*,iostat=io) y%m_cutoff
    call r_eread( io,arg1 )
  case('pbc')
    read(arg1,*,iostat=io) y%m_pbcs
    call r_eread( io, arg1)
  case('nsteps')
    read(arg1,*,iostat=io) y%m_nsteps
    call r_eread( io,arg1 )
  case('translation')
    read(arg1,*,iostat=io) y%m_transl
    call r_eread( io,arg1 )
  case('rotation')
    read(arg1,*,iostat=io) y%m_rotats
    call r_eread( io,arg1 ) 
  case('maxvol')
    read(arg1,*,iostat=io) y%m_maxvol
    call r_eread( io,arg1 )
  case('frqvol')
    read(arg1,*,iostat=io) y%m_ifvol
    call r_eread( io,arg1 )
  case('correlation')
    read(arg1,*,iostat=io) y%m_corr
    call r_eread( io,arg1 )
  case('average')
    read(arg1,*,iostat=io) y%m_averag
    call r_eread( io,arg1 )
  case('solute')
    read(arg1,*,iostat=io) y%m_solute
    call r_eread( io,arg1 )
  case('exclude')
    read(arg1,*,iostat=io) y%m_excfile
    call r_request( inexcl, y%m_excfile ) 
  case('cutsol')
    read(arg1,*,iostat=io) y%m_cutsol
    call r_eread( io,arg1 )
  case('drho') 
    read(arg1,*,iostat=io) y%m_drho
    call r_eread( io,arg1 )
  case('lambda')
    read(arg1,*,iostat=io) y%m_yukawa
    call r_eread( io,arg1 )
  case('plot')
    read(arg1,*,iostat=io) y%m_plot
    call r_eread( io,arg1 )
  case('orca')
    read(arg1,*,iostat=io) y%m_orca
    call r_eread( io,arg1 )
  case('type')
    read(arg1,*,iostat=io) y%m_torca
    call r_eread( io,arg1 )
  case('constraints')
    read(arg1,*,iostat=io) y%m_iorfile
    call r_request(inforca,y%m_iorfile) 
  case('internal')
    read(arg1,*,iostat=io) y%m_intfile
    call r_request( intcor,y%m_intfile )
  case default
    write(*,*) 'Warning: Invalid argument in line ', arg2
    write(*,*) 'this line will be ignored. Check the input file.'
  end select
  end if ! io
  end do ! io

   write(*,*)
   write(*,*)"--------------------------------------------------------------------------------------------------"
   write(*,'("                  Conditions for simulation                ")')
   write(*,*) y%m_titles 

  select case(y%m_ensbls)
  case('nvt')
    write(*,'(" Canonical ensemble [NvT] ")')
  case('npt')
       write(*,'(" Isothermal-isobaric ensemble [NpT]")')
  case('mvt')
    write(*,'(" Grand canonical ensemble [mVT]")')
  case default
    write(*,'(" Undefined ensemble type")')
  end select

  select case(y%m_potens)
  case('m_ljc')
    write(*,*) 'Using Lennard-Jones + Coulombic Potentials'
  case('m_bc')
    write(*,*) 'Using Buckinham + Coulombic potentials'
  case('m_ljy')
    write(*,*) 'Using Lennard-Jones + Yukawa potentials'
  case default
    write(*,*) 'Undefined potential type. Check the input file.'
    stop
  end select

  write(rkval,'(f8.1)') y%m_temper
  write(*,'(" Temperature                    : ",a8," ℃ ")') adjustl(rkval)
  write(rkval,'(f8.1)') y%m_press
  write(*,'(" Pressure                       : ",a8," atm ")') adjustl(rkval)
  write(rkval,'(f8.1)') y%m_cutoff
  write(*,'(" Cut-off radius                 : ",a8," Å ")') adjustl(rkval)
  write(ival,'(i8)') y%m_nsteps
  write(*,'(" Nsteps                         : ",a8," steps ")') adjustl(ival)
  write(rkval,'(f8.1)') y%m_transl
  write(*,'(" Maximum translation value      : ",a8," Å ")') adjustl(rkval)
  write(rkval,'(f8.1)') y%m_rotats
  write(*,'(" Maximum rotation value         : ",a8," degree")') adjustl(rkval)

  if(y%m_ensbls == 'npt') then
    write(rkval,'(f8.1)') y%m_maxvol
    write(*,'(" Displacement of volume       : ",a8," Å ")') adjustl(rkval)
    write(ival,'(i8)') y%m_ifvol
    write(*,'(" Frequency for move volume    : ",a8," steps ")') adjustl(ival)
  end if

  if(y%m_solute) then
    write(*,*)
    write(*,'(" ---- Using oversize rigid body ----")')
  write(rkval,'(f8.1)') y%m_cutsol
    write(*,'(" Cut-off for o. rbody         : ",a8," Å ")') adjustl(rkval)
    write(*,'(" ----                           ----")')
  end if
  write(*,*)"------------------------------------------------------------------------------"

  end subroutine r_sim
!
!bop
! 
! !Iroutine: r_box
!
! !Description: In this routine the array variables are allocated.
!\\
!\\
! !Interface:
!
  subroutine r_box( y, x )
!
  implicit none
!
! !Input parameters:
!
  type(simulation), intent(inout) :: y
  type(box), pointer :: x
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
! !Remarks:
! Verify herarchy of pointers of pointers: each level pointers must be
! allocated before being used: allocate or nullify safely.
!
!eop
!----------------------------------------------------------------------------------
!   Local variables
  integer, parameter :: min_comp = 1
  integer :: i, j, k, n, error, iskip
  integer :: imol, isite
  real(rkind) :: s, e, q
  character(len=8)  :: ival, jval, kval
  character(len=11), dimension(4) :: labels



  integer :: ff_maximo
  integer :: bin
  integer :: iterator
!  
! The variable types: c_<name> are linked to pointer respective.
  character(len=12) :: c_x,&
                     & c_nmol, c_nsite, c_mass, c_edge, c_hedge,&
                     & c_ndiv, c_cell, c_celli, c_nedge, c_cells, c_cellsi,&
                     & c_cmass, c_site, c_ns, c_symbol, c_idpar,&
                     & c_sigma, c_epsilon, c_charge, c_fn, c_fns,&
                     & c_ncell, c_ncelsa, c_ncelsb,&
                     & c_head, c_list, c_hedsa, c_lista, c_hedsb, c_listb,&
                     & c_rot, c_molname, c_idtype, c_idlabel,&
                     & c_svcom, c_extype, c_rcom, c_rsite,&
                     & c_bucka, c_buckr, c_buckc, c_yuka, c_tconstr, c_sconstr,&
                     & c_rconstr, c_misymbol, c_miparam, c_mitypename, c_miidtype,&
                     & c_minsite, c_misigma, c_miepsilon, c_micharge,&
                     & c_s3ii, c_s9ii, c_s3ij, c_s9ij, c_lrcii, c_lrcij,&
                     & c_minmol,&
                     & c_sehpy, c_serho, c_sevol, c_sete

  labels =(/' name_type ','  id_type  ','   nmol    ','  nsites   '/)                 
  iterator = 0
!----- combinatorial label types

!-------------
  read(inbox,*) y%m_ntf
   write(*,'("                Box simulation properties               ")') 
   write(ival,'(i8)') y%m_ntf
   write(*,'(" Total number of components     : ",a8)') , adjustl(ival)

!boc
  nullify(x)
  if(.not.associated(x)) then
    allocate(x,stat=error)
    call r_ealloc( c_x, error )
  end if
!eoc
!----------
!----------
! Statistical error
  
  allocate( x%m_sehpy(kacum), stat=error ) ! Enthalpy
  call r_ealloc( c_sehpy, error )
  allocate( x%m_serho(kacum), stat=error ) ! Density
  call r_ealloc( c_serho, error )
  allocate( x%m_sevol(kacum), stat=error ) ! Volume
  call r_ealloc( c_sevol, error )
  allocate( x%m_sete(kacum), stat=error )  ! Total energy
  call r_ealloc( c_sete, error )
!-------------------------------------------------------------------------
 read(ininfo,*) y%mi_ntf
!-------------------------------------------------------------------------
 bin = (y%mi_ntf * (y%mi_ntf - 1))/2
!
  if(y%m_ntf >= min_comp) then
!---------------------------------------------------------------
   allocate(x%mi_typename(y%mi_ntf),stat=error) ! mi_typename
   call r_ealloc( c_mitypename, error )
   allocate(x%mi_idtype(y%mi_ntf),stat=error)   ! mi_idtype
   call r_ealloc( c_miidtype, error )
   allocate(x%mi_nsite(y%mi_ntf),stat=error)    ! mi_nsite
   call r_ealloc( c_minsite, error )
   allocate(x%mi_nmol(y%mi_ntf),stat=error)     ! mi_nmol
   call r_ealloc( c_minmol, error )

   read(ininfo,*) (x%mi_typename(i), i = 1, y%mi_ntf)
   read(ininfo,*) (x%mi_idtype(i), i = 1, y%mi_ntf)
   read(ininfo,*) (x%mi_nmol(i), i = 1, y%mi_ntf)
   read(ininfo,*) (x%mi_nsite(i), i = 1, y%mi_ntf)

   ff_maximo = x%mi_nsite(1)

   do i = 2, y%mi_ntf
    if( x%mi_nsite(i) > ff_maximo ) ff_maximo = x%mi_nsite(i) 
   end do

   ff_maximo = ff_maximo

!  Long-ranage correction pointers   

   allocate(x%mi_symbol(ff_maximo,y%mi_ntf),stat=error)
   call r_ealloc( c_misymbol, error )
   allocate(x%mi_param(ff_maximo, y%mi_ntf),stat=error)
   call r_ealloc( c_miparam, error)
   allocate(x%mi_sigma(ff_maximo,y%mi_ntf),stat=error)
   call r_ealloc( c_misigma, error )
   allocate(x%mi_epsilon(ff_maximo,y%mi_ntf),stat=error)
   call r_ealloc( c_miepsilon, error )
   allocate(x%mi_charge(ff_maximo,y%mi_ntf),stat=error)
   call r_ealloc( c_micharge, error )
   allocate(x%s3_ii(y%mi_ntf),stat=error)
   call r_ealloc( c_s3ii, error )
   allocate(x%s9_ii(y%mi_ntf), stat=error)
   call r_ealloc( c_s9ii, error )
   allocate(x%s3_ij(bin),stat=error)
   call r_ealloc( c_s3ij, error )
   allocate(x%s9_ij(bin),stat=error)
   call r_ealloc( c_s9ij, error )
   allocate(x%lrc_ii(y%mi_ntf),stat=error)
   call r_ealloc( c_lrcii, error )
   allocate(x%lrc_ij(bin),stat=error)
   call r_ealloc( c_lrcij, error)
!-----
   do i = 1, y%mi_ntf
    do j = 1, x%mi_nsite(i)

     read(ininfo,*) x%mi_symbol(j,i), x%mi_param(j,i)
  
      iskip = x%mi_param(j,i)
      rewind(inparm)

      do n = 1, iskip - 1
      read(inparm,*) s, e, q
      end do

      read(inparm,*) s, e, q
      x%mi_sigma(j,i) = s
      x%mi_epsilon(j,i) = e
      x%mi_charge(j,i) = q      

    end do
  end do
   
!--------------------------------------------------------------
    allocate(x%m_nmol(y%m_ntf),stat=error)  ! m_nmol
    call r_ealloc( c_nmol,error )
    allocate(x%m_nsite(y%m_ntf),stat=error) ! m_nsite
    call r_ealloc( c_nsite,error )
    allocate(x%m_mass(y%m_ntf),stat=error)  ! m_mass
    call r_ealloc( c_mass,error )
    allocate(x%m_edge(ndim),stat=error)     ! m_edge
    call r_ealloc( c_edge,error )
    allocate(x%m_hedge(ndim),stat=error)    ! m_hedge
    call r_ealloc( c_hedge,error )
    allocate(x%m_ndiv(ndim),stat=error)     ! m_ndiv
    call r_ealloc( c_ndiv,error )
    allocate(x%m_cell(ndim),stat=error)     ! m_cell
    call r_ealloc( c_cell,error )
    allocate(x%m_celli(ndim),stat=error)    ! m_celli 
    call r_ealloc( c_celli,error )
    allocate(x%m_nedge(ndim),stat=error)    !m_nedge
    call r_ealloc( c_nedge,error )
    allocate(x%m_cells(ndim),stat=error)    ! m_cells
    call r_ealloc( c_cells,error )
    allocate(x%m_cellsi(ndim),stat=error)   ! m_cellsi
    call r_ealloc( c_cellsi,error )

    read(inbox,*)(x%m_nmol(i),i=1,y%m_ntf),(x%m_nsite(i),i=1,y%m_ntf),&
               & (x%m_mass(i),i=1,y%m_ntf),(x%m_edge(i),i=1,ndim)
    read(inbox,*) m_seed

!----------------------------------------------------------------
! Set default seed in RANLUX
    if( .not. y%m_averag ) then
        m_seed = 314159265
    end if    
!----------------------------------------------------------------

! Half edges
   x%m_hedge(1) = x%m_edge(1) / 2.0_rkind
   x%m_hedge(2) = x%m_edge(2) / 2.0_rkind
   x%m_hedge(3) = x%m_edge(3) / 2.0_rkind

   y%m_mxmol = f_suma(y,x)   
   y%m_mxns  = f_maximo(y,x)  
   y%m_mxatms = 0
   imol =  0
   isite = 0

   write(ival,'(i8)') y%m_mxmol
   write(*,'(" Total number of molecules      : ",a8)') adjustl(ival)
   write(*,*)"---------------------------------------------"
   write(*,'(1x,4a11)') labels
   write(*,*)"---------------------------------------------"
   do i = 1, y%mi_ntf
   write(kval,'(i8)') x%mi_idtype(i)
   write(ival,'(i8)') x%mi_nmol(i)
   write(jval,'(i8)') x%mi_nsite(i)
   write(*,'(1x,a6,7x,a8,5x,a8,5x,a8)') adjustl(x%mi_typename(i)),adjustl(kval),adjustl(ival), adjustl(jval)
   end do
   write(*,*)"------------------------------"
   write(*,'(" Lenght side box(xyz)           : ",3f10.1," [Å]")')x%m_edge(:)
   
   write(*,*)"------------------------------"
   if( y%m_pbcs ) then
    write(*,'(" Using periodical boundary conditions ")')
   else
    write(*,'(" Without periodical boundary conditions ")')
   end if  
   write(*,*)"------------------------------"

   do i = 1, y%m_ntf
   y%m_mxatms = y%m_mxatms + x%m_nmol(i) * x%m_nsite(i) 
   end do 

   write(ival,'(i8)') y%m_mxatms
   write(*,'(" Total number of particles      : ",a8)') adjustl(ival)
   write(*,'("---------------------------------------------------------------------------------")')

   allocate(x%m_cmass(ndim,y%m_mxmol),stat=error)         ! m_cmass
   call r_ealloc( c_cmass,error )
   allocate(x%m_site(ndim,y%m_mxns,y%m_mxmol),stat=error) ! m_site
   call r_ealloc( c_site,error )

   allocate(x%m_rcom(ndim,y%m_mxmol),stat=error)          ! Restart com
   call r_ealloc( c_rcom,error )
   allocate(x%m_rsite(ndim,y%m_mxns),stat=error)          ! Restart site
   call r_ealloc( c_rsite,error )
   
   allocate(x%m_ns(y%m_mxmol),stat=error)                 ! m_ns
   call r_ealloc( c_ns,error )
   allocate(x%m_symbol(y%m_mxns,y%m_mxmol),stat=error)    ! m_symbol
   call r_ealloc( c_symbol,error )
   allocate(x%m_idpar(y%m_mxns,y%m_mxmol),stat=error)     ! m_idpar 
   call r_ealloc( c_idpar,error )

   allocate(x%m_sigma(y%m_mxns,y%m_mxmol),stat=error)     ! m_sigma
   call r_ealloc( c_sigma,error )                          
   allocate(x%m_epsilon(y%m_mxns,y%m_mxmol),stat=error)   ! m_epsilon
   call r_ealloc( c_epsilon,error )
   allocate(x%m_charge(y%m_mxns,y%m_mxmol),stat=error)    ! m_charge
   call r_ealloc( c_charge,error )
   allocate(x%m_bucka(y%m_mxns,y%m_mxmol),stat=error)     ! m_bucka
   call r_ealloc( c_bucka,error )
   allocate(x%m_buckr(y%m_mxns,y%m_mxmol),stat=error)     ! m_buckr
   call r_ealloc( c_buckr,error )
   allocate(x%m_buckc(y%m_mxns,y%m_mxmol),stat=error)     ! m_buckc
   call r_ealloc( c_buckc,error )
   allocate(x%m_yuka(y%m_mxns,y%m_mxmol),stat=error)     ! m_yuka
   call r_ealloc( c_yuka,error )

   allocate(x%m_fn(y%m_ntf),stat=error)                   ! m_fn
   call r_ealloc( c_fn,error )
   allocate(x%m_fns(y%m_ntf),stat=error)                  ! m_fns
   call r_ealloc( c_fns,error )

   allocate(x%m_ncell(y%m_mxmol),stat=error)              ! m_ncell
   call r_ealloc( c_ncell,error )
   allocate(x%m_ncelsa(y%m_mxmol),stat=error)             ! m_ncelsa
   call r_ealloc( c_ncelsa,error )
   allocate(x%m_ncelsb(y%m_mxmol),stat=error)             ! m_ncelsb
   call r_ealloc( c_ncelsb,error )
   allocate(x%m_head(0:y%m_mxmol),stat=error)             ! m_head
   call r_ealloc( c_head,error )
   allocate(x%m_list(y%m_mxmol),stat=error)               ! m_list
   call r_ealloc( c_list,error )
   allocate(x%m_hedsa(0:y%m_mxmol),stat=error)            ! m_hedsa
   call r_ealloc( c_hedsa,error )
   allocate(x%m_lista(y%m_mxmol),stat=error)              ! m_lista
   call r_ealloc( c_lista,error )
   allocate(x%m_hedsb(0:y%m_mxmol),stat=error)            ! m_hedsb
   call r_ealloc( c_hedsb,error )
   allocate(x%m_listb(y%m_mxmol),stat=error)              ! m_listb
   call r_ealloc( c_listb,error )
   allocate(x%m_rot(ndim,y%m_mxns),stat=error)            ! m_rot   
   call r_ealloc( c_rot,error )
   allocate(x%m_molname(y%m_mxmol),stat=error)            ! m_molname
   call r_ealloc( c_molname,error )
   allocate(x%m_idtype(y%m_mxmol),stat=error)             ! m_idtype
   call r_ealloc( c_idtype,error )  
   allocate(x%m_idlabel(y%m_mxmol),stat=error)             ! m_idtype
   call r_ealloc( c_idlabel,error )  


   if(y%m_solute) then
     read(inexcl,'(a)')
     read(inexcl,*) y%m_iexcl
     write(*,*)'---------------------------------------------------------'
     write(*,*)'              Oversize Rigid Body in use                 '
     write(ival,'(i8)') y%m_iexcl
     write(*,'(" Excluded components         : ",a8)') adjustl(ival)
     
     allocate(x%m_extype(y%m_iexcl),stat=error)
     call r_ealloc( c_extype,error )
     
     read(inexcl,'(a)')
     read(inexcl,*) (x%m_extype(k), k = 1, y%m_iexcl)    
     
     do i = 1, y%m_iexcl
     write(ival,'(i8)') x%m_extype(i)
     write(*,'(" Excluded types              : ",a8)') adjustl(ival)
     end do
     write(*,*)'---------------------------------------------------------'
   end if ! active

! Junk ORCA
   if(y%m_orca) then
     write(*,*)'------------------------------------------------------------'
     write(*,*)' Using Quantum program ...                                     '
     select case(y%m_torca)
     case('single')
       write(*,'("  The Single Point Calculation")')
     case('opt')
       write(*,'("  Optimization Geometry Calculation")')
     case('HF-3c')
       write(*,'("  HF Single Point + Grimm Corrections")')
     case('cxyz')
       write(*,'("  Constrained Cartesian Coordinates")')
     case('cfrag')
       write(*,'("  Constrained Cartesian Fragments")')
     case default
       write(*,'("  Warning: Undefined calculation type in ORCA")')
     end select   

     read(inforca,'(a)')
     read(inforca,*) y%m_carga, y%m_multi

     read(inforca,'(a)')
     read(inforca,*) y%m_icons  


      if(y%m_icons == 0 .and. y%m_torca == "opt") then       
       write(*,'("  Unconstraint Calculation")')
      end if  


     if(y%m_icons >=1) then
      allocate(x%m_tconstr(y%m_icons),stat=error)
       call r_ealloc(c_tconstr,error)
       read(inforca,'(a)')
       read(inforca,*) (x%m_tconstr(k), k = 1, y%m_icons)
       allocate(x%m_sconstr(y%m_mxns,y%m_icons),stat=error)
       call r_ealloc(c_sconstr,error)
       allocate(x%m_rconstr(y%m_icons),stat=error)
       call r_ealloc(c_rconstr,error)

       read(inforca,'(a)')
       read(inforca,*)(x%m_rconstr(k), k = 1, y%m_icons) 
  
     write(*,*)'----------------------------------------------------' 
     write(*,'(" Constrained type components                        ")')
     write(*,*) (x%m_tconstr(k), k = 1, y%m_icons)
     write(*,'(" Constrained Cartesian coordinates                  ")')
  
     read(inforca,'(a)')
     do k = 1, y%m_icons
     read(inforca,*) (x%m_sconstr(j,k), j = 1, x%m_rconstr(k))
     write(*,*) (x%m_sconstr(j,k), j = 1, x%m_rconstr(k))
     end do ! k
     write(*,*)'----------------------------------------------------'
    end if

   end if ! Junk ORCA

   allocate(x%m_svcom(ndim,y%m_mxmol), stat=error)
    call r_ealloc(c_svcom,error)   
!-----

! Reading the simulation box    

   do k = 1, y%m_ntf      ! Total components
     do j = 1, x%m_nmol(k) ! Total number of molecules
      imol = imol + 1
      read(inbox,*) x%m_molname(imol), x%m_idtype(imol), x%m_cmass(1,imol), x%m_cmass(2,imol), x%m_cmass(3,imol)
!-----
      do i = 1, x%m_nsite(k) ! Total number of sites in molecules
       isite = isite + 1
       x%m_ns(imol) = x%m_nsite(k)
       read(inbox,*)x%m_symbol(isite,imol), x%m_idpar(isite,imol), x%m_site(1,isite,imol),&
                    & x%m_site(2,isite,imol), x%m_site(3,isite,imol)

       iskip = x%m_idpar(isite,imol)
       rewind(inparm) ! Read parameters file

        do n = 1, iskip - 1
        read(inparm,*) s, e, q    
        end do ! n

        read(inparm,*) s, e, q
        x%m_sigma(isite,imol)   = s
        x%m_epsilon(isite,imol) = e
        x%m_charge(isite,imol)  = q

      end do ! i
     isite = 0
     end do ! j
   end do ! k
 else


   nullify(x%mi_typename)
   nullify(x%mi_idtype)
   nullify(x%mi_nsite)
   nullify(x%mi_sigma)
   nullify(x%mi_epsilon)
   nullify(x%mi_charge)
   nullify(x%mi_nmol)
   nullify(x%mi_symbol)
   nullify(x%mi_param)
   nullify(x%s3_ii)
   nullify(x%s9_ii)
   nullify(x%s3_ij)
   nullify(x%s9_ij)
   nullify(x%lrc_ii)
   nullify(x%lrc_ij)
   nullify(x%m_nmol)
   nullify(x%m_nsite)
   nullify(x%m_mass)
   nullify(x%m_edge)
   nullify(x%m_hedge)
   nullify(x%m_ndiv)
   nullify(x%m_cell)
   nullify(x%m_celli)
   nullify(x%m_cells)
   nullify(x%m_cellsi)
   nullify(x%m_cmass)
   nullify(x%m_site)
   nullify(x%m_rcom)
   nullify(x%m_rsite)
   nullify(x%m_ns)
   nullify(x%m_symbol)
   nullify(x%m_idpar)
   nullify(x%m_sigma)
   nullify(x%m_epsilon)
   nullify(x%m_charge)
   nullify(x%m_bucka)
   nullify(x%m_buckr)
   nullify(x%m_buckc)
   nullify(x%m_yuka)
   nullify(x%m_fn)
   nullify(x%m_fns)
   nullify(x%m_ncell)
   nullify(x%m_ncelsa)
   nullify(x%m_ncelsb)
   nullify(x%m_head)
   nullify(x%m_list)
   nullify(x%m_hedsa)
   nullify(x%m_lista)
   nullify(x%m_hedsb)
   nullify(x%m_listb)
   nullify(x%m_rot)
   nullify(x%m_nedge)
   nullify(x%m_molname)
   nullify(x%m_idtype)
   nullify(x%m_idlabel)
   nullify(x%m_extype)
   nullify(x%m_tconstr)
   nullify(x%m_sconstr)
   nullify(x%m_rconstr)
   nullify(x%m_svcom)
   nullify(x%m_sehpy)
   nullify(x%m_sevol)
   nullify(x%m_sete)
   nullify(x%m_serho)
 
 end if

  end subroutine r_box
!
!bop
!
! !Iroutine: f_suma
!
! !Description: This function pass a derived type pointer as argument 
!and calculate total value in lot of data.
!\\
!\\
! !Interface:
!
  function f_suma( y, x )
!
  implicit none
!
! !Input parameters:
!
  type(simulation), intent(inout) :: y
  type(box), pointer :: x
!
! !Output parameters:
!
  integer :: f_suma
!eop
!--------------------------------------------------------------------
!  Local variables
  integer :: i

  f_suma = 0
  
  do i = 1, y%m_ntf
  f_suma = f_suma + x%m_nmol(i)
  end do

  f_suma = f_suma

 end function f_suma
!
!bop
!
! !Iroutine:  f_maximo
!
! !Description: This function pass a derived type pointer as argument 
!and calculate maximum value in lot of data.
!\\
!\\
! !Interface:
!
  function f_maximo( y, x )
!    
  implicit none
!
! !Input parameters:
!
  type(simulation), intent(inout) :: y
  type(box), pointer :: x
!  
! !Output parameters:
!
  integer :: f_maximo
!eop
!---------------------------------------------------------------------
!   Local variables
  integer :: i

  f_maximo = x%m_nsite(1)

  do i = 2, y%m_ntf
  if(x%m_nsite(i) > f_maximo) f_maximo = x%m_nsite(i)
  end do

  f_maximo = f_maximo
 end function f_maximo

end module m_read

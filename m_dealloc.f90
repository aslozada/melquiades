!---------------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                       !
!---------------------------------------------------------------------------
!bop
!
!
! !Module:  m_dealloc
!
!Desription: This module contains a routine for deallocate memory space.
!This is execute to ending calculation
!\\
!\\
! !Interface:
module m_dealloc     
! !Uses:
!
   use m_boxtype  
   implicit none 
!      
! !Public member functions:
!
  private

  public :: r_dealloc
!
! !Revision history
! 07Nov 2015 Asdrubal Lozada
!
!eop
!---------------------------------------------------------------------------
  contains
!
!bop
! !Iroutine: r_dealloc
!\\
! !Interface
  subroutine r_dealloc( x )
  implicit none
  type(box), pointer :: x
  integer :: error

   

   deallocate(x%mi_typename,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating ")')
   end if    
   deallocate(x%mi_idtype,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating ")')
   end if    
   deallocate(x%mi_nsite,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating ")')
   end if    
   deallocate(x%mi_sigma,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%mi_epsilon,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%mi_charge,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    

   deallocate(x%m_nmol,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_nsite,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_mass,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_edge,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_hedge,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_ndiv,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_cell,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_celli,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_cells,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_cellsi,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_cmass,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_site,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_ns,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_symbol,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_idpar,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_sigma,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_epsilon,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_charge,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_bucka,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_buckr,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_buckc,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_yuka,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_fn,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_fns,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_ncell,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_ncelsa,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_ncelsb,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_head,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_list,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_hedsa,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_lista,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_hedsb,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_listb,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_rot,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_nedge,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
   deallocate(x%m_molname,stat=error)
   if(error /=0) then
       write(*,'(" Error deallocating")')
   end if    
! Nullify
  

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

  end subroutine r_dealloc
end module m_dealloc      

!----------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program
!----------------------------------------------------------------------
!bop
!
! !Module:  melquiades
!
! !Description: This is a main program in MELQUIADES. This block of program
!call the routines for inicialization and run the calculation energy.
!MELQUIADES is a Metropolis Monte Carlo Program for simulation of multicomponents
!systems using arbitrary potentials models
!-------------------------------------------------------------------------
! Main program : MELQUIADES
! This program use the Metropolis Monte Carlo Algorithm
! for calculate of thermodynamical properties in system 
! complex or multicomponents and multiple potential types.
! This program is stand-alone, except the use library xtc.
! Some interfaces call system have implemented with the ab initio
! program orca, the semiempirical program mopac and the plotting 
! program gnuplot.
! The current version write in the "Laboratorio de Química Teórica"
! at the "universidade federal de são carlos" by Asdrubal Lozada Blanco.
! The principal motive is that program can be adapted to multiples study 
! case, only without modification the core of program and with adition of
! other potential type or analysis tools
!
!e-mail: aslozada@gmail.com
!e-mail: alozada@ufscar.br
!---------------------------------------------------------------------------  
!
! !INTERFACE:

  program melquiades
!
! !DESCRIPTION:
!
! This is the main program of MELQUIADES.
! 
!\\
!\\
! !USES:
  use m_kind
  use m_head
  use m_simtype
  use m_boxtype
  use m_read
  use m_references
  use m_init
  use m_zeros
  use m_configuration
  use m_metropolis
  use m_dealloc
  use m_error
  use m_unit
  use m_setup
!  use nmetropolis
!  use orca
!!!*use m_references
!
! !Revision history:
! 20Aug 2015 Asdrubal Lozada
! 24Sep 2015 add explicit parallelization Asdrubal Lozada
!
!eop
!------------------------
  implicit none
 
  type(simulation) :: y
  type(box), pointer :: x
  type(temporary) :: t
!--------------------------------------------------------
  real(rkind), dimension(3) :: edge ! Static edges
! Timing
  real(rkind) :: start_time, end_time
!--------------------------------------------------------
! Transfer pointer edges to static edges
!--------------------------------------------------------
 call cpu_time(start_time)
!------------------------------------------------------
  call r_info
  call r_comline
  call r_reference
  call r_sim(y)
!---
  call r_box(y, x)
  call r_thermo(edge, t, y, x)
!----------------------------   
  call r_setup( edge, t , y, x) 
! call pipe_orca(y,x)

!------
!------------------------------------------------------
 call cpu_time(end_time) 
 call r_dealloc( x )
!----------------------------------------------------------------------------
  write(*,1000) (end_time - start_time)
  write(*,*)'--------------------------------------------------------------------'    
  call r_foot
1000 format(' Total run time: ',e20.10,' seconds')
!------------------------------------------------------

 end  program melquiades

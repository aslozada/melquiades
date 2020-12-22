!------------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                    !
!------------------------------------------------------------------------
!bop
!
! !Module:  m_unit
!
! !Description: This module assigns numerical values to unit files.
!\\
!\\
! !Interface:
!
module m_unit
!
  implicit none
!
! !Public data members:
!
  integer, public, parameter :: infile   = 10 ! Initial input file
  integer, public, parameter :: inbox    = 20 ! Box simulation    
  integer, public, parameter :: inparm   = 30 ! L-J parameters file
  integer, public, parameter :: influt   = 34 ! Store flutuations
  integer, public, parameter :: ininfo   = 35 ! Info 
  integer, public, parameter :: incorr   = 36 ! Correlation
  integer, public, parameter :: inaver   = 37 ! Store averages
  integer, public, parameter :: incomp   = 38 ! Store component division
  integer, public, parameter :: intcor   = 40 ! Internal coordinates
  integer, public, parameter :: inexcl   = 50 ! Exclusion type 
  integer, public, parameter :: inforca  = 60 ! Orca's simulation
  integer, public, parameter :: ioutorca = 70 ! Output Orca
  integer, public, parameter :: ioutplot = 80 ! Input plot program
  integer, public, parameter :: ixtc     = 90 ! xtc binary file
  integer, public, parameter :: igro     = 95 ! Auxiliar file xtc 
  integer, public, parameter :: intemp   = 11 ! Scratch file
  integer, public, parameter :: insep    = 15 ! Scratch file
  integer, public, parameter :: inplot   = 21 ! Scratch file
  integer, public, parameter :: iscript  = 22 ! Bash script to gnuplot
  integer, public, parameter :: inorca   = 31 ! Scratch file
  integer, public, parameter :: ixyz     = 41 ! xyz unit
  character(len=40), public  :: m_inpfile ! Linked to unit 10
  character(len=40), public  :: m_oorfile ! Linked to unit 70
  character(len=40), public  :: m_optfile ! Linked to unit 80
  character(len=40), public  :: m_xtcfile ! Linked to unit 90
  character(len=40), public  :: m_grofile ! Linked to unit 95 
  character(len=40), public  :: m_xyz     ! Linked to unit 41
!--------------------------------------------------------------------------
! !Revision history:
! 05Aug 2015 Asdrubal Lozada
!
! !Remarks:
! Unit numbers cannot be negative, and the range 1-99 is allowed.
!For more details on this topic, check for the Fortran90 use.
!eop
!-------------------------------------------------------------------------------------
end module m_unit

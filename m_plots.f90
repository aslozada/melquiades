!----------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                  !
!----------------------------------------------------------------------
!bop
!
! !Module: m_plots
!
! !Description: This module contains a routine for writting formatted input
!files to graphing software, e.g., GNUPLOT or XMGRACE.
!\\
!\\
! !Interface:
!
module m_plots
!  
! !Uses:
!
  use m_kind
  use m_simtype
  use m_unit, only : inplot, ioutplot, iscript
  implicit none
!
! !Public member functions:
  public :: r_plot
!
! !Revision history:
! 07Aug2015 Asdrubal Lozada
!
!eop
!----------------------------------------------------------------------
  contains
!
!bop
!
! !Iroutine:  r_plot
!
! !Description: This routine write a input file to graphing software GNUPLOT.
!The call to program is made through the internal function: system().
!The use this routine is optional and only at the end of calculation.
!A graph with a flutuating average energy value and standard deviation are obtained.
!\\
!\\
! !Interface:
!
  subroutine r_plot( y )
  implicit none
!
! !Revision history:
! 07Aug 2015 Asdrubal Lozada
!
! !See also:
! See GNUPLOT references for more details about graphic properties.
!
!eop
!----------------------------------------------------------------------
!   Local variables
  type(simulation), intent(inout) :: y
  integer, parameter :: nmax = 10000
  integer :: i, ii
  real(rkind) :: avg_ene, std_ene
  character(len=100) :: argument
  character(len=40) :: script

  script = trim(y%m_titles)//'.gpl'
 
  open(ioutplot, file='energy.dat', status='unknown')
  open(iscript, file=script, status='unknown')  

  rewind(inplot) 

  write(iscript,*)'#!/usr/bin/gnuplot'
  write(iscript,*)'#'
  write(iscript,*)'# Plot the energy per molecule value.'
  write(iscript,*)'# by Asdrubal Lozada Blanco'
  write(iscript,*)'# Laboratorio de Química Teórica - UFSCar' 
  write(iscript,*)'# 2015'
  write(iscript,*)'reset'
  write(iscript,*)'# color definitions'
  write(iscript,*)'set style line 1 lc rgb "#000000" lt 1 lw 0.5'
  write(iscript,*)'set style line 2 lc rgb "#949599" lt 1 lw 0.5'
  write(iscript,*)'set multiplot layout 2,1'
  write(iscript,*)'set ylabel "⟨E⟩ (kcal/mol)"'
  write(iscript,*)'set xrange[0:10000] '
  ! write(iscript,*)'unset key'
  write(iscript,*)'set tmargin 0.5'
  write(iscript,*)'set bmargin 0 '
  write(iscript,*)'set xtics format""'
  write(iscript,*)'#'
  write(iscript,*)'plot "energy.dat" u 1:2 w l ls 1 title "⟨E⟩"'
  write(iscript,*)'#'
  write(iscript,*)'set tmargin 0'
  write(iscript,*)'set bmargin at screen 0.15'
  write(iscript,*)'set ylabel "σ²(E) (kcal/mol)"'
  write(iscript,*)'set xlabel "N•10³"'
  write(iscript,*)'set xrange[0:10000] '
  write(iscript,*)'set xtics autofreq'
  write(iscript,*)'set xtics format"%g"'
! write(iscript,*)'unset key'
  write(iscript,*)'#'
  write(iscript,*)'plot "energy.dat" u 1:2 w l ls 2 title "σ²(E)"'
  write(iscript,*)'unset multiplot'

  do i = 1, nmax
  read(inplot,*) ii, avg_ene, std_ene
  write(ioutplot,*) ii, avg_ene, std_ene 
  end do
!
!boc
  argument ='gnuplot -persist '//trim(script)
  call system(argument)  
!eoc

  end subroutine r_plot
end module m_plots 

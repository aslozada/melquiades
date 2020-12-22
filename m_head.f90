!-------------------------------------------------------------------------------------
!                   MELQUIADES: Metropolis Monte Carlo Program                       !
!-------------------------------------------------------------------------------------
!bop
!
! !Module:   m_head
!
! !Description: This module shows the MELQUIADES's info.
!\\
!\\
! !Interface:
!
module m_head
!
! !Use:
 use m_unit
!
  implicit none
!
! !Public member functions:
!
  private
  public :: r_info, r_foot
!
! !Revision history:
! 05Aug 2015 Asdrubal Lozada
!
!eop
!-------------------------------------------------------------------------------------
  contains
!
!bop
!
! !Iroutine: r_info
!\\
!\\
! !Interface:
!
  subroutine r_info
!
  !Revision history:
! 05Aug 2015 Asdrubal Lozada
!
!eop
!-------------------------------------------------------------------------------------
  write(*,*) '--------------------------------------------------------------------------------------------------'
  write(*,*) '                                        M E L Q U I A D E S                                       '
  write(*,*) '--------------------------------------------------------------------------------------------------'
  write(*,*) '                   Metropolis Monte Carlo Program for interaction energy calculation              '
  write(*,*) '--------------------------------------------------------------------------------------------------'
  write(*,*) '                                                by                                                '
  write(*,*) '                                       Asdrubal Lozada Blanco                                     '
  write(*,*) '                                 Laboratorio de Química Teorica - LQT                             '
  write(*,*) '                                      Departamento de Química                                     '
  write(*,*) '                                 Universidade Federal de São Carlos                               '
  write(*,*) '                                              Brazil                                              '
  write(*,*) '                                               2015                                               '
  write(*,*) '--------------------------------------------------------------------------------------------------'
  write(*,*) '--------------------------------------------------------------------------------------------------'
  write(*,*)
!--------------------------------------------------------------------------
  end subroutine r_info

  subroutine r_foot
!
  !Revision history:
! 05Aug 2015 Asdrubal Lozada
!
!eop
!-------------------------------------------------------------------------------------
  implicit none

  character(8)  :: date
  character(10) :: time
  character(5)  :: zone
  integer,dimension(8) :: values

  call date_and_time(date,time,zone,values)
  call date_and_time(DATE=date,ZONE=zone)
  call date_and_time(TIME=time)
  call date_and_time(VALUES=values)

  write(*,'(" End calculation: ", a,2x,a,2x,a)') date, time, zone
  write(*,'(8i5)') values
  write(*,*) '--------------------------------------------------------------------------------------------------'
  write(*,*) 'M E L Q U I A D E S                                                                               '
  write(*,*) 'Laboratorio de Química Teorica - LQT                                                              '
  write(*,*) 'Departamento de Química                                                                           '
  write(*,*) 'Universidade Federal de São Carlos                                                                '
  write(*,*) 'Brazil                                                                                            '
  write(*,*) 'alozada@ufscar.br    2015                                                                         '
  write(*,*) '--------------------------------------------------------------------------------------------------'
!---------------------------------------------------------------------------
 end subroutine r_foot
end module m_head

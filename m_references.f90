! Write by Asdirubal Lozada
! Laboratorio de Química Teórica - LQT
! Universidade Federal de São Carlos
! 2015
 module m_references
 
  private
 
  public :: r_bibtex, r_reference

 contains 
! Structure like-bibtex 
 subroutine r_bibtex(journal,volume,pages,year,author)
 implicit none
 character(len=100), intent(in) :: journal,volume,pages,year,author
! ---------------------------
! cat
  write(*,*)trim(author)//' '//trim(journal)//' '//trim(volume)//'('//trim(year)//')'//' '//trim(pages)

 end subroutine r_bibtex

!--------------------
  subroutine r_reference
  implicit none
!-------------
! Local variables
  character(len=100) ::journal, volume, pages, year, author, doi  
!-------------------
! External subroutines
  write(*,*)'This program uses:'
  write(*,*)'High-quality pseudorandom number generator of Luscher: "RANLUX"'
   journal = 'Comp. Phys. Comm.'
   volume  = '79'
   pages   = '111'
   year    = '1994'
   author  = 'F.James'
   doi     = ''
   call r_bibtex(journal,volume,pages,year,author)
  write(*,*)
!------------------
  write(*,*)'Reference papers and books:'
!--------------------
! Format bibtex-like
!*******************
   journal = 'J.Chem.Phys.'
   volume  = '21'
   pages   = '1087'
   year    = '1953'
   author  = 'N.Metropolis, A.W.Rosenbluth, M.N.Rosenbluth, A.H.Teller, E.Teller, ' 
   doi     = ''
   call r_bibtex(journal,volume,pages,year,author)
   journal = 'J.Chem.Phys.'
   volume  = '91'
   pages   = '461'
   year    = '1989'
   author  = 'H.Flyvbjerg, H.G.Petersen, '
   doi     = ''
   call r_bibtex(journal,volume,pages,year,author)
   journal = 'J.Bra.Chem. Soc.'
   volume  = '20'
   pages   = '1541'
   year    = '2009'
   author  = 'L.C.G.Freitas, '
   doi     = ''
   call r_bibtex(journal,volume,pages,year,author)
   journal = 'Computer Simulation of Liquids'
   volume  = '1' 
   pages   = '149'
   year    = '1989'
   author  = 'M.P.Allen, D.J.Tildesley, '
   doi     = ''
   call r_bibtex(journal,volume,pages,year,author)
   journal = 'Understanding Molecular Simulation'
   volume  = '1'
   pages   = '550'
   year    = '2002'
   author  = 'D.Frenkel, B.Smit, '
   doi     = ''
   call r_bibtex(journal,volume,pages,year,author)
  write(*,*)'--------------------------------------------------------------------------------------------------'
        
  end subroutine r_reference

end module m_references 

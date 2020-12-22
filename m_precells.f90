!---------------------------------------------------------------------------
!         MELQUIADES: Metropolis Monte Carlo Program                       !
!---------------------------------------------------------------------------
!bop
!
! !Module:  m_precells
!
! !Description: This module contains a group of routines needed to build
!the Linked-Cells List
!\\
!\\
! !Interface:
!
module m_precells
!
! !Uses:
!
  use m_kind
  use m_simtype
  use m_boxtype
  use m_constants, only : ndim
  use m_cutoffs
  implicit none
!
! !Public member functions:
!
  public :: r_intcell
  public :: f_idcell
  public :: r_idneighs
  public :: r_headinit
  public :: r_headinit2
  public :: r_headinit3
  public :: r_linkscell
  public :: r_linkscell2
  public :: r_linkscell3
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
!eop
!----------------------------------------------------------------------
  contains
!
!bop
!
! !Iroutine: r_intcell
!
! !Description: This routine calculates the total
!number of cells on which the box simulation will be divided.
!\\
!\\
! !Interface:
!
  subroutine r_intcell( y, x )
!
  implicit none
!
! !Input parameters:
  type(simulation), intent(inout) :: y
  type(box), pointer :: x
!
! !Revision history:
! 06Aug 2015 Asdubal Lozada
!
!eop
!----------------------------------------------------------------------
!   Local variables
  character(len=8) :: ival


  y%m_ncelln = 0 ! Default value
  
  if( .not.y%m_solute ) then
    x%m_ndiv(1) = int(x%m_edge(1))
    x%m_ndiv(2) = int(x%m_edge(2))
    x%m_ndiv(3) = int(x%m_edge(3))
    do while((x%m_edge(1)/x%m_ndiv(1)) < y%m_cutoff .and. (x%m_edge(2)/x%m_ndiv(2)) < y%m_cutoff&
           & .and. (x%m_edge(3)/x%m_ndiv(3)) < y%m_cutoff)

     x%m_ndiv(1) = x%m_ndiv(1) - 1
     x%m_ndiv(2) = x%m_ndiv(2) - 1
     x%m_ndiv(3) = x%m_ndiv(3) - 1

     x%m_cell(1) = x%m_ndiv(1)
     x%m_cell(2) = x%m_ndiv(2)
     x%m_cell(3) = x%m_ndiv(3)

    end do

  if(x%m_cell(1) >= 1 .and. x%m_cell(2) >= 1 .and. x%m_cell(3) >= 1) then
    y%m_ncellc = x%m_cell(1) * x%m_cell(2) * x%m_cell(3)
  else
    write(*,'("Error cells number too small [<1]")')
    write(*,'("Check cut-off")')
    stop
  end if
  
  x%m_celli(1) = x%m_cell(1)/x%m_edge(1)
  x%m_celli(2) = x%m_cell(2)/x%m_edge(2)
  x%m_celli(3) = x%m_cell(3)/x%m_edge(3)

  write(*,'("---------------------------------------------------------")')
  write(*,'("          Linked-Cells method specfications")')
  write(ival,'(i8)') y%m_ncellc
  write(*,'(" Total number of cells          : ",a8)') adjustl(ival)
  write(*,'(" Size of subcells               : ",3f8.1," [Å]")') x%m_edge(:)/x%m_cell(:)
  write(*,'("---------------------------------------------------------")')

 else
   x%m_ndiv(1) = int(x%m_edge(1))
   x%m_ndiv(2) = int(x%m_edge(2))
   x%m_ndiv(3) = int(x%m_edge(3))

   do while((x%m_edge(1)/x%m_ndiv(1)) < y%m_cutoff .and. (x%m_edge(2)/x%m_ndiv(2)) < y%m_cutoff&
           & .and. (x%m_edge(3)/x%m_ndiv(3)) < y%m_cutoff)

   x%m_ndiv(1) = x%m_ndiv(1) - 1
   x%m_ndiv(2) = x%m_ndiv(2) - 1
   x%m_ndiv(3) = x%m_ndiv(3) - 1

   x%m_cell(1) = x%m_ndiv(1)
   x%m_cell(2) = x%m_ndiv(2)
   x%m_cell(3) = x%m_ndiv(3)

   end do

   if(x%m_cell(1) >= 1 .and. x%m_cell(2) >= 1 .and. x%m_cell(3) >= 1) then
     y%m_ncellc = x%m_cell(1) * x%m_cell(2) * x%m_cell(3)
   else
  
     write(*,'("Error cells number too small [<1]")')
     write(*,'("Check cut-off")')
     stop
   end if

   x%m_celli(1) = x%m_cell(1)/x%m_edge(1)
   x%m_celli(2) = x%m_cell(2)/x%m_edge(2)
   x%m_celli(3) = x%m_cell(3)/x%m_edge(3)

   write(*,'("---------------------------------------------------------")')
   write(*,'("            Linked-Cells method specfications")')
   write(*,'("           corresponding to components not excluded      ")') 
   write(ival,'(i8)') y%m_ncellc
   write(*,'(" Total number of cells         : ",a8)') adjustl(ival)
   write(*,'(" Size of subcells              : ",3f8.1," [Å]")') x%m_edge(:)/x%m_cell(:)
   write(*,'("---------------------------------------------------------")')

   x%m_ndiv(1) = int(x%m_edge(1))
   x%m_ndiv(2) = int(x%m_edge(2))
   x%m_ndiv(3) = int(x%m_edge(3))
 
   do while((x%m_edge(1)/x%m_ndiv(1)) < y%m_cutsol .and. (x%m_edge(2)/x%m_ndiv(2)) < y%m_cutsol&
           & .and. (x%m_edge(3)/x%m_ndiv(3)) < y%m_cutsol)

   x%m_ndiv(1) = x%m_ndiv(1) - 1
   x%m_ndiv(2) = x%m_ndiv(2) - 1
   x%m_ndiv(3) = x%m_ndiv(3) - 1

   x%m_cells(1) = x%m_ndiv(1)
   x%m_cells(2) = x%m_ndiv(2)
   x%m_cells(3) = x%m_ndiv(3)

   end do

   if(x%m_cells(1) >= 1 .and. x%m_cells(2) >= 1 .and. x%m_cells(3) >= 1) then
      y%m_ncelln = x%m_cells(1) * x%m_cells(2) * x%m_cells(3)
   else
    write(*,'("Error cells number too small [<1]")')
    write(*,'("Check cut-off")')
    stop
   end if

   x%m_cellsi(1) = x%m_cells(1)/x%m_edge(1)
   x%m_cellsi(2) = x%m_cells(2)/x%m_edge(2)
   x%m_cellsi(3) = x%m_cells(3)/x%m_edge(3)

   write(*,*)'---------------------------------------------------------'
   write(*,*)'            Linked-Cell method specfications'
   write(*,*)'           corresponding to oversize compounds           ' 
   write(ival,'(i8)') y%m_ncelln
   write(*,'(" Total number of cells         : ",a8)') adjustl(ival)
   write(*,'(" Size of subcells              : ",3f8.1," [Å]")') x%m_edge(:)/x%m_cells(:)
   write(*,*)'---------------------------------------------------------'
   
  end if ! End active

 end subroutine r_intcell
!
!bop
!
! !Iroutine: idcell
!
! !Description: This function calculates the index cell.
!\\
!\\
! !Interface:
!
  function f_idcell( com, y, x )
  implicit none
!
! !Input parameters:
  type(simulation), intent(inout) :: y
  type(box), pointer :: x
  real(rkind), dimension(:), intent(in) :: com
!  
! !Output parameters:  
  integer :: f_idcell
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
!eop
!------------------------------------------------
!
!   Local variables
  real(rkind) :: xx, yy, zz
  integer :: ix, iy, iz

  xx = com(1) + x%m_hedge(1)
  yy = com(2) + x%m_hedge(2)
  zz = com(3) + x%m_hedge(3)

  if( .not.y%m_solute ) then 
    ix = int(xx * x%m_celli(1))
    iy = int(yy * x%m_celli(2))
    iz = int(zz * x%m_celli(3))
!boc
  f_idcell = ix + iy * x%m_cell(1) + iz * x%m_cell(1) * x%m_cell(2)
!eoc
  else

  ix = int(xx * x%m_cellsi(1))
  iy = int(yy * x%m_cellsi(2))
  iz = int(zz * x%m_cellsi(3))

  f_idcell = ix + iy * x%m_cells(1) + iz * x%m_cells(1) * x%m_cells(2)

  end if ! Active

  end function f_idcell
!
!bop
!
! !Iroutine: r_idneighs
!
! !Description: This routine specifies and saves the list of neighbors.
!\\
!\\
! !Interface:
!
  subroutine r_idneighs( idcell, y, x)
  implicit none
! 
! !Input parameters:
!
  type(simulation), intent(inout) :: y
  type(box), pointer :: x
  integer,intent(in) :: idcell  
!
! !Revision history:
! 06Aug 2015 Asdrubal Lozada
!
!eop
!-----------------------------------------------------------
!   Local variables
  integer :: jd_cell, iter, k
  integer :: cell_ix, cell_iy, cell_iz
  integer :: cell_jx, cell_jy, cell_jz

  integer, dimension(ndim*ndim*ndim) :: nx, ny, nz
  data nz/-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1/
  data ny/-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1/
  data nx/-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1/
    
  iter = 0

  if( .not.y%m_solute ) then
    cell_iz = idcell / (x%m_cell(1)*x%m_cell(2))
    cell_iy = mod((idcell/x%m_cell(1)),x%m_cell(2))
    cell_ix = mod(idcell,x%m_cell(1))

   do k = 1, y%m_viz
   
    cell_jz = cell_iz + nz(k)
    if( cell_jz < 0 ) then
      cell_jz = cell_jz + x%m_cell(3)
    else if(cell_jz >= x%m_cell(3)) then
      cell_jz = cell_jz - x%m_cell(3)
    end if

    cell_jy = cell_iy + ny(k)
    if( cell_jy < 0 ) then
      cell_jy = cell_jy + x%m_cell(2)
    else if(cell_jy >= x%m_cell(2)) then
      cell_jy = cell_jy - x%m_cell(2)
    end if

    cell_jx = cell_ix + nx(k)
    if( cell_jx < 0 ) then
      cell_jx = cell_jx + x%m_cell(1)
    else if( cell_jx >= x%m_cell(1) ) then
      cell_jx = cell_jx - x%m_cell(1)
    end if

   jd_cell = cell_jx + cell_jy * x%m_cell(1) + cell_jz * x%m_cell(1) * x%m_cell(2)

   iter = iter + 1

   x%m_ncell(iter) = jd_cell ! Saves the list of neighbors
  
  end do !k
  else
   
!boc   
   cell_iz = idcell / (x%m_cells(1)*x%m_cells(2))
   cell_iy = mod((idcell/x%m_cells(1)),x%m_cells(2))
   cell_ix = mod(idcell,x%m_cells(1))
!eoc   

  do k = 1, y%m_viz
  cell_jz = cell_iz + nz(k)
    if( cell_jz < 0 ) then
      cell_jz = cell_jz + x%m_cells(3)
    else if( cell_jz >= x%m_cells(3) ) then
      cell_jz = cell_jz - x%m_cells(3)
    end if

  cell_jy = cell_iy + ny(k)
    if( cell_jy < 0 ) then
    cell_jy = cell_jy + x%m_cells(2)
    else if( cell_jy >= x%m_cells(2) ) then
    cell_jy = cell_jy - x%m_cells(2)
    end if

    cell_jx = cell_ix + nx(k)
    if(cell_jx < 0) then
    cell_jx = cell_jx + x%m_cells(1)
    else if(cell_jx >= x%m_cells(1)) then
    cell_jx = cell_jx - x%m_cells(1)
    end if

    jd_cell = cell_jx + cell_jy * x%m_cells(1) + cell_jz * x%m_cells(1) * x%m_cells(2)

   iter = iter + 1

   x%m_ncelsa(iter) = jd_cell ! Store neighbor list common

  end do !k

  end if ! Active

end subroutine r_idneighs
!
!bop
!
! !Iroutine: r_linkscell
!
! !Description: In this routine the Linked-Cells List are builded. 
!Head is at the initial position in the list. The last position have zero value.
!\\
!\\
! !Interface:
!
 subroutine r_linkscell( y, x )
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
!eop
!-----------------------------------------------------------
!   Local variables
  integer :: idcell
  integer :: j
  real(rkind), dimension(3) :: com

  do j = 1, y%m_mxmol 
   com(1) = x%m_cmass(1,j)
   com(2) = x%m_cmass(2,j)
   com(3) = x%m_cmass(3,j)
!boc
   idcell = f_idcell( com, y, x )

   x%m_list(j) = x%m_head(idcell)
   x%m_head(idcell) = j
!eoc
  end do ! j

 end subroutine r_linkscell

  subroutine r_linkscell2(y,x)
  implicit none
! Dummy arguments
  type(simulation), intent(inout) :: y
  type(box), pointer :: x
! Local variables
  integer :: idcell, j

  real(rkind), dimension(3) :: com

  do j = 1, y%m_mxmol
 
   com(1) = x%m_cmass(1,j)
   com(2) = x%m_cmass(2,j)
   com(3) = x%m_cmass(3,j)

   idcell = f_idcell( com, y, x )

   x%m_list(j) = x%m_head(idcell)
   x%m_head(idcell) = j

   end do ! j
 end subroutine r_linkscell2
  
 subroutine r_linkscell3(y,x)  

  implicit none
! Dummy arguments
  type(simulation), intent(inout) :: y
  type(box), pointer :: x
! Local variables
  integer :: idcell, j

  real(rkind), dimension(3) :: com

   do j = 1, y%m_mxmol  

    com(1) = x%m_cmass(1,j)
    com(2) = x%m_cmass(2,j)
    com(3) = x%m_cmass(3,j)

    idcell = f_idcell( com, y, x )

    x%m_listb(j) = x%m_hedsb(idcell)
    x%m_hedsb(idcell) = j 
  
  end do ! j

 end subroutine r_linkscell3
!
!bop
!
! !Iroutine: r_headinit
!
! !Description: This routine initializes the head values
!in Linked-Cell List.
!\\
!\\
! !Interface:
!
  subroutine r_headinit( y, x )
  implicit none
!  
! !Input parameters:
  type(simulation), intent(inout) :: y
  type(box), pointer :: x
!
! !Revision history: 
! 06Aug 2015 Asdrubal Lozada
!
!eop
!------------------------------------------------------
!  Local variables
  integer :: i
!boc
  do i = 0, y%m_ncellc-1   
   x%m_head(i) = 0  ! Check
  end do
!eoc 

 end subroutine r_headinit

 subroutine r_headinit2( y, x )
 implicit none
! Dummy arguments
  type(simulation), intent(inout) :: y
  type(box), pointer :: x
! Local variables
  integer :: i

   do i = 0, y%m_ncellc-1
    x%m_head(i) = 0
   end do ! i     

  end subroutine r_headinit2

 subroutine r_headinit3( y, x )
 implicit none
! Dummy arguments
  type(simulation), intent(inout) :: y
  type(box), pointer :: x
! Local variables
  integer :: i
  
  if(y%m_ncelln >= 1) then
   do i = 0, y%m_ncelln-1
    x%m_hedsb(i) = 0
   end do ! i    
  end if  ! ncells

  end subroutine r_headinit3

 end module m_precells

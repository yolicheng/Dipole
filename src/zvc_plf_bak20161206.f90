!***********************************************************************
!     ****   zvc_plf   ****
! 2016-12-06 08:59:54 discard! Very bad option.
! For Planar normal LF problem, we compute the zero velocity curves in the
! first quadrant with x>0,y>0, since H is dependent only on x and y 
! and use the symmetry to obtain the full curve. 

! so for a given energy level, we only need to discretize x in a grid, and compute y 
! as a function of H and x, 
!  y = sqrt( 4/(h-3x^2) - x^2 )

!       Input Variables 
!  h0           the energy level to be studies 
!  xf           the  upper bound of x-axis 

!  Finally Revised by Yu -- 20160
!***********************************************************************
subroutine zvc_plf( h0, xf)

use dp_mod
implicit none
 
! Input  and Output Declaration   
real(kind=dp), intent(in)      ::  h0,  xf 
!real(kind=dp), intent(out)     ::    
 
! Local Variable
integer     ::  i, nx, isy2 
real(kind=dp)  :: dx, x, y, x2, y2
  
  dx = 1.d-4
  nx = int( xf / dx)
  
  open(333, file='./dat/zvc_plf.dat')
  write(333, *) '# H0 = ', h0
  
  x = 0.d0
  dx = 1.d-4 
  isy2  = 0
  
  ! this is for sure not a good way to do zero velocity curve...
  do 
    if(x > .6214d0 .and. x < 0.622d0) then 
      dx = 1.d-7
      print*; read* !ck
    else 
      dx = 1.d-4
    endif   
    
    if(x > xf) exit 
    
    x2 = x**2 
    y2 = 4.d0 / (h0 - 3.d0*x2)**2 - x2 
    x = x + dx
    
    if(y2 < 0.d0) then 
      isy2 = isy2 + 1
      if(isy2 == 1)  write(333,*)
      print*, 'y2<0', y2; print*; read*; !ck
      cycle
    endif    
    
    y = sqrt(y2)
    write(333, *) x, y
    
    write(*, *) x, y; ! print*; read* !ck
    
    
  enddo 
  close(333)  
  return  
end subroutine zvc_plf


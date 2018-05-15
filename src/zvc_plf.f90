!***********************************************************************
!     ****   zvc_plf   ****
! For Planar normal LF problem, produce the isoline data in y for contour plot  

!       Input Variables 
!  xlim(3)     x0:dx:xf, the lower bound, increment and upper bound for x
!  ylim(3)     y0:dy:yf, the lower bound, increment and upper bound for y

!  Output file: ./dat/zvc_plf/zvc_plf.dat

!  Finally Revised by Yu -- 20161206 -- CKD!            
!***********************************************************************
subroutine zvc_plf( xlim, ylim)

use dp_mod
use plf_mod
implicit none
 
! Input  and Output Declaration   
real(kind=dp), intent(in)     ::  xlim(3), ylim(3) 

! Local Variable
integer :: i, j, nx, ny
real(kind=dp)  :: x0, dx, xf, y0, dy, yf,  x, y, pv(4), cj
  
  x0 = xlim(1); dx = xlim(2); xf = xlim(3)
  nx = ceiling( (xf - x0) / dx )
  y0 = ylim(1); dy = ylim(2); yf = ylim(3)
  ny = ceiling( (yf - y0) / dy )
  
!  print*, 'CHECK zvc_plf: xlim, ylim, nx, ny'
!  print*, xlim, ylim, nx, ny; print*; read*
  
  open(333, file='./dat/zvc_plf/zvc_plf.dat', status='replace', access='append')
  write(333, *) '#  x       y       H'
  pv = 0.d0
  
  do i = 1, ny, 1
    y = y0 + (i-1) * dy
    
!    if(dabs(y) < 1.d-1 ) cycle
    
    do j = 1, nx, 1
      x = x0 + (j-1) * dx
!      if(dabs(x) < 1.d-1 ) cycle
      
      pv(1:2) = (/x, y/)
      
      call gr_cjplf( pv, cj)
      write(333, *) x, y, cj 
             
    end do
    
    write(333, *) ! add a blank line
  end do  
  
  close(333)  
  return  
end subroutine zvc_plf


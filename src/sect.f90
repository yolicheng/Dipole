subroutine sect(y,ind, y0, g,dg)

! the surface of section, defined by g =  y(ind) - y0 

! 	Input parameters:
! y(*)      point
! ind  		the index of the component of y to be used as poincare section  
! y0 		the poincare section to be specified as y(ind) - y0

!  	Output parameters:
! g 		funciton that equated to 0 gives the surface of section
! dg(*) 	gradient of function g

! 20160119 - modify to general form, poincare section could be any general plane difined by 
!            g =  y(ind) - y0, previously it is only y=0 plane



implicit none 
integer, parameter :: dp = kind(1.d0)

integer, intent(in) :: ind
real(kind=dp), intent(in) :: y(42), y0
real(kind=dp), intent(out) :: g, dg(42)

!  Local variables
integer :: i


g =  y(ind) - y0  

!print*, 'ck-- inside sect' 
!print*, 'g,ind, y0,y(ind),  y(1:6)', g,ind, y0,y(ind),  y(1:6)

do  i = 1,42
  dg(i) = 0
enddo 

dg(ind) = 1.d0

return       
end

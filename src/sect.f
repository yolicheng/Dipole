! the surface of section, defined by g =  y(ind) - y0
! 	Input parameters:
! y(*)      point
! ind  		the index of the component of y to be used as poincare section  
! y0 		the poincare section to be specified as y(ind) - y0

!  	Output parameters:
! g 		funciton that equated to 0 gives the surface of section
! dg(*) 	gradient of function g

subroutine sect(y,ind, y0, g,dg)

implicit none 
integer, parameter :: dp = kind(1.d0)

integer, intent(in) :: ind
real(kind=dp), intent(in) :: y(42), y0
real(kind=dp), intent(out) :: g, dg(42)

integer :: i

g = y(ind) - y0
       
do  i = 1,42
  dg(i) = 0
enddo 

dg(ind) = 1.d0

return       
end

subroutine zerovc(x0,xf,,ind,y0, cj, fzvc)
! compute the planar velocity curve, with ind-th component specified to y0
! The idea is like ezplot function in Matlab

! Provided that vx=vy=vz=0, for a given energy, we need to fix one component for example z=0
!  and use do loop to compute for a given range of value in x, the value of y that gives us the given energy.


! 	Input parameters:
! x0-xf		the range of value for x
! ind  		the index of the component of y to be used as poincare section  
! y0 		the poincare section to be specified as y(ind) - y0
! cj 		the given energy
! fzvc		the data file of x-y in the domain [x0,xf] on zero velocity curve, for a given energy on plane 

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

! use include to split this module into different files, with gr_rk78, deriv, gr_lf in external fiels
!***********************************************************************
!  for the case of lorentz force, compute the energy integral
!  c = 3x²-z² - sgn* 2(y²+z²)/r³ - (dotx² + doty² + dotz²), sgn:1 if q/m>0; -1 if q/m<0

!	input
!  y(6) 	the initial state
!  sgn   	1 if q/m>0; -1 if q/m<0
!***********************************************************************
subroutine gr_cj(y, cj)
        
use lf_mod, only : sgn !, dp
implicit none 

integer, parameter:: dp = kind(1.d0)
real(kind=dp), intent(in) :: y(6)
real(kind=dp), intent(out) :: cj 
 
real(kind=dp) :: y1, y2, y3, r, r3, dv2   

print*, 'test if this is called?'
read(*,*) 
y1=y(1)**2
y2=y(2)**2
y3=y(3)**2
        
r = dsqrt(y1 + y2 + y3) 
r3 = r**3
dv2 = y(4)**2 + y(5)**2 + y(6)**2 

 cj = 3*y1 - y3 - sgn * 2 * (y2+y3)/r3 - dv2 

print*, 'cj', cj 
read(*,*)
return
end subroutine gr_cj

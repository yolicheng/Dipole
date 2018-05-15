subroutine cjlf_vy(y, ind, cj)
!  for the case of lorentz force, given the energy level, with the initial state as 
!  (x,0,z, 0,vy,0), keep the value of x and z, change the value of vy to get the required cj

!  c = 3x²-z² - sgn* 2(y²+z²)/r³ - (dotx² + doty² + dotz²), sgn:1 if q/m>0; -1 if q/m<0

!	Input
!  y(6) 	the initial state and also the updated one with energy as cj
!  cj 		the desired energy
!  ind 		the index of component of y to be modified, all the other components remain the same


! 	Global variable 
!  sgn   	1 if q/m>0; -1 if q/m<0 for the module lf_mod

        
use lf_mod, only : sgn, dp

implicit none 

real(kind=dp), intent(inout) :: y(6)
real(kind=dp), intent(in) :: cj
integer, intent(in) :: ind

real(kind=dp) :: y1, y2, y3, r, r3, vy, aaa
        
y1=y(1)**2
y2=y(2)**2
y3=y(3)**2
        
r = dsqrt(y1 + y2 + y3) 
r3 = r**3
        
!       dv2 = y(4)**2 + y(5)**2 + y(6)**2 
!       cj = 3*y1 - y3 - sgn * 2 * (y2+y3)/r3 - dv2 

vy = dsqrt( 3*y1 - y3 - sgn*2*(y2+y3)/r3 - cj - y(4)**2 - y(6)**2 ) 

!print*, 'check vy', vy
!read(*,*) aaa
if (y(5) .ge. 0.d0) then   ! keep the sign of vy
  y(5) = vy
else 
  y(5) = -vy  
endif

return
end


subroutine lf2cj(y, ind, cj, iscj)
!  for the case of lorentz force, given the energy level, with the initial state as 
!  (x,0,z, 0,vy,0), keep the value of x and z, change the value of vy to get the required cj
!  always keep the positive square root for v(ind-3)

!  c = 3x²-z² - sgn* 2(y²+z²)/r³ - (dotx² + doty² + dotz²), sgn:1 if q/m>0; -1 if q/m<0

!	Input
!  y(6) 	the initial state and also the updated one with energy as cj
!  cj 		the desired energy
!  ind 		the index of component of y to be modified, all the other components remain the same
!  	 	for my problem, only vx,vy,vz are allowed to be modified, range of value is [4,5,6]
! 		the initial value of y(ind) is 0


! 	Global variable 
!  sgn   	1 if q/m>0; -1 if q/m<0 for the module lf_mod

        
use lf_mod, only : sgn, dp

implicit none 

real(kind=dp), intent(inout) :: y(6)
real(kind=dp), intent(in) :: cj
integer, intent(in) :: ind
integer, intent(out) :: iscj

real(kind=dp) :: y1, y2, y3, r, r3, dv2, v2, aaa
 
iscj = 1        
y1=y(1)**2
y2=y(2)**2
y3=y(3)**2
        
r = dsqrt(y1 + y2 + y3) 
r3 = r**3
        
dv2 = y(4)**2 + y(5)**2 + y(6)**2 
!       cj = 3*y1 - y3 - sgn * 2 * (y2+y3)/r3 - dv2 

!print*, 'lf2cj,  ind, y(ind), y', ind, y(ind), y

v2 = 3*y1 - y3 -  sgn*2*(y2+y3)/r3 - cj - dv2

if( v2 < 0 )  then 
  print*, 'v^2<0', v2
  iscj = 0
!  read(*,*)  
endif  

y(ind)  = dsqrt( 3*y1 - y3 - sgn*2*(y2+y3)/r3 - cj - dv2 ) 

!print*, 'post-lf2cj,  ind, y(ind), y', ind, y(ind), y

!print*, 'check vy', vy
!read(*,*)  
! 
!if (y(ind) .ge. 0.d0) then   ! keep the sign of vy
!  y(ind) = v 
!else 
!  y(ind) = -vy  
!endif

return
end


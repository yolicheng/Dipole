!***********************************************************************
!     ****   gamma   ****
! the inverse of map PI, \gamma: (y, vy)_0 --> (x,y,vx,vy)_0
! given the initial y0 and vy0, compute the full state (x,y,vx,vy)_0
! for a prescribed energy level, in princile we need to compute the velocity in the direction of ind-th component in pvin 

!       Input Variables 
!  pvin     dimension n, the knowns components of the state 
!  n        the size of pvin 
!  h0       prescribed energy level 
! ind,p0    pv(ind) = p0, the fixed value of ind-th component    


!       Output Variables 
!  isv     flag, 1: succeed to get pv, 0: fail, v**2 < 0
!  pv      dimension n+2, the  full state, used as initial state for the integration of vector field. 

!  Routine Used:
!     cj2v : compute the missing velocity 

!  Finally Revised by Yu -- 20160831
!***********************************************************************
subroutine gamm( pvin, n, ind_fun, h0, ind, p0, pv, isv, cj2v)

use dp_mod
implicit none
 
! Input  and Output Declaration   
integer, intent(in)      ::  n, ind_fun(n), ind
integer, intent(out)     ::  isv  
real(kind=dp), intent(in)      :: pvin(n), h0, p0   
real(kind=dp), intent(out)     :: pv(n+2)   
 
external :: cj2v  
! Local Variable
integer :: ind_vel 
!real(kind=dp)  :: cj 


  ! the initial position and velocity
  pv = 0.d0 
  pv(ind_fun) = pvin  ! array operation 

  ! the component of x and vx 
  pv(ind) = p0 
  
  
  ! given the energy and other components, compute the velocity (vx ) in ind-th dirction (x)

  ind_vel  =  ind + n/2 + 1
  call cj2v(pv, h0, ind_vel, isv) ! TODO: debug 
  
!  if(isv == 0) read*  ! ckd
  
!  print*, 'pv: ', pv 
!  print*; read*  

  return  
end subroutine gamm  

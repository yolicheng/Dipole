!***********************************************************************
!     ****   gamma   ****
! the inverse of map gamma, with 
!    \gamma:     R^(n-2) --> R^n, i.e.,  (x,vx, y, vy)_0 --> (x,y,z, vx,vy, vz)_0
!    \gamm_inv:  R^n --> R^(n-2), i.e.,  (x,y, z, vx,vy, vz)_0 --> (x, vx, y, vy)_0 

! Given the full state, return its independent components apart from the Poincare map 
! for Fourier analysis 

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
subroutine gamm_inv( pvin, n, ind_fun, pv)

use dp_mod
implicit none
 
! Input  and Output Declaration   
integer, intent(in)            :: n, ind_fun(n)
real(kind=dp), intent(in)      :: pvin(n+2)   
real(kind=dp), intent(out)     :: pv(n)   
 
  pv = pvin(ind_fun)
  
  print*, 'check Gamma^-1: '
  print*, 'Input: ', pvin 
  print*, 'Output:', pv 
  print*; read*
  return  
end subroutine gamm_inv 

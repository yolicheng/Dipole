!***********************************************************************
!     ****   std_map   ****
! This is the standard  map, which is area-preserved, used to check the computation of invariant curve

! Todo: use this map the final decide, what is the necessary input and output of a map 
!       for Poincare map, which is special case, use module to pass the relevent parameters 

!  \hat p = p + K * sinx
!  \hat x = x + \hat p 

! we take K =  -0.5, and (0, 0) is an elliptic point, 

!       Input Variables 
!  ptin     the initial point (x, p)
!  isdp     1: compute the differential, 0: only compute the map 

!       Output Variables 
!  ptf      the image of ptin under the standard map 
!  dmap     the differential of the standard map w.r.t. (x, p)       

!  Routine Used:
!     

!  Finally Revised by Yu -- 20160
!***********************************************************************
!subroutine poinc_prtbp( pvin, ind, p0, h, n, ndim, ind_fun, dir, imax, tmax, &
!                        tf, pvf, dp_dx)

subroutine std_map( ptin, ptf,  isdp, dpdx)

! to keep the same form as Poincare map, we put the second parameter as k 
!   while in Poincare map, we have the return time in the second parameter 

use dp_mod
use pi_mod
implicit none
 
! Input  and Output Declaration   
integer, intent(in)     ::  isdp
real(kind=dp), intent(in)      ::  ptin(2) !, k   
real(kind=dp), intent(out)     ::  ptf(2), dpdx(2, 2)   
 
! Local Variable
real(kind=dp)  :: k, x0, p0, x1, p1
  
  k = -0.5d0 
  
  x0 = ptin(1)
  p0 = ptin(2)
  
!  p0 = dmod(p0, pi2)
!  x0 = dmod(x0, pi2)
  
  !  \hat p = p + K * sinx
  !  \hat x = x + \hat p 
  p1 = p0 + k*dsin(x0)  ! \hat p 

!  p1 = dmod(p1, pi2)
  
  x1 = x0 + p1     ! \hat x
!  x1 = dmod(x1, pi2)
  
  ptf = (/x1, p1/)
  
!  print*, 'check std_map: (x0, p0) = ', ptin   
!  print*, '(x1, p1) = ', ptf
!  print*; read* 
   
  ! dpf / dx = k*cos(x);   dpf / dp  = 1.d0 
  
  ! dxf / dx = 1.d0 + d pf / dx = 1.d0 + k*cos(x)
  ! dxf / dp = d pf / dp = 1.d0 
  if(isdp == 1) then
    dpdx(2, 1) = k*dcos(x0)  ! dpf / dx
    dpdx(2, 2) = 1.d0        ! dpf / dp 
  
    dpdx(1, 1) = 1 + dpdx(2, 1)  ! dxf / dx  = 1 + dpf / dx
    dpdx(1, 2) = dpdx(2, 2)      ! dxf / dp  = dpf / dp
  endif 
  
  return  
end subroutine std_map


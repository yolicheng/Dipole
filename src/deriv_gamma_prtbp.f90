!***********************************************************************
!     ****   deriv_gamma   ****
! Compute the differential of \gamma( (x,y,vx,vy)_0 ) = (y, vy)_0 
! w.r.t.  (x,y,vx,vy)_0 


!       Input Variables 
!  a            

!       Output Variables 
!  a            

!  Routine Used:
!     

!  Finally Revised by Yu -- 20160
!***********************************************************************
subroutine deriv_gamma_prtbp( pv, ind, dvind_dx, dvdx)

use dp_mod
implicit none
 
! Input  and Output Declaration   
real(kind=dp), intent(in)      ::  pv(4)  
real(kind=dp), intent(out)     ::  dv_dx(2)  
 
! Local Variable
real(kind=dp)  ::
 
external :: dvind_dx 
  
  ! TODO: we need also compute the differential of vx_0  w.r.t. (x,y,vx,vy)_0 
  call dvind_dx(pv0, ind, dv_dx)
  
  
  ! compute the differential of the map \gamma:  (y, vy)_0 --> (x,y,vx, vy)_0 
  dgamm_dx = 0.d0
  dgamm_dx(ind_vel, :) = dv_dx 
 
  k = 1
  do i = 1, ndim, 1
    if(i == ind .or. i == ind_vel ) cycle ! zero
    dgamm_dx(i, k) = 1.d0 
    k = k+1
  end do
  
  ! check dgamm_dx 
  print*,'dgamm_dx'
  do i = 1, ndim, 1
    print*, dgamm_dx(i, :)  
  end do
  print*; !read*  
  return  
end subroutine deriv_gamma


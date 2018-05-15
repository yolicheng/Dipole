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
subroutine deriv_gamma( pv, ndim, ind, dgamm_dx, dvind_dx)

use dp_mod
implicit none
 
! Input  and Output Declaration  
integer, intent(in)     ::  ndim, ind 
real(kind=dp), intent(in)      ::  pv(ndim)  
real(kind=dp), intent(out)     ::  dgamm_dx(ndim, ndim-2)  
 
! Local Variable
integer :: k, i , ind_vel
real(kind=dp)  :: dvdx1(ndim-2)
 
external :: dvind_dx  ! compute d vx_0 / d {y, vy}
  
  ind_vel = ind + ndim/2
  
  ! TODO: we need also compute the differential of vx_0  w.r.t. (x,y,vx,vy)_0 
  call dvind_dx(pv, ind, dvdx1)
  
  
  
  ! compute the differential of the map \gamma:  (y, vy)_0 --> (x,y,vx, vy)_0 
  dgamm_dx = 0.d0      ! initialize to zeros
  dgamm_dx(ind_vel, :) = dvdx1  ! int-th row 
 
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


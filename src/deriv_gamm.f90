!***********************************************************************
!     ****   deriv_gamma   ****
! Compute the differential of the projection map: R^n --> R^(n-2)
! \gamma( (x,y,vx,vy)_0 ) = (x, vx)_0     for ndim = 4 planar case 

! For spacial case, \gamma( (x,y,z,  vx,vy,vz)_0 ) = (x,z, vx,vz)_0 

!       Input Variables 
!  pv        the current state 
!  ndim      dimension of the phase space 
!  ind       index of the compont for Poincare section 


!       Output Variables 
!  dgamm_dx  the differential of  \gamma             

!  Routine Used:
!   dvind_dx  the partial derivative of x(ind) w.r.t. the initial state

!  Finally Revised by Yu -- 20170519
!***********************************************************************
subroutine deriv_gamm( pv, ndim, ind, dgamm_dx, dvind_dx)

use dp_mod
implicit none
 
! Input  and Output Declaration  
integer, intent(in)          ::  ndim, ind 
real(kind=dp), intent(in)    ::  pv(ndim)  
real(kind=dp), intent(out)   ::  dgamm_dx(ndim, ndim-2)  
 
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
  
  ! check dgamm_dx  !  -- ckd
!  print*,'dgamm_dx'
!  do i = 1, ndim, 1
!    print*, dgamm_dx(i, :)  
!  end do
!  print*; !read*  
  
  return  
end subroutine deriv_gamm


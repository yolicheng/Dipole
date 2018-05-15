!***********************************************************************    
!  this subroutine is to integrate n equilibrium points simultaneously, to make sure the sampling points are all isochronous with each other. 
!  time interval:  [0, tf],   always save the energy or jacobi constant in the eight-th column to check the numerical error.

!  
! Note: 
!    1. the subroutine to compute the vector field shoube be deriv_n (instead of deriv for only one point)
!    2. for backward in time integration, the stepsize h should be negative, and because h=1.d-3 is assigned inside the subroutine
!       we introduce  tdir to control the integration direction    :  1:unstable; -1:stable  
!       h = tdir * h   
   
!  better to do the write loop outside the subroutine--- 
!  in case to change the name of the data file, and write during the do loop

!	Input variables  
!  y(*)	        inital state of trajectory, must be 1-dimensional row array. instead of multi-dimensional (np, 6) matrix
!  t0		start time for the integration
!  tf		end time for the integration 
!  n		the dimension of state   
!  fob		file tag of orbit (every 4 rows: t, p1, p2, p3, p4)
!  fdrv		file tag of relative position and velocity (every 6 rows: t, dr(3), dv(3))
!  fh 		file tag of angular momentum(t, h(3), <r14, h> ) 
!  tdir		control the integration direction :  1: forward;  -1: backward
!  deriv_n	external subroutine to compute the vector field for  n points at the same time
!  gr_cj	external subroutine to compute the conservative quantity: energy or Jacobi Constant

  
!  ROUTINE USED:  GR_RK78  gr_lf_n 

!  revised by Yu  -- 20160428
! --------------------------------------------------------------------------
subroutine plob_eq_n(y0, np, t0,tf, tdir, fob, fdrv, fh, deriv_n, gr_cj) 

use dp_mod
use pi_mod
implicit none 


integer, intent(in)       :: np, tdir, fob, fdrv, fh  
real(kind=dp), intent(in) :: y0(np*6), t0, tf

! Local Variables
real(kind=dp) ::   hmin, hmax, e, h, t, y(np*6), r(13, np*6), b(np*6), f(np*6), & ! rk78 
                   yi(np, 6), rm, vm, r_rel(np*(np-1)/2, 3), v_rel(np*(np-1)/2, 3),  &   !relative position and velocity
                   r12(3), r13(3), r14(3), hpl(3), rh14  ! check plane 

integer :: i, j, irel, debug  


external :: deriv_n, gr_cj ! vector field&& energy or Jacobi constant
real(kind=dp), external ::  dnrm2



debug = 0

! specify the error control 
! But after test, we conclude that we have to use fixed stepsize, otherwise, we are not comparing at the same epoch
h =  tdir*1.d-3 

! adaptive step control 
hmin = 1.d-11
hmax = 1.d0
e    = 1.d-13

! initial time and state
t =  t0  
y =  y0
  

! if the evaluation condition is (t+h .lt. tf), only consider the positive step size, integate forward
! ( dabs(t+h-t0) .lt. dabs(tf-t0) ), both positive and negative cases are considered.


do while( dabs(t-t0) .lt. dabs(tf-t0) ) 
  ! assignment for each point, 6-dimension(position+velocity)
  do i = 1, np 
    yi(i, :) = y( 6*(i-1)+1 : 6*i )
    
    ! the state of each point
    write(fob, '(7e24.14)')  t, yi(i, :)
!    if(debug ==1)  
     write(*, '(7e24.14)')  t, yi(i, :)
  enddo 
  
  if (debug == 1) print*, 't,h', t, h !; read*
  
  ! Save the relative position and velocity in the order: r12, r13, r14, r23, r24, r34 
  irel = 0 ! dimension of r_rel = np * (np-1) / 2 
  do i = 1, np-1  ! j > i 
    do j = i+1, np
      irel = irel + 1 
      
      if(debug == 1) print*, 'i,j,irel', i,j,irel !ck 
      
      r_rel(irel, :) = yi(j, 1:3) - yi(i, 1:3)  ! j-th point - i-th point 
      v_rel(irel, :) = yi(j, 4:6) - yi(i, 4:6)  ! j-th point - i-th point  
    enddo 
  enddo   
  
  if(debug == 1)  then 
    print*, 'irel=', irel ; read*
  endif 
  
  ! the relative position w.r.t the first point 
  if( np .ge. 2) then 
    r12 = r_rel(1, :) 
  endif 
  
  if(np == 4) then    
    r13 = r_rel(2, :)
    r14 = r_rel(3, :)
  
  !to unity normal vector of the first 3 points  
    call cross_product(r12, r13, hpl)
    hpl = hpl / dnrm2( 3, hpl, 1) ! normalize
  
  ! check if the 4th point is on this plane?  r14h= dcos<r14,h> - pi/2 
  ! the dot product with the angular momentum should be 0, the angular with h should be 90 degree 
    rh14 = dot_product(r14, hpl) / dnrm2(3, r14, 1)
    rh14 = dacos(rh14)-pi/2.d0 ! in radian
!  ******* save the useful information in data file to plot later **********    
    write(fh, '(5e24.14)') t, hpl, rh14  ! angular momentum 
  endif 
  
    
!  if( debug == 1 ) then 
!    print*; print*, 'angular momentum: r12 X r13'
!    print*, 'r12', r12
!    print*, 'r13', r13
!    print*, 'r14', r14
!    print*, 'h    ', hpl
!  
!    print*;  print*, 'check the 4th point on the plane defined by the first three points: r41 X h'
!    print*, '=0?  ', rh14;  print*
!    read*
!  endif 

  
  if(debug ==1) print*, 'Relative state check!'
  
  ! the relative position and velocity
  do i = 1, irel 
  
!  no need to save addition data as vm and rm, we can do this data process in gnuplot
!    rm = dnrm2 (3, r_rel(i, :), 1) 
!    vm = dnrm2 (3, v_rel(i, :), 1) 

    write(fdrv, '(7e24.14)') t, r_rel(i, :), v_rel(i, :) 
    
    if(debug == 1) write(*, '(7e24.14)') t, r_rel(i, :), v_rel(i, :)
  enddo 
  
  if(debug == 1) print*  ! for sepearation on screen
  
  ! write to file and then do the integration, in this case, the first row is the equilibrium point
!  print*, 'before gr_rk78'
  call gr_rk78(t, y, np*6, h, hmin,hmax,e, r,b,f, deriv_n)
!  print*, 'after gr_rk78'; read*
  
enddo
  
  
end subroutine plob_eq_n


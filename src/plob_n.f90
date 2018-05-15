!***********************************************************************    
!  this subroutine is to integrate n points simultaneously,  for two purposes:
!  1 -- For station keeping:  to make sure the sampling points are all isochronous with each other.
!  2 -- Compute the Time-T map for a curve (i.e. a periodic orbit)

!  time interval:  [0, tf],   always save the energy or jacobi constant in the n+2-th column to check the numerical error.

!  ** NOTE ** 
! ----- tdir is real, not integer  ---------


!    1. for this routine, we do not care about the variational matrix, so in principle, always save the orbit 
!    2. the subroutine to compute the vector field shoube be deriv_n (instead of deriv for only one point)
!    3. for backward in time integration, the stepsize h should be negative, and because h=1.d-3 is assigned inside the subroutine
!       we introduce  tdir to control the integration direction    :  1:unstable; -1:stable  
!       h = tdir * h   
   
!  better to do the write loop outside the subroutine--- 
!  in case to change the name of the data file, and write during the do loop

!	Input variables  
!  y(*)           inital state of trajectory, must be 1-dimensional row array. 
!                 transfrom the input multi-dimensional (np, 6) matrix into a row array, handle with the reshape inside this routine is more convenient
!  ndim           the dimension of state 
!  np             number of points
!  t0             start time for the integration
!  tf             end time for the integration 
!  tdir		control the integration direction :  1: forward;  -1: backward

!  fob            file tag of orbit (every 4 rows: t, p1, p2, p3, p4)
!  deriv_n        external subroutine to compute the vector field for  n points at the same time
!  gr_cj          external subroutine to compute the conservative quantity: energy or Jacobi Constant

  
!  ROUTINE USED:  GR_RK78  gr_lf_n 

!  revised by Yu  -- 20161018
! --------------------------------------------------------------------------
subroutine plob_n(y0, ndim, np, t0,tf, tdir, ispl, fob, yf, deriv, gr_cj) 

use dp_mod
use pi_mod
implicit none 

integer, intent(in)           :: ndim, np, ispl, fob 
real(kind=dp), intent(in)     :: y0(np, ndim), t0, tf, tdir
real(kind=dp), intent(out)    :: yf(np, ndim)

! Local Variables
real(kind=dp) ::   hmin, hmax, e, h, t, y(np*ndim), r(13, np*ndim), b(np*ndim), f(np*ndim), & ! rk78 
                   yi(ndim), cj, t0copy, tfcopy 

integer :: i, debug  

character(len=70) ::  fmt1

external :: deriv, gr_cj      ! vector field   && energy or Jacobi constant
real(kind=dp), external ::  dnrm2

  write(fmt1,  fmt='(a,i0,a)')  '(', ndim+2, 'e24.14)'
  
  debug = 1

  ! specify the stepsize and error control 
  ! But after test, we conclude that we have to use fixed stepsize, otherwise, we are not comparing at the same epoch
  h =  tdir * dmin1(1.d-3,  tf / 1.d2)
  t0copy = tdir * t0 
  tfcopy = tdir * tf 
  
  ! Adaptive step control 
  hmin = dmin1(1.d-6, tf / 1.d2)
  hmax = 1.d-1
  e    = 1.d-14

  ! initial time and state
  t =  t0copy 
  y =  reshape( transpose(y0),  (/np*ndim/) ) ! column-wise, work on np*ndim array
  
  if(debug == 1) then 
    print*, 'Check the state of the initial 5  points: '
    do i = 1, 5, 1
      print*, y0(i, :)
    end do
    print*; read*
  endif 
  
  ! if the evaluation condition is (t+h .lt. tf), only consider the positive step size, integate forward
  ! ( dabs(t+h-t0) .lt. dabs(tf-t0) ), both positive and negative cases are considered.

  do while( dabs(t-t0copy) .lt. dabs(tfcopy-t0copy) ) 
    ! assignment for each point, ndim-dimension(position+velocity)
 
    do i = 1, np 
      yi   = y( ndim*(i-1)+1 : ndim*i )
      call gr_cj(yi, cj) 
      
      ! the state and energy  at each point
      if(ispl == 1)  write(fob, fmt = fmt1)  t, yi, cj 
      
      if(debug ==1 .and. i < 5)  then 
!        if(i == 5) read* !the reshapce ckd!
        write(*,   fmt = fmt1)  t, yi, cj 
      endif     
    enddo 
  
    if (debug == 1) then 
      print*; print*, 't,h', t, h ; print*; read*; !ck
    endif 
   
    if(ispl == 1) write(fob, *)  ! add a blank line between epoch 
    
    ! Save with the initial epoch
    call gr_rk78(t, y, np*ndim, h, hmin, hmax,e, r,b,f, deriv)
  
  enddo

  ! --- Only 1 more step to reach the final time tf -----
  if(dabs(t-tfcopy) .gt. 1.d-14) then
    
    h =  tfcopy - t ! for the stable orbit, the final time should be -tf
    
    ! to make sure the next step can be executed within allowable interval [hmin hmax]
    ! hmin should be dabs(h)
    call gr_rk78(t,y, np*ndim, h, dabs(h), hmax, e,r,b,f, deriv)
    
    
    do i = 1, np 
      yi = y( ndim*(i-1)+1 : ndim*i )
      call gr_cj(yi, cj) 
      
      ! the state and energy  at final point
      if(ispl == 1)  write(fob, fmt = fmt1)  t, yi, cj 
      
      if(debug ==1)  write(*,   fmt = fmt1)  t, yi, cj 
    enddo 
  
   
!    if(ispl == 1) write(fob, *)  ! add a blank line between epoch 
        
  endif 
   
  ! The final state in np X ndim array 
  yf = transpose( reshape(y,  (/ndim, np/)) )
  
end subroutine plob_n


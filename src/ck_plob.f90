!***********************************************************************    
!  this subroutine is to plot a segment of trajectory, with the intial state y0(6), 
!  and time span for integration is [0, tf], always save the energy or jacobi constant in the eight-th column to check the numerical error.
!  And save all the date in file ftag, with 2 blank lines between 2 different orbits

! Note: for backward in time integration, the stepsize h should be negative, and because h=1.d-3 is assigned inside the subroutine
!       we introduce  tdir to control the integration direction    :  1:unstable; -1:stable  
!       h = tdir * h   
   
!  better to do the write loop outside the subroutine--- 
!  in case to change the name of the data file, and write during the do loop

!  Input variables  
!    y0(nvar) :  inital state of trajectory 
!    ndim     :  dimension of the phase space
!    nvar     :  dimension of state, ndim: position + velocity; ndim*(ndim+1):   6(4)+the variational matrix 
!    t0    :  start time for the integration
!    tf    :  end time for the integration 

!    tdir  :  control the integration direction :  1: forward;  -1: backward
!    ispl  :  flag to decide if we save the integration(1) or not(0).  
!             Sometimes we only need the integration without needing the data 
!    ftag  :  data file to save the result of integration, t -- state -- cj 
!    deriv :  external subroutine to compute the vector field
!   gr_cj  :  external subroutine to compute the conservative quantity: energy or Jacobi Constant

! Output Variables
!    y(nvar)  :  finial state of trajectory 
  
!  ROUTINE USED: GR_RK78,  DERIV,  GR_CJ

!  revised by Yu  -- 20160219
! --------------------------------------------------------------------------
subroutine ck_plob(y0, t0, tf, ndim, nvar, tdir, ispl, ftag, deriv, gr_cj,  y) 
use dp_mod
implicit none 
    
integer, intent(in)        :: ndim, nvar, ftag, tdir, ispl     
real(kind=dp), intent(in)  :: y0(nvar),t0, tf
real(kind=dp), intent(out) :: y(nvar) ! As the final state
  
! External declaration
external :: deriv, gr_cj ! vector field  && energy or Jacobi constant

! Local Variables
real(kind=dp) :: r(13,nvar), b(nvar), f(nvar), t, h, hmin, hmax, e,  & ! rk78 
                 cj, dh, yst(ndim)  ! energy 
    
                 
integer           ::  debug  ! to output the orbit on the screen to check 
character(len=70) :: fmt1

  debug = 0
  
  write(fmt1,  fmt='(a,i0,a)')  '(', ndim+2,'e24.14)'
!  print*, fmt1; read* 
  
  
  ! copy the initial condition and time 
  y = y0    
  t = t0    
  
  
  ! 
  ! maybe small value is better, for a p.o. with small period, the orbit could be quite coarse 
  h = dmin1(1.d-3,  tf/1.d2)
  
  h =  tdir*h

! specify the error control 
  hmin = dmin1(1.d-6,  tf/1.d2)
  hmax = 1.d-1
  
!  hmax = 1.d-1
  e    = 1.d-14

    if(debug == 1) then
      print*, 'before plob: h,hmin, hmax', h, hmin, hmax 
      print*, 't0, tf, y0', t0, tf, y0(1:ndim)
      read*
    endif 

! if the evaluation condition is (t+h .lt. tf), only consider the positive step size, integate forward
! ( dabs(t+h-t0) .lt. dabs(tf-t0) ), both positive and negative cases are considered.

  
  if(ispl == 1) then 
    ! better to save as a block to use index for plot, put the two lines at the begining of each block, makes more sense
    write(ftag,*); write(ftag,*)  
  endif 
  

! integrate for the required time interval  
  do while( dabs(t-t0) .lt. dabs(tf-t0) ) 
  
    if(ispl == 1) then 
      yst = y(1:ndim)
      call gr_cj(yst, cj)
      write(ftag,  fmt = fmt1)  t, yst, cj
    endif 
    
    if(debug == 1) then 
      yst = y(1:ndim);   call gr_cj(yst, cj)
      write(*, fmt = fmt1)  t, y(1:ndim), cj !ck
    endif  
    call gr_rk78(t,y, nvar, h,hmin,hmax,e, r,b,f, deriv )
    
  enddo
  
  ! --- for only 1 step to get the exact final time -----
  if(dabs(t-tf) .gt. 1.d-14) then
    
    dh =  tdir*tf - t ! for the stable orbit, the final time should be -tf

    ! to make sure the next step can be executed within allowable interval [hmin hmax]
    ! hmin should be dabs(h)
    call gr_rk78(t, y, nvar, dh, dabs(dh), hmax, e,r,b,f, deriv)
    
  endif 
  
  ! -- the last point --- 
  if(ispl == 1)  then 
    yst = y(1:ndim)
    call gr_cj(yst, cj) 
    write(ftag, fmt = fmt1) t, yst, cj 
    
   if(debug == 1)  write(*, fmt = fmt1)  t, y(1:ndim), cj 
  endif 
     
  if(debug == 1) then 
    yst = y(1:ndim);   call gr_cj(yst, cj)
    write(*, fmt = fmt1)  t, y(1:ndim), cj !ck
    print*, 'Finish plob, t0, t, tf, yf', t0, t, tf, y(1:ndim) 
  endif
  
  
  return 

end subroutine ck_plob

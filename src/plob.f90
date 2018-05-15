!***********************************************************************    
!  Integrate a trajectory  with the intial state y0(ndim) for time span [t0, tf]
!  t0 and tf are positive numbers, we take into consideration of the time sense
!  
!  and contol the integration time by the modulus of the elapsed time 
!  where we take 1:unstable (forwards); -1:stable(backwards, time starts with 0.d0 and goes to negative) 
!  The energy or Jacobi constant in the ndim+2-th column to check the numerical error.

!  Save all the date in file ftag, with 2 blank lines at the beginning to seperate 2 different orbits

! ** Note ** :

! ----- tdir is real, not integer  ---------


!  for backward integration, the stepsize h should be negative, we take h = tdir * h
   

!  Input variables  
!    y0(nvar) :  inital state of trajectory 
!    ndim     :  dimension of the phase space
!    nvar     :  dimension of state, ndim: position + velocity; ndim*(ndim+1):   6(4)+the variational matrix 
!    t0    :  start time for the integration
!    tf    :  end time for the integration 

!    tdir  :  control the integration direction :  1: forward for stable manifold;  -1: backward for unstale manifold 
!    ispl  :  flag to decide if we save the integration(1) or not(0).  
!             Sometimes we only need the integration without needing the data 
!    ftag  :  data file to save the result of integration, t -- state -- cj 
!    deriv :  external subroutine to compute the vector field
!   gr_cj  :  external subroutine to compute the conservative quantity: energy or Jacobi Constant

! Output Variables
!    y(nvar)  :  finial state of trajectory 
  
!  ROUTINE USED: GR_RK78,  DERIV,  GR_CJ

! ----------  by Yu  -- 20160219
! Version 2.0 by Yu  -- 20161021  --remove tmax?? 

! --------------------------------------------------------------------------
subroutine plob(y0, t0, tf, ndim, nvar, tdir, ispl, ftag, deriv, gr_cj,  y) 

use dp_mod
implicit none 
    
integer, intent(in)        :: ndim, nvar, ftag,  ispl     
real(kind=dp), intent(in)  :: y0(nvar),t0, tf, tdir 
real(kind=dp), intent(out) :: y(nvar) ! As the final state
  
! External declaration
external :: deriv, gr_cj ! vector field  && energy or Jacobi constant

! Local Variables
real(kind=dp) :: r(13,nvar), b(nvar), f(nvar), t, h, hmin, hmax, e,  & ! rk78 
                 cj, dh, yst(ndim), t0copy, tfcopy   ! energy 
    
                 
integer           ::  debug  ! to output the orbit on the screen to check 
character(len=70) :: fmt1

  debug = 0
  
  ! Take into account of the time sense
  t0copy = tdir * dabs(t0) 
  tfcopy = tdir * dabs(tf)
  
  ! save time + state + energy (ndim+2) numbers
  write(fmt1,  fmt='(a,i0,a)')  '(', ndim+2,'e24.14)' 
  
  ! Copy the initial condition and time 
  y = y0    
  t = t0copy   
 
  
  ! 2016-12-05 18:00:08 
  ! We have checked the final state with different values of h, hmin, hmax,
  ! for the same time interval 
  ! Finnaly suggested value: 
  !     -- hmin = 1.d-6 (not too small);  
  !     -- hmax = 1.d-1 (not too big);
  !     -- e    = 1.d-14
  !     -- h    = 1.d-3, but to make sure the integration proceed, we take the minimum 
  !                      between 1.d-3 and tf/1.d1   
  
  
  ! maybe small value is better, for a p.o. with small period, the orbit could be quite coarse 
  h = dmin1(1.d-3,  tf/1.d2)
  h = tdir*h

  if(debug == 1) then
    print*,'h, tdir,tf', h, tdir, tf
  endif 
  
    
  ! Stepsize and error control 
  
  ! for hmin, 1.d-10 is illegally too small, the accumulation of round error will give us unreliable results 
  ! for hmax, 1.d-1 is reasonable, and 1.d0 is a little bit too much 
  hmin = dmin1(1.d-6,  tf/1.d2)  
  hmax = 1.d-1
  
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
    write(ftag,*); 
    write(ftag,*)  
  endif 
  

! integrate for the required time interval  
  do while( dabs(t-t0copy) .lt. dabs(tfcopy - t0copy) ) 
  
    if(ispl == 1) then 
      yst = y(1:ndim);   call gr_cj(yst, cj)
      write(ftag,  fmt = fmt1)  t, yst, cj
    endif 
    
    if(debug == 1) then 
      yst = y(1:ndim);   call gr_cj(yst, cj)
      write(*, fmt = fmt1)  t, y(1:ndim), cj !ck
    endif  
    call gr_rk78(t,y, nvar, h,hmin,hmax,e, r,b,f, deriv )
    
    if(hmin > 1.d-5) then 
      print*, 'hmin, t, y,h ', hmin, t, y(1:ndim), h
      print*; read*; !ck 
    endif 
  enddo
  
  ! --- for only 1 step to get the exact final time -----
  if(dabs(t-tfcopy) .gt. 1.d-14) then
    
    dh =  tfcopy - t ! for the stable orbit, the final time should be -tf
    if(debug == 1) then 
      print*, 'dh, h ', dh, h; print*; read*; !ck
    endif  
     
    ! to make sure the next step can be executed within allowable interval [hmin hmax]
    ! hmin should be dabs(h)
    call gr_rk78(t, y, nvar, dh, dmin1(dabs(dh), hmin), hmax, e,r,b,f, deriv)
    
  endif 
  
  ! -- the last point --- 
  if(ispl == 1)  then 
    yst = y(1:ndim);   call gr_cj(yst, cj) 
    write(ftag, fmt = fmt1)  t, yst, cj 
    
   if(debug == 1)  write(*, fmt = fmt1)  t, y(1:ndim), cj 
  endif 
     
  if(debug == 1) then 
    yst = y(1:ndim);   call gr_cj(yst, cj)
    write(*, fmt = fmt1)   t, y(1:ndim), cj !ck
    print*, 'Finish plob, t0, t, tf, yf', t0, t, tf, y(1:ndim) 
  endif
  
  return 

end subroutine plob

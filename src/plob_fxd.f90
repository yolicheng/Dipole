subroutine plob_fxd(y0, n, t0, h0, np, xmax, ftag, yall,  isescp, deriv, gr_cj) 
!********************************************************************
!  this subroutine is to plot a segment of trajectory, the intial state y0(n), 
!  time span for integration is (0, tf]  to keep in consistent with Furier Analysis

! TODO: be careful with the time span, (0, tf] or [0, tf)?? 

! And save all the data in file ftag, with 2 blank lines between 2 different orbits

! The idea of fixed step integration: 
!   --1. Assume h0 is the fixed stepsize.  Use hmin = h0/10, and hmax = h0
!   --2. we integrate until the t > i * h0, i is integer, save all the data needed for integration at epoch t :  t, h, y(n), ... 
!   --3. integrate back with a very small step, dh = i*h0 - t, remember to modify temprorily the value of  hmin... 

!  --discard : if we forcely fix hmin=hmax=h0, it will violates the error control during integration  
! 

! 20160128  
!  add the escape constraint, stop integration once we are outside this domain. 
   
!  better to do the write loop outside the subroutine--- 
!  in case to change the name of the data file, and write during the do loop

!  Insmut variables  
!    y0(n) :  inital state of trajectory 
!    n     :  number of dimension 
!    t0    :  start time for the integration
!    h     :  the fixed step size
!    np    :  number of sampling points  
!   xmax   :  the maximum modulus of any component, to detect escape
!             **NOTE** have tried 1 and 2, in most cases, 1 is enough 
!    ftag  :  file tag to save the orbit 
 
!   deriv  :  the vector field, use as input and declare external, 
!             so we can deal with different problems then 

! Output Variables
!  yall(np):  the state of all the sampling points
!  isescp  :  1: the orbit escape, one of the position component exceeds the bound specified by xmax
  
!  ROUTINE USED: GR_RK78     DERIV(A,B,N,F)
!********************************************************* ***********

implicit none 
integer, parameter:: dp = kind(1.d0) 

integer, intent(in)           :: n, np, ftag 
integer, intent(out)          :: isescp  
real(kind=dp), intent(in)     :: y0(n), t0, h0, xmax
real(kind=dp), intent(out)    :: yall(np, n) ! As the final state
external deriv

! 	Local variables
integer :: ip, debug   
real(kind=dp) ::   r(13,n), b(n), f(n), h, hmin, hmax, e1, t,  y(n), cj, &  !rk78
                   tc, yc(n), dh, ti

debug = 0

! by default, the orbit will not escape
isescp = 0 
! Initial value of stepsize, has to be small enough 
h = dmin1(h0 /1.d1, 1.d-3) 
 
! use the same value to make it as fixed step size -- discard
hmin = 1.d-6
hmax = 1.d-1

!hmax = h0
  
e1   = 1.d-13
 
if(debug == 1) then
  print*, 'before plob_fxd: h,hmin, hmax', h, hmin, hmax 
  read*
endif   

! initial point --- not saved, start from 1*h0 
ip = 0   !  index of sampling point
t  = t0  !  0.d0
y  = y0  !  initial state    


! discard since we want to have time span as (0, tf] --- 2016-06-14 22:13:56 
!yall(1, :) = y0
!call gr_cj(y0(1:n), cj)
!write(ftag,'(8e24.14)')  t0, y0,  cj

! Integrate the periodic orbits by gr_rk78, write the x, y in txt file with ftag

do while( ip .lt. np ) 

  ip = ip + 1  ! index of samples 
  ti = t0 + ip * h0  ! the epoch of ip-th sample 

! ---- integrate unitl the next epoch, ti = t0 +  ip * h0 -------------

! -- 1:  we reach to the first epoch that t > t0 + ip * h0
  do while( t < ti ) 
    call gr_rk78(t,y,n,  h,hmin,hmax,e1, r,b,f,  deriv)
    
!    if (debug == 1)   write(*,'(8e24.14)')  t, h, y 
  enddo 
  
  ! always use t, y, and h as the main line integration  
  ! so we use tc, yc and hc for the backward integration to obtain data at each epoch ti
  tc = t 
  yc = y 
  
! -- 2:  integrate backwards with a very small step size dh = ti - tc  
  dh = ti - tc 
  
  if(debug == 1)  print*, 'check dh=', dh, '   at ', ip, '-th sampling point';! print*;! read*
  
  call gr_rk78(tc, yc, n, dh, dabs(dh), hmax,e1,r,b,f,  deriv)
  
  if(debug == 1)  print*, 'After a new step backwards, t=', tc, '== ti?', ti ;! print*;! read*
!-------------------------- 
  
  if( maxval( dabs(yc(1:3)) ) > xmax  ) then ! once escape, return
    write(ftag, *) ;      print*, 'Escape!';     
    print*, tc, yc(1:n);  read*
    isescp = 1
    return
  endif  
  
  ! save ip-th sampling point
  yall(ip, :) = yc
  call gr_cj( yc(1:n), cj )
  
  write(ftag, '(8e24.14)')  tc, yc, cj ! write to the file 
  
  
  if(debug == 1)  then 
    print* ;  print*, ip, '-th sampling point' 
    write(*, '(8e24.14)')  tc, yc, cj  ! print to screen
    read*
  endif   
  
enddo

  if(debug == 1) then 
    write(*,*) 'after plob_fxd, t0, tf, yf'
    print*,  t0, tc,  yc; print*; read*
  endif 

write(ftag,*); write(ftag,*)   ! add two blank line to seperate orbit

 close(ftag)
return 
end subroutine plob_fxd



subroutine eqprop_num( p0, np, tf, xmax, ramp, vamp, fds, fdsmax, deriv_n, iscap, dsmax, ymax)

! here we try the first approach ---- purely numerical 

! 1.  numerical : 
!   integrate samples from a ball around the equilibrium point,  for a given time interval [0:100] days 
!   see how this ball evolves, to which value of the radius will the initial points from this ball escape during this time interval .... 

!   for the position and velocity, we use 2-sigam Gaussian random numbers with deviation in position and velocity as ramp and vamp, respectively

! during the integration, at each epoch we only keep the maximum relative distance from the reference point P0,  
 
!  2. analytical:  --- discard this approach 
!    use the linear solution of the variational equation as the deviation from the equilibrium point, 
!    and we have three modes related to the three pairs of pure imaginary eigenvalues of the Jacobi Matrix of the differential of the vector field.
 
!    and for each mode, we study the stability by incease the value of the magnitude of this mode, and then add the state of the equilibrium point, 
!    use this as the initial condition and integrate with the nonlinear equations of motion also for 100 days, to see up to which value will each mode remain in the neibougher 


! . 2-sigma Gaussian random deviations are added as the injection error on the initial state. 

!      SUBROUTINE DLARNV( IDIST, ISEED, N, X ) from LAPACK 
!      -- provided 4-dimension initial seed iseed for dlarnv(), with iseed(4) odd
!      -- idist specifies the distribution of the random numbers:
!           = 1:  uniform (0,1)   
!           = 2:  uniform (-1,1)
!           = 3:  normal  (0,1)   -- we take this one


! 3. the cost to move the orbit back to  the initial point (the cost of station keeping), when a control is necessary to be applied
!    check the magnitude of the position deviation... 



!       Input Variables 
!  p0           the initial point 
!  np           number of samples 
!  tf           time interval for integration is [0 : tf]
!  ramp, vamp   the 2-sigma deviation of the position and velocity: 
!               ramp can be seen as the radius of the sampling ball around the initial point p0 
!  fds       the maximum distance at all the epoches


!  Routine Used:
!     dlarnv, dnrm2, gr_rk78, gr_lf_n

! Finally Revised by Yu -- 20160525
!----------------------------------------

use lf_mod, only : runit, vunit 

implicit none
integer, parameter  ::  dp = kind(1.d0) 

 
!  Input  and Output Declaration   
integer, intent(in)     ::  np,   fds, fdsmax
integer, intent(out)     ::  iscap 
real(kind=dp), intent(in)       ::  p0(6), tf, xmax,  ramp, vamp   
real(kind=dp), intent(out)      ::  dsmax,  ymax(6) 
 
! Local Variables
real(kind=dp) ::  devn(np*6), devi(np,6), mean, sd, &  ! deviation 
                  y(np*6), t, h, hmin, hmax, e, r(13, np*6), b(np*6), f(np*6), &  ! rk78
                  yi(np, 6), dsmax_t, ymax_t(6), ds, dr(3)   ! dsmax 
                   
integer :: i, j, iseed(4), idist, iserr    

!	subroutine from packages declared external here
real(kind=dp), external :: dnrm2 ! from BLAS-LEVEL1
external :: dlarnv ! from LAPACK to generate normal distribution random data 
external :: deriv_n

! try to save the generated random data devi, but failed, because array with automatic size cannot have save attributes 
! but we can save the updated iseed to produce new random data 
data iserr /0/
save iserr, iseed    

! only use the initial seed iseed at the first call 
if (iserr  ==  0)  then 

! use Lapack Driver routine DLARNV to generate normal distribution random data
!  SUBROUTINE DLARNV( IDIST, ISEED, N, X )
! different seed generate different random data, but the same seed always produce the same one...

  !  iseed = (/29, 18, 3, 11/)
  iseed = (/0, 0, 0, 1/) !  use this one, or should we use another one?
   
!  iserr = 1 ! use the same random data, the only difference is the radius of the ball around P0
endif 

 
! ----------- random data in normal distribution --------------- 
idist = 3 ! normal distribution 
call dlarnv( idist, iseed, np*6,  devn)
  
! Check mean and standard deviation
mean = SUM(devn)/ np/ 6
sd = DSQRT(SUM((devn - mean)**2)/ np / 6)

print*, 'check the normal distribution deviation'  
WRITE(*, "(A,F18.10)") "Mean = ", mean
WRITE(*, "(A,F18.10)") "Standard Deviation = ", sd
  
! write the deviation in real unit into file as reference
open(88,   file='./dat/dev.dat', access ='append',status='replace') 
write(88,*) '#The deviation in real unit: (x,y,z) km	(vx,vy,vx) km/s'

do i =  1, np 
  write(88,  '(6f20.10)')  ramp * devi(i, 1:3)*runit, vamp * devi(i, 4:6)*vunit 
enddo

 close(88) 
 

 ! mark that we already have this np*6 initial random data   
devi = reshape( devn, (/np, 6/))

do i = 1, np 

  ! add the deviation  -- for position 
  yi(i, 1:3)  = p0(1:3)  + ramp * devi(i, 1:3) 

  ! -- for velocity
  yi(i, 4:6)  = p0(4:6)  + vamp * devi(i, 4:6)  

enddo 


! ---   initialize  the initial epoch ----------------  
t = 0.d0 
dsmax = 0.d0 ! for all ds 
dsmax_t  = 0.d0 ! for each epoch 

iscap = 0

! construct one-dimensional array of the initial condition for integration
do i = 1, np 
  y((i-1)*6+1 : i*6) = yi(i, :) 

! compute the distance from the fixed point, and record only the maximum
  dr =  yi(i,1:3)-p0(1:3)
  ds = dnrm2( 2, dr, 1 ) 
  
  if( ds > dsmax_t) then ! at the initial epoch, dsmax =  dsmax_t = 0
    dsmax_t = ds 
    ymax_t  = yi(i,:)
    
    dsmax = ds
    ymax  = ymax_t
  endif
    
  ! save the state of each point -- generally, we do not save the orbit data, it is not necessary. 
  
enddo     

write(fds, '(2e24.14)')  t, dsmax_t


! ---------------------------------------------


! -------------- start integration -------------------------
! specify the error control for rk78
! But after test, we conclude that we have to use fixed stepsize, otherwise, we are not comparing at the same epoch
h  =  1.d-3 

! adaptive step control 
hmin = 1.d-10
hmax = 1.d0
e    = 1.d-14


! question 1:   save all the points at the initial epoch? do we need to save all the points??? 
! or it is only necessary to save  the maximum distance from the fixed point P0 


do while( t .lt. tf ) 

! for every epoch, reset dsmax 
  dsmax_t = 0.d0 
  
  call gr_rk78(t, y,  np*6, h, hmin,hmax,e,  r,b,f, deriv_n)


! assignment for each point, 6-dimension(position+velocity)
  do i = 1, np 
  
    yi(i, :) = y( 6*(i-1)+1 : 6*i )
    
    ! compute the distance from the fixed point, and record only the maximum 
    dr =  yi(i,1:3)-p0(1:3)
    ds = dnrm2( 3, dr, 1 ) 
  
    if( ds > dsmax_t) then 
      dsmax_t = ds
      ymax_t = yi(i, :)
    endif
     
  enddo 
  
 ! for each epoch t, check if the orbit escapes     
  if( dsmax_t > xmax ) then 
    iscap = 1
    exit 
  endif 

  ! save only the maximum distance for each epoch
  write(fds, '(2e24.14)') t, dsmax_t


  ! obtain the maximal ds for all dsmax_t     
  if( dsmax_t > dsmax ) then 
    dsmax = dsmax_t
    ymax  = ymax_t
  endif
  
enddo 

! save the deviation a, maximum ds and the corresponding I.C. 
write(fdsmax, '(9e24.14, i5)') ramp, vamp, dsmax, ymax, iscap 

print*, 'ramp, vamp, dsmax, ymax(6)'
write(*, '(9e24.14)') ramp, vamp, dsmax,  ymax
 
write(fds, *); write(fds, *) ! two blank lines to seperate different values of ramp-vamp


return  
end subroutine eqprop_num
 





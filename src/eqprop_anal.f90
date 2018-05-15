subroutine eqprop_anal( p0, tf, xmax, ramp,  isob, fob, fdsmax, deriv_n, iscap)

! here we try the second approach ---- analytical + numerical  

! 1.  numerical : 
!   integrate samples from a ball around the equilibrium point,  for a given time interval [0:100] days 
!   see how this ball evolves, to which value of the radius will the initial points from this ball escape during this time interval .... 

!   for the position and velocity, we use 2-sigam Gaussian random numbers with deviation in position and velocity as ramp and vamp, respectively

! during the integration, at each epoch we only keep the maximum relative distance from the reference point P0,  
 
!  2. analytical: 
!    use the linear solution of the variational equation as the deviation from the equilibrium point, 
!    and we have three modes related to the three pairs of pure imaginary eigenvalues of the Jacobi Matrix of the differential of the vector field.
 
!    and for each mode, we study the stability by incease the value of the magnitude of this mode, and then add the state of the equilibrium point, 
!    use this as the initial condition and integrate with the nonlinear equations of motion also for 100 days, to see up to which value will each mode remain in the neibougher 

!       Input Variables 
!  p0           the initial point 
!  np           number of samples 
!  tf           time interval for integration is [0 : tf]

!  isob         save the orbit data or not 
!  fob          file tag of the propagation of all the points 
!  fdsmax       the maximum distance at all the epoches


!  Routine Used:
!     dlarnv, dnrm2, gr_rk78, gr_lf_n

! Finally Revised by Yu -- 20160526
!----------------------------------------

use lf_mod, only : runit, vunit 

implicit none
integer, parameter  ::  dp = kind(1.d0) 

 
!  Input  and Output Declaration   
integer, intent(in)     ::  np,  isob, fob, fdsmax
integer, intent(out)     ::  iscap 
real(kind=dp), intent(in)         ::  p0(6), tf, xmax,  ramp, vamp   

 
! Local Variables
real(kind=dp) ::  devn(np*6), devi(np,6), mean, sd, &  ! deviation 
                  y(np*6), t, h, hmin, hmax, e, r(13, np*6), b(np*6), f(np*6), &  ! rk78
                  yi(np, 6), dsmax, ds, dr(3) ! eq_num
                   
                  
integer :: i, j, iseed(4), idist   

!	subroutine from packages declared external here
real(kind=dp), external :: dnrm2 ! from BLAS-LEVEL1
external :: dlarnv ! from LAPACK to generate normal distribution random data 
external :: deriv_n

  

! use Lapack Driver routine DLARNV to generate normal distribution random data
!  SUBROUTINE DLARNV( IDIST, ISEED, N, X )
! different seed generate different random data, but the same seed always produce the same one...

iscap = 0

! The initial deviation by linearized equation, three modes, each related to one pair of pure imaginary eigenvalues
ksi1 =  

do i = 1, np 

! print*, devi(i, :)
  
  ! add the deviation  -- for position 
  yi(i, 1:3)  = p0(1:3)  + ramp * devi(i, 1:3) 

  ! -- for velocity
  yi(i, 4:6)  = p0(4:6)  + vamp * devi(i, 4:6)  

enddo 

!print* ;   read*
  
! Check mean and standard deviation
mean = SUM(devn)/ np/ 6
sd = DSQRT(SUM((devn - mean)**2)/ np / 6)
  
WRITE(*, "(A,F18.10)") "Mean = ", mean
WRITE(*, "(A,F18.10)") "Standard Deviation = ", sd
  

! write the deviation in real unit into file as reference
open(88,file='./dat/dev.dat', access ='append',status='replace') 
write(88,*) '#The deviation in real unit: (x,y,z) km	(vx,vy,vx) km/s'

do i =  1, np 
  write(88,  '(6f20.10)')  ramp * devi(i, 1:3)*runit, vamp * devi(i, 4:6)*vunit 
enddo

 close(88)   

! initialize 
dsmax = 0.d0

! ---- the initial epoch 
t = 0.d0 

! construct one-dimensional array for integration
do i = 1, np 
  y((i-1)*6+1 : i*6) = yi(i, :) 

! compute the distance from the fixed point, and record only the maximum
  dr =  yi(i,1:3)-p0(1:3)
  ds = dnrm2( 2, dr, 1 ) 
  
!  print*, 'ds=', ds 
!  read* 
  
  if( ds > dsmax) dsmax = ds 
  
  ! save the state of each point
  if(isob == 1) write(fob, '(7e24.14)')  t, yi(i, :)
  
enddo     

write(fdsmax, '(2e24.14)') t, dsmax
if(isob == 1) write(fob, *)  ! a blank line between each epoch  
! ---------------------------------------------


! specify the error control for rk78
! But after test, we conclude that we have to use fixed stepsize, otherwise, we are not comparing at the same epoch
h  =  1.d-4 

! adaptive step control 
hmin = 1.d-10
hmax = 1.d0
e    = 1.d-13


! question 1:   save all the points at the initial epoch? do we need to save all the points??? 
! or it is only necessary to save  the maximum distance from the fixed point P0 


do while( t .lt. tf ) 
! for every epoch, reset dsmax 
  dsmax = 0.d0 
  
  call gr_rk78(t, y,  np*6, h, hmin,hmax,e,  r,b,f, deriv_n)
! 
!  if( dabs( y(1) ) > xmax  .or. dabs( y(2) ) > xmax  .or. dabs( y(3) ) > xmax ) then 
!    
!  endif 
  
  ! assignment for each point, 6-dimension(position+velocity)
  do i = 1, np 
  
    yi(i, :) = y( 6*(i-1)+1 : 6*i )
    
    ! compute the distance from the fixed point, and record only the maximum 
    dr =  yi(i,1:3)-p0(1:3)
    ds = dnrm2( 3, dr, 1 ) 
  
    if( ds > dsmax) dsmax = ds
    
  ! save the state of each point
    if(isob == 1) write(fob, '(7e24.14)')  t, yi(i, :)
  
  enddo 
  
  ! for each epoch t 
  write(fdsmax, '(2e24.14)') t, dsmax
!  write(*, '(2e24.14)') t, dsmax
  if(isob == 1) write(fob, *)  ! a blank line between each epoch   
  
  if( dsmax > xmax ) then 
    iscap = 1
    exit 
  endif 
    
enddo 

write(fdsmax, *)
write(fdsmax, *)

return  
end subroutine eqprop_anal
 





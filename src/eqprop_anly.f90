subroutine eqprop_anly( p0, nfam, tf, xmax, ramp,  isob, fob, fds, deriv_n, iscap)

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
!  nfam         number of families 
!  tf           time interval for integration is [0 : tf]

!  isob         save the orbit data or not 
!  fob          file tag of the propagation of all the points 
!  fds       the maximum distance at all the epoches


!  Routine Used:
!    dnrm2, gr_rk78, gr_lf_n

! Finally Revised by Yu -- 20160526
!----------------------------------------
use lf_mod 

implicit none

!integer, parameter  ::  dp = kind(1.d0) 

 
!  Infamut  and Output Declaration   
integer, intent(in)      ::  nfam,  isob, fob, fds
integer, intent(out)     ::  iscap 
real(kind=dp), intent(in)         ::  p0(6), tf, xmax,  ramp  

 
! Local Variables
real(kind=dp) ::  dlf(6,6), wr(6), wi(6), vr(6,6),  vrchs(6), &    ! anly :  eigenspace of variational matrix   
                  y(nfam*6), t, h, hmin, hmax, e, r(13, nfam*6), b(nfam*6), f(nfam*6), &  ! rk78
                  yi(nfam, 6), ds(nfam),  dr(3) ! eq_num
                   
                  
integer :: i, j    

!	subroutine from packages declared external here
real(kind=dp), external :: dnrm2 ! from BLAS-LEVEL1
external :: deriv_n

  

iscap = 0

! The initial deviation by linearized equation, three modes, each related to one pair of pure imaginary eigenvalues

! Jacobi matrix of the lorentz force with respective to the state 
!call dflrtz(p0, dlf)
!do i = 1, n
! write(*,'(6f8.4)') dlf(i,:) 
!enddo

! compute the eigenvalues and eigenvectors  of dlf
!call eigrg(dlf, 6 ,1, wr,wi, vr)
!print*; read*
 
!do i = 1, nfam 
!  
!  print*, 'check', 2*i,'-th column of vr to use', vr(:, 2*i)
!  vrchs = vr(:, i*2)
!  vrchs = vrchs/dnrm2(3,vrchs(1:3), 1) ! the deviation in position is ramp, so | vrchs(1:3) | = 1
!  
!  print*, dnrm2(6,vrchs, 1), vrchs
!  
!    ! add the deviation  
!  yi(i, :)  = p0  + ramp * vrchs 

!enddo 


! ---- the initial epoch 
t = 0.d0 

! construct one-dimensional array for integration
do i = 1, nfam 
  
  y((i-1)*6+1 : i*6) = yi(i, :) 

! compute the distance from the fixed point, and record only the maximum
  dr =  yi(i,1:3) - p0(1:3)
  ds(i) = dnrm2( 2, dr, 1 ) 
  
  ! save the state of each point
  if(isob == 1) write(fob, '(7e24.14)')  t, yi(i, :)
  
enddo     

write(fds,  *)  t, ds 

if(isob == 1) write(fob, *)  ! a blank line between each epoch  
! ---------------------------------------------


! specify the error control for rk78
! But after test, we conclude that we have to use fixed stepsize, otherwise, we are not comparing at the same epoch
h  =  1.d-4 

! adaptive step control 
hmin = 1.d-10
hmax = 1.d0
e    = 1.d-14


do while( t .lt. tf ) 
  
  call gr_rk78(t, y,  nfam*6, h, hmin,hmax,e,  r,b,f, deriv_n)
  
  ! assignment for each point, 6-dimension(position+velocity)
  do i = 1, nfam 
  
    yi(i, :) = y( 6*(i-1)+1 : 6*i )
    
    ! compute the distance from the fixed point, and record only the maximum 
    dr =  yi(i,1:3)-p0(1:3)
    ds(i) = dnrm2( 3, dr, 1 ) 
    
  ! save the state of each point
    if(isob == 1) write(fob, '(7e24.14)')  t, yi(i, :)
  
  enddo 
  
  ! for each epoch t 
  write(fds, *) t, ds 
  
  if(isob == 1) write(fob, *)  ! a blank line between each epoch   
  
  if( maxval( ds ) > xmax ) then 
    iscap = 1
    exit 
  endif 
    
enddo 

write(fds, *)
write(fds, *)

return  
end subroutine eqprop_anly
 





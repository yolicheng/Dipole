! ******************************************************
! For PRTBP, with a given Jacobi Constant and the three values of (x,y, vx(or vy)),  
! compute the differential of the missing velocity '+'vx  w.r.t. (y, vy)
! for the purpose to compute the Jacobi Matrix of the linear equation for Newton method

!   vx   = f(h0, x0, y, vy) = sqrt( 2*omega - h0 - vy**2 )
!   d vx / d x =   Omega_x / vx
!   d vx / d y =  - vy / vx 
 
! To make it general, introduce ind_vel = ind + 2 as the index of ind-th velocity 
! and 5-ind is the index of the missing velocity  

!     Input Varaibles
!  x    the state vector, of dimension ndim (4 for PRTBP) (x,y,vx,vy)
! ind   the index of component fixed for Poincare section, x(ind) = p0 
!        the velocity in this direction is a function of vx = f(h0, x0, y, vy)

!     Output Varaibles
! dvind_dx   the differential of v(ind) w.r.t. to the free components in state x

! TODO: haven't checked for ind=2, if necessary, check every index carefully

! Finally revised by Yu  2016-08-30 18:17:14 
!******************************************************
subroutine dvind_dx_prtbp(x, ind, dvdx)
use dp_mod
use rtbpconst_mod, only : mu 
implicit none 

!  Input and Output
integer, intent(in)     ::  ind
real(kind=dp), intent(in)       :: x(4) 
real(kind=dp), intent(out)      :: dvdx(2) 

!  Local Varaibles
integer ::   ind_vel
real(kind=dp) ::  y1, y12, y22,  r1, r13, r15, r2, r23, r25, p1, p2, q, & ! velocity + acceleration
                  dwdx !, cj  

  ind_vel = ind + 2 

  ! acceleration
  ! Omega = 1/2(x**2 + y**2) + (1-mu)/r1 + mu / r2 +  MU*(1-MU) /2 
  ! r1**2 = (x-mu)**2 + y**2 + z**2
  ! r2**2 = (x-mu+1)**2 + y**2 + z**2

  ! Omega_x = x - 1/r1**3 * (1-mu) * (x-mu) - 1/r2**3 * mu * (x-mu+1) 
  ! Omega_y = y - 1/r1**3 * (1-mu) * y      - 1/r2**3 * mu * y
  ! Omega_z =    -1/r1**3 * (1-mu) * z      - 1/r2**3 * mu * z
     
  y1  = x(1) - mu
  y12 = y1*y1
  y22 = x(2)**2

  r1  = y12 + y22 ! the square of the distance to the big primary
  r13 = r1 * dsqrt(r1)
  r15 = r13 * r1

  r2  = (y1+1.d0)**2 + y22 
  r23 = r2 * dsqrt(r2)
  r25 = r23*r2

  p1 = (1.d0-mu) / r13
  p2 = mu/r23
  q  = -1.d0*(p1 + p2)

  ! to check if this routine is with no error, check the energy 
  !  call gr_cjprtbp(x, cj) ! ckd, cj = cj0  
  !  print*, 'cj = ', cj 
  !  read* 

  ! original one from gr_rtbp.f90 
  !   dwdx = x(1) - y1 * p1 - (y1 + 1.d0) * p2
  !   dwdy = x(2) * (1.d0 + q)

  ! vx = sqrt( 2*omega - h0 - X(3)*X(3) - X(4)*X(4) )
  ! the input vx is already computed by the above equation before the call this routine

  ! if ind = 1, x=p0 is fixed, we need compute d vx / d {y, vy}
  !   d vx / d y  =  Omega_y / vx
  !   d vx / d vy =   -vy    / vx

  ! similarly, if ind == 2, y=p0 is fixed, 
  ! we need compute d vy / d {x, vx}

  ! d vy / d x  =  Omega_x / vy
  ! d vy / d vx =   -vx    / vy

  ! vx and vy can be specified by ind_vel and 5-ind 
  !   if we take x=p0(ind = 1), the free velocity is vy, ind_vy = 4 = 5-ind 
  !   if we take y=p0(ind = 2), the free velocity is vx, ind_vx = 3 = 5-ind 

  if(ind == 1)   dwdx = x(2) * (1.d0 + q)                 ! Omega_y
  if(ind == 2)   dwdx = x(1) - y1 * p1 - (y1 + 1.d0) * p2 ! Omega_x 

  ! d Vy / d x =  Omega / vx
  dvdx(1) =  dwdx / x(ind_vel)
   
  ! d Vy/ d V x = - vy / vx
  dvdx(2) = -x(5-ind) / x(ind_vel)
  
  
!  print*, 'vx=', x(3), 'Omega_x = ', dwdx, 'Omega_x/vx = ', dwdx/x(3); read*  
!  print*, 'vy=', x(4), '-vy/vx=', x(4)/x(3)
!  print*, 'dvdx =', dvdx; print*;  read*
!  
!  print*, 'check the index, ind=', ind, 'ind_vel =', ind_vel, & 
!          'the free velocity: ', 5-ind;  read* !ckd
!  print*, 'Omega_x = ', dwdx, 'x(2) * (1.d0 + q)', x(2) * (1.d0 + q); read*  !ckd

!  print*, 'v(ind) = ', x(ind_vel); read*  !ckd
!  print*, 'dvind / dx', dvind_dx ; read*  !ckd
 
return 
end 
  

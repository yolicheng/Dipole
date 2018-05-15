! For the planar normal case problem, where we have invariance in z-direction normal case,   
! Put the related common parameters and  subroutines here 
! N = [0 0 1]   --- the most important one 

! for the safe side, declare the variables that are only used in this module private attribute, 
! and use subroutine init_lf_mod to do the assignment
!
! --- Public -----

! beta = n/w_c, the most important parameter,  where
! 	   n     is angular velocity of the mean orbital rate of the chief’s circular reference orbit
!        w_c   is the magnitude of the dipole’s rotational angular rate

! eq 	 the chosen equilibrium 

! -- 20160731, for beta, ieq,  read the value from the sreen, to make more general routines...
!               
module plf_mod


use  dp_mod
implicit none
save

integer, parameter, private   :: n = 4   ! 4 dimensional state 
integer, private    ::  debug
                                
!  ----  Public ---- 
! the most common used ones, do not need to declare in the main routine 
real(kind=dp) :: eq(n),  beta, runit, vunit,  tunit
integer       :: sgn  

external :: lfunit 

contains


!*******************************************************************************
  !  Initialize  the coordinates of the only equilibirum point 
  
  !  Parameters to be initialized:  beta,  sgn, ieq
  !         and units of length, time, velocity:    runit, tunit, vunit
!*******************************************************************************
  
subroutine init_plf 
  implicit none   
  debug = 0               
  ! keep only case 1 N = [0 0 1]; 

  ! initialization
  eq = 0.d0 
  
  ! the first kind of equilibrium point, q/m <  0
  ! x=\pm ( 1/3 )^(1/3), y=0, z= 0;  
  
  sgn   = -1
  eq(1) = ( 1.d0/3.d0 ) ** (1.d0/3.d0)


! the state of the equilibrium point, positon + velocity 
  print*, 'The equilibirum point: (x  -  y  - z)'
  write(*,'(3f16.10)')   eq(1:3)
    
  print*, 'check init_plf:   sgn, eq',  sgn, eq !ckd
  read*
    
  
end subroutine init_plf
  
  
subroutine init_beta(beta0)
  
  real(kind=dp),intent(in) :: beta0  
  beta = beta0
  
  ! Compute the unit of distance, time and velocity, which are all function of beta 
  call lfunit(beta, runit, tunit, vunit)
 
end 
  
  
!***********************************************************************
! The quantity energy 

! For the conservative quantity, which we refer to as Energy or energy level
! it consists of two parts:
!  1. the common part that is independent on the different cases:   3x²-z²
!  2. Lorentz force-related part, only part that is case sensitive: (x_i+1 ^ 2 + x_i+2 ^ 2)
!     Normal case 1:          2(x²+y²)/r³    -- N = [ 0 0 1]
   
!  3. the sum of the velocity-squared:   - (dotx² + doty²)

!  so we have: 
!  1.  c = 3x²-z² - sgn* 2(x²+y²)/r³ - (dotx² + doty²)
 
! where sgn = :1   if q/m>0;   :-1   if q/m<0

!	Input Parameter
!  pv(4) 	the state vector(position + velocity)

!	Output Parameter
!  cj             the energy

!	Module-based Parameter
!  sgn      :1   if q/m>0;   :-1   if q/m<0

!------------------------------------------------------------------------------
subroutine gr_cjplf(pv, cj)

implicit none 

real(kind=dp), intent(in)  :: pv(n)
real(kind=dp), intent(out) :: cj 
 
! Local Variables 
real(kind=dp) :: x2, y2,  r, r3, dv2 

  x2 = pv(1)**2  ! x^2
  y2 = pv(2)**2  ! y^2
        
  r   = dsqrt( x2 + y2 ) 
  r3  = r**3
  dv2 = pv(3)**2 + pv(4)**2 ! vx^2 + vy^2 

  cj = 3*x2 - dv2  - sgn * 2.d0 * (x2 + y2) / r3 

  return 
end subroutine gr_cjplf


!***********************************************************************
!  *********       deriv_cjplf    ******************
! This subroutine computes the derivative of energy w.r.t the state (position+velocity)
! use: to get a presribed energy for either p.o. or invariant tori 

! According to the subroutine gr_cjlf, we write cj with two parts (common + case_sensitive):  

!  cj = cj_common  - 2*sgn * cj_case,  sgn:1   if q/m>0; -1   if q/m<0

!  -- where cj_common =  3x²-z² - (dotx² + doty² + dotz²),

!  -- cj_case   =  (x_i+1 ^ 2 + x_i+2 ^ 2)
!      Normal case 1:          2(x²+y²)/r³    -- N = [ 0 0 1]

! Common part.... 
!  d cj_common / d x =  6*x 
!  d cj_common / d y =  0
!  d cj_common / d z =  -2*z
 
!  d cj / d vx =  - 2 vx 
!  d cj / d vy =  - 2 vy 
!  d cj / d vz =  - 2 vz 

! Comment: d (1/r³) / dx = -3x / r^5

! --- 1. Normal case:  cj_case = (x²+y²)/r³  ------- 
!  d cj_case / dx =   -x * (x² + y² ) / r^5
!  d cj_case / dy =   -y * (x² + y² ) / r^5 

!	Input Varaibles
!  pv(4) 	the state vector(position + velocity)

!	Output Varaibles
!  dcj(4) 	the derivative of the energy w.r.t the state

!     Module base Varaibles
!  sgn      1 if q/m>0; -1 if q/m<0

!  Finally revised by Yu -- 20160919

!------------------------------------------------------------------------------
subroutine deriv_cjplf(pv, dcj)

implicit none 

real(kind=dp), intent(in)  :: pv(n)
real(kind=dp), intent(out) :: dcj(n) 
 
! Local Variables 
real(kind=dp) :: x, y, x2, y2,  r,  r5,  & 
                 coef_case, aux  

  x  = pv(1);  y  = pv(2)  
  x2 = x*x;    y2 = y*y;   
        
  r  = dsqrt(x2 + y2 ) 
  r5 = r**5 

  ! --- common part ----- 
  dcj(1) = 6*x    ! 6*x 
  dcj(2) = 0.d0   ! 0 

  dcj(3:4) = -2*pv(3:4)    ! -2* [ vx, vy ] 

  !  ----- case sensitive part -------------
  coef_case = -2.d0*sgn 

  ! --- 1. Normal case:  cj_case = (x²+y²)/r³  ------- 
  !  d cj_case / dx =   -x * (x² + y² ) / r^5
  !  d cj_case / dy =   -y * (x² + y² ) / r^5 
  aux = ( x2 + y2 ) / r5
  dcj(1) = dcj(1) - coef_case * (-x) * aux  
  dcj(2) = dcj(2) - coef_case * (-y) * aux 
 

!  print*, 'check dcj'  !ck 
!  print*, 'y0 = ', y
!  print*, 'dcj=', dcj 
!  read* 

  return 
end subroutine deriv_cjplf


!***********************************************************************
! for the case of planar normal case with lorentz force, given the energy level,
! with the initial state as  (x0,y0, vx0,vy0 )  

! keep the position coordinates(x0, y0) and 1 components of velocities,  
! compute the value of X(indv) to achieve the required energy cj0

!  always keep the positive square root for X(indv)


!	Input Varaibles
!  X(n)     the initial state and also the updated one with energy cj
!  cj0      the desired energy
!  indv     the index of component of v to be modified, all the other components remain the same
!           only vx and vy   are allowed to be modified, indv \in [3,4]
!           the initial value of X(indv) is 0


! 	Global variable from the module lf_mod 
!  sgn   	1 if q/m>0; -1 if q/m<0 for the module lf_mod
!
! TODO:
! first try to pass cs as input, then include this file into lf_mod 

! Finally revised by Yu 2016-09-19 18:41:24 
!------------------------------------------------------------------------
subroutine cj2v_plf(X, cj0, indv, isv)

implicit none 
real(kind=dp), intent(inout)  :: X(n)
real(kind=dp), intent(in)  :: cj0
integer, intent(in)  :: indv
integer, intent(out) :: isv

! Local Variables 
real(kind=dp) :: x2, y2, r, r3, sv2, v2, coef_case, cj  

! set the default value to be 1, succeed! 
  isv = 1 

  X(indv) = 0.d0 !put as zero 
 
  x2 = X(1)**2  ! x^2
  y2 = X(2)**2  ! y^2
        
  r   = dsqrt(x2 + y2 ) 
  r3  = r**3
  
  ! vx^2 + vy^2, since X(indv)=0, so sv2 is only the sum of square root of the nonzero components
  sv2 = X(3)**2 + X(4)**2   

  coef_case = sgn * 2.d0 

  ! the common part which include v**2
  ! cj_common = 3*x2 -sv2 
  ! cj = cj_common  - coef_case * (x2 + y2) / r3
  
  ! the real sv2 =  3*x2 - cj_common, substract the nonzero terms, we will have the square of X(indv)
  ! v2 = 3*x2 - cj_common - sv2 
  
  v2 = 3*x2 - ( cj0 + coef_case * (x2 + y2) / r3 ) - sv2

  if( v2 < 0 )  then 
    print*, 'v^2 < 0 !', v2
    isv = 0
    return
!   read(*,*)  
  endif  

  ! always the positive velocity
  X(indv)  = dsqrt(v2) 

  call gr_cjplf(x, cj)
!  print*, 'cj = cj0? ', cj, cj0
!  print*, 'cj2v_lf,  indv, X', indv, X
!  read*
 
  return
end subroutine cj2v_plf


! **********************  dvind_dx_plf **************************
! Compute the differential of the missing velocity '+'PV(indv)  w.r.t. 
! the 2 remaining free coordinates among PV = (x,y, vx, vy), except PV(ind)=p0 and PV(indv)

! For LF, with a given Jacobi Constant and the five values among (x,y, vx, vy), where one of the velocity PV(indv) 
! will be expressed as a function of the energy and the other coordinates

! In order to compute the differential of the linear system for Newton method
! indv = ind + 2  is the index of ind-th velocity 


!     Input Varaibles
! pv    the state vector, of dimension ndim (4  for PLF) (x,y, vx,vy)
! ind   the index of component fixed for Poincare section, pv(ind) = p0 
!       the velocity in this direction is a function of  f(h0, pv)

!     Output Varaibles
! dvdx   the differential of pv(ind) w.r.t. to the free components in state x


! Finally revised by Yu  2016-09-19 18:56:14  ! ckd
!***********************************************************************
subroutine dvind_dx_plf(pv, ind, dvdx)

implicit none 

!  Input and Output
integer, intent(in)     ::  ind
real(kind=dp), intent(in)       :: pv(n) 
real(kind=dp), intent(out)      :: dvdx(n-2) 

! Local Variables 
integer       :: indv, i, k
real(kind=dp) :: x, y, x2, y2, r,  r5,  coef_case, dv2dx(n), aux, cj  
  
  ! the index of ind-th velocity in the full state vector  
  indv = ind + n/2 
  
  x  = pv(1);  y  = pv(2); 
  x2 = x*x;    y2 = y*y;   
    
!  print*, 'Initial state: PV: ', pv 
!  print*, 'x, y', x, y; print*; read*
        
  r  = dsqrt(x2 + y2 ) 
  r5 = r**5 
  
  !  cj =  3x²-z² - (dotx² + doty² + dotz²)  - 2*sgn * cj_case

  ! v(ind)**2 = v2 = 3*x2  - ( cj0 + 2*sgn*cj_case ) - (vx**2 + vy**2  )  
  ! v(ind)    = sqrt(v2) 
  ! **NOTE**: here v(ind) is already computed by cj2v_lf, so it is not zero
  
  ! By chain rule, we have 
  !  d v(ind) / d pv  = d v(ind) / d v2 *  d v2 / d pv 
  !                   = .5 / sqrt(v2) * d v2/ d X = 0.5 / v(ind) * d v2 / d X
  
  ! d v2 / d x = 6x - 2*sgn * d cj_case / dx
  ! d v2 / d y =    - 2*sgn * d cj_case / dy
  
  ! For differential w.r.t. the velocity, ind-th is zero...
  ! d v2 / d vx = - 2 vx
  ! d v2 / d vy = - 2 vy
  
  ! we need to skip  ind-th and  indv-th components from dvdx_all(1:6)  

  ! --- common part ----- 
  dv2dx(1) = 6.d0*x    ! 6*x 
  dv2dx(2) = 0.d0      ! 0 
  dv2dx(3:4) = -2.d0 * pv(3:4)    ! -2* [ vx, vy] 
 
!   to check if this routine is with no error, check the energy 
  call gr_cjplf(pv, cj) ! ckd, cj = cj0  
 
  if(debug == 1) then 
    print*, 'cj = ', cj ;   read* 
  endif   
  
  !  ----- case sensitive part -------------
  coef_case = 2.d0*sgn 
  
!  print*, 'coef_case = ', coef_case, 'sgn =',sgn; print*; read*

  ! --- 1. Normal case:  cj_case = (x²+y²)/r³  ------- 
  !  d cj_case / dx =   -x * (x² + y² - 2z²) / r^5
  !  d cj_case / dy =   -y * (x² + y² - 2z²) / r^5 
  aux = ( x2 + y2 ) / r5
  dv2dx(1) = dv2dx(1) - coef_case * (-x) * aux  
  dv2dx(2) = dv2dx(2) - coef_case * (-y) * aux 

!  print*, 'dv2dx :', dv2dx; print*; read*
 ! skip ind-th and indv-th components 
 !  d v(ind) / d pv = 0.5 / v(ind) * d v2 / d pv
 
  k = 1
  do i = 1, n, 1
  
    if(i == ind .or. i == indv) cycle
    
    dvdx(k) = .5d0 / pv(indv) * dv2dx(i)
    
    k = k+1
    
  end do
  
  !-- ckd 
!  print*, 'check dv2dx = ', dv2dx 
!  print*, 'pv(indv) = ', pv(indv)
!  print*, 'dvdx = ', dvdx ; read* 
 
  return 
end subroutine dvind_dx_plf


!************************************************************************
!  computation of the vector field +  variational equations to be integrated for lorentz force problem 
!  include this model subroutine into the respective module, so we can directly use the system-based paramters
!  without adding to the input parameter list, keep the form deriv(t, y, n, f) to be called in gr_rk78

!     Input parameters:
!  x        time epoch
!  y        state vecto (neq = ndim*(ndim+1) + Variational matrix)
!  neq      neqnumber of equations,  if it equals ndim the variational equation are skipped

!     Output parameters:
!  f(*)     vector field
!           in y(*) and f(*) the first ndim components correspond to the position and velocity
!           the rest are the variational equations
!
!  lf_mod Variable 
!    sgn, beta,  n !  inside the module, no need for use declaration

! Subroutine use: 
!     dflrtz2 

! Finaly revised by Yu -- 20160921 !ckd 
! -----------------------------------------------
subroutine gr_plf(t, y, neq, f)
implicit none 
 
! Input and Output
integer, intent(in) :: neq  ! dimension of the state
real(kind=dp), intent(in)  :: t, y(neq)  
real(kind=dp), intent(out) :: f(neq)

! Local Variable
integer :: i, debug 
real(kind=dp) :: dlf(n,n), x1,x2, dx1,dx2, &
                 b, c, cbt, ccst, d, dbt, dcst, & ! dlfrtz
                 r2, r, r5,  tmp, &
                 dphi(n, n), phi(n,n) ! variatioal matrix
                 

debug = 0 !ckd

!print*, 'For vector field, n has to be 4, n=', n; print*; read*

! x1 -- x;  x2 -- y
! dx, dy 
f(1) = y(3)
f(2) = y(4)
 
x1  = y(1)
x2  = y(2)

dx1 = y(3)
dx2 = y(4)

r2 = x1*x1 + x2*x2
r = dsqrt(r2)
r5 = r**5 
tmp = sgn / r5! the commom first term


! **************** fl_i+1 = fl_x ****************************
! fl_i+1 = sgn*1/R^ 5 * c
! b = x_i+1^ 2 + x_i+2^ 2  
! c = -beta * b * dx_i+2  + x_i+1 * b

b = x1 * x1 + x2 * x2  

 cbt = -b * dx2   ! coefficient of beta
 ccst = x1 * b   ! constant part
 c = cbt * beta + ccst
 
 f(3) = tmp * c
 
 ! **************** fl_i+2  = fl_y ******************************

! fl_i+2 = sgn*1/R^ 5 * d
! b = x_i+1^ 2 + x_i+2 ^ 2 
! d = beta * b * dx_i+1   + x_i+2 * b

!b = x1 * x1 + x2 * x2  ! same as the one in fl_i+1

dbt = b * dx1  ! coefficient of beta
dcst = x2 * b   ! constant part
d = dbt * beta + dcst
f(4) = tmp * d 
 
! the previous just got the lorentz force, what follows is the acceleration 
! adding the left hand side: 
f(3) = f(3) + 2 * y(4) + 3 * y(1) ! fx + 2vy + 3x
f(4) = f(4) - 2 * y(3)            ! fy - 2vx
 
if(neq.eq. 4)  return

! ------------------- Position + Velocity -------------------------------- 


! ----------------Variational Matrix --------------------------------------
! this part need to be check really carefully!

!              \dot phi = Dlf * phi

!print*, 'n , neq ', n, neq ; read*
  
! so here, we need to do a matrix multiplication.
phi = reshape( y(n+1:neq), (/n, n/) )

! To compute Jacobi Matrix of the vector field, only y(1:n) is used, where n is the dimension of state vector, specified in lf_mod 
call dflrtz2(y(1:n), dlf) 

dphi = matmul(dlf, phi)

f(n+1:neq) =  reshape(dphi, (/n*n/)) ! checked, discard the use of equivalence

if (debug == 1) then ! checked   fine!
  print*, 'check the parameter of the system, beta, sgn'
  print*, beta, sgn
  print*, 'check the variational matrix'
  do i = 1, n
    write(*,'(4e18.8)') dlf(i,:) 
  enddo
  print*; read*
endif 

return
end



!**********************************************************************
! Compute the differential of the vector field 

!     Input
!  x0       the state of point of which to compute

!     Output
!  dlf      dimension n X n,  the differential of vector field 


! 	Module Global Variable 
! beta, n, sgn   variable precision, and the demension of the system

! Subroutine use: None
! Finaly revised by Yu -- 20160921 !--ckd 
! -----------------------------------------------
subroutine dflrtz2(x0, dlf2)
implicit none

real(kind=dp), intent(in) :: x0(n)
real(kind=dp), intent(out) :: dlf2(n,n)

! local variables
real(kind=dp) :: x1, x2, dx1, dx2, r2, r, r5, &
                 tmp, tmp1, tmp2, &
                 b, cbt, ccst, dbt, dcst, &
                 dlf(n, n), dlfsub(n/2, n)  

! rule of     x_i, x_i+1, x_i+2
! x -> y  -> z  -> x -> y  
! case 1: N =  [0 0 1]  z  -> x  -> y 
 
! position
x1 = x0(1)
x2 = x0(2)

! velocity
dx1 = x0(3)
dx2 = x0(4)

dlf = 0.d0 ! initialize as zero matrix

r2 = x1*x1 + x2*x2
r = dsqrt(r2)
r5 = r**5 

! fl_z ...... remove this part 


! **************** fl_i+1 = x ******************************
! fl_i+1 = sgn*1/R⁵ * c
! b = x²_i+1 + x²_i+2 - 2*x²_i
! c = -beta * b * dx_i+2 - 3*beta* x_i+2 * x_i * dx_i + x_i+1 * b

b = x1 * x1 + x2 * x2 

 cbt = -b * dx2  ! coefficient of beta
 ccst = x1 * b   ! constant part

!   df_i+1/dx_i+1, df_i+1/dx_i+2,  df_i+1/d dx_i+1, df_i+1/d dx_i+2

tmp = sgn / r5 ! the commom first term


! 1 column df_i+1/dx_i+1
! 1 row, coefficient of beta; 2 row, constant part
tmp1 =  5*x1 / r2  ! independent on beta
dlf(1,1) = tmp *( -2 * x1 * dx2  - tmp1 * cbt )  !beta
dlf(2,1) = tmp *( 3 * x1 * x1 + x2 * x2   - tmp1 * ccst ) !cst


! 2 column df_i+1/dx_i+1
! 1 row, coefficient of beta; 2 row, constant part
tmp2 =  5*x2 / r2  ! independent on beta
dlf(1,2) = tmp *( -2 * x2 * dx2  - tmp2 * cbt )  !beta
dlf(2,2) = tmp *( 2 * x1 * x2 - tmp2 * ccst ) !cst


! 0 !cst

! 3 column df_i+1 / d dx_i+1
! 1 row, coefficient of beta; 2 row, constant part
! 0 
! 0 


! 4 column df_i+1 / d dx_i+2
! 1 row, coefficient of beta; 2 row, constant part

dlf(1,4) =  - sgn * b / r5 !beta 
! 0 !cst


! **************** fl_i+2 =  y *****************************

! fl_i+2 = sgn*1/R^5 * d
! b = x²_i+1 + x²_i+2  
! d = beta * b * dx_i+1   + x_i+2 * b

!b = x1 * x1 + x2 * x2   ! same as the one in fl_î+1

 dbt = b * dx1  ! coefficient of beta
 dcst = x2 * b   ! constant part

!tmp = sgn / r5 ! the commom first term, same with f_i+1

! df_i+2/dx_i+1, df_i+2/dx_i+2,  df_i+2/d dx_i+1, df_i+2/d dx_i+2


! 1 column df_i+2/dx_i+1 
! 3 row, coefficient of beta; 4 row, constant part
!tmp1 =  5*x1 / r2  ! independent on beta
dlf(3,1) = tmp *( 2 * x1 * dx1 - tmp1 * dbt )  !beta
dlf(4,1) = tmp *( 2 * x1 * x2  - tmp1 * dcst ) !cst


! 2 column df_i+2/dx_i+1
! 3 row, coefficient of beta; 4 row, constant part
!tmp2 =  5*x2 / r2  ! same as in fl_i+1, independent on beta
dlf(3,2) = tmp *( 2 * x2 * dx1   - tmp2 * dbt )  !beta
dlf(4,2) = tmp *( 3 * x2 * x2 + x1 * x1 - tmp2 * dcst ) !cst


! 3 column df_i+2 /d dx_i+1  -- dvy -- to check 
! 3 row, coefficient of beta; 4 row, constant part

!b = x1 * x1 + x2 * x2  

dlf(3, 3) =  sgn * b / r5 !beta 
! 0 !cst

! 4 column df_i+2/d dx_i+2
! 3 row, coefficient of beta; 4 row, constant part
! 0 
! 0 

! --- check dlf with the coefficient of beta and the constant part seperately
!print*; print*, 'DF_lz ( coef of beta \\ const)'
!do i = 1, 4, 1
!  print*, dlf(i,:)
!end do
!print*

! initialize dlf2- the 2 dimensional Jacobi Matrix  
dlf2 = 0.d0 
dlf2(1,3) = 1.d0
dlf2(2,4) = 1.d0

! multiply by beta for d v / d x 
dlfsub(1,:) = dlf(1,:) * beta + dlf(2,:)
dlfsub(2,:) = dlf(3,:) * beta + dlf(4,:)

dlf2(3:4, :) = dlfsub                  

! with the normal order, x -> y  
! add the constant matrix 
 
dlf2(3,1) = dlf2(3,1) + 3.d0 ! d vx / d x
dlf2(3,4) = dlf2(3,4) + 2.d0 ! d vx / d vy
dlf2(4,3) = dlf2(4,3) - 2.d0 ! d vy / d vx 
 
return

end subroutine dflrtz2 


!***********************************************************************
!  Compute the vector field of n points simultaneously in  planar lorentz force problem 
!  for the purpose to integrate several points simultaneously or obtain the time-t map for a curve 

! ** NOT ** in order to do integration using rk78, we have to pass y as an  array with only 1 row.... 

!  so we have 
!  1         - n(=4),    the first point 
!  n+1       - 2n,       the second point 
! ......
!  n*(i-1)+1 - 6*i       for i-th point
! where n=4 for plf, is the dimension of the vector field 

!  include this subroutine into the respective module, so we can directly use the system-based paramters
!  without adding to the input parameter list,  keep the form deriv(t, y, n, f) to be called in gr_rk78

!     Input parameters:
!  t              independent variable (typically time)
!  y(*,*)         state of point, of dimension nrow X 4
!  nvar             number of variables (np/6 points) 

!  Output parameters:
!     f(*)       vector field
!                in y(*) and f(*) the first 6 components correspond to the position and velocity
!
! Subroutine use:  gr_lf 

! Finaly revised by Yu -- 20160426
! -----------------------------------------------
subroutine gr_plf_n(t, y, nvar, f)
implicit none 
 
!     Input and Output
!integer,  parameter ::  n = 4
integer, intent(in) ::  nvar  
real(kind=dp), intent(in)  :: t, y(nvar)  
real(kind=dp), intent(out) :: f(nvar)

!     Local Variable
integer       :: i, nrow, debug  
real(kind=dp) :: yi(n), fi(n) 

  debug = 0
!     Resource 
  nrow = nvar / n

  if(debug == 1) then 
    print*, 'nvar=', nvar, 'nrow=', nrow; read*
  endif   

 ! make sure the reshape  of the array is done in the same way, which are are column-wise in fortran 
  do i = 1, nrow 

  ! the index of i-th point saved in y
    yi = y( n*(i-1)+1 : n*i )
  
    call  gr_plf(t, yi, n, fi)
  
    f( n*(i-1)+1 : n*i ) = fi
    
    if(debug == 1 ) then 
      print*, i,'-th point:'
      print*,  yi
      print*,  fi
    endif
     
  enddo 

  if(debug == 1) then   !ckd
    print*; read*
  endif   


  return
end subroutine gr_plf_n



end module plf_mod


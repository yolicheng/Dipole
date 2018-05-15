module lf_mod
!  Put all the commom variables into this module for lorentz force problem, to avoid passing a lot of shared arguments down to all  subroutines that use them

!  for the safe side, declare the variables that are only used in this module private attribute, and use subroutine init_lf_mod to do the assignment

!  use multidimensional array for the equilibirum points,  more convenient for the assignment 

!   --------------------  Index of cases --------------------------
!  case 1 --  normal case,   N = [0 0 1]   --- the most important one 
!  case 2 --  radial case,   N = [1 0 0]   --- studied by Chao Peng, more familes of p.o.s 
!  case 3 --  normal case,   N = [0 1 0]   --- the least interesting one (my point ) 

! TODO - 20151223 -finished for planar symmetry, and x-axis symmetry --20160302
! define the type of symmetry and the corresponding index of the componts to 
!   be as control and target variables 
! at this moment, we only deal with symmetric wrt xz plane, and initial point is 
!   on the y=0 plane
!
! --- Common variables for all 3 cases -----
! beta 	 = n/w_c, the most important parameter,  where
! 	 n   is angular velocity of the mean orbital rate of the chief’s circular reference orbit
!        w_c is the magnitude of the dipole’s rotational angular rate
! eq 	 the chosen equilibrium 

implicit none
save
integer, parameter            :: dp = kind(1.d0),  n = 6 ! 6 dimensional state
integer, parameter, private   :: debug = 0 

! Public -- the most common used ones, do not need to declare in the main routine 
real(kind=dp) :: eq(n), runit, vunit, tunit 
integer       :: sgn !, ieq  
!real(kind=dp) :: beta   

! -- 20160731,  read the value from the sreen, more general routines...
!   beta, ieq,   
 
! Private -- accessible only by the module subroutines and functions
integer, parameter, private  :: neq1 = 3, neq2 = 3,   neq3 = 3
integer, private             :: cs, sgnl(3, 3), ieq      
real(kind=dp),  private      :: x0, eql(3, 3, n) , beta  
                               

contains

! system based subroutines -- in seperate files to avoid a single file being too long  
  include 'dflrtz.f90'       ! put the first three into dlfrtz 

! vector feild   
  include 'gr_lf.f90'        ! for one orbit
  include 'gr_lf_n.f90'      ! for the integration of mutliple orbits simultaneously
 

!  include 'lfunit.f90' ! the function to compute the unit of distance, time and velocity
!  --20160428  added to this file, delete the stand-alone file
  
! ******************************** Subroutines  ********************************  
!1.	subroutine init_lf_mod    		 private constants initialization of basic module parameter
!2.	subroutine init_lf(beta0, cs0, ieq)	
!3.	subroutine gr_cjlf(y, cj)		
!4	subroutine gr_lf(t, y, neq, f)
!5	subroutine gr_lf_n(t, y, np*6, f)	for np points to compute simultaneously

!subroutine dflrtz(x0, dlf)


! ************************** Initialization ************************************
  subroutine init_lf_mod  
  
! According to the index of case, specify the equilibirum points, sgn and beta 
! use the multi-dimensional array for all 3 cases, this is more convenient for our problems, that always have 3 kinds of equilibirum points for each case.

! for the non-zero component, we always take the positive value, and if we need the other ones, handle it in the main routine
  
! finish only case 1, with eq1, sgn1 -20160302
 ! cs : case ---  1: N=[ 0 0 1]; 2: N=[ 1 0 0]; 3: N=[ 0 1 0]; 
  
  ! initialization
  eql = 0.d0 
   
!***************************** 1: N=[0 0 1 ]  *******************************  
! -- 20160519 -- finish case 3 

!  It seems that this case is of the greatest potential,  because of the full symmetry
!  and also because the characteristic polynomials are the same despite the sign of the
!  equilibrium point, so for four equilibirum points, we only need to study one, and the 
!  other three images can be obtained using symmetries.
!  
!     
  ! the first equilibrium point, q/m <  0
  sgnl(1, 1) = -1
  ! x=\pm (  sqrt(1/3))^(1/3), y=0, z= 0;   
  eql(1, 1,  1) = (dsqrt(1.d0/3.d0) ) ** (1.d0/3.d0)
    
  ! the second equilibrium point, q/m > 0 
  ! x 0, y² = 2z², z = = \pm (2 / (3*sqrt(3) )^(1/3) 
  sgnl(1, 2) =  1
  x0 = ( 2.d0 / (3.d0*dsqrt(3.d0)) ) ** (1.d0/3.d0)
 
  eql(1, 2, 3) = x0               ! z
  eql(1, 2, 2) = dsqrt(2.d0) * x0 ! y 

  ! the third equilibrium point, q/m  > 0 
  ! we are taking z = x in this case. 
  sgnl(1, 3) =  1
  ! x = \pm (1/12/sqrt(6) ^(1/3), y = 0,  z² = 5 * x²
  x0 = ( 1.d0/ (12.d0*dsqrt(6.d0)) ) ** (1.d0/3.d0)
  eql(1, 3, 1) = x0 ! x
  eql(1, 3, 3) =  dsqrt(5.d0) * x0 ! z 

 
 !*****************************  2: N=[ 1 0 0]  ********************************
    
  ! the first equilibrium point, q/m > 0
  sgnl(2, 1) = 1
  ! x=0,y=0,z= \pm 1
  eql(2, 1, 3) = 1.d0
    
  ! the second equilibrium point, q/m < 0
  ! x = \pm (2 sqrt(9))^(1/3), y² = 2x²,z = 0
  sgnl(2, 2) = -1
  x0 = ( 2.d0 / (9.d0*dsqrt(3.d0)) ) ** (1.d0/3.d0)
 
  eql(2, 2,1) = x0
  eql(2, 2,2) = dsqrt(2.d0) * x0

  ! the third equilibrium point, q/m < 0, this is what we choose to study, at beta=10, has 3 famililes of center manifold
  ! we are taking z = x in this case. 
  sgnl(2, 3) = -1
  ! x = \pm (1/4/sqrt(2)^(1/3), y = 0,z² = x²
  x0 = ( 1.d0/ (4.d0*dsqrt(2.d0)) ) ** (1.d0/3.d0)
  eql(2, 3,1) = x0
  eql(2, 3,3) = x0 !z0 
  !-----------------------------------------------------------------------------
  

!***************************** 3: N=[0 1 0 ]  ******************************* - TODO
  
!-----------------------------------------------------------------------------


  
  end subroutine init_lf_mod

!*******************************************************************************
  !  Initialize relative parameters for the specified case to study
  
  !  Parameters to be initialized:   beta, cs, sgn, ieq
  !         and units of length, time, velocity:    runit, tunit, vunit
!*******************************************************************************
  
  subroutine init_lf(cs0)

!    real(kind=dp),intent(in) :: beta0  
    integer, intent(in) :: cs0 
    
    ! Local Variables 
    integer :: i !, isgn 

    cs = cs0
! the state of the equilibrium point, positon + velocity 
    print*, 'The equilibirum points: (x  -  y  - z)'
    do i = 1, 3, 1
      write(*,'(3f16.10)')   eql(cs, i, 1:3)
    end do
    
    print*,   'Please input the index of equilibirum point to study: '
    read(*,*) ieq 
      
    eq  = eql(cs,  ieq, :)
    sgn = sgnl(cs, ieq)
    
    print*, 'check init_lf: ieq, sgn, eq', ieq, sgn, eq !ckd
    read*
    
!!  Instead of pass the value of beta from routine... try read from the screen, which is more flexible    
!   if(cs == 1) then 
!     ! if(ieq==1)  beta = 10.d0 ! TODO 
!     if(ieq == 2)  beta = 1.d0
!     if(ieq == 3)  beta = 10.d0 
!  
!   elseif(cs == 3) then 
!     if(ieq == 2) beta = 1.d0
!   endif  
!   
!!   beta = beta0  ! still use this one, to save time...
!    
!!    print*, 'Please inpute beta:'
!!    read(*,*) beta 
!     
!    print*, 'check, cs, ieq,beta', cs, ieq, beta 
!    read* 
!    
!    
!! also compute the unit of distance, time and velocity, which are all function of beta 
!    call lfunit(beta, runit, tunit, vunit)
  
  end subroutine init_lf
  
  
  subroutine init_beta(beta0)
    real(kind=dp),intent(in) :: beta0  
    beta = beta0
  end 
  
  
!***********************************************************************
! The quantity energy is case sensitive  

! For the conservative quantity, which we refer to as Energy or energy level
! it consists of two parts:
!  1. the common part that is independent on the different cases:   3x²-z²
!  2. Lorentz force-related part, only part that is case sensitive: (x_i+1 ^ 2 + x_i+2 ^ 2)
!     Normal case 1:          2(x²+y²)/r³    -- N = [ 0 0 1]
!     Radial case 2:          2(y²+z²)/r³    -- N = [ 1 0 0] 
!     Tangential tcase 3:     2(z²+x²)/r³    -- N = [ 0 1 0]  
   
!  3. the sum of the velocity-squared:   - (dotx² + doty² + dotz²)

!  so we have: 
!  1.  c = 3x²-z² - sgn* 2(x²+y²)/r³ - (dotx² + doty² + dotz²)
!  2.  c = 3x²-z² - sgn* 2(y²+z²)/r³ - (dotx² + doty² + dotz²)
!  3.  c = 3x²-z² - sgn* 2(z²+x²)/r³ - (dotx² + doty² + dotz²)
 
! where sgn = :1   if q/m>0;   :-1   if q/m<0

!	Input Parameter
!  pv(6) 	the state vector(position + velocity)

!	Output Parameter
!  cj 		the energy

!	Module-based Parameter
!  sgn   	:1   if q/m>0;   :-1   if q/m<0
!  cs       case 
!------------------------------------------------------------------------------
subroutine gr_cjlf(pv, cj)

implicit none 

real(kind=dp), intent(in)  :: pv(6)
real(kind=dp), intent(out) :: cj 
 
! Local Variables 
real(kind=dp) :: x2, y2, z2, r, r3, dv2, coef_case, cj_common   

x2 = pv(1)**2  ! x^2
y2 = pv(2)**2  ! y^2
z2 = pv(3)**2  ! z^2
        
r = dsqrt(x2 + y2 + z2) 
r3 = r**3
dv2 = pv(4)**2 + pv(5)**2 + pv(6)**2 

coef_case = sgn * 2.d0 

! the common part 
 cj_common = 3*x2 - z2 - dv2 
 
! case-sensitive part, related to lorentz force 
if ( cs == 1)  then ! x^2 + y^2 

  cj = cj_common  - coef_case * (x2 + y2) / r3 
  
elseif( cs == 2  ) then  ! y^2 + z^2 

  cj = cj_common  - coef_case * (y2 + z2) / r3 
  
elseif( cs == 3 ) then ! z^2 + x^2 

  cj = cj_common  - coef_case * (z2 + x2) / r3 
  
end if

return 
end subroutine gr_cjlf


!***********************************************************************
!  *********       deriv_cjlf    ******************
! This subroutine computes the derivative of energy w.r.t the state (position+velocity)
! use: to get a presribed energy for either p.o. or invariant tori 

! According to the subroutine gr_cjlf, we write cj with two parts (common + case_sensitive):  

!  cj = cj_common  - 2*sgn * cj_case,  sgn:1   if q/m>0; -1   if q/m<0

!  -- where cj_common =  3x²-z² - (dotx² + doty² + dotz²),

!  -- cj_case   =  (x_i+1 ^ 2 + x_i+2 ^ 2)
!      Normal case 1:          2(x²+y²)/r³    -- N = [ 0 0 1]
!      Radial case 2:          2(y²+z²)/r³    -- N = [ 1 0 0] 
!      Tangential tcase 3:     2(z²+x²)/r³    -- N = [ 0 1 0]     

! Common part.... 
!  d cj_common / d x =  6*x 
!  d cj_common / d y =  0
!  d cj_common / d z =  -2*z
 
!  d cj / d vx =  - 2 vx 
!  d cj / d vy =  - 2 vy 
!  d cj / d vz =  - 2 vz 

! Comment: d (1/r³) / dx = -3x / r^5

! --- 1. Normal case:  cj_case = (x²+y²)/r³  ------- 
!  d cj_case / dx =   -x * (x² + y² - 2z²) / r^5
!  d cj_case / dy =   -y * (x² + y² - 2z²) / r^5 
!  d cj_case / dz =   ( -3z / r^5 ) *   (x²+y²)


! --- 2. Radial case:  cj_case = (y²+z²)/r³  ------- 
 
!  d cj_case / dx =   ( -3x / r^5 ) *   (y²+z²)
!  d cj_case / dy =   -y * (y² + z² - 2x²) / r^5
!  d cj_case / dz =   -z * (y² + z² - 2x²) / r^5 


! --- 3. Tangential case:  cj_case = (z²+x²)/r³  ------- 
!  d cj_case / dx =   -x * (z² + x² - 2y²) / r^5  
!  d cj_case / dy =   ( -3y / r^5 ) *   (z²+x²)
!  d cj_case / dz =   -z * (z² + x² - 2y²) / r^5

!	Input Varaibles
!  pv(6) 	the state vector(position + velocity)

!	Output Varaibles
!  dcj(6) 	the derivative of the energy w.r.t the state

!     Module base Varaibles
!  cs, 
!  sgn   	1 if q/m>0; -1 if q/m<0

!  Finally revised by Yu -- 20160731 -- TODO  unchecked.... 
!------------------------------------------------------------------------------
subroutine deriv_cjlf(pv, dcj)

implicit none 

real(kind=dp), intent(in)  :: pv(6)
real(kind=dp), intent(out) :: dcj(6) 
 
! Local Variables 
real(kind=dp) :: x, y, z, x2, y2, z2, r, r3, r5,  & 
                 coef_case, cj_case, cj_common, aux  

x  = pv(1);  y  = pv(2);  z  = pv(3)
x2 = x*x;    y2 = y*y;    z2 = z*z
        
r  = dsqrt(x2 + y2 + z2) 
r5 = r**5 

! --- common part ----- 
dcj(1) = 6*x    ! 6*x 
dcj(2) = 0.d0   ! 0 
dcj(3) = -2*z   ! -2*z 

dcj(4:6) = -2*pv(4:6)    ! -2* [ vx, vy, vz ] 

!  ----- case sensitive part -------------
coef_case = -2*sgn 

if( cs == 1 ) then !

! --- 1. Normal case:  cj_case = (x²+y²)/r³  ------- 
!  d cj_case / dx =   -x * (x² + y² - 2z²) / r^5
!  d cj_case / dy =   -y * (x² + y² - 2z²) / r^5 
!  d cj_case / dz =   ( -3z / r^5 ) *   (x²+y²)
  aux = ( x2 + y2- 2 * z2 ) / r5
  dcj(1) = dcj(1) - coef_case * (-x) * aux  
  dcj(2) = dcj(2) - coef_case * (-y) * aux 
  dcj(3) = dcj(3) - coef_case * (-3*z) / r5 * (x2 + y2)
  
else if( cs == 2  ) then

! --- 2. Radial case:  cj_case = (y²+z²)/r³  ------- 
 
!  d cj_case / dx =   ( -3x / r^5 ) *   (y²+z²)
!  d cj_case / dy =   -y * (y² + z² - 2x²) / r^5
!  d cj_case / dz =   -z * (y² + z² - 2x²) / r^5 

  aux = (y2 + z2- 2 * x2) / r5
  dcj(1) = dcj(1) - coef_case * (-3*x) / r5 * (y2 + z2)
  dcj(2) = dcj(2) - coef_case * (-y) * aux   
  dcj(3) = dcj(3) - coef_case * (-z) * aux 


else if( cs == 3  ) then

! --- 3. Tangential case:  cj_case = (z²+x²)/r³  ------- 
!  d cj_case / dx =   -x * (z² + x² - 2y²) / r^5  
!  d cj_case / dy =   ( -3y / r^5 ) *   (z²+x²)
!  d cj_case / dz =   -z * (z² + x² - 2y²) / r^5

  aux = ( z2 + x2- 2 * y2 ) / r5 
  dcj(1) = dcj(1) - coef_case * (-x) * aux  
  dcj(2) = dcj(2) - coef_case * (-3*y) / r5 * (z2 + x2)
  dcj(3) = dcj(3) - coef_case * (-z) * aux  

end if

print*, 'check dcj'  !ck 
print*, 'y0 = ', y
print*, 'dcj=', dcj 
read* 

return 
end subroutine deriv_cjlf
 
!*******************************************************************************  
! Compute the unit of distance for given beta, the parameters B0 and q2m are specified already
! runit = (abs( B0/n * q2m /beta ) ) **(1./3.) 

! this file saves all the rescaling unit, including distance, time, velocity for Lorentz force problem
! Note: The unit of distance for HCW equations is 1 

!  	Input variables  
!   beta   :    n / wc

! 	Output Variables
!   runit : unit of distance  : km
!   tunit : unit of time      : s
!   vunit : unit of velocity  : km/s
  
!  ROUTINE USED:  none
!  Finally revised by Yu  -20160428
!*************************************************************************
subroutine lfunit(beta, runit, tunit, vunit)

implicit none
integer, parameter:: dp = kind(1.d0) 

real(kind=dp), intent(in):: beta  
real(kind=dp), intent(out):: runit, tunit, vunit 
 
! 	Local variables
integer :: ni, ismp, nsmp 
real(kind=dp) ::  muE, rE, pi, B0, alt, q2m, rc, n  

! a = (B0/n * q/m * 1/beta)^(1/3)  ---  the unit of distance
!  where n is the mean orbital motion of the chief around the Earth, which is in a keplerian,circular orbit 

muE  = 3.986d5 	 ! u_earth -- 		km^3/s^2
rE   = 6371d0    ! radius of Earth --	km

pi =  3.141592653589793 ! common used parameter

!  magnetic moment, take the value of Geomagnetic moment at the surface of the Earth's equator 
!  Peck's paper is 8e15 Wb-m, is the same value in different unit
!  T(tesla) = 1Wb/m^2
B0 = 8.d6; ! Tkm^3

! -----  variables needed to be specified ------------------------
alt = 20200d0  	 ! altitude -- km  --- Medium Earth Orbit

! q2m = [1e-6 1e-5 1e-4 1e-3] -- take the first one to get the minimum a 
q2m = 1.d-6

rc  = rE + alt 	 ! the radius of the chief from the Earth's center
n   = dsqrt(muE/rc**3)   ! rad/s

! -- time 
tunit = 1.d0 / n 

!tunit = 6860.30148727185	! s    -- time  !ckd
!print*, 'tunit=', tunit; ! read*
 

! v = n*rc   ! test with r_geo=42,164(alt=35,786), proves right 

! r_star  = (B0/n * q/m * 1/beta)^(1/3)  --- test with the function assignment, fine!
runit = (dabs( B0/n * q2m / beta ) ) ** (1.d0/3.d0) 

! -- distance  
! runit = 38.0024033332631 	! km -- beta = 1
! runit = 30.1625275142708 	! km -- beta = 2
! runit = 17.6391530962123   	! km -- beta = 10  

vunit = runit / tunit 		! km/s -- velocity - pay attention here, because there 

print*, 'runit(km), tunit(s), vunit(km/s)',  runit, tunit, vunit  !ckd
print* !read*

return
end subroutine lfunit


end module lf_mod


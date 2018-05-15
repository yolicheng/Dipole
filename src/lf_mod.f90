module lf_mod
!  Put all the commom variables into this module for lorentz force problem, to avoid passing a lot of shared arguments down to all  subroutines that use them

!  for the safe side, declare the variables that are only used in this module private attribute, and use subroutine init_lf_mod to do the assignment

!  use multidimensional array for the equilibirum points,  more convenient for the assignment 

!   --------------------  Index of cases --------------------------
!  case 1 --  normal case,   N = [0 0 1]   --- the most important one 
!  case 2 --  radial case,   N = [1 0 0]   --- studied by Chao Peng, more familes of p.o.s 
!  case 3 --  normal case,   N = [0 1 0]   --- the least interesting one (my point ) --TODO 

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
integer, parameter, private   :: dp = kind(1.d0) 
integer, parameter, private   :: n = 6, &  ! 6 dimensional state 
                                 debug = 0 

! Public -- the most common used ones, do not need to declare in the main routine 
real(kind=dp) :: eq(n), runit, vunit,  tunit 
integer       :: sgn, cs, ieq  
real(kind=dp) :: beta   

! -- 20160731,  read the value from the sreen, more general routines...
!   beta, ieq,   
 
! Private -- accessible only by the module subroutines and functions
integer, parameter, private  :: neq1 = 3,   neq2 = 3,  neq3 = 3

! TODO: if we claim cs to be private, it cannot be used by other routines that are not included in this module 

integer, private             :: sgnl(3, 3) !, ieq   ! cs,   
real(kind=dp),  private      :: x0, eql(3, 3, n) !, beta  
                               

contains

! system based subroutines -- in seperate files to avoid a single file being too long  
  include 'dflrtz.f90'       ! put the first three into dlfrtz 

! vector feild   
  include 'gr_lf.f90'        ! for one orbit
  include 'gr_lf_n.f90'      ! for the integration of mutliple orbits simultaneously
 

!  include 'lfunit.f90' ! the function to compute the unit of distance, time and velocity
!  --20160428  added to this file, delete the stand-alone file
  
! ******************************** Subroutines  ********************************  
!1.	subroutine init_lf  
!2.	subroutine init_lf 
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
! TODO: check the coordinates of the equilibria
!     
! update the order of the three kinds of equilibria, q/m. +,+,- 
!  2017-03-03 18:07:00 
if( cs == 1) then 
  
  ! the first equilibrium point, q/m > 0 
  ! x 0, y² = 2z², z = = \pm (2 / (3*sqrt(3) )^(1/3) 
  sgnl(1, 1) =  1
  x0 = ( 2.d0 / (3.d0*dsqrt(3.d0)) ) ** (1.d0/3.d0)
 
  eql(1, 1, 3) = x0               ! z
  eql(1, 1, 2) = dsqrt(2.d0) * x0 ! y 
 
  ! the second equilibrium point, q/m  > 0 
  ! we are taking z = x in this case. 
  sgnl(1, 2) =  1
  ! x = \pm (1/12/sqrt(6) ^(1/3), y = 0,  z² = 5 * x²
  x0 = ( 1.d0/ (12.d0*dsqrt(6.d0)) ) ** (1.d0/3.d0)
  eql(1, 2, 1) = x0               ! x
  eql(1, 2, 3) = dsqrt(5.d0) * x0 ! z 
  
  ! the third equilibrium point, q/m <  0
  sgnl(1, 3) = -1
  ! x=\pm ( 1/3 )^(1/3), y=0, z= 0;    
  eql(1, 3,  1) = ( 1.d0/3.d0 ) ** (1.d0/3.d0)
  

  return


else if (cs == 2)  then 
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
  
  return 


else 
!***************************** 3: N=[0 1 0 ]  ******************************* - TODO
 ! the first equilibrium point, q/m > 0
  sgnl(3, 1) = 1
  
  ! x=0,y=0,z= \pm 1
  eql(3, 1, 3) = 1.d0
    
  ! the second equilibrium point, q/m < 0
  ! x = \pm (2 sqrt(9))^(1/3), y² = 2x²,z = 0
  sgnl(3, 2) = -1
  x0 = ( 1.d0 / 3.d0 ) ** (1.d0/3.d0)
 
  eql(3, 2,1) = x0


  ! the third equilibrium point is quite special 
  ! Regardless of the sign of charge, we have the whole y-axis 
  ! expect the origin as equilibria  
  print*, 'The whole y-axis except (y=0) is equilibiria, Input the value of y'
  read(*,*) x0
  ! x = 0, y = Input, z =0 
  eql(3, 3,2) = x0 ! y
!-----------------------------------------------------------------------------
endif 
  
  end subroutine init_lf_mod
  
  !***********************************************************************
  !     **** eq_dat   ****
  ! Save the coordinates of the equilibria by case in the same file, but
  ! different blocks for plot 
  
  !  Finally Revised by Yu -- 20170228
  !***********************************************************************
  subroutine  eq_dat( ftag)
  implicit none
   
  ! Input  and Output Declaration   
  integer, intent(in)      ::    ftag ! the file tag for equilibria
   
  ! Local Variable
  integer :: ieq, nzero, id1, id2  
  real(kind=dp)  :: pos(3), pos2(3)
  
  
    
    do cs = 1, 3, 1
    
      call init_lf_mod
      ! the state of the equilibrium point, positon + velocity 
      print*, 'The equilibirum points: (x  -  y  - z)'
      
      do ieq = 1, 3, 1
        pos =  eql(cs, ieq, 1:3)
        write(*,'(3f16.10)')   pos
        
        print*, 'Input the number and index of non-zero components'
        read(*,*) nzero;     read(*,*) id1 
        pos2 = pos;   pos2(id1) = -pos(id1)
        write(ftag, '(3f16.10)') pos 
        write(ftag, '(3f16.10)') pos2
        
        if(nzero > 1)  then 
          read(*,*) id2 
          pos2(id2) = -pos(id2) 
          write(ftag, '(3f16.10)') pos2
          pos2(id1) = -pos2(id1)  
          write(ftag, '(3f16.10)') pos2
        endif 
        
        write(ftag, *) ! one blank line 
      end do
      
      write(ftag, *);  write(ftag, *)  
      
    end do
  
    return  
  end subroutine  eq_dat
  
  
  

!*******************************************************************************
  !  Initialize relative parameters for the specified case to study
  
  !  Parameters to be initialized:   cs, ieq, beta, sgn
  !         and units of length, time, velocity:    runit, tunit, vunit
!*******************************************************************************
  
  subroutine init_lf 
    
    ! Local Variables 
    integer :: i  

    print*, 'Select the case: 1:  N=[0 0 1];    2: N=[1 0 0];    3:  N=[0 1 0]' 
    
    read(*,*) cs
    
    call init_lf_mod

   

   ! the state of the equilibrium point, positon + velocity 
    print*, 'The equilibirum points: (x  -  y  - z)'
    do i = 1, 3, 1
      write(*,'(3f16.10)')   eql(cs, i, 1:3)
    end do
    
    print*,   'Please input the index of equilibirum point to study: '
    read(*,*) ieq 
      
    eq  = eql(cs,  ieq, :)
    sgn = sgnl(cs, ieq)
    
    print*, 'check init_lf: ieq, sgn, eq' 
    print*,  ieq, sgn, eq(1:3) !ckd
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
        
r   = dsqrt(x2 + y2 + z2) 
r3  = r**3
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
!  cs       1: normal; 2: radial; 3:tangential  
!  sgn      1 if q/m>0; -1 if q/m<0

!  Finally revised by Yu -- 20160919
! -- TODO  unchecked.... 
!------------------------------------------------------------------------------
subroutine deriv_cjlf(pv, dcj)

implicit none 

real(kind=dp), intent(in)  :: pv(6)
real(kind=dp), intent(out) :: dcj(6) 
 
! Local Variables 
real(kind=dp) :: x, y, z, x2, y2, z2, r,  r5,  & 
                 coef_case,  aux  

! --- check dcj_dx by centre difference
!integer :: i
!real(kind=dp) :: dcj2(6), temp, pv2p(6), pv2m(6), cjp, cjm  

!temp = 1.d-4 
!do i = 1, 6
!  pv2p = pv;  pv2m = pv  
!  pv2p(i) = pv(i) + temp;  
!  pv2m(i) = pv(i) - temp
!  call gr_cjlf(pv2p, cjp)
!  call gr_cjlf(pv2m, cjm)
!  
!  dcj2(i) = (cjp-cjm)/temp/2.d0
!enddo 
!print*, 'check dcj_dx by routine centre difference:'  !ck 
!print*, dcj2 ; print*


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
coef_case =  2*sgn 

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


! --- check dcj_dx ------ 

!print*, 'check dcj_dx by routine deriv_cjlf:'  !ck 
!print*, dcj ; print*

!read* 

return 
end subroutine deriv_cjlf


subroutine cj2v_lf(X, cj0, indv, isv)
! Given the energy level, and the I.C.  (x0,y0,z0, vx0,vy0,vz0),  
! Fix the position coordinates(x0,y0,z0) and two components of velocities, 
! change the value of X(indv) to achieve the required energy cj0

!  always keep the positive square root for X(indv)


!	Input Varaibles
!  X(6)     the initial state and also the updated one with energy cj
!  cj0      the desired energy
!  indv      the index of component of v to be modified, all the other components remain the same
!           only vx,vy,vz are allowed to be modified, indv \in [4,5,6]
!           the initial value of X(indv) is 0


! 	Global variable from the module lf_mod 
!  sgn   	1 if q/m>0; -1 if q/m<0 for the module lf_mod
!
! TODO:
! first try to pass cs as input, then include this file into lf_mod 

! Finally revised by Yu : 2017-05-19 18:03:00 
!------------------------------------------------------------------------

implicit none 

real(kind=dp), intent(inout)  :: X(6)
real(kind=dp), intent(in)     :: cj0
integer, intent(in)  :: indv
integer, intent(out) :: isv

! Local Variables 
real(kind=dp) :: x2, y2, z2, r, r3, dv2, v2, coef_case, cj  

! set the default value to be 1, succeed! 
!print*, 'check cj2v_lf:  cj0, X0, indv?, cs, ieq'
!print*, cj0, X, indv, cs, ieq
!print*; read*


isv = 1 

X(indv) = 0.d0 !put as zero 
 
x2 = X(1)**2  ! x^2
y2 = X(2)**2  ! y^2
z2 = X(3)**2  ! z^2
        
r   = dsqrt(x2 + y2 + z2) 
r3  = r**3
dv2 = X(4)**2 + X(5)**2 + X(6)**2 

coef_case = sgn * 2.d0 

! the common part which include v**2
! cj_common = 3*x2 - z2 - dv2 
! v2 =   3*x2 - z2  - cj_common
 
! case-sensitive part, related to lorentz force 
if ( cs == 1)  then ! x^2 + y^2 

!  cj = cj_common  - coef_case * (x2 + y2) / r3
 
  v2 = 3*x2 - z2 - ( cj0 + coef_case * (x2 + y2) / r3 ) - dv2
  
  
elseif( cs == 2  ) then  ! y^2 + z^2 

!  cj = cj_common  - coef_case * (y2 + z2) / r3 
  
  v2 = 3*x2 - z2 - ( cj0 + coef_case * (y2 + z2) / r3 ) - dv2
  
elseif( cs == 3 ) then ! z^2 + x^2 

!  cj = cj_common  - coef_case * (z2 + x2) / r3 
  
  v2 = 3*x2 - z2 - ( cj0 + coef_case * (z2 + x2) / r3 ) - dv2
  
end if


if( v2 < 0 )  then 
  print*, 'v^2 < 0 !', v2
  isv = 0
  return
  read(*,*)  
endif  

! always the positive velocity
X(indv)  = dsqrt(v2) 

! check the energy 
call gr_cjlf(x, cj)

!print*, 'cj = cj0? ', cj, cj0   ! ckd  2017-05-21 17:05:19 
!print*, 'cj2v_lf,  indv, X', indv, X
!read*
 
return
end subroutine cj2v_lf


! ******************************************************
! Compute the differential of the missing velocity '+'PV(indv)  w.r.t. the four remaining free coordinates among PV=(x,y,z, vx, vy, vz), except PV(ind)=p0 and PV(indv)

! For LF, with a given Jacobi Constant and the five values among (x,y,z, vx, vy, vz), where one of the velocity PV(indv) is expressed as a function of the energy and the other coordinates

! In order to compute the differential of the linear system for Newton method
! indv = ind + 3  is the index of ind-th velocity 


!     Input Varaibles
! pv    the state vector, of dimension ndim (6  for LF) (x,y,z, vx,vy,vz)
! ind   the index of component fixed for Poincare section, pv(ind) = p0 
!       the velocity in this direction is a function of  f(h0, pv)

!     Output Varaibles
! dvdx   the differential of pv(ind) w.r.t. to the free components in state x


! Finally revised by Yu  2016-09-19 18:56:14 
!******************************************************
subroutine dvind_dx_lf(pv, ind, dvdx)
implicit none 

!  Input and Output
integer, intent(in)     ::  ind
real(kind=dp), intent(in)       :: pv(6) 
real(kind=dp), intent(out)      :: dvdx(4) 

! Local Variables 
integer    ::  indv, i, k, debug 
real(kind=dp) :: x, y, z, x2, y2, z2, r,  r5,  coef_case, & 
                 dv2dx(6), aux, cj  
  
  debug = 0
  ! the index of ind-th velocity in the full state vector  
  indv = ind + 3 
  
  x  = pv(1);  y  = pv(2);  z  = pv(3)
  x2 = x*x;    y2 = y*y;    z2 = z*z
        
  r  = dsqrt(x2 + y2 + z2) 
  r5 = r**5 
  
  ! v(ind)**2 = v2 = 3*x2 - z2 - ( cj0 + coef_case* cj_case ) - (vx**2 + vy**2 + vz**2)  
  ! v(ind)    = sqrt(v2) 
  ! **NOTE**: here v(ind) is already computed by cj2v_lf, so it is not zero
  
  ! By chain rule, we have 
  !  d v(ind) / d pv  = d v(ind) / d v2 *  d v2 / d pv 
  !                   = .5 / sqrt(v2) * d v2/ d X = 0.5 / v(ind) * d v2 / d X
  
  ! d v2 / d x = 6x - d cj_case / dx
  ! d v2 / d y =  - d cj_case / dy
  ! d v2 / d z = -2z - d cj_case / dz
  
  ! For differential w.r.t. the velocity, ind-th is zero...
  ! d v2 / d vx = - 2 vx
  ! d v2 / d vy = - 2 vy
  ! d v2 / d vz = - 2 vz
  
  ! we need to skip  ind-th and  indv-th components from dvdx_all(1:6)  

  ! --- common part ----- 
  dv2dx(1) = 6.d0*x    ! 6*x 
  dv2dx(2) = 0.d0      ! 0 
  dv2dx(3) = -2.d0*z   ! -2*z 
  dv2dx(4:6) = -2*pv(4:6)    ! -2* [ vx, vy, vz ] 
 
!   to check if this routine is with no error, check the energy 
  if(debug == 1) then 
    call gr_cjlf(pv, cj) ! ckd, cj = cj0  
    print*, 'cj = ', cj 
    read* 
  endif 
  
  !  ----- case sensitive part -------------
  coef_case =  2.d0*sgn 

  if( cs == 1 ) then !

  ! --- 1. Normal case:  cj_case = (x²+y²)/r³  ------- 
  !  d cj_case / dx =  -x * (x² + y² - 2z²) / r^5
  !  d cj_case / dy =  -y * (x² + y² - 2z²) / r^5 
  !  d cj_case / dz =  ( -3z / r^5 ) * (x²+y²)
    aux = ( x2 + y2- 2 * z2 ) / r5
    dv2dx(1) = dv2dx(1) - coef_case * (-x) * aux  
    dv2dx(2) = dv2dx(2) - coef_case * (-y) * aux 
    dv2dx(3) = dv2dx(3) - coef_case * (-3*z) / r5 * (x2 + y2)
  
  else if( cs == 2  ) then

  ! --- 2. Radial case:  cj_case = (y²+z²)/r³  ------- 
 
  !  d cj_case / dx =  ( -3x / r^5 ) *   (y²+z²)
  !  d cj_case / dy =  -y * (y² + z² - 2x²) / r^5
  !  d cj_case / dz =  -z * (y² + z² - 2x²) / r^5 

  aux = (y2 + z2- 2 * x2) / r5
  dv2dx(1) = dv2dx(1) - coef_case * (-3*x) / r5 * (y2 + z2)
  dv2dx(2) = dv2dx(2) - coef_case * (-y) * aux   
  dv2dx(3) = dv2dx(3) - coef_case * (-z) * aux 


  else if( cs == 3  ) then

  ! --- 3. Tangential case:  cj_case = (z²+x²)/r³  ------- 
  !  d cj_case / dx =  -x * (z² + x² - 2y²) / r^5  
  !  d cj_case / dy =  ( -3y / r^5 ) *   (z²+x²)
  !  d cj_case / dz =  -z * (z² + x² - 2y²) / r^5

    aux = ( z2 + x2- 2 * y2 ) / r5 
    dv2dx(1) = dv2dx(1) - coef_case * (-x) * aux  
    dv2dx(2) = dv2dx(2) - coef_case * (-3*y) / r5 * (z2 + x2)
    dv2dx(3) = dv2dx(3) - coef_case * (-z) * aux  

  end if

  ! skip ind-th and indv-th components 
  ! d v(ind) / d pv = 0.5 / v(ind) * d v2 / d pv

  k = 1
  do i = 1, 6, 1
  
    if(i == ind .or. i == indv) cycle
    
    dvdx(k) = .5d0 / pv(indv) * dv2dx(i)
!    print*, 'i, k, dvdx:',i, k, dvdx ! ckd
    k = k+1
  enddo
  
  if(debug == 1) then 
    print*, 'check dv2dx = ', dv2dx  !ckd  2017-05-21 16:53:05 
    print*, 'indv, pv(indv) = ', indv,    pv(indv)
    print*, 'dvdx = ',        dvdx 
  endif 
  
return 
end subroutine dvind_dx_lf
  
  


end module lf_mod


program main_tori 

!  this program starts with a purely numerical invariant tori. and try to use the approach by Josep Maria to parameterize.
!  Reference: 
!        1. Advanced Course on Long Time Integrations, course note, 2007
!        2. The dynamics around the collinear equilibrium points of the RTBP, 2001

!   Notation:  
!     flow: \phi,    parameterization of invariant tori: \psi,    invariant curve: \varphi   (use numerically computed poincare map \gamma to compare) 
!     \psi is a map goes from R2 to R6:  (theta1, theta2)  -- map to -->  \psi(theta1, theta2) (the position and velocity)

!  We start with a numerically known invariant tori, so this is easier than the real goal... 
!    To get the initial guess of two frequencies, we do the FFT for the curve \gamma 
!   --0.1  Compute the invariant curve \gamma numerically, just integrare and compute the poincare map at the section z=0.56 (initial guess), 
!          could be the one of the equilibrium point
! 
!   --0.2  Do Fourier Analysis on \gamma to compute (omega1, omega2)  --->  (rou, t2)
!          the problem is the points we obtained are not equally-spaced in time, so FFT will not be suitable, try Gerard's routine

!   --0.3  Compute the Fourier series representation of \gamma  


!   --1.  Compute the invariant curve   \varphi(xi) = \psi(xi,0)
!         invariance condition:  \varphi ( xi + rou) = \phi_{t2} ( \varphi(xi) )  

! try to finish the first one....
!        instead of using omega_1 and omega_2 as the unknowns, we use the rotation number rou and time t2 of one revolution as the unknows. 
!        1.1  -- First we parameterize the space \xi  into  i*2pi/ (1 + 2*Nf), i = 1, ..., 2*Nf

!        1.2  -- slove \varphi as a truncated Fourier series.
!                \varphi (xi) = A_0 + Sigma_ k from 1 to Nf ( C_k cos(k*xi) + S_k cos(k*xi))

!        1.3  -- for each xi_i, the unknowns coefficients A's and B's are the same, because they are supposed to be on the same curve. 
!                write the equations carefully, 

!     initial gusee for Nf:  5, as suggested by Gerard, this is enough.
! 

!  Note1:  To avoid indetermination in invariant curve, the initial section , theta2 = 0 (could any value in [0:2pi]), which we take to obtain on the invariant curve
!     we fix a value for a coordinate of A0
!     first approach:  A0_x  in x component 
!
!  Note2:  To avoid indetermination in phase shift, any xi_0 such that xi-xi0 will satisfies the invariant conditon 
!         we fix a value for a coordinate of A1 to be zero 
!         first approach:  A1_x  in x component 

!  Note3:  we fix the energy level

! 

!           -- 1.1   integration for one revolution,  

!   --2.  Globalize the invariant tori    \psi(theta1 + omega1 * t, theta2 + omega2 * t) = \phi_t ( theta1, theta2)
!  TODO: 
!         invariance condition:  \varphi ( xi + rou) = \phi_{t2} ( \varphi(xi) )  

! 


! take eq3, beta=10,   4 elliptic equilibrium points in this case  --- radial case, just for test, finnaly we will go to normal case(the most important one)
! this is the most stable case 



! module based parameters
 !n =6,  runit, vunit, tunit, initialized along with beta
use fft_mod 
use lf_mod ! the system related parameters of lorentz force problem

implicit none

! 	Module lf_mod based Parameter
!  dp, n, 		! the dimension of the problem is defined in lf_mod, n=6 
!  runit, vunit, tunit  ! unit of the fundamental quantities


integer, parameter  ::   np    = 200, &   ! number of sample points in the perturbed ball 
                         namp  = 10 ! number of values of deviation taken into consideration 
                         
real(kind=dp)       :: pi = 4.d0*datan(1.d0),  day = 24.d0 * 60.d0 * 60.d0  

integer ::        cs, ieq,  debug 
real(kind=dp) ::  beta 

! Local Variables
real(kind=dp) :: rsgm, vsgm, ramp, vamp, & ! Num: error in normal distribution
                 tic, toc, tf , xmax, dsmax, ymax(6), ymaxl(8,6), & ! dsmax 
                 tol, tolres, rmm, vmm, dsmm !fft
                  
integer ::  i, j,  idprop,  foball, fdsall, fdsmax, isob, iscap, nfam, & 
            nfmax, nsm, fob, ffc, ffcbas,  nitmax 
            
character(len=70) :: fnoball, fndsall, fndsmax,  & ! for all the np points
                     fnob, fnfc, fnfcbas,  fndr   

!	subroutine from packages declared external here
real(kind=dp), external :: dnrm2 ! from BLAS-LEVEL1

! radius of outer bound for escape, nondimensinal unit 1 is big enough....
xmax =1.d0 

! for sampling, 2**12=4096, in general, it is enough  
nsm  = 2**14


nfmax = 10 

!  time interval:  1.d4 ~= 800 days, 
!tf = 1.d2 / 2.d0  ! day-- discard 
!print*, 'day=', day 
!tf = tf * day / tunit  ! 1.d3 

!  but it is better to take adimensional unit as the power of 10  
! integration time: adimensional
tf = 1.d3 ! for 3 digits in frequency

! debug or not 
debug =  0

! Initialize the private variables eq1, sgn1 for  case 1
call init_lf_mod

! case 1: N=[1 0 0]; 2:  N=[0 1 0]; 3:  N=[0 0 1]; 
 cs = 1 ! use 2 to test the swap rule is coded correctly--ckd!

! the index of the equilibrium points
! 1:  q/m > 0,  x=0,y=0,z= \pm 1
! 2:  q/m < 0,  x = \pm (2 sqrt(9))^(1/3), y² = 2x²,z = 0
! 3:  q/m < 0   x = \pm (1/4/sqrt(2)^(1/3), y = 0,z² = x²  ! this the the case we study now 

ieq = 3 ! 3, x,0,z is the case that we study currently

! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy
beta = 10.d0   !  6 dimension center manifold 


! Provided the case, beta and ieq,  beta, cs, sgn and eq are initialized by subroutine init_lf from module lf_mod
call init_lf(beta, cs, ieq) 

print*, 'check the assignment with module' !--ckd
print*, 'beta, cs, sgn, eq', beta, cs, sgn, eq
! -- the real units-----
print*, 'runit(km), tunit(s), vuint(km/s): ', runit, tunit, vunit
!read*


! assignment of the 4 symmetric equilibrium points
!do i = 1, np 
!  eqi(i, : ) = eq
!enddo 

! 1st eq: x = x_eq,  z = x 
! the second eq, x=-x_eq, z = x 
!eqi(2, 1:3:2)  = -eq(1:3:2)

!! the third:  x = x_eq,  z = -x 
!eqi(3, 3) = -eq(3)

!! the fourth eq, x= -x_eq, z = -x 
!eqi(4, 1) = -eq(1)


! ---- integrate the np points together, to obtain position and velocity at the same epoch --------------
! the filename to save the useful data, the integration, relative position+velocity, angular momentum
fdsmax = 21;  fdsall = 22

! use new name, to avoid overwritten the good ones, and the type of approach selected is used as suffix
fndsall    =  './dat/tori/dsall.dat'
fndsmax    =  './dat/tori/dsmax.dat'

print*,   fndsmax, fndsall; print*   !ck

!open(fdsall, file=fndsall, access ='append',status='replace')
!write(fdsall,*) '# The maximum distance at each epoch  : t       ds'

!open(fdsmax, file=fndsmax, access ='append',status='replace')
!write(fdsmax,*) '# The maximum distance : ramp     vamp    dsmax   (x y z vx vy vz)_0'

print*, 'tf = ', tf;  !read* 


! ---- initialize the deviation in position and velocity -----------------------
!  to obtain 2-sigma deviation, we set the standard deviation to be the maximum vaule/2
! but we are supposed to check the maximum value of ramp and vamp such that the error make the orbit escape...
rsgm =  1.d-1 / 2.d0 !  km     	100 m
vsgm =  3.d-6 / 2.d0 !  km/s 	3 mm/s

! try smoother increase, 100 m each step 
!amp = (/1.d0, 2.d0, 5.d0, 8.d1, 1.d1, 2.d1/)


isob = 0 ! do not save the propagation data for all the points, not necessary 

!do i = 1, namp, 1  

!  call cpu_time(tic)
!  
!!  ramp = rampl(i) ! the size of rampl is fixed, not flexible enough...
!!  vamp = vampl(i)

!  ramp =  i * rsgm / runit  
!  vamp =  i * vsgm / vunit   

!  call  eqprop_num( eq, np, tf, xmax, ramp, vamp, fdsall, fdsmax,  gr_lf_n, iscap, dsmax, ymax)
!  
!  call cpu_time(toc)
!  print*, 'Elapsed time for sampling: ', toc-tic; print*
!  
!  if (iscap == 1) then 
!    print*, 'Escape ! ramp=', ramp*runit, '   ,vamp=', vamp*vunit 
!    read*
!    exit
!  endif 
!  
!enddo 
 
! close(fdsmax)

!! -------- Fourier Analysis ------------ 
write(*,*) ' Maximum amplitude of the residual //  Tolerance in correction for refinement? (suggestion: 1.d-4, 1.d-6)'
read(*,*)  tolres,  tol 
        
!  tol = 1.d-8; tolres= 1.d-5;  ! too much, fail with 1.d-5... I don't know why
nitmax = 15  
        
! save the data for Fourier analysis of the orbit which has the largest distance from the equilibrium
fob   = 25;    ffc   = 26  ; ffcbas  = 27
fnob      =  './dat/tori/ob.dat'
fnfc      = './dat/tori/fc.dat'
fnfcbas   = './dat/tori/fcbas.dat'

print*, fnob, fnfc, fnfcbas;  print*   !ck

open(fob  ,file=fnob  , access ='append',status='replace')
write(fob  ,*)  '# Each block for one ramp:  t	(x y z vx vy vz)'

open(ffc  ,    file=fnfc  ,   access ='append',status='replace')
open(ffcbas  , file=fnfcbas,  access ='append',status='replace')


open(fdsmax, file=fndsmax) 
read(fdsmax,*) ! first line is comment 

! i is the number of amplitudes from above 
i = 1

tf = 4.d2 

do j = 1, i, 1  

!write(fdsmax, '(9e24.14)') ramp, vamp, dsmax, ymax 
  read(fdsmax, '(9e24.14, i5)') rmm, vmm, dsmm, ymax, iscap
  write (*,*) ymax
  
   if(iscap == 1) exit 
  ! we need to handle specially for the escape case TODO 
  xmax =  2.d0            
  if(j==8) call fft_ob(ymax, tf, nsm, tol, tolres, nitmax, nfmax, &
              xmax, fob, ffc, ffcbas, gr_lf,  gr_cjlf)  
  read*
enddo   
  
stop
end 





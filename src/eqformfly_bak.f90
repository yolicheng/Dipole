program eqformfly 

! take eq3, beta=10, 4 elliptic equilibrium points in this case  
! this is the most stable case 


! ******************************** two options ********************************* 
! 1.  numerical : 
!   integrate samples from a ball around the equilibrium point,  for a given time interval [0:100] days 
!   see how this ball evolves, to which value of the radius will the initial points from this ball escape during this time interval .... 
 
!  2. analytical: 
!    use the linear solution of the variational equation as the deviation from the equilibrium point, 
!    and we have three modes related to the three pairs of pure imaginary eigenvalues of the Jacobi Matrix of the differential of the vector field.
 
!    and for each mode, we study the stability by incease the value of the magnitude of this mode, and then add the state of the equilibrium point, 
!    use this as the initial condition and integrate with the nonlinear equations of motion also for 100 days, to see up to which value will each mode remain in the neibougher 


! here we try the first approach ---- purely numerical 


!! 1. The relative position and velocity evolution of the 4(or 3) symmetric equilibrium points. Does the plane(the angular momentum) that contains any three points remains?

! 2. 2-sigma Gaussian random deviations are added as the tracking error on the initial state. 

!      SUBROUTINE DLARNV( IDIST, ISEED, N, X ) from LAPACK 
!      -- provided 4-dimension initial seed iseed for dlarnv(), with iseed(4) odd
!      -- idist specifies the distribution of the random numbers:
!           = 1:  uniform (0,1)   
!           = 2:  uniform (-1,1)
!           = 3:  normal  (0,1)   -- we take this one


! 3. the cost to move the orbit back to  the initial point (the cost of station keeping), when a control is necessary to be applied
!    check the magnitude of the position deviation... 

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
foball = 20;  fdsmax = 21;  fdsall = 22

! use new name, to avoid overwritten the good ones, and the type of approach selected is used as suffix
fnoball    =  './dat/eq/oball.dat'  ! in fact, doesn't use this one...
fndsall    =  './dat/eq/dsall.dat'
fndsmax    =  './dat/eq/dsmax.dat'

print*, fnoball, fndsmax, fndsall; print*   !ck
!read(*,*)

! in fact, this one is not used 
open(foball,file=fnoball, access ='append',status='replace') 
write(foball,*)  '# Every 200 rows for 1 epoch:  t	(x y z vx vy vz)'

open(fdsall, file=fndsall, access ='append',status='replace')
write(fdsall,*) '# The maximum distance at each epoch  : t       ds'

open(fdsmax, file=fndsmax, access ='append',status='replace')
write(fdsmax,*) '# The maximum distance : ramp     vamp    dsmax   (x y z vx vy vz)_0'

print*, 'tf = ', tf;  !read* 

write(fdsall,'(a,f8.2,a)') '# Time interval: (0: ', tf*tunit/day, '] days'


! ---- initialize the deviation in position and velocity -----------------------
!  to obtain 2-sigma deviation, we set the standard deviation to be the maximum vaule/2
! but we are supposed to check the maximum value of ramp and vamp such that the error make the orbit escape...
rsgm =  1.d-1 / 2.d0 !  km     	100 m
vsgm =  3.d-6 / 2.d0 !  km/s 	3 mm/s


! we test the 6 values  ! 1.d1 is already too large, the maximum distance escapes after less than 20 days
!amp = (/1.d0, 2.d0, 5.d0, 1.d1, 2.d1, 5.d1/)

! try smoother increase, 100 m each step 
!amp = (/1.d0, 2.d0, 5.d0, 8.d1, 1.d1, 2.d1/)

!! in  adimensional unit
!rampl =  amp * rsgm / runit  
!vampl =  amp * vsgm / vunit   

 ! the result from eqprop_num
! ymaxl(1,:) = (/0.57495730308022d0, 0.27466593181008d-1, 0.58670598733340d0, &
!            -0.42129195807698d-2, 0.64065383550117d-2, -0.43714800803090d-2/)
! ymaxl(2,:) = (/ 0.59390093276441d0, 0.84833507039811d-1, 0.62180471971656d0, &
!                -0.10416616116962d-1, 0.15091427402557d-1, -0.16060114064986d-1/)
! ymaxl(3,:) = (/0.60437839476031d0, -0.10961977545897d0, 0.64293222376783d0, &
!                 0.12907794668116d-1,  0.21674386316663d-1, 0.22218488179684d-1/)
! ymaxl(4,:) = (/0.61631497079643d0, 0.13710213002118d0, 0.67009617556981d0, &
!                -0.15457745064804d-1, 0.28958023443637d-1, -0.28526480727998d-1/)
! ymaxl(5,:) = (/0.62501368846202d0, 0.18983878030009d0, 0.68440234711280d0, &
!                -0.23314650388083d-1, 0.33512191020045d-1, -0.38900086569518d-1/)
! ymaxl(6,:) = (/0.63901166698081d0, -0.27637804434103d0, 0.71223372389798d0, &
!                0.35684239302859d-1,  0.41709022341248d-1, 0.56418315737744d-1/)
! ymaxl(7,:) = (/0.63115721940913d0, -0.29058618573802d0, 0.68457239331459d0, &
!                 0.37923221160710d-1, 0.32809074843776d-1, 0.54762477753652d-1/)
! ymaxl(8,:) = (/0.62139211385632d0, 0.43019333662142d0, 0.66634514666535d0, &
!                -0.50813367195214d-1, 0.19891230584956d-1, -0.55631192634020d-1/)
!    x0=  0.61687701626196212       0.43811945228551907       0.68750882575066419      -0.10510648053364617        3.3449333705340804E-002  -6.9270976139865312E-002
!             

isob = 0 ! do not save the propagation data for all the points, not necessary 

do i = 1, namp, 1  

  call cpu_time(tic)
  
!  ramp = rampl(i) ! the size of rampl is fixed, not flexible enough...
!  vamp = vampl(i)

  ramp =  i * rsgm / runit  
  vamp =  i * vsgm / vunit   

  call  eqprop_num( eq, np, tf, xmax, ramp, vamp, isob, foball, fdsall, fdsmax,  gr_lf_n, iscap, dsmax, ymax)
  
  call cpu_time(toc)
  print*, 'Elapsed time for sampling: ', toc-tic; print*
  
  if (iscap == 1) then 
    print*, 'Escape ! ramp=', ramp*runit, '   ,vamp=', vamp*vunit 
    read*
    exit
  endif 
  
enddo 
 
 close(fdsmax)

! -------- Fourier Analysis ------------ 
write(*,*) ' Maximum amplitude of the residual //  Tolerance in correction for refinement? (suggestion: 1.d-4, 1.d-6)'
read(*,*)  tolres,  tol 
      
        
!  tol = 1.d-8; tolres= 1.d-5;  ! too much, fail with 1.d-5... I don't know why
nitmax = 15  
        
! save the data for Fourier analysis of the orbit which has the largest distance from the equilibrium
fob   = 25;    ffc   = 26  ; ffcbas  = 27
fnob      =  './dat/eq/ob.dat'
fnfc      = './dat/eq/fc.dat'
fnfcbas   = './dat/eq/fcbas.dat'

print*, fnob, fnfc, fnfcbas;  print*   !ck

open(fob  ,file=fnob  , access ='append',status='replace')
write(fob  ,*)  '# Each block for one ramp:  t	(x y z vx vy vz)'

open(ffc  ,    file=fnfc  ,   access ='append',status='replace')
open(ffcbas  , file=fnfcbas,  access ='append',status='replace')

 
!0.99793136994836d0 0.12825255985786d1 0.28087189607075d0 0.11910819065319d1 0.10451573176095d1 0.10144395476379d1
! ymax = ymaxl(i,:)
! 0.22676825685350E-01    0.46670958292742E-02    0.45933779953807E+00   0.61687701626196E+00    0.43811945228552E+00    0.68750882575066E+00   -0.10510648053365E+00    0.33449333705341E-01   -0.69270976139865E-01
!   ymax = (/ 0.61687701626196d0,    0.43811945228552d0,    0.68750882575066d0,  &
!            -0.10510648053365d0,     0.33449333705341d-1,   -0.69270976139865d-1/)
  
!  tf   = 0.34858334376319d1 
!  ymax = (/ 0.d0,  0.95416926423625d0,  0.54917304800717d0,  0.34667273261016d0, 0.d0, 0.d0 /) 
  
!   tf = 0.42841573866374d1
!   ymax =  (/0.55905421061878d0,  0.d0,  0.55534599921232d0, 0.d0, -0.10440740360204d-2,  0.d0/)  
!          
 ! TODO -- dsmax, ymax, will be used here for fft analysis, use fixed-step rk78 integration 
!  subroutine fft_ob(x0, tf, nfmax, np, fob, ffc  , deriv, gr_cj) 

! iamp == 1, the first quasi-periodic orbit
! ymax = (/0.57495730308022d0,    0.27466593181008d-1,    0.58670598733340d0 , &
!        -0.42129195807698d-2,    0.64065383550117d-2,   -0.43714800803090d-2/)
   
! subroutine fft_ob(x0, lt, nsm, tol, tolres, nitmax, nfmax, &
!                  xmax, fob, ffc, deriv, gr_cj) 
! read from the file  fdsmax 
open(fdsmax, file=fndsmax) 
read(fdsmax,*) ! first line is comment 

! i is the number of amplitudes from above 

do j = 1, i, 1  

!write(fdsmax, '(9e24.14)') ramp, vamp, dsmax, ymax 
  read(fdsmax, '(9e24.14, i5)') rmm, vmm, dsmm, ymax, iscap
  write (*,*) ymax
  
   if(iscap == 1) exit 
  ! we need to handle specially for the escape case TODO 
              
  call fft_ob(ymax, tf, nsm, tol, tolres, nitmax, nfmax, &
              xmax, fob, ffc, ffcbas, gr_lf,  gr_cjlf)  
  read*
enddo   
  
stop
end 





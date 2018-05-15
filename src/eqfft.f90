!***********************************************************************
!     ****      ****
!

!       Input Variables 
!  a            

!       Output Variables 
!  a            

!  Routine Used:
!     

!  Finally Revised by Yu -- 20160
!***********************************************************************

subroutine ( in, out  optional)
implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration   
real(kind=dp), intent(in)      ::  in  
real(kind=dp), intent(out)      ::  out 
 
! Local Variable
real(kind=dp)  ::
  
    
  source

  return  
end subroutine 


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
fndsall    =  './dat/eq/dsall.dat'
fndsmax    =  './dat/eq/dsmax.dat'

print*,   fndsmax, fndsall; print*   !ck

open(fdsall, file=fndsall, access ='append',status='replace')
write(fdsall,*) '# The maximum distance at each epoch  : t       ds'

open(fdsmax, file=fndsmax, access ='append',status='replace')
write(fdsmax,*) '# The maximum distance : ramp     vamp    dsmax   (x y z vx vy vz)_0'

print*, 'tf = ', tf;  !read* 


! ---- initialize the deviation in position and velocity -----------------------
!  to obtain 2-sigma deviation, we set the standard deviation to be the maximum vaule/2
! but we are supposed to check the maximum value of ramp and vamp such that the error make the orbit escape...
rsgm =  1.d-1 / 2.d0 !  km     	100 m
vsgm =  3.d-6 / 2.d0 !  km/s 	3 mm/s

! try smoother increase, 100 m each step 
!amp = (/1.d0, 2.d0, 5.d0, 8.d1, 1.d1, 2.d1/)


isob = 0 ! do not save the propagation data for all the points, not necessary 

do i = 1, namp, 1  

  call cpu_time(tic)
  
!  ramp = rampl(i) ! the size of rampl is fixed, not flexible enough...
!  vamp = vampl(i)

  ramp =  i * rsgm / runit  
  vamp =  i * vsgm / vunit   

  call  eqprop_num( eq, np, tf, xmax, ramp, vamp, fdsall, fdsmax,  gr_lf_n, iscap, dsmax, ymax)
  
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
  
end 





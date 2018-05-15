program curve
!  numerically compute the invariant curve, data saved in ./dat/curve/ob.dat(w2)   ./dat/curve/ob.dat(w1), where w1 and w2 means along this frequency
!  and the time step is the period of another frequency ( 1./ w2 (w1) )

! 
! take eq3, beta=10, 4 elliptic equilibrium points in this case  
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
                 tol, tolres, rmm, vmm, dsmm, & !fft
                 h,  xl(np,6), w1, w2 ! curve
                  
integer ::  i, j,  idprop,  foball, fdsall, fdsmax, isob, iscap, nfam, & 
            nfmax, nsm, fob, ffc, ffcbas,  nitmax 
            
character(len=70) :: fnoball, fndsall, fndsmax,  & ! for all the np points
                     fnob, fnfc, fnfcbas,  fndr   

!	subroutine from packages declared external here
real(kind=dp), external :: dnrm2 ! from BLAS-LEVEL1

! radius of outer bound for escape, nondimensinal unit 1 is big enough....
xmax = 1.d0 

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
 cs = 2 ! use 2 to test the swap rule is coded correctly--ckd!

! the index of the equilibrium points
! 1:  q/m > 0,  x=0,y=0,z= \pm 1
! 2:  q/m < 0,  x = \pm (2 sqrt(9))^(1/3), y² = 2x²,z = 0
! 3:  q/m < 0   x = \pm (1/4/sqrt(2)^(1/3), y = 0,z² = x²  ! this the the case we study now 

ieq = 3 ! 3, x,0,z is the case that we study currently

! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy
!beta = 10.d0   !  6 dimension center manifold 


! Intialize the state vector of the equilibirum points 
call init_lf(cs) 

!  Instead of pass the value of beta from routine... try read from the screen, which is more flexible    
if(cs == 2) then 
   ! if(ieq==1)  beta = 10.d0 ! TODO 
   if(ieq == 2)  beta = 1.d0
   if(ieq == 3)  beta = 10.d0 
  
 elseif(cs == 1) then 
   if(ieq == 2) beta = 1.d0
 endif  
 
call init_beta(beta)   
!   beta = beta0  ! still use this one, to save time...
    
!    print*, 'Please inpute beta:'
!    read(*,*) beta 
     
print*, 'check, cs, ieq,beta', cs, ieq, beta 
read* 
    
    
! also compute the unit of distance, time and velocity, which are all function of beta 
 call lfunit(beta, runit, tunit, vunit)
    
! -- the real units-----
print*, 'runit(km), tunit(s), vuint(km/s): ', runit, tunit, vunit
!read*

! ---- integrate the np points together, to obtain position and velocity at the same epoch --------------
! the filename to save the useful data, the integration, relative position+velocity, angular momentum
fdsmax = 21;  !fdsall = 22

! use new name, to avoid overwritten the good ones, and the type of approach selected is used as suffix
fndsmax    =  './dat/tori/dsmax.dat'

!! -------- Fourier Analysis ------------ 
!write(*,*) ' Maximum amplitude of the residual //  Tolerance in correction for refinement? (suggestion: 1.d-4, 1.d-6)'
!read(*,*)  tolres,  tol 
        
!  tol = 1.d-8; tolres= 1.d-5;  ! too much, fail with 1.d-5... I don't know why
nitmax = 15  
        
! save the data for Fourier analysis of the orbit which has the largest distance from the equilibrium
fob   = 25;    ffc   = 26  ; ffcbas  = 27
!fnob      =  './dat/curve/ob1.dat'
!fnob      =  './dat/curve/tori_fxd.dat'
fnob      =  './dat/curve/ob2.dat'

fnfc      = './dat/curve/fc.dat'
fnfcbas   = './dat/curve/fcbas.dat'

print*, fnob, fnfc, fnfcbas;  print*   !ck

open(fob  ,file=fnob  , access ='append',status='replace')
write(fob  ,*)  '# Each block for one ramp:  t	(x y z vx vy vz)'

open(ffc  ,    file=fnfc  ,   access ='append',status='replace')
open(ffcbas  , file=fnfcbas,  access ='append',status='replace')


open(fdsmax, file=fndsmax) 
read(fdsmax,*) ! first line is comment 

! i is the number of amplitudes from above 
i = 1
!tf = 4.d2 
w1 = 0.233347
w2 = 0.021967

!h = 1.d0/w1

h = 1.d0/w1 

!np = w1/w2 ! np is declared as parameter, cannot change value
i = 1

do j = 1, i, 1  

!write(fdsmax, '(9e24.14)') ramp, vamp, dsmax, ymax 
  read(fdsmax, '(9e24.14, i5)') rmm, vmm, dsmm, ymax, iscap
  write (*,*) ymax
  
!   if(iscap == 1) exit 
  ! we need to handle specially for the escape case TODO 
  xmax =  2.d0            
!  call fft_ob(ymax, tf, nsm, tol, tolres, nitmax, nfmax, &
!              xmax, fob, ffc, ffcbas, gr_lf,  gr_cjlf)  

 
   print*, 'Before plob_fxd, x0='
   print*,  ymax ; print*

!  subroutine plob_fxd(y0,t0, h0, np, xmax, ftag, yall, deriv, gr_cj)
   call plob_fxd(ymax, 0.d0, h, np, xmax, fob, xl, gr_lf,  gr_cjlf)
enddo   
  
stop
end 


 




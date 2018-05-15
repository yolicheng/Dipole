program curve_poinc
!  Use Poincare map method to numerically compute the invariant curve, 
!  data saved in ./dat/curve_poinc/ob.dat 

! Advantage of this approach, no information of the frequencies and rotation number needed in advance 

! Still, start with an easy example, 2d tori, where we have eq3, beta=10, 4 elliptic equilibrium points in this case  
! this is the most stable case 

! module based parameters
 !n =6,  runit, vunit, tunit, initialized along with beta
 
use fft_mod 
use lf_mod ! the system related parameters of lorentz force problem

implicit none

! 	Module lf_mod based Parameter
!  dp, n                ! the dimension of the problem is defined in lf_mod, n=6 
!  runit, vunit, tunit  ! unit of the fundamental quantities


integer, parameter  ::  npoinc =  1000 !20   ! number of Poincare maps, to start test, 11 is enough?
                         
real(kind=dp)       ::  pi = 4.d0*datan(1.d0),  day = 24.d0 * 60.d0 * 60.d0  

! For the case to study 
integer ::        cs, ieq,  debug 
real(kind=dp) ::  beta 

! Local Variables
real(kind=dp) :: tol, tolres, rmm, vmm, dsmm, & !fft
                 p0, tmax, xmax, xi(6), tf, xf(6), tic, toc  ! curve

                  
integer ::  i, j, fdsmax, iscap, fob, ffc, ffcbas, & 
            ind, dir, imax, ncrs  ! poinc 
            
character(len=70) :: fndsmax, fnob, fnfc, fnfcbas   

!	subroutine from packages declared external here
real(kind=dp), external :: dnrm2 ! from BLAS-LEVEL1


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
        
        
! save the data for Fourier analysis of the orbit which has the largest distance from the equilibrium
fob   = 25;    ffc   = 26; ffcbas  = 27
fnob      =  './dat/curve_poinc/ob.dat'
fnfc      = './dat/curve_poinc/fc.dat'
fnfcbas   = './dat/curve_poinc/fcbas.dat'

print*, fnob, fnfc, fnfcbas;  print*   !ck

open(fob  ,file=fnob  , access ='append',status='replace')
write(fob  ,*)  ' # t      (x y z vx vy vz)   cj'

open(ffc  ,    file=fnfc  ,   access ='append',status='replace')
write(ffc  ,*)  ' # freq      ' ! TODO 

open(ffcbas  , file=fnfcbas,  access ='append',status='replace')


open(fdsmax, file=fndsmax) 
read(fdsmax,*) ! first line is comment 

! index of initial point, first one corresponds to a nice tori  
i = 1

! outer bound for escape, nondimensinal unit 1 is big enough....
xmax =  2.d0 
tmax =  1.d2 !a rough guess?? or not...

ind  = 3           ! z component 
p0   = 0.56d0      ! rough guess, by observing directly at the plot of the tori  
dir  = 1           ! take the positive velocity 
imax = 1           ! consider every intersection, and choose by the direction of velocity


do j = 1, i, 1  

  read(fdsmax, '(9e24.14, 1i5)') rmm, vmm, dsmm,  xi, iscap
  write (*,*) xi
  
!   if(iscap == 1) exit 
  ! we need to handle specially for the escape case TODO 
  
!  subroutine poinc_n(xi, ind, p0, dir, fpc, npoinc, imax, tmax, xmax, tf,xf, ncrs, deriv, gr_cj)          
  call poinc_n(xi, ind, p0, dir, fob, npoinc, imax, tmax, xmax, tf,xf, ncrs, gr_lf, gr_cjlf) 
  
enddo   
  
! -- compute the rotation number -- 
!subroutine rotnum( pt, np, rou)

!call  rotnum( pt, np, rou)


stop

end 


 




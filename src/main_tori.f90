program main_tori 

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
                         
real(kind=dp)       :: pi = 4.d0*datan(1.d0) !,  day = 24.d0 * 60.d0 * 60.d0  

integer ::        cs, ieq,  debug 
real(kind=dp) ::  beta 

! Local Variables
real(kind=dp) ::  tf , xmax,  ymax(6),  & ! dsmax 
                 tol, tolres, rmm, vmm, dsmm !fft
                  
integer ::  i, j,  fdsmax, iscap,  & 
            nfmax, nsm, fob, ffc, ffcbas,  nitmax 
            
character(len=70) :: fndsmax,  & ! for all the np points
                     fnob, fnfc, fnfcbas  

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

! case 1:  N=[0 0 1];      2: N=[1 0 0];         3:  N=[0 1 0]; 
 cs = 2  

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
write(*,*) ' Maximum amplitude of the residual //  Tolerance in correction for refinement? (suggestion: 1.d-4, 1.d-6)'
read(*,*)  tolres,  tol 
        
!  tol = 1.d-8; tolres= 1.d-5;  ! too much, fail with 1.d-5... I don't know why
nitmax = 15  
        
! save the data for Fourier analysis of the orbit which has the largest distance from the equilibrium
fob   = 25;    ffc   = 26  ; ffcbas  = 27
fnob      =  './dat/tori/ob.dat'
fnfc      = './dat/tori/fc.dat'
fnfcbas   = './dat/tori/fcbas.dat'

end 




 

program refncur
!  refine the invariant curve \varphi obtained by time T2 map, using as the initial guess the coefficients computed by ffcur.exe  ( main_fourcur.f90 )
!  
!  The approximated functions (x,y,z) are saved in ./dat/curve/fun.dat, probably we need also vx, vy, vz 

!  and the Fourier coefficients, of dimension (0:nf), and save in ./dat/curve/fcs.dat
!  the first line is comment, need to skip 
!  the second line is omega1, omega2 
!  the rest ones are  Fourier coefficientscs:
!    fx(k), sifx(k), csfy(k), sify(k), csfz(k), sifz(k), k = 0, ..., nf 

!  the values from fcs.dat are taken as the initial guess for the refinement of invariant curve \varphi by 'refine_curve'

!   
use fft_mod 
use lf_mod ! the system related parameters of lorentz force problem
use tori_mod 

implicit none

! 	Module lf_mod based Parameter
!  dp, n, 		! the dimension of the problem is defined in lf_mod, n=6 
!  runit, vunit, tunit  ! unit of the fundamental quantities


integer, parameter        :: nf    = 5    ! number of Fourie modes
real(kind=dp), parameter  :: pi    = 4.d0*datan(1.d0)  

! Global Variables
integer       ::  cs, ieq 
real(kind=dp) ::  beta 


! Local Variables
real(kind=dp) ::  omega1, omega2, t,  dt,  cj  ! curve
                 
real(kind=dp), dimension(6, 0:11) ::  CSF6, SIF6 !  approximate fourier coef of  all six compoments
  
! -- refinement -- 
real(kind=dp) ::  h,  rou, t2, c0(6),  ck(6,nf), sk(6,nf), fun6(6)  ! fun_refine
              
integer ::  i,  k,  &      ! do loop counter
            ffcs, ffun, ffcs_refn, ffun_refn,&  ! file tag
            isref ! refinement 
            
            
character(len=70) ::    fnfcs, fnfun, fnfcs_refn, fnfun_refn !  file name 
                      

!	subroutine from packages declared external here
real(kind=dp), external :: dnrm2 ! from BLAS-LEVEL1
double precision fun 


! Initialize the private variables eq1, sgn1 for  case 1
call init_lf_mod

! case 1:  N=[0 0 1];  2: N=[1 0 0]; 3:  N=[0 1 0];  
 cs = 2 ! TODO: use 1 to test the swap rule is coded correctly--ckd!

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

if(cs == 1) then 
! normal case
  if(ieq == 2) beta = 1.d0
 
  
else if(cs == 2) then 
! radial case 

 ! if(ieq==1)  beta = 10.d0 ! TODO 
  if(ieq == 2)  beta = 1.d0
  if(ieq == 3)  beta = 10.d0 
  
! else if(cs == 3) then
! TODO: tangential case -- at this moment hasn't explored

endif  
 
call init_beta(beta)  
 
!  beta = beta0  ! still use this one, to save time...
    
!  print*, 'Please inpute beta:'
!  read(*,*) beta 
     
print*, 'check, cs, ieq,beta', cs, ieq, beta 
read* 
    
    
! Compute the unit of distance, time and velocity, which are all function of beta 
! TODO: use nondimensinal unit, is more useful... 

call lfunit(beta, runit, tunit, vunit)
    
! -- the real units, related to the value of beta -----
print*, 'runit(km), tunit(s), vuint(km/s): ', runit, tunit, vunit
!read*

! Use new name for each exploration, to avoid overwritten the good ones, 
! and the type of approach selected is used as suffix, such ***_refn 
        
! Read refinement of Fourier coefficients + rou + t2 in ./dat/curve/fsc_refn.dat 
ffcs = 111; fnfcs = './dat/curve/fcs.dat'
open(ffcs, file = fnfcs) 

ffcs_refn = 112; fnfcs_refn = './dat/curve/fcs_refn.dat'
open(ffcs_refn, file = fnfcs_refn) 

! The approximated and refined Fourier function  
ffun = 222; fnfun = './dat/curve/fun.dat'
open(ffun,  file = fnfun) 

ffun_refn = 223; fnfun_refn = './dat/curve/fun_refn.dat'
open(ffun_refn,  file = fnfun_refn) 

! TODO: all the six components + energy, now we only have [x,y,z] 
!write(ffun  ,*)  '# Each block for one fitted curve:  t	(x y z vx vy vz)  cj'

! the two dominant frequencies obtained by fft_mod FFT, better to read from fcs.dat, 
! just to keep coherent with the numerical data 
  
!omega1 = 0.233347
!omega2 = 0.021967

read(ffcs, *)  ! skip the first comment line 
read(ffcs, *) omega1, omega2 
print*, 'omega1, omega2', omega1, omega2

! Read the coefficients C0, C's and S's from the file ./curve/fcs.dat 
do k = 0, nf, 1
!  write(ffcs, *)  csf6(k), sif6(k), csfy(k), sify(k), csfz(k), sifz(k)
   read(ffcs, *)  csf6(1, k), sif6(1,k), csf6(2,k), sif6(2,k), csf6(3,k), sif6(3,k), &
                  csf6(4, k), sif6(4,k), csf6(5,k), sif6(5,k), csf6(6,k), sif6(6,k)
enddo
  
 c0 = csf6(:, 0)
 ck = csf6(:, 1:nf)  
 sk = sif6(:, 1:nf)  
 
! ** Note:  t2 doesn't have to be 1./omega2, could be 1./omega1, we use the same notation
!           for both cases 
!    ** Always keep in mind the unit of frequecy we obtained from fftsc is cycle per unit time
!       so we DO NOT need to multiply 2*pi to compute t2 

t2  = 1.d0 / omega1  ! period along omega1 
rou = t2 * omega2    ! rotation number 

! If we want to obtain the curve along frequency omega2, just fixed the time interval to be 
! one period along the other frequency omega1

! h  = 1.d0/omega1 
! np = omega1/omega2  ! the minimum points need for 1 period in omega2.... TODO: not sure 


! Time interval for discretisize the invariant curve with in one period along frequency omega2 
dt = 2.d0*pi/ omega2 /100.d0

!subroutine refine_curve( h, rou, t2, c0, ck, sk, isref )

call gr_cjlf(c0, h)
call refine_curve( h, rou, t2, c0, ck, sk, isref )

do i = 1, 100  
  t = (i-1) * dt
  do k = 1, 6, 1
    fun6(k) = fun( t, nf, csf6(k,:), sif6(k,:), omega2)
  end do

  call gr_cjlf(fun6, cj)
  print*, 'check the energy: h = ', h, 'refined one:', cj 
  read* 
   
  write(ffun_refn, *) t, fun6, cj 
    
enddo 
  
 ! only need to check ..... cool!!!
stop
end 






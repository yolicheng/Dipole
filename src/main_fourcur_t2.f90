program fourcur_t2      
!  compute the Fourier representation of the invariant curve(numerical computed), the data is saved in ./dat/curve/ob_t2.dat

!  which is the obtained by a time t2 flow, with t2 = 1./w1 (w2), if we want to have a curve along w1, just do sections after 1 period of the other frequencty 

! Then we use the general routine 'GR_FOUN' that,  given the frequecies and number of Fourier modes nf, 
! computes the Fourier coefficients, of dimension (0:nf), and save in ./dat/curve/fcs_t2.dat
!  csfx(k), sifx(k), csfy(k), sify(k), csfz(k), sifz(k), k = 0, ..., nf 

! The approximated functions (x,y,z) are saved in ./dat/curve/fun_t2.dat
  

use fft_mod 
use lf_mod ! the system related parameters of lorentz force problem

implicit none

! 	Module lf_mod based Parameter
!  dp, n, 		      ! the dimension of the problem is defined in lf_mod, n=6 
!  runit, vunit, tunit  ! unit of the fundamental quantities


integer, parameter        :: nf = 5   ! number of  Fourier modes
                             np = 11  ! number of points for one revolution to do gr_foun for the initial guess of C's and S's
                             
real(kind=dp), parameter  :: pi = 4.d0*datan(1.d0)  

! Global Variables
integer       ::  cs, ieq 
real(kind=dp) ::  beta 


! Local Variables
! we take the dimension to be 11, because w1/w2 ~= 11, the least points to finish 1 revolution 
real(kind=dp) ::  w1, w2, t, ymax(6), cj,  pv(6, 11),  &  ! curve
                  h, dt, fun6(6), pv_in(6) ! fun
                 
real(kind=dp), dimension(6, 0:11) ::  CSF6, SIF6 !  fourier coef of  all six compoments
real(kind=dp), dimension(0:11)    ::  CSF, SIF  ! fourier  coef of 1 component 
  
integer ::  i, j, k, l, j_st,  &      ! do loop counter
            fpoincdt, ffcs, ffun  ! file tag
            
character(len=70) ::   fnpoincdt, fnfcs, fnfun !  file name 
                      

!	subroutine from packages declared external here
real(kind=dp), external :: dnrm2 ! from BLAS-LEVEL1
double precision fun 


! Initialize the private variables eq1, sgn1 for  case 1
call init_lf_mod

! case 1:  N=[0 0 1];   2: N=[1 0 0];   3:  N=[0 1 0];  
 cs = 2 !  radial case 

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

! the filename to save the useful data, the integration, relative position+velocity, angular momentum
! use new name, to avoid overwritten the good ones, and the type of approach selected is used as suffix
fpoincdt = 21;   fnpoincdt    =  './dat/curve/ob_t2.dat'

!call system("head -n 10 ./dat/curve/ob.dat > tmp.dat" )

! the two dominant frequencies obtained by fft_mod FFT with the initial torus ./dat/tori/ob.dat
w1 = 0.233347
w2 = 0.021967


! If we want to obtain the curve along frequency w2, just fixed the time interval to be 1 period along the other frequency w1
!  first, we try the smaller period with a bigger frequecy 
!h = 1.d0/w2
h = 1.d0/w1 

!np = w1/w2  ! the minimum points need for 1 period in w2.... TODO: not sure...

! the coefficients computed by general Fourier analysis:  fourier.f + fun.f 
ffcs = 111; fnfcs = './dat/curve/fcs_t2.dat'
open(ffcs, file= fnfcs) 
write(ffcs,*) '# The first line: omega1, omega2. The rest: ck(i), sk(i), i=0,...,nf' 
write(ffcs, *) w1, w2


! the save the approximated Fourier function, just used as reference
ffun = 222; fnfun = './dat/curve/fun.dat'
open(ffun,  file = fnfun) 
write(ffun, *) ' # t    (x,y,z, vx,vy,vz)   cj'

! Time interval for discretisize the invariant curve with in one period along frequency w2 
! NOTE: be careful if we should multiply 2pi or not 
!dt = 2*pi/w2 
dt = 1.d0 / w2 

do i = 1, 5 

! open the file that contains the poincare section in time, and skip the first commented line
  open(fpoincdt, file=fnpoincdt) 
  read(fpoincdt,*) ! first line is comment 

  j_st = 1+ 5*(i-1)
  
  ! every time we take 11 points, with 5 overlaped with the previous set
  do j = 1, j_st+10, 1  

    read(fpoincdt, '(8e24.14)')  t, ymax, cj
!    write (*,*) ymax
  
    if( j >= j_st)  then
      pv(:, j-j_st+1) = ymax  ! all the 6 components
    endif 
  
  enddo  

  ! general Fourier analysis for nf Fourier modes, with 11 sample points 
  ! deal with the six components one by one 
  do k = 1, 6, 1
    pv_in =  pv(k,:) ! one component 
    call gr_foun( pv_in, 11, nf, csf, sif  )
    csf6(k, :) = csf
    sif6(k, :) = sif
  end do
  
!  call  GR_FOUN(x, 11, nf, CSFx,SIFx) 

  
  ! write the coefficients C's and S's into file fcs.dat  cfx-sfx -- cfy-sfy -- ......
  do k = 0, nf, 1
!    write(ffcs, *)  csfx(k), sifx(k), csfy(k), sify(k), csfz(k), sifz(k)
     write(ffcs, *)  csf6(1, k), sif6(1,k), csf6(2,k), sif6(2,k), csf6(3,k), sif6(3,k), &
                     csf6(4, k), sif6(4,k), csf6(5,k), sif6(5,k), csf6(6,k), sif6(6,k)
  enddo

  do l = 1, 100  
    t = (l-1) * dt
!    funx = fun( t, nf, csfx, sifx, w2)
!    funy = fun( t, nf, csfy, sify, w2)   
!    funz = fun( t, nf, csfz, sifz, w2)
!    write(ffun, *) t, funx, funy, funz
    do k = 1, 6, 1
      fun6(k) = fun( t, nf, csf6(k,:), sif6(k,:),  w2)
    end do
    call gr_cjlf(fun6, cj)
    write(ffun, *) t, fun6, cj 
    
  enddo 
  
  write(ffcs, *)
  write(ffun, *) 
    
  close(fpoincdt) ! in order to read different lines(not consecutively, but overlap 5 line each time with the previous one), we have to close and reopen for each loop 
enddo 

 
stop
end 






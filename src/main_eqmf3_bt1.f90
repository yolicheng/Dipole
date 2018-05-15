program main_eqmf 

! 20160216  --- epsl= 1.d-6 for both eigenvectors, then the dominant one moves much faster, which can not produce a surface(will move along the dominant eigenvectors) 
! save in main_eqmf_dsmf.90 

! try different distance for different eigenvectors , 1.d-4 for dominant, 1.d-3 for smaller one-- the first one is the dominant one


! 20160214  
! compute the stable(with eigenvalues of positive real part) and stable(respective with eigenvalues of negative real part)manifold of the equilibrium point 

!  consider the 3rd equilibrium point, beta = 1,  we have 2 pairs of real, 1 pair of pure imaginary eigenvalues 
!  By Alex, we choose the initial state by adding a perturbation on the equilibrium along the combination of the eigenvectors we are interested in.
!   stable (corresponding to eigenvalue with positive real part) and unstable (corresponding to eigenvalue with negative real part)
!   for 2d stable(unstable) manifold, dx0 = a*v1 + b*v2, where v1 and v2 are normalized eigenvectors, a=cos(t) and b=sin(t), we test t=[0:2*pi]

use lf_mod

implicit none

integer, parameter ::  neq = 42  
real(kind=dp), parameter :: pi = 4*datan(1.d0) 

! specify the case, equilibrium point and parameter beta that we are going to study
integer ::  cs0, ieq 
real(kind=dp) ::  beta0,  dlf(n,n)

! local variables
integer :: i, j, imf, imax, ncrs, feqmf, feqpoinc, feqcj, &
	   iv1, iv2, it, nt, stb  ! for eqmf

real(kind=dp) ::  wr(n),wi(n), vr(n,n), y0(neq), &   ! stability of equilibrium points
                  x0(6), dx0(6), epsl, cj_tar, & !manifold of eq 
                  xf(6), tf, xpc(6), tmax, tfpl, & !  poinc_n of eq
                  cj_eq, cj, cj2, aaa, & ! debug auxillary variables
                  a, b, dsa, dsb, v1(6), v2(6), dt  !eqmf

 character(len=70) :: fnmf      

real(kind=dp), external :: dnrm2 ! from package NAG Fortran Library Routine

! use this to initialize all the equilibrium points and the respective sign of (q/m), 
! which are the private variable in module lf_mod
call init_lf_mod

! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy
beta0 = 1.d0

! the branch of manifold  -- compute only one branch once 
stb = 1 !unstable  iv1 = 1; iv2 = 4; 
!stb = -1! stable  iv1 = 2; iv2 = 3;

! case 1: N=[1 0 0]; 2:  N=[0 1 0]; 3:  N=[0 0 1]; 
 cs0 = 1 ! use 2 to test the swap rule is coded correctly--ckd!

! 3:  q/m < 0   x = \pm (1/4/sqrt(2)^(1/3), y = 0,z² = x²  ! this the the case we study now 
ieq = 3 ! 3, x,0,z is the case that we study currently

! beta, cs, sgn, eq  are initialized by subroutine init_lf in  module lf_mod
! provided the case, beta, and ieq 
call init_lf(beta0, cs0, ieq) 

print*, 'check the assignment with module' !--ckd
print*, 'beta, cs, sgn, eq', beta, cs, sgn, eq

! the differential matrix of the lorentz force with respective to state 
!subroutine dflrtz( x0, dlf)
call dflrtz( eq, dlf)

! check the energy 
!C1 = -3   
!C2 = 2.2894  
!C3 = 1.8899 

call gr_cjlf(eq, cj_eq)
print*,'check energy of equilibrium point!, cj, ieq,', cj_eq, ieq, eq
read(*,*)  
read(*,*)

! set the target energy to be   ej_eq
 cj_tar = cj_eq  ! 
 
 
do i = 1, n
 write(*,'(6f8.4)') dlf(i,:) 
enddo

! compute the eigenvalues and eigenvectors 
!subroutine eigrg(a,n,isv, wr,wi,vr)
call eigrg(dlf,n,1, wr,wi,vr)

read(*,*)  
read(*,*)  


! for eq3 bt1 
! for 2d stable(unstable) manifold, dx0 = a*v1 + b*v2, where v1 and v2 are normalized eigenvectors, a=cos(t) and b=sin(t), we test t=[0:2*pi] 
!   2.62631284  -2.62631284  -0.94880137   0.94880137 (    0.00000000,    2.40785033) (    0.00000000,   -2.40785033)
! eigenvectors

!  -0.3292  -0.3292   0.0691  -0.0691 (  0.0000,  0.0658) (  0.0000, -0.0658)
!  -0.1155   0.1155  -0.7077  -0.7077 (  0.1467,  0.0000) (  0.1467, -0.0000)
!  -0.0699  -0.0699   0.1435  -0.1435 (  0.0000, -0.3482) (  0.0000,  0.3482)
!  -0.8647   0.8647  -0.0656  -0.0656 ( -0.1585,  0.0000) ( -0.1585, -0.0000)
!  -0.3033  -0.3033   0.6715  -0.6715 (  0.0000,  0.3532) (  0.0000, -0.3532)
!  -0.1837   0.1837  -0.1362  -0.1362 (  0.8385,  0.0000) (  0.8385, -0.0000)
 
! For stable: we use  1-4 eigenvectors. Unstable, we use 2-3 eigenvectors


!  ***************** compute the stable and unstable manifold ***************************
feqmf = 20; 

! only works for ieq = 3 ... because the eigenvectors are specified
! remember to rename the data file for different families of po when beta = 10 
if(stb == 1) then 
  iv1 = 1; iv2 = 4; 
  write(fnmf, fmt='(a,i0,a)') './dat/eqmf3u_bt', idint(beta), '.dat' 
else 
  iv1 = 2; iv2 = 3;
  write(fnmf, fmt='(a,i0,a)') './dat/eqmf3s_bt', idint(beta), '.dat' 
endif 

print*, fnmf 
read*
open(feqmf,file=fnmf,access ='append',status='replace') 


dsa = 1.d-7; dsb = 1.d-9

tf = 1.d2 ! integration time for manifold 

nt = 20 ! number of steps for initial point to compute manifold 
dt = pi/(nt-1) ! step size

v1 = vr(:,iv1) / dnrm2(6, vr(:,iv1), 1 )  ! unit vector
v2 = vr(:,iv2) / dnrm2(6, vr(:,iv2), 1 )  ! unit vector

print*, 'check two vectors we used'
print*, vr(:,iv1), vr(:,iv2), dnrm2(6, v1, 1), dnrm2(6, v2, 1)
read*
  
do it = 1, nt
  a = dsin(dt*(it-1))
  b = dcos(dt*(it-1))
  
  dx0 = a * dsa * v1 + b * dsb  * v2

  print*,'dx0m, a,b, dx', dnrm2(6, dx0, 1), a, b, dx0
!  read(*,*)
    
  x0 =  eq + dx0  ! the linear combination of the 2 eigenvectors
  
  call gr_cjlf(x0, cj)
  print*, 'x0, cj', x0, cj
  
!  subroutine plob(y0,t0,tf,stb,ftag, deriv, gr_cjlf, y) 
  call plob(x0,0.d0, tf, stb, feqmf, gr_lf, gr_cjlf, xf) 
  write(feqmf,*); write(feqmf,*) !add two blank lines
  
  call gr_cjlf(xf, cj2)
  
  write(*, *) 'dcj, cj, cj2', cj2-cj, cj, cj2  
!  read(*,*)  
  
enddo  

 close(feqcj)
stop
end program main_eqmf



















  

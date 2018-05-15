program main_eqmf 
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
	   iv1, iv2, it, nt ! for eqmf

real(kind=dp) ::  wr(n),wi(n), vr(n,n), y0(neq), &   ! stability of equilibrium points
                  x0(6), vr_usd(6,3), dx0(6), epsl, cj_tar, & !manifold of eq 
                  xf(6), tf, xpc(6), tmax, tfpl, & !  poinc_n of eq
                  cj_eq, cj, cj2, aaa, & ! debug auxillary variables
                  a, b, dt, stb   !eqmf
                  
external :: gr_lf ! the differential of the vector field and the variational matrix


! use this to initialize all the equilibrium points and the respective sign of (q/m), 
! which are the private variable in module lf_mod
call init_lf_mod


! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy
beta0 = 1.d0

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
!subroutine dflrtz( x0, beta, cs, sgn, dlf)
call dflrtz( eq, beta, cs, sgn, dlf)

! check the energy 
!C1 = -3   
!C2 = 2.2894  
!C3 = 1.8899 

call gr_cjlf(eq, cj_eq)
print*,'check energy of equilibrium point!, cj, ieq,', cj_eq, ieq, eq
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
! save the data of  eqmf
!feqpoinc = 21; feqcj = 22 
!open(feqpoinc,file='./dat/eqpoinc.dat',access ='append',status='replace')  
!open(feqcj,file='./dat/eqcj3s.dat',access ='append',status='replace')  


feqmf = 20; 
!open(feqmf,file='./dat/eqmf3u_bt1.dat',access ='append',status='replace') 
!iv1 = 1; iv2 = 4; stb = 1 !unstable

open(feqmf,file='./dat/eqmf3s_bt1.dat',access ='append',status='replace') 
iv1 = 2; iv2 = 3; stb = -1! stable 

epsl   = 1.d-6 ! the magnitude of the variation 
tf = 1.d2 ! integration time for manifold 

nt = 10 ! test 3 
dt = pi/(nt-1) ! step size


do it = 1, nt
  a = dsin(dt*(it-1))
  b = dcos(dt*(it-1))
  
  
  dx0 = a*vr(:,iv1) + b*vr(:,iv2) 
  print*,'check unit vector,v1,v2', norm2(vr(:,iv1)), norm2(vr(:,iv2))
  print*,'dx,dx0m, a,b', dx0, norm2(dx0), a, b
!  read(*,*)
    
  x0 =  eq +  epsl * dx0/norm2(dx0) ! the ifam-th column of vr corresponds to the first eigenvalue
  
  call gr_cjlf(x0, cj)
  print*, 'x0, cj', x0, cj
  
!  subroutine plob(y0,t0,tf,stb,ftag, deriv, gr_cjlf, y) 
  call plob(x0,0.d0, tf, stb, feqmf, gr_lf, gr_cjlf, xf) 
  write(feqmf,*); write(feqmf,*) !add two blank lines
  
  call gr_cjlf(xf, cj2)
  
  write(*, *) tf, cj_eq, cj, cj2, cj-cj_eq, cj2-cj
!  read(*,*)  
  
enddo  

 close(feqcj)
stop
end program main_eqmf



















  

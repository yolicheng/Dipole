program main_eqmf  ! eq3, bt2, 4 complex eigenvalues, 2 pure imaginary ones
! 20160214  
! compute the stable(with eigenvalues of positive real part) and stable(respective with eigenvalues of negative real part)manifold of the equilibrium point 
 
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
integer :: i, j, feqmf(2), iv(2), stb(2), dir(2)  ! for eqmf

real(kind=dp) ::  wr(n),wi(n), vr(n,n), y0(neq), &   ! stability of equilibrium points
                  x0(6), vr_usd(6,3), dx0(6), epsl, cj_tar, & !manifold of eq 
                  xf(6), tf, xpc(6), tmax, tfpl, & !  poinc_n of eq
                  cj_eq, cj, cj2   ! debug auxillary variables
                
! use this to initialize all the equilibrium points and the respective sign of (q/m), 
! which are the private variable in module lf_mod
call init_lf_mod


! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy
beta0 = 2.d0

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
!subroutine dflrtz( x0,  dlf)
call dflrtz( eq,  dlf)

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

! for eq3 bt2
! (    0.00000000,    4.12389242) (    0.00000000,   -4.12389242) (   -0.98948989,    0.68981583) (   -0.98948989,   -0.68981583) (    0.98948989,    0.68981583) (    0.98948989,   -0.68981583)
! eigenvectors

! (  0.1177,  0.0000) (  0.1177, -0.0000) (  0.2622, -0.0645) (  0.2622,  0.0645) (  0.2622,  0.0645) (  0.2622, -0.0645)
! (  0.0000, -0.1637) (  0.0000,  0.1637) ( -0.4031, -0.2810) ( -0.4031,  0.2810) (  0.4031, -0.2810) (  0.4031,  0.2810)
! ( -0.1221,  0.0000) ( -0.1221, -0.0000) (  0.2986,  0.0621) (  0.2986, -0.0621) (  0.2986, -0.0621) (  0.2986,  0.0621)
! (  0.0000,  0.4852) (  0.0000, -0.4852) ( -0.2150,  0.2447) ( -0.2150, -0.2447) (  0.2150,  0.2447) (  0.2150, -0.2447)
! (  0.6750,  0.0000) (  0.6750, -0.0000) (  0.5927,  0.0000) (  0.5927, -0.0000) (  0.5927,  0.0000) (  0.5927, -0.0000)
! (  0.0000, -0.5034) (  0.0000,  0.5034) ( -0.3383,  0.1445) ( -0.3383, -0.1445) (  0.3383,  0.1445) (  0.3383, -0.1445)

! we only have 1d stable(unstable manifold), because the focuse have quadrupole(\pm a + \pm b*i type eigenvalues)
 
! For stable: we use  3 eigenvectors. Unstable, we use 5 eigenvectors, checked! 


!  ***************** compute the stable and unstable manifold ***************************
! save the data of  eqmf
feqmf(1)  = 20; feqmf(2) = 21
open(feqmf(1),file='./dat/eqmf3u_bt2.dat',access ='append',status='replace') 
open(feqmf(2),file='./dat/eqmf3s_bt2.dat',access ='append',status='replace') 

iv(1) = 5;  stb(1) = 1 !unstable
iv(2) = 3;  stb(2) = -1! stable 

print*, 'two vector2,3,5', vr(:,iv(1)), vr(:,iv(2))
print*,'check unit vector,v1,v2', norm2(vr(:,iv(1))), norm2(vr(:,iv(2)))

epsl   = 1.d-6 ! the magnitude of the variation 
tf = 1.d2! integration time for manifold 

dir = (/1,-1/)

do i = 1, 2 ! compute stable and unstable in a loop, save in two index 
  dx0 = vr(:,iv(i))
  
  do j = 1, 2 ! we need to consider both negative direction and positive direction on the eigenvector 
  
    x0 =  eq +  dir(j) * epsl * dx0/norm2(dx0) ! the ifam-th column of vr corresponds to the first eigenvalue
    print*,'dx', dx0 
    read(*,*)
  
    call gr_cjlf(x0, cj)
    print*, 'x0, cj', x0, cj
  
    print*, 'ck i, stb, feqmf', i, stb(i), feqmf(i)
!  subroutine plob(y0,t0,tf,stb,ftag, deriv, gr_cjlf, y) 
    call plob(x0,0.d0, tf, stb(i), feqmf(i), gr_lf, gr_cjlf, xf) 
  
    write(feqmf(i),*); write(feqmf(i),*) !add two blank lines
  
    call gr_cjlf(xf, cj2)
  
    write(*, *) tf, cj_eq, cj, cj2, cj-cj_eq, cj2-cj
    read(*,*)  
  enddo ! for two branches
enddo  


 close(feqmf(1)); close(feqmf(2))
stop
end program main_eqmf



















  

program main_eqmf  !! eq2 bt 1 :   2 real eigenvalues, 4 pure imaginary ones

! 20160308   -- remember to change to filename of the data file, not to overwrite the old ones...

! compute the stable(with eigenvalues of positive real part) and stable(respective with eigenvalues of negative real part)manifold of the equilibrium point 
 
 
use lf_mod

implicit none

integer, parameter ::  neq = 42  
real(kind=dp), parameter :: pi = 4*datan(1.d0) 

! specify the case, equilibrium point and parameter beta that we are going to study
integer ::  cs0, ieq 
real(kind=dp) ::  beta0,  dlf(n,n)

! local variables
integer :: i,j, feqmf(2), iv(2), stb(2), dir(2)  ! for eqmf

real(kind=dp) ::  wr(n),wi(n), vr(n,n), y0(neq), &   ! stability of equilibrium points
                  x0(6), vr_usd(6,3), dx0(6), epsl, cj_tar, & !manifold of eq 
                  xf(6), tf, xpc(6), tmax, tfpl, & !  poinc_n of eq
                  cj_eq, cj, cj2   ! debug auxillary variables
                
! use this to initialize all the equilibrium points and the respective sign of (q/m), 
! which are the private variable in module lf_mod
call init_lf_mod


! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy
beta0 = 1.d0

! save the data of  eqmf
feqmf(1)  = 20; feqmf(2) = 21
open(feqmf(1),file='./dat/eqmf2u_bt1.dat',access ='append',status='replace') 
open(feqmf(2),file='./dat/eqmf2s_bt1.dat',access ='append',status='replace') 


! case 1: N=[1 0 0]; 2:  N=[0 1 0]; 3:  N=[0 0 1]; 
 cs0 = 1 ! use 2 to test the swap rule is coded correctly--ckd!

! 3:  q/m < 0   x = \pm (1/4/sqrt(2)^(1/3), y = 0,z² = x²  ! this the the case we study now 
ieq = 2 ! 3, x,0,z is the case that we study currently

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
 
! Eigenvalues
!    0.00000000      2.69451269    0.00000000     -2.69451269    1.41421356    -1.41421356     0.00000000      1.11337386    0.00000000     -1.11337386
! Eigenvectors
!  0.1730 -0.0908    0.1730  0.0908    0.4082   0.2887  -0.0609  0.0000   -0.0609 -0.0000  
!  0.0000  0.2260    0.0000 -0.2260    0.0000   0.4082  -0.2264  0.1782   -0.2264 -0.1782  
!  0.0829  0.1579    0.0829 -0.1579   -0.4082   0.2887   0.0000 -0.5998    0.0000  0.5998  
!  0.2446  0.4661    0.2446 -0.4661    0.5774  -0.4082   0.0000 -0.0678    0.0000  0.0678  
! -0.6090  0.0000   -0.6090 -0.0000    0.0000  -0.5774  -0.1984 -0.2520   -0.1984  0.2520  
! -0.4255  0.2233   -0.4255 -0.2233   -0.5774  -0.4082   0.6679  0.0000    0.6679 -0.0000  
 
! we only have 1d stable(unstable manifold), 4d center manifold
 
! For stable: we use 4 -th eigenvectors. Unstable, we use 3 -th eigenvectors, checked! 


!  ***************** compute the stable and unstable manifold ***************************

iv(1) = 3;  stb(1) = 1  !unstable
iv(2) = 4;  stb(2) = -1 ! stable 

print*, 'two vectors', iv, vr(:,iv(1)), vr(:,iv(2))
print*,'check unit vector,v1,v2', norm2(vr(:,iv(1))), norm2(vr(:,iv(2)))

epsl   = 1.d-7 ! the magnitude of the variation 
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



















  

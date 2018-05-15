program main_eqmf 
! 20160118  
!  consider the third equilibrium point
!  take beta = 10, we have three pairs of pure imaginary eigenvalues

use lf_mod

implicit none

integer, parameter ::  neq = 42  

! specify the case, equilibrium point and parameter beta that we are going to study
integer ::  cs0, ieq 
real(kind=dp) ::  beta0,  dlf(n,n)

! local variables
integer :: i, j, imf, imax, ncrs, feqmf, feqpoinc, feqcj

real(kind=dp) ::  wr(n),wi(n), vr(n,n), y0(neq), &   ! stability of equilibrium points
                  x0(6), vr_usd(6,3), dx0(6), epsl, cj_tar, & !manifold of eq 
                  xf(6), tf, xpc(6), tmax, tfpl, & !  poinc_n of eq
                  cj_eq, cj, cj2, aaa ! debug auxillary variables
                  
external :: gr_lf ! the differential of the vector field and the variational matrix


! use this to initialize all the equilibrium points and the respective sign of (q/m), 
! which are the private variable in module lf_mod
call init_lf_mod


! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy
beta0 = 10.d0

! case 1: N=[1 0 0]; 2:  N=[0 1 0]; 3:  N=[0 0 1]; 
 cs0 = 1 ! use 2 to test the swap rule is coded correctly--ckd!

! the index of the equilibrium points
! 1:  q/m > 0,  x=0,y=0,z= \pm 1
! 2:  q/m < 0,  x = \pm (2 sqrt(9))^(1/3), y² = 2x²,z = 0
! 3:  q/m < 0   x = \pm (1/4/sqrt(2)^(1/3), y = 0,z² = x²  ! this the the case we study now 
  
ieq = 3 ! 3, x,0,z is the case that we study currently

! the state of the equilibrium point, positon+velocity
!eqst = (/1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0/)! test

! beta, cs, sgn, eq  are initialized by subroutine init_lf in  module lf_mod
! provided the case, beta, and ieq 

! subroutine init_lf(beta0, cs0, ieq)
call init_lf(beta0, cs0, ieq) 

print*, 'check the assignment with module' !--ckd
!print*, 'beta, cs, sgn, eq', beta, cs, sgn, eq

! the differential matrix of the lorentz force with respective to state 
!subroutine dflrtz( x0, beta, cs, sgn, dlf)
call dflrtz( eq, beta, cs, sgn, dlf)

! check the energy 
!C1 = -3   
!C2 = 2.2894  
!C3 = 1.8899 

call gr_cjlf(eq, cj_eq)
print*,'check energy of equilibrium point!, cj, ieq,', cj_eq, ieq, eq
read(*,*) aaa

! set the target energy to be a little (1.d-4) less than ej_eq
 cj_tar = cj_eq - 1.d-4
 
 
do i = 1, n
 write(*,'(6f8.4)') dlf(i,:) 
enddo

! compute the eigenvalues and eigenvectors 
!subroutine eigrg(a,n,isv, wr,wi,vr)
call eigrg(dlf,n,1, wr,wi,vr)

read(*,*) aaa
! Add a small perturbation on the equilibrium point,  using poincare section to see 
! is the motion around the equilibrium point is stable or not. 
! how to decide the perturbation?? 

! play with the pertubation to observe the behavior of the dynamic around the equilibrium point 

! 1. should be small enough, with magnitude of 1.d-3 
! 2. first along one center eigenspace, the first family
!    then 1*1 + 1*2
!    1*1 + 1*3
!    1*1 + 1*2 + 1*3
!    1*2 + 1*3

! poincare section: y=0 plane, if we fix the energy, we will have 4 dimension  manifold
! it is difficult to display, 
 
!  the three eigenvectors that belong to three families of po

! 1, 4, 6 -th column of the eigenvector matrix, corresponds to 1:real part, 2,3, imaginary part
vr_usd(:,1) = vr(:,1)
vr_usd(:,2) = vr(:,4)
vr_usd(:,3) = vr(:,6)

!vr_usd = vr(:, (/1, 4, 6/)
 
!  ***************** compute the 3 families of PO ***************************
! save the data of po and  its initial state 
feqmf = 20;  feqpoinc = 21; feqcj = 22
open(feqmf,file='./dat/eqmf.dat',access ='append',status='replace')  
open(feqpoinc,file='./dat/eqpoinc.dat',access ='append',status='replace')  
open(feqcj,file='./dat/eqcj.dat',access ='append',status='replace')  

epsl   = 1.d-2 ! the magnitude of the variation 

!!number of crossing needed to be complete
imax =  500
tmax = 1.d2

!do imf = 7,7 
!do imf = 6,6 
!do imf = 5,5 
!do imf = 4,4 
!do imf = 3,3
!do imf = 2,2
!do imf = 1,1

do imf = 1, 16
! this is the orbit with different initial condition to be tested, the perturbation:
!    1
!    2
!    3
!    1*1 + 1*2
!    2*1 + 1*2 
!    1*1 + 2*2 

!    1*1 + 1*3
!    2*1 + 1*3 
!    1*1 + 2*3 

!    1*2 + 1*3
!    2*2 + 1*3
!    1*2 + 2*3

!    1*1 + 1*2 + 1*3
!    2*1 + 1*2 + 1*3
!    1*1 + 2*2 + 1*3
!    1*1 + 1*2 + 2*3

  if(imf .le. 3) then 
    dx0 =  vr_usd(:, imf)
    
  elseif(imf == 4) then    ! 1*1 + 1*2
    dx0 =  vr_usd(:, 1) + vr_usd(:,2) 
    
  elseif(imf == 5) then   !    2*1 + 1*2 
    dx0 =  2*vr_usd(:, 1) + vr_usd(:,2)  
  
  elseif(imf == 6) then  !    1*1 + 2*2
    dx0 =  vr_usd(:, 1) + 2*vr_usd(:,2) 

!    1*1 + 1*3
!    2*1 + 1*3 
!    1*1 + 2*3 

  elseif(imf == 7) then   
    dx0 =  vr_usd(:, 1) +  vr_usd(:,3) 

  elseif(imf == 8) then   
    dx0 = 2*vr_usd(:, 1) + vr_usd(:,3) 

  elseif(imf == 9) then   
    dx0 =  vr_usd(:, 1) + 2*vr_usd(:,3) 

!    1*2 + 1*3
!    2*2 + 1*3
!    1*2 + 2*3

  elseif(imf == 10) then   
    dx0 =  vr_usd(:, 2) + vr_usd(:,3) 

  elseif(imf == 11) then  
    dx0 =  2*vr_usd(:, 2) + vr_usd(:,3) 

  elseif(imf == 12) then  
    dx0 =  vr_usd(:, 2) + 2*vr_usd(:,3) 


!    1*1 + 1*2 + 1*3
!    2*1 + 1*2 + 1*3
!    1*1 + 2*2 + 1*3
!    1*1 + 1*2 + 2*3
 
    
  elseif(imf == 13) then   
    dx0 =  vr_usd(:, 1) + vr_usd(:,2)  + vr_usd(:,3)
  
  elseif(imf == 14) then   
    dx0 =  2*vr_usd(:, 1) + vr_usd(:,2)  + vr_usd(:,3)

  elseif(imf == 15) then   
    dx0 =  vr_usd(:, 1) + 2*vr_usd(:,2)  + vr_usd(:,3)

  elseif(imf == 16) then   
    dx0 =  vr_usd(:, 1) + vr_usd(:,2)  + 2*vr_usd(:,3)
  
        
  endif
  

  x0 =  eq +  epsl * dx0/norm2(dx0) ! the ifam-th column of vr corresponds to the first eigenvalue
   
  ! use the poincare section y=0  to study the poincare map of the manifold 
  ! pay attention that the energy level is fixed here, how to realize this? 
  ! for the given initial state, which are all of form (x,0,z, 0,vy,0), fix x and z, modify vy

! check the energy here, and then decide which energy level to take  
  call gr_cjlf(x0, cj)
  print*, imf,'-th perturbation, x0, cj'
  print*, x0, cj
  
!  print*, 'x0,cj, cj_eq', x0, cj,  cj_eq
  
  ! modify vz to select the same energy level
  ! check if we use cj_eq, the energy of the equilibrium point, it will almost stay on x-z plane
 
!  call cjlf_vy(x0, cj_tar)
!!  call gr_cjlf(x0, cj2)
!!  print*,'new state and energy', x0, cj2, cj_eq
!!  read(*,*) aaa
!  

!  print*, x0, cj_tar
!  read(*,*) aaa
  
  
!  tf = 14.83839473d0 + 1.d0
!  call plob(x0,0.d0, tf, feqmf, gr_lf, xf)
! the stop condition: integration time or number of crossing through y=0 plane?
!                    choose the crossing, imax = 1000, to have enough point
!                    and use the stop time to integrate the orbit    

!                     at last, choose the maximum time as the stop conditon,
!                     because some may not cross the y=0 plane due to the combination of initial condition 

! subroutine poinc_n(y0, tmax, fpc, yf, tf, ncrs)    
  call poinc_n(x0, imax, tmax, feqpoinc, tf, xpc, ncrs)  
!  print*, ncrs,' crossing of poincare section'
  
  write(feqpoinc,*); write(feqpoinc,*) !add two blank lines

!  subroutine plob(y0,t0,tf,ftag, deriv, y) 
  call plob(x0,0.d0, tf, feqmf, gr_lf, xf) 
  write(feqmf,*); write(feqmf,*) !add two blank lines
  
! check if xf and xpc is the same one? 
  print*, 'check the final state in poinc_n'
  print*, xpc, xf
  
  call gr_cjlf(xf, cj2)
  
!  write(feqcj,'(3f12.8)') tf, cj2-cj, cj, cj2
  write(feqcj, *) tf,cj, cj2 , cj2-cj
  
  write(*, *) tf,cj, cj2 , cj2-cj
!  read(*,*) aaa
  
enddo  

 close(feqcj)
stop
end program main_eqmf



















  

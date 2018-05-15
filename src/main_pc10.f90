program main_pc10 !eq = 3, beta = 10 
! 20160310  


! 20160127
!   for each energy level, take a point that escape, which explain the domain in which the orbit is stable. 
!   and will also be used for the fourier analysis

! 20160118  
!  consider the third equilibrium point
!  take beta = 10, we have three pairs of pure imaginary eigenvalues, do the poincare section as Gerard explained

!  With the equilibrium point, as (x0,y0,z0, xdot=0, ydot=0, zdot =0 )
!  1. fix the energy level H = H0 (the equilibrium point)
!  2. fix also z = z0, this is chosen to be the poincare section 
!  3. For this energy level, select different initial point nearby the equilibrium point 
!      5 * dx   by 5 * dy 
!  4. integrate for 200 crossings, only take the one with zdot > 0, to avoid 2 intersections in one loop
!  5. plot the poincare map in (x,y), (xdot,ydot), (x,xdot),(y,ydot) - planar 

! 6. for the change of energy level, H = H0 + dH, keep all the initial value, only modify vz


use lf_mod

implicit none

integer, parameter ::  neq = 42  

! specify the case, equilibrium point and parameter beta that we are going to study
integer ::  cs0, ieq 
real(kind=dp) ::  beta0,  dlf(n,n)

! local variables
integer :: i, j, ix, iy, icj, iscj, imf, imax, ncrs, feqmf, feqpc, feqcj

real(kind=dp) ::  wr(n),wi(n), vr(n,n), y0(neq), &   ! stability of equilibrium points
                  x0(6), vr_usd(6,3), dx0(6), epsl, epsl_dx(8), cj_tar, & !manifold of eq 
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
read(*,*) aaa

! set the target energy to be cj_eq
cj_tar = cj_eq  

 
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


!  With the equilibrium point, as (x0,y0,z0, xdot=0, ydot=0, zdot =0 )
!  1. fix the energy level H = H0 (the equilibrium point)
!  2. fix also z = z0, this is chosen to be the poincare section 
!  3. For this energy level, select different initial point nearby the equilibrium point 
!      5 * dx   by 5 * dy 
!  4. integrate for 200 crossings, only take the one with zdot > 0, to avoid 2 intersections in one loop
!  5. plot the poincare map in (x,y), (xdot,ydot), (x,xdot),(y,ydot) - planar 

! 6. enlarge the domain to be explored, change the initial coordinate along dx = dy line

!  ***************** compute the 3 families of PO ***************************
! save the data of po and  its initial state 
feqmf = 20;  feqpc = 21; feqcj = 22
open(feqmf,file='./dat/mf10.dat',access ='append',status='replace')  
open(feqpc,file='./dat/pc10.dat',access ='append',status='replace')  
!open(feqcj,file='./dat/eqcj.dat',access ='append',status='replace')  


!!number of crossing needed to be complete
!imax =  1 ! 10! 800
imax = 800 ! 400
!epsl = 2.d-3
tmax = 1.d2  ! 100 as the maximum integration time for one orbit

! 5.d-2 is the allowable biggest value, 5.5d-2 begin to escape
!epsl_dx = (/5.d-3,  7.5d-3, 1.d-2, 2.5d-2, 5.d-2, 5.5d-2, 6.d-2, 1.d-1/)  ! the second is too close 

!epsl_dx = (/5.d-3, 1.5d-2, 3.d-2, 4.d-2, 5.d-2, 6.d-2, 5.5d-2, 1.d-1/)   !
!epsl_dx = (/5.d-3,  1.d-2,  5.d-2, 1.d-1, 2.d-1, 4.d-1, 8.d-1, 1.6d0/) ! 1.d-1 is too much
 

epsl_dx = (/5.d-3, 1.5d-2, 3.d-2, 4.d-2, 5.d-2, 5.5d-2, 6.d-2,  1.d-1/)   !this is the one used to plot 
 
!epsl_dx = (/5.d-3,  2.d-2,  4.d-2, 6.d-2, 6.1d-2, 1.d-1, 2.d-1, 1.6d0/)  
  
do icj =  0, 0 !-1,0 ! 0,0! -1, 1

! modify vy to select the same energy level
  
  if (icj == -2) then 
    cj_tar = cj_eq - 1.d-1

  elseif (icj == -1) then 
    cj_tar = cj_eq - 1.d-2
        
  elseif (icj == 0) then 
      
    cj_tar = cj_eq
        
  elseif (icj == 1) then 
     
    cj_tar = cj_eq + 2.d-4
  endif   


! increase the domain to be explored, instead of the first five, we take the first 6 to get escape orbit for every energy level 
!epsl_dx = (/5.d-3, 1.5d-2, 3.d-2, 4.d-2, 5.d-2, 5.5d-2, 6.d-2, 1.d-1/)   !  
     
  do ix = 1, 5 ! 3, 4! 8 ! 1, 2  ! 1, 1   !2
      
      x0(1) =  eq(1) +  epsl_dx(ix) ! x 
      x0(2) =  eq(2) +  epsl_dx(ix) ! y 
      x0(3) =  eq(3) ! z
    
  ! choose xdot and ydot randomly, and modify the value of vz to obtain specified energy   
!      x0(4) = 0.d0
      x0(4) = 1.d-4 ! vx
      x0(5) = -1.d-5 ! vy
     
      x0(6) = 0.d0 ! vz
   
      call gr_cjlf(x0, cj)
   
    
 ! check the energy level before changing the value of vz 
      print*, 'before vz, cj, cj_tar, x0', cj, cj_tar, x0
!       read(*,*) aaa

      call lf2cj(x0, 6, cj_tar, iscj) ! 6: vz to be modified
      
      if(iscj ==0) then 
        print*, 'Vz*Vz < 0 in the energy level', cj_tar,'! icj=',icj,' ix=',ix,' iy=',iy
        read(*,*) aaa
        cycle
      endif
      
      call gr_cjlf(x0, cj)
   
    ! check the energy here, and then decide which energy level to take  
      print*, 'initial point x0 on energy level cj_tar', x0, cj, cj_tar
      read(*,*) aaa
    
 
  ! pay attention that the energy level is fixed here, how to realize this? 
  ! for the given initial state, which are all of form (x,0,z, 0,vy,0), fix x and z, modify vy
 
  ! check if we use cj_eq, the energy of the equilibrium point, it will almost stay on x-z plane
 

! the stop condition:  number of crossing through poincare section
!                    set imax = 800, to have enough point,  and  then use the stop time to integrate the orbit    


!subroutine poinc_n(y0, imax, tmax, ind, z0, fpc, tf, yf, ncrs)  
      call poinc_n(x0, imax, tmax, 3, eq(3), feqpc, tf, xpc, ncrs)  
!  print*, ncrs,' crossing of poincare section'
      write(feqpc,*);  ! to seperate block for different orbit

!  subroutine plob(y0,t0,tf,ftag, deriv, y) 
      call plob(x0,0.d0, tf, feqmf, gr_lf, xf)  
  
! check if xf and xpc is the same one? 
      print*, 'check the final state in poinc_n'
      print*, xpc, xf
  
      call gr_cjlf(xf, cj2)
  
!  write(feqcj,'(3f12.8)') tf, cj2-cj, cj, cj2
       write(feqcj, *) ncrs, tf,cj, cj2 , cj2-cj
  
      write(*, *) ncrs, tf,cj, cj2 , cj2-cj
       read(*,*) aaa
  
!    enddo  
  enddo
 
 ! for different energy, add  two blank lines, plot with index
  write(feqmf,*); write(feqmf,*) !  
  write(feqpc,*)    
  write(feqcj, *); write(feqcj, *) 
enddo 
 close(feqcj)
stop
end program main_pc10



















  

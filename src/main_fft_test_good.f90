program main_fft

!  consider the third equilibrium point
!  take beta = 10, we have three pairs of pure imaginary eigenvalues, do the poincare section
!  as Gerard explained

! 20160118
!  With the equilibrium point, as (x0,y0,z0, xdot=0, ydot=0, zdot =0 )
!  1. fix the energy level H = H0 (the equilibrium point)
!  2. fix also z = z0, this is chosen to be the poincare section 
!  3. For this energy level, select different initial point nearby the equilibrium point 
!      5 * dx   by 5 * dy 
!  4. integrate for 200 crossings, only take the one with zdot > 0, to avoid 2 intersections in one loop
!  5. plot the poincare map in (x,y), (xdot,ydot), (x,xdot),(y,ydot) - planar 

! 6. for the change of energy level, H = H0 + dH, keep all the initial value, only modify vz

! 20160121 
!   1. use fixed stepsize integrator, possible with gr_rk78?, test with the same h, hmin, hmax 
!      to generate equispaced sample points 
!   2. use fft from numerical recipes to do the Fast Fourier Transformation
!   3. also start from H0, and then test the other energy

! 20160122 -- test fft 
!   1. test the function in Gerard's book, to see if we get the right magnitude and frequency

!    the function is f(t) = cos(2*pi*0.13*t) - 1/2*sin(2*pi*0.27*t) + 3/4 *sin(2*pi*0.41*t)


use lf_mod
!use, intrinsic :: iso_c_binding
!include 'fft3.f03'

implicit none

! number of sampling points on the orbit: 2**32, 2**24, 2**16 points should be good
integer, parameter ::   nsm = 2**18, nsm2 = 2*nsm, &  !number of sampling points 
                        smpl = 1, np = nsm*smpl,  & !  !  2**10 = 1024, steps of integration  
                        
			neq = 42 ! number of dimensions for monodrony matrix computation 


real(kind=dp), parameter :: pi = 4*datan(1.d0), lt = 2.d2 * 2 * pi !
!, & ! time length of sampling 
!			    fs = np / lt /2/pi  !  sampling frequency: per cycle

! specify the case, equilibrium point and parameter beta that we are going to study
integer ::  cs0, ieq 
real(kind=dp) ::  beta0,  dlf(n,n)

! local variables
integer :: i, j, ix, iy, icj, iscj,  imf, imax, ncrs, feqmf, feqpc, feqcj, &
	   isign, f_fft1, f_fft2, f_fft3, f_fft4, f_fft !, nsm2 !--fft

real(kind=dp) ::  wr(n),wi(n), vr(n,n), y0(neq), &   ! stability of equilibrium points
                  x0(6), vr_usd(6,3), dx0(6), epsl, epsl_dx(8), cj_tar, & !manifold of eq 
                  xf(6), tf, xpc(6), tmax, tfpl, & !  poinc_n of eq
                  cj_eq, cj, cj2, aaa, & ! debug auxillary variables
                  xl(np,6), data1(nsm),data2(nsm),  h, t, & !, h fft  
                  f, amp1, amp2,  amp3, amp4, &! write to file  
                  tst, tend ! evaluate time  
      
DOUBLE COMPLEX fft1(nsm),fft2(nsm), fft3(nsm),fft4(nsm) 
                 
external :: gr_lf ! the differential of the vector field and the variational matrix


! use this to initialize all the equilibrium points and the respective sign of (q/m), 
! which are the private variable in module lf_mod
call init_lf_mod

! the fixed stepsize for the integration of the orbit
h = lt/np ! * 2 * pi

! dt = 2*pi * 1./fs; % time interval  

!h = 8.d-3 
     
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

!do i = 1, n
! write(*,'(6f8.4)') dlf(i,:) 
!enddo

! check the energy 
!C1 = -3   
!C2 = 2.2894  
!C3 = 1.8899 

call gr_cjlf(eq, cj_eq)
print*,'check energy of equilibrium point!, cj, ieq,', cj_eq, ieq, eq
!read(*,*) aaa

! set the target energy to be cj_eq
 cj_tar = cj_eq  


! compute the eigenvalues and eigenvectors 
!subroutine eigrg(a,n,isv, wr,wi,vr)
call eigrg(dlf,n,1, wr,wi,vr)

!read(*,*) aaa

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
open(feqmf,file='./dat/mf_fft.dat',access ='append',status='replace')  
open(feqpc,file='./dat/pc_fft.dat',access ='append',status='replace')  
!open(feqcj,file='./dat/eqcj.dat',access ='append',status='replace')  

f_fft1 = 30; f_fft2 = 31;  f_fft3 = 32;  f_fft4 = 33; f_fft = 35
open(f_fft1, file='./dat/fft1.dat',access ='append',status='replace')  
open(f_fft2, file='./dat/fft2.dat',access ='append',status='replace')  

open(f_fft3, file='./dat/fft3.dat',access ='append',status='replace')  
open(f_fft4, file='./dat/fft4.dat',access ='append',status='replace')  

open(f_fft, file='./dat/fft.dat',access ='append',status='replace')  

! test the same function, and compare with the result of matlab
! test f(t) = cos(2*pi*0.13*t) - 1/2*sin(2*pi*0.27*t) + 3/4 *sin(2*pi*0.41*t)
call cpu_time(tst)    
do i = 1, np ! tested!
  t =  (i-1) * h  !* 2 * pi ! we have 2 cycles
  
  data1(i)  = dcos(0.13*t) - 1./2*dsin(0.27*t) + 3./4. *dsin(0.41*t);
  data2(i)  = dcos(1.3*t) + 1./4*dsin(2.7*t) - 1./2. *dsin(4.1*t);


!  data1(i) = 3.d0*dsin(8.d0*t) 
!  data2(i) = dcos(13.d0*t) - 1.d0/2.d0*dsin(3.d0*t) + 3.d0/4.d0 *dsin(5.d0*t)
end do
 
print*, 'data sampling finished!'
call cpu_time(tend)
print*, 'Elapsed time for a and b sampling is ', tend-tst


call cpu_time(tst)    
call dtwofft(data1,data2,fft1,fft2, nsm) 
    
call cpu_time(tend)
print*, 'Elapsed time for twofft is ', tend-tst


do i = 1, np/2
! if we use real time
   f  = (i-1) / lt  * pi * 2! ! for twofft, amp = amp/(nsm/2),  but if t is the real time, f = f*2*pi  
  
! cycle per second
!   f  = i / lt   ! for twofft, amp = amp/(nsm/2), f = f , but if t is the real time, f = f*2*pi/2 = f*pi
     
  amp1 = cdabs( fft1(i) )/ nsm * 2
  amp2 = cdabs( fft2(i) )/ nsm * 2
 
  if(i == 1) then 
    amp1 = amp1 / 2 
    amp2 = amp2 / 2 
  endif
  
  write(f_fft1, '(2f20.10)')  f, amp1 
  write(f_fft2, '(2f20.10)')  f, amp2 
  
  if (dabs(amp1) > 1.d-2 ) then 
    print*,'f1:', f, amp1
    write(f_fft, *) 'f1:', f, amp1
  endif
  
  if (dabs(amp2) > 1.d-2 ) then 
   print*, 'f2:', f, amp2
   write(f_fft, *) 'f2:', f, amp2
  endif
enddo 


 
!read(*,*)  
! close(f_fft1); close(f_fft2)
! stop



!!number of crossing needed to be complete
!imax =  1 ! 10! 800
imax = 800 ! 400
tmax = 1.d2  ! 100 as the maximum integration time for one orbit



! 5.d-2 is the allowable biggest value, 5.5d-2 begin to escape
!epsl_dx = (/5.d-3,  7.5d-3, 1.d-2, 2.5d-2, 5.d-2, 5.5d-2, 6.d-2, 1.d-1/)  ! the second is too close 

epsl_dx = (/5.d-3, 1.5d-2, 3.d-2, 4.d-2, 5.d-2, 5.5d-2, 6.d-2, 1.d-1/)   !
 
!epsl_dx = (/5.d-3,  1.d-2,  5.d-2, 6.d-2, 2.d-1, 4.d-1, 8.d-1, 1.6d0/)  
 
  
! test H0  
do icj =  0,0  !-2, 1 !-1,0 ! 0,0! -1, 1

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
  

  ! try with the first point   
  do ix = 1, 1 ! 1,5 ! 3, 4! 8 ! 1, 2     !2
      
    x0(1) =  eq(1) +  epsl_dx(ix) ! x 
    x0(2) =  eq(2) +  epsl_dx(ix) ! y 
    x0(3) =  eq(3) ! z
    
! choose xdot and ydot randomly, and modify the value of vz to obtain specified energy   
    x0(4) = 1.d-4 ! vx
    x0(5) = -1.d-5 ! vy
     
    x0(6) = 0.d0 ! vz
   
    call gr_cjlf(x0, cj)
   
    
 ! check the energy level before changing the value of vz 
!    print*, 'before vz, cj, cj_tar, x0', cj, cj_tar, x0
!     read(*,*) aaa

    call lf2cj(x0, 6, cj_tar, iscj) ! 6: vz to be modified
      
    if(iscj ==0) then 
      print*, 'Vz*Vz < 0 in the energy level', cj_tar,'! icj=',icj,' ix=',ix,' iy=',iy
      read(*,*) aaa
      cycle
    endif
      
    call gr_cjlf(x0, cj)
   
    ! check the energy here, and then decide which energy level to take  
!    print*, 'initial point x0 on energy level cj_tar', x0, cj, cj_tar
!    read(*,*) aaa
    
 ! integrate for enough long time? how long? fixed step size
 !  and save the data 
!    tf = 1.d2 
!    call plob(x0,0.d0, tf, feqmf, gr_lf, xf) 
!    print*,' fininsh plob'
!    read(*,*)
    
    call cpu_time(tst)    
 
!subroutine plob_fxd(y0, t0,  h, np, ftag, deriv, yall) 
    call plob_fxd(x0, 0.d0, h/smpl, np, feqmf, gr_lf, xl) 
    
    call cpu_time(tend)
    print*, 'Elapsed time for orbit sampling is ', tend-tst
    
    
!! use the components of the state vector from the integration with fixed stepsize
    print * ,'Integration of orbit finished!'
    
    data1 = xl(1:np:smpl, 1) ! x
    data2 = xl(1:np:smpl, 2) ! y 
 
    call cpu_time(tst)    
    call dtwofft(data1,data2,fft1,fft2, nsm)
    call cpu_time(tend)

    print*, 'Elapsed time for twofft is ', tend-tst
    
do i = 1, np/2
! if we use real time
  f  = (i-1) / lt  * pi * 2 ! ! for twofft, amp = amp/(nsm/2), f = f/2, but if t is the real time, f = f*2*pi/2 = f*pi
  
! cycle per second
!  f  = i / lt  ! ! for twofft, amp = amp/(nsm/2)  
    
  amp1 = cdabs( fft1(i) )/ nsm * 2
  amp2 = cdabs( fft2(i) )/ nsm * 2
  
  if(i == 1) then 
    amp1 = amp1 / 2 
    amp2 = amp2 / 2 
  endif
    
  write(f_fft3, '(2e20.10)')  f, amp1 
  write(f_fft4, '(2e20.10)')  f, amp2 
  
  if (dabs(amp1) > 1.d-3 ) then 
    print*,'f3:', f, amp1
    write(f_fft, *) 'f3:', f, amp1
  endif
  
  if (dabs(amp2) > 1.d-3 ) then 
    print*, 'f4:', f, amp2
    write(f_fft, *) 'f4:', f, amp2
  endif
  
enddo 
    
  
! check the final state and stop time? 
!    print*, 'check the final state'
!    print*, np, xl(np,:)
  
    call gr_cjlf(xl(np,:), cj2)
  
    write(*, *)  cj, cj2 , cj2-cj
    read(*,*) aaa
  enddo
  
!  write(feqmf,*); write(feqmf,*) !add two blank lines
  
enddo 


 close(feqcj)
stop
end program main_fft



















  

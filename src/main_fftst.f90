program main_fft
!  consider the third equilibrium point ! 3:  q/m < 0   x = \pm (1/4/sqrt(2)^(1/3), y = 0,z² = x²  ! this the the case we study now 
!  take beta = 10, we have three pairs of pure imaginary eigenvalues, do the poincare section

! 20160203
! To make sure the analysis the frequency for escape orbit is convincable, we will do more specific Fourier Transformation, 
!   keep the same step size, and seperate the total sampling time into different overlapped part, for example, 
!   now we have almost T=2500, take T1=500, T2=[250:750], T3=[500:1000], and see if the frequency is consistent with each other.


! 20160129
! 1.  look at the transition part before the orbit escape.
!     epsl_dx = (/5.d-3, 1.5d-2, 2.5d-2, 3.5d-2, 4.5d-2, 5.5d-2, 6.d-2,  6.2d-2/)   !this is the one used to plot --- 6.2 is large 
!     the range of value of dx is choosen to be 6.d-2, 6.2d-2, 6.22d-2,  6.24d-2 
! 2. the data files are named with suffix ts


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

! 20160122 -- test fft - finished! now we need to do the Fourier analysis for x-y-z at the same time!
!   1. test the function in Gerard's book, to see if we get the right magnitude and frequency
!    the function is f(t) = cos(2*pi*0.13*t) - 1/2*sin(2*pi*0.27*t) + 3/4 *sin(2*pi*0.41*t)
!  it is better to take the step size in time as 2*pi/2^n, not 1/2^n

! 20160126
!  analyze more orbit to see, to see when it escapes, how the orbit behave in frequency domain
!  at this moment, forget about fftw, because the domimant time is orbit integration, not fft

! by a lot of testing, now we will use nsm = 2^21, lt = 2.d3*2*pi, which seperate the fundamental frequency for x and z perfectly, but for y there is still problem

! 20160128 -- a little trick to speed up
! 1  because the most time spent is writting the sampling point into mf_fft.dat, no the intergration rk78,
!    but 2^18 points are enough to do fft. However, with h=lt/nsm, the tolerance during intergration cannot be satisfied.
!    so, we introduce smpl to set the fixed stepsize for integration as h=lt/nsm/smpl, in this way, we can save the writting time and 
!    keep enough sampling points

! 2  todo: plot frequency vesus delta x, to see how the dominant frequency evolves with respect to the distance from the equilibrium points
!    keep the first 4 frequency that has the maximum magnitude  

! 3. nsm = 2*18, lt = 400*2*pi, is enough for fft. to have a good step size for rk78, we take smpl = 2**3.

use lf_mod
use po_mod 

implicit none

!np = 2^21, lt = 2.d3 *2*pi is the best match, we have x-z clearly extracted.
! fs = nsm / lt /2/pi  !  sampling frequency: per cycle


! number of sampling points on the orbit: 2**32, 2**24, 2**16 points should be good
integer, parameter ::   nsm = 2**18,  &  !number of sampling points 
                        smpl = 2**5, np = nsm*smpl,  & ! steps of integration  
			neq = 42, & ! number of dimensions for monodrony matrix computation 
			nfq = 4, & ! keep the first nfq frequencies with the biggest amplitude
			nob = 10, & !number of initial points to explore
			nf = 4, dnp = nsm/nf !number of sampling points to do FFT seperately to check the result makes sense or not,
			            !give us dt=500/8*2pi=62.5*2*pi
			
			 
real(kind=dp), parameter :: pi = 4*datan(1.d0), lt = 5.d2 * 2 * pi, & !! time length of sampling 
 			    xmax = 1.d0 !  when x > xmax, escape,  stop the integration


! specify the case, equilibrium point and parameter beta that we are going to study
integer ::  cs0, ieq 
real(kind=dp) ::  beta0,  dlf(n,n)

! local variables
integer :: i, j, ix, iy, icj, iscj,  imf, imax, ncrs, feqmf, feqpc, feqcj, &
	   isign, f_fft1, f_fft2, f_fft3, f_fft4, f_fft, obi, obf, cji, cjf, &!  --fft
	   dt, fmx_xy, fmx_z, p0, pf ! save the first 4 domimant frequencies, x-y together(planar component), z out-of-plane component

real(kind=dp) ::  wr(n),wi(n), vr(n,n), y0(neq), &   ! stability of equilibrium points
                  x0(6), vr_usd(6,3), dx0(6), epsl, epsl_dx(3,nob), dx, cj_tar, & !manifold of eq 
                  xf(6), tf, xpc(6), tmax, tfpl, & !  poinc_n of eq
                  cj_eq, cj, cj2, aaa, & ! debug auxillary variables
!                  xl(nsm,6), x(nsm),y(nsm), z(nsm), a(nsm), f(nsm/2),  h, t, & !  fft for nsm points 
                  
                  xl(nsm,6), a(nsm),x(dnp),y(dnp), z(dnp), afi(dnp), f(dnp/2),  h, t, & !  fft_sp  
                  amp_min, iamp, ramp, &! write to file  
                  tst, tend   ! evaluate time  
      
!DOUBLE COMPLEX fft1(nsm),fft2(nsm), fft3(nsm),fft4(nsm)  ! fft for nsm points

DOUBLE COMPLEX fft1(dnp),fft2(dnp), fft3(dnp),fft4(dnp) 
                 
external :: gr_lf ! the differential of the vector field and the variational matrix

! use this to initialize all the equilibrium points and the respective sign of (q/m), 
! which are the private variable in module lf_mod
call init_lf_mod

! system constant to be specified, case 1, eq 3
 cs0 = 1 ! use 2 to test the swap rule is coded correctly--ckd!
ieq = 3 ! 3, x,0,z is the case that we study currently

! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy
beta0 = 10.d0

! ***************************  Variables Initialization ************************ 
! 1: H0, 2: H1-lower, 3: H2-higher
!6.7 is too big! already escape
!epsl_dx = (/5.d-3, 1.5d-2, 2.5d-2, 3.5d-2, 4.5d-2, 5.5d-2, 6.d-2,  6.2d-2/)   !this is the one used to plot --- 6.2 is large 
!epsl_dx = (/6.d-2,  6.24d-2, 6.25d-2, 6.255d-2, 6.26d-2, 6.29d-2, 6.7d-2, 6.9d-2, 7.d-2, 0.d0/)   !for test, look at the one around 6.2d-2

!epsl_dx = (/5.d-3, 2.d-2, 3.5d-2, 5.d-2, 6.d-2, 6.2d-2, 6.24d-2, 6.25d-2, 6.255d-2, 6.26d-2/)   !for test, look at the one around 6.2d-2, the poincare map projection in x-y plane is too crowded.

!epsl_dx = (/5.d-3, 2.d-2, 4.d-2,  6.d-2, 6.2d-2, 6.24d-2, 6.25d-2, 6.255d-2, 6.26d-2/)   !for test, look at the one around 6.2d-2, wrong

!************** Energy level: H1 = H0, interval 5.9 ~ 6
! First, we look at the one around 6.2d-2,  save the result in ob9 subfolder, wrong!!!!
! 5.9995 and 5.999 are almost the same in frequencies and manifold, so use one only
! cji = 1; cjf = 1
!epsl_dx(1,:) = (/5.d-3, 2.d-2, 4.d-2, 6.0d-2, 5.99d-2, 5.992d-2, 5.994d-2,5.99601d-2, 5.998d-2, 5.999d-2/) ! discard this 

! instead of looking into a neighbourhood of 1.d-5 unit distance around the equilibrium, we take equispaced 20 points between 4.d-2 and 6.d-2 
! to look at the frequencies
!epsl_dx(1,:) = (/5.d-3, 2.d-2, 4.d-2, 6.0d-2, 5.99d-2, 5.992d-2, 5.994d-2,5.99601d-2, 5.998d-2, 5.999d-2/) ! discard this 

epsl_dx(1,:) = (/5.d-3, 2.d-2, 4.d-2, 6.0d-2/) ! use 4 to show the trend to escape, and get the coarse domain .

!For the orbits between the last and the second last orbit, take 20 orbit in between.

! instead of assigning one by one, use for loop to specify the increment
!obi = 4; obf = 4 ! test the fourth one, which escapes

!!************** Lower Energy level: H2 = H0-0.1, interva  5.1 ~ 5.2
! cji = 2; cjf = 2
!epsl_dx(2,:) = (/5.d-3, 2.d-2, 3.5d-2, 5.2d-2, 5.1d-2, 5.12d-2, 5.14d-2,5.1601d-2, 5.18d-2, 5.19d-2/)!for test, look at the one around 6.2d-2, wrong
!obi = 4; obf = 4 !9

!!!************** Higher Energy level: H3 = H0+0.0002, interva 6.015...
! cji = 3; cjf = 3  
!epsl_dx(3,:) = (/0.5d-2, 2d-2, 4d-2, 6.015d-2, 6.d-2, 6.004d-2, 6.008d-2, 6.012d-2, 6.013d-2, 6.014d-2 /)
!obi = 1; obf = 10


! -------------------------------------------------------------

!!number of crossing needed to be complete
!imax =  1 ! 10! 800
imax = 800 ! 400
tmax = 1.d2  ! 100 as the maximum integration time for one revolution


! the fixed stepsize for the integration of the orbit
h = lt/np   !

print*, 'stepsize for integration', h

amp_min = 1.d-3 ! the minimum value of magnitude we are going to consider 

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
!read(*,*)  

! set the target energy to be cj_eq
 cj_tar = cj_eq  


! compute the eigenvalues and eigenvectors 
!subroutine eigrg(a,n,isv, wr,wi,vr)
call eigrg(dlf,n,1, wr,wi,vr)

!read(*,*)  

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
! save the data of po and  its initial state, better to write the data struction of the data file in the first line
feqmf = 20;  feqpc = 21; feqcj = 22
open(feqmf,file='./dat/mf_fftst.dat',access ='append',status='replace')
write(feqmf,*) '# t	 x	 y	 z	 vx	 vy	vz'  

open(feqpc,file='./dat/pc_fftst.dat',access ='append',status='replace')  
write(feqpc,*) '# tpc, ypc(1:6), ncrs'  

open(feqcj,file='./dat/eqcj_fftst.dat',access ='append',status='replace')  
write(feqcj, *) '# ncrs, tf, cj, cj2, cj2-cj'

f_fft1 = 30; f_fft2 = 31;  f_fft3 = 32;  f_fft4 = 33;  f_fft = 35;! save x-y-z and test function a

open(f_fft1, file='./dat/fft1st.dat',access ='append',status='replace')  
open(f_fft2, file='./dat/fft2st.dat',access ='append',status='replace')  
open(f_fft3, file='./dat/fft3st.dat',access ='append',status='replace')  
open(f_fft4, file='./dat/fft4st.dat',access ='append',status='replace')  

open(f_fft, file='./dat/fftst.dat',access ='append',status='replace')
write(f_fft,*) '# f, amp, ramp, iamp, the frequency which has amplitude greater than 1.d-3'

fmx_xy = 36; fmx_z = 37 ! save the available magnitude and frequency   
open(fmx_xy, file='./dat/fqmxesc_xy1.dat',access ='append',status='replace')  
write(fmx_xy,*) '# For x-y: dx	f	 Am' 

open(fmx_z, file='./dat/fqmxesc_z1.dat',access ='append',status='replace')  
write(fmx_z,*) '# For z: dx	f	 Am' 
 
! ---------------- test function sampling-------------- 
! test f(t) = cos(2*pi*0.13*t) - 1/2*sin(2*pi*0.27*t) + 3/4 *sin(2*pi*0.41*t)
call cpu_time(tst)    
do i = 1, nsm ! --tested! result is fine!
  t =  (i-1) * h * smpl !* 2 * pi ! we have 2 cycles
  a(i)  = dcos(0.13*t) - 1./2*dsin(0.27*t) + 3./4. *dsin(0.41*t);
end do
 
print*, 'data sampling for a finished!'
call cpu_time(tend)
print*, 'Elapsed time for a and b sampling is ', tend-tst
 
do icj =  cji, cjf  !-2, 1 !-1,0 ! 0,0! -1, 1

! modify vy to select the same energy level
  
  if (icj == 1) then 
    cj_tar = cj_eq 

!  elseif (icj == 2) then 
!    cj_tar = cj_eq - 1.d-2
        
  elseif (icj == 2) then 
    cj_tar = cj_eq - 1.d-1
        
  elseif (icj == 3) then 
    cj_tar = cj_eq + 2.d-4! doesn't make too much sense with this energy level
!     cj_tar = cj_eq + 1.d-2
  endif   

  ! try with the first point   
  do ix = obi, obf ! 1,5 ! 3, 4! 8 ! 1, 2     !2
    dx = epsl_dx(icj,ix)
    print* , 'dx=',   epsl_dx(icj,ix)
!    read(*,*)
    x0(1) =  eq(1) +  dx ! x 
    x0(2) =  eq(2) +  dx ! y 
    x0(3) =  eq(3) ! z
    
! choose xdot and ydot randomly, and modify the value of vz to obtain specified energy   
    x0(4) = 1.d-6 	! vx
    x0(5) = -1.d-7	! vy
    x0(6) = 0.d0 	! vz
    
 ! check the energy level before changing the value of vz 
    call gr_cjlf(x0, cj)
    print*, 'before vz, cj, cj_tar, x0', cj, cj_tar, x0
!    read(*,*)  
  
    call lf2cj(x0, 6, cj_tar, iscj) ! 6: vz to be modified
      
    if(iscj ==0) then 
      print*, 'Vz*Vz < 0 in the energy level', cj_tar,'! icj=',icj,' ix=',ix,' iy=',iy
      read(*,*) 
      cycle
    endif
      
    call gr_cjlf(x0, cj)
   
 ! check the energy here, and then decide which energy level to take  
    print*, 'initial point x0 on energy level cj_tar', cj, cj_tar, x0 
     read(*,*)  
    
!    cycle
    
 ! integrate for enough long time? how long? fixed step size
 !  and save the data 
!    tf = 1.d2 
!    call plob(x0,0.d0, tf, feqmf, gr_lf, xf) 
!    print*,' fininsh plob'
!    read(*,*)


! ****************************** Poincare maps ******************************    
!subroutine poinc_n(y0, imax, tmax, ind, z0, fpc, tf, yf, ncrs)  
!    call poinc_n(x0, imax, tmax, 3, eq(3), feqpc, tf, xpc, ncrs)
!    write(feqpc, *)   
!    call gr_cjlf(xpc, cj2)
!    write(feqcj, *) ncrs, tf,cj, cj2 , cj2-cj
!    write(*, *) ncrs, tf,cj, cj2 , cj2-cj
!     read(*,*)
! ****************************** orbit sampling ******************************    
    call cpu_time(tst)    
!subroutine plob_fxd(y0,t0, h, nsm, smpl, xmax, ftag, deriv, yall) 
    xl = 0.d0 ! to avoid mistake, for every orbit, reinitialize the state vector
    
    call plob_fxd(x0, 0.d0, h, nsm, smpl, xmax, feqmf, gr_lf, xl) 
    
    call cpu_time(tend)
    print*, 'Elapsed time for orbit sampling is ', tend-tst
    print * ,'Integration of orbit finished!'

! check the final state and stop time? 
    print*, 'check the final state'
    print*, nsm, xl(nsm,:)
!   read(*,*)
    
    call gr_cjlf(xl(nsm,:), cj2)
    write(*, *) lt, cj, cj2 , cj2-cj
!    read(*,*)  
    
!    cycle

! for every 2**15 points, we do a fft, and compare the result
!! use the components of the state vector from the integration with fixed stepsize 
  do j = 1, 2*3-1 ! the available time length is almost 2400, 6 times of our new time interval 62.5*2*pi
    
    p0 = (j-1)*dnp/2 + 1
    
    pf = p0 + dnp -1
    print*, 'p0,pf', p0,pf
    read(*,*)   
    
    x = xl(p0:pf, 1) ! x
    y = xl(p0:pf, 2) ! y 
    z = xl(p0:pf, 3) ! z 
    
    afi = a(p0:pf) !a
    
    call cpu_time(tst)    
    call dtwofft(x, y,fft1,fft2, dnp) ! for x and y
    
    call dtwofft(z, afi,fft3,fft4, dnp) ! for z and a
    call cpu_time(tend)

    print*, 'Elapsed time for  2 * twofft is ', tend-tst

!      read(*,*)
! Extract real magnitude and frequency from fft1 returned by twofft, save them in seperate files   
! output the available frequency and magnitude  

    print*,'#********* for x **********' 
    write(f_fft, *) '# ********* x *********'
 
!subroutine fqext(fft, nsm, lt, f_fft1, f_fft, fmx, nfq, dx)
    call fqext(fft1, dnp, lt/nf, f_fft1, f_fft, fmx_xy, nfq, dx) 
!    print*, 'finshi fft for x'
!    read(*,*)
    
    write(*, *) '# ********* y *********'
    write(f_fft, *) '# ********* y *********'
    call fqext(fft2, dnp, lt/nf, f_fft2, f_fft, fmx_xy, nfq, dx) 

    write(*, *) '# ********* z *********'  
    write(f_fft, *) '# ********* z *********'
    call fqext(fft3, dnp, lt/nf, f_fft3, f_fft, fmx_z, nfq, dx) 
  
    write(f_fft, *) '# ********* a *********'
    call fqext(fft4, dnp, lt/nf, f_fft4, f_fft, 6, nfq, dx) 

  enddo ! for one time interval
  
! for different orbit, there is one blank line to seperate orbit as block 

  enddo ! this is for one orbit
 
! for different energy, add 2 blank lines to use index to plot 
  write(f_fft1, *); write(f_fft2, *);  write(f_fft3, *);  write(f_fft4, *); write(f_fft, *); write(fmx_xy, *); write(fmx_z, *)
  write(f_fft1, *); write(f_fft2, *);  write(f_fft3, *);  write(f_fft4, *); write(f_fft, *); write(fmx_xy, *); write(fmx_z, *)
  write(feqpc, *);  write(feqcj, *) 
  write(feqpc, *);  write(feqcj, *)  
enddo  ! one energy level

 close(f_fft1); close(f_fft2);  close(f_fft3); close(f_fft4); close(f_fft); close(fmx_xy); close(fmx_z)
 
stop
end program main_fft



















  

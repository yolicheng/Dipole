program main_pc ! eq = 3, beta = 10
!  we have three pairs of pure imaginary eigenvalues, do the poincare section as Gerard explained

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

implicit none

!np = 2^21, lt = 2.d3 *2*pi is the best match, we have x-z clearly extracted.
! fs = nsm / lt /2/pi  !  sampling frequency: per cycle


! number of sampling points on the orbit: 2**32, 2**24, 2**16 points should be good
integer, parameter ::   nsm = 2**18,  &  !number of sampling points 
                        smpl = 2**5, np = nsm*smpl,  & ! steps of integration  
			neq = 42, & ! number of dimensions for monodrony matrix computation 
			nfq = 4 ! keep the first 4 frequencies with the biggest amplitude

real(kind=dp), parameter :: pi = 4*datan(1.d0), lt = 5.d2 * 2 * pi, & !! time length of sampling 
 			    xmax = 1.d0 !  when x > xmax, escape,  stop the integration


! specify the case, equilibrium point and parameter beta that we are going to study
integer ::  cs0, ieq 
real(kind=dp) ::  beta0,  dlf(n,n)

! local variables
integer :: i, j, ix, iy, icj, iscj,  imf, imax, ncrs,  fpc, fmfpc, & !poinc map
	   isign, fmf_fft, fx_fft, fy_fft, fz_fft, fa_fft, ffqmx, nfrg !  --fft

real(kind=dp) ::  runit, tunit, vuint, & ! the unit of dimensionless scaling
		  wr(n),wi(n), vr(n,n), y0(neq), &   ! stability of equilibrium points
                  x0(6), vr_usd(6,3), dx0(6), epsl, epsl_dx(4), incdx, dx, cj_tar, & !manifold of eq 
                  xf(6), tf, xpc(6), tmax, tfpl, & !  poinc_n of eq
                  cj_eq, cj, cj2, aaa, & ! debug auxillary variables
                  xl(nsm,6), x(nsm),y(nsm), z(nsm), a(nsm), h, t, & !, h fft  
                  f(nsm/2), amp1, amp2,  amp3, amp4, iamp, ramp, &! write to file  
                  tst, tend, &  ! evaluate time  
                  fqmx_x(nfq,2), fqmx_y(nfq,2), fqmx_z(nfq,2), fqmx_a(nfq,2) !the domimant frequency which has the first nfq maximum amplitude
      
DOUBLE COMPLEX fft1(nsm),fft2(nsm), fft3(nsm),fft4(nsm) 
 character(len=70) ::  fnmf, fnx, fny, fnz, fna, fnmx, fnpc, fnmfpc  


! use this to initialize all the equilibrium points and the respective sign of (q/m), 
! which are the private variable in module lf_mod
call init_lf_mod

! the fixed stepsize for the integration of the orbit
h = lt/np   !
nfrg = 20 ! 20 ! compute more point in between the escape and the bounded one next to it

print*, 'stepsize for integration', h

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

call lfunit(beta, runit, tunit, vuint)
print*, 'runit, tunit, vuint', runit, tunit, vuint

! the differential matrix of the lorentz force with respective to state 
!subroutine dflrtz( x0, beta, cs, sgn, dlf)
call dflrtz( eq, dlf)

!do i = 1, n
! write(*,'(6f8.4)') dlf(i,:) 
!enddo

! check the energy 
!C1 = -3   
!C2 = 2.2894  
!C3 = 1.8899 

call gr_cjlf(eq, cj_eq)
print*,'check energy of equilibrium point!, cj, ieq,', cj_eq, ieq, eq
!read*

! set the target energy to be cj_eq
 cj_tar = cj_eq  


! compute the eigenvalues and eigenvectors 
!subroutine eigrg(a,n,isv, wr,wi,vr)
call eigrg(dlf,n,1, wr,wi,vr)

!read*

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

!!number of crossing needed to be complete
!imax =  1 ! 10! 800
!imax = 800 ! 400
!tmax = 1.d2  ! 100 as the maximum integration time for one orbit

! 5.d-2 is the allowable biggest value, 5.5d-2 begin to escape
 
epsl_dx = (/5.d-3, 2.d-2, 4.d-2, 6.0d-2/)    !this is the one used to plot --- 6.2 is large 
incdx = (epsl_dx(4) - epsl_dx(3) )/ (nfrg +1) ! select 20 points in between, in total  22 points


! *************   select the energy level ************** 
 cj_tar = cj_eq; icj = 0 ! h0
! cj_tar = cj_eq - 1.d-1; icj = 1 ! h1
! cj_tar = cj_eq - 1.d-2; icj = 2 ! h2 
! cj_tar = cj_eq + 2.d-4; icj = 3 ! h3

! the idea to write to mulitiple files, 
! ! build filename -- i.dat
!write(fn,fmt='(a,i0,a)') filenum, '.dat'


!  feqpc = 21; feqcj = 22
!open(feqpc,file='./dat/pc_fft.dat',access ='append',status='replace')  


fmf_fft = 20; fx_fft = 30; fy_fft = 31;  fz_fft = 32;  fa_fft = 33;  ! save x-y-z and test function a
ffqmx = 36 ! save the dominant magnitude and frequency

write(fnmf, fmt='(a,i0,a)') './dat/mf_fft', icj, '.dat' 
write(fnx, fmt='(a,i0,a)') './dat/x_fft', icj, '.dat' 
write(fny, fmt='(a,i0,a)') './dat/y_fft', icj, '.dat' 
write(fnz, fmt='(a,i0,a)') './dat/z_fft', icj, '.dat' 
write(fna, fmt='(a,i0,a)') './dat/a_fft', icj, '.dat' 
write(fnmx, fmt='(a,i0,a)') './dat/fqmx', icj, '.dat' 

!! open an empty file using the fname --- not to overwrite to good and completed result
!open(fmf_fft, file= fnmf, access ='append',status='replace')  
!open(fx_fft, file= fnx, access ='append',status='replace')  
!open(fy_fft, file= fny, access ='append',status='replace')
!open(fz_fft, file= fnz, access ='append',status='replace')
!open(fa_fft, file= fna, access ='append',status='replace')
!open(ffqmx,  file= fnmx, access ='append',status='replace')  
!write(ffqmx,*) '# dx 	fx	ax 	fy 	ay 	fz 	az'


write(fnpc, fmt='(a,i0,a)') './dat/pc', icj, '.dat'     ! poincare maps
write(fnmfpc, fmt='(a,i0,a)') './dat/mfpc', icj, '.dat'  ! the orbit used to get the poincare maps
open(fpc,  file=fnpc,  access ='append',status='replace')  
open(fmfpc,file=fnmfpc,access ='append',status='replace')  
nfrg = 0 ! only consider the first four orbits
print* fpc, fmfpc 
read*

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
 
nfrg = 0
! instead of using epsl_dx to specify the dx, we use a do loop and assign the value of dx there.   
do ix = 1,  nfrg + 4  !   the first four are 
  if (ix .le. 4) then  
    dx = epsl_dx(ix)
  else 
    dx = epsl_dx(3) + (ix-4) * incdx
  endif     
  
!  ! to check if dx is okay  -- ckd
!  print*, 'dx', dx 
!  cycle 
  
  
  x0(1) =  eq(1) + dx ! x 
  x0(2) =  eq(2) + dx ! y 
  x0(3) =  eq(3) ! z
    
! choose xdot and ydot randomly, and modify the value of vz to obtain specified energy 
! give a very small velocity  
  x0(4) = 1.d-6 	! vx
  x0(5) = -1.d-7	! vy
  x0(6) = 0.d0 	! vz
    
 ! check the energy level before changing the value of vz 
  call gr_cjlf(x0, cj)
  print*, 'before vz, cj, cj_tar, x0', cj, cj_tar, x0
!  read*

  call lf2cj(x0, 6, cj_tar, iscj) ! 6: vz to be modified
      
  if(iscj ==0) then 
    print*, 'Vz*Vz < 0 in the energy level', cj_tar,'! icj=',icj,' ix=',ix,' iy=',iy
  ! read* 
    cycle
  endif
     
  call gr_cjlf(x0, cj)
   
 ! check the energy here, and then decide which energy level to take  
  print*, 'initial point x0 on energy level cj_tar', x0, cj, cj_tar
!  read*
 
 ! ******************************poincare maps ****************************** 
 
! the stop condition:  number of crossing through poincare section
!                    set imax = 800, to have enough point,  and  then use the stop time to integrate the orbit    

! if we want to use poinc in pomod we need to initialize before 

! Step size and error control for rk78 integration 
hmin = 1.d-10
hmax = 1.d0
e    = 1.d-14! 1.d-13 is from Gerard's book 

symt = 2
sec0 = 0.d0 

! subroutine init_po(symt, sec, hmin, hmax, e, tmax, tol, prsc) - for pomod
call  init_po(symt, sec0, hmin, hmax, e, tmax, tol, prsc)

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
  
  
   cycle  ! at this moment, regardless of the fft analysis, just do the poincare section
    
    
! ****************************** orbit sampling ******************************    
  call cpu_time(tst)    
!subroutine plob_fxd(y0,t0, h, nsm, smpl, xmax, ftag, deriv, yall) 
  xl = 0.d0 ! to avoid mistake, for every orbit, reinitialize the state vector
 
! subroutine plob_fxd(y0,t0, h0, nsm, smpl, xmax, ftag, yall,  deriv) 
  call plob_fxd(x0, 0.d0, h, nsm, smpl, xmax, fmf_fft, xl, gr_lf) 
    
  call cpu_time(tend)
  print*, 'Elapsed time for orbit sampling is ', tend-tst
  print * ,'Integration of orbit finished!'

! check the final state and stop time? 
  print*, 'check the final state'
  print*, nsm, xl(nsm,:)
!  read*
    
  call gr_cjlf(xl(nsm,:), cj2)
  write(*, *)  cj, cj2 , cj2-cj
!  read*  
    
!! use the components of the state vector from the integration with fixed stepsize    
  x = xl(:, 1) ! x
  y = xl(:, 2) ! y 
  z = xl(:, 3) ! z 
    
  call cpu_time(tst)    
  call dtwofft(x, y,fft1,fft2, nsm) ! for x and y
    
  call dtwofft(z, a,fft3,fft4, nsm) ! for z and a
  call cpu_time(tend)

  print*, 'Elapsed time for  2 * twofft is ', tend-tst
!   ! read*
! Extract real magnitude and frequency from fft1 returned by twofft, save them in seperate files   
  print*,'#********* for x **********'  
  write(*, *) '# ********* x *********'
 
!subroutine fqext(fft, nsm, lt, f_fft, nfq, fqmx)
  call fqext(fft1, nsm, lt, fx_fft, nfq, fqmx_x) 
  
  write(*, *) '# ********* y *********'
  call fqext(fft2, nsm, lt, fy_fft, nfq, fqmx_y) 
  
  write(*, *) '# ********* z *********'
  call fqext(fft3, nsm, lt, fz_fft, nfq, fqmx_z) 
  
  write(*, *) '# ********* a *********'
  call fqext(fft4, nsm, lt, fa_fft, nfq, fqmx_a)  

! Pick the first nfq=4 dominant frequencies and then save them in the form 
! dx - fx- ax - fy - ay - fz - az, keep 4 digits after the decimer, for convenient use in latex tabular

  print*, 'fqmx: dx - fx- ax - fy - ay - fz - az ' 
  do j = 1, nfq
    write(ffqmx,'(7f10.4)') dx, fqmx_x(j, :), fqmx_y(j, :), fqmx_z(j, :)
    write(*, '(7f10.4)')    dx, fqmx_x(j, :), fqmx_y(j, :), fqmx_z(j, :) 
    write(*, '(2f10.4)'), fqmx_a(j,:)   
  enddo    
  
  write(ffqmx,*) ! add a blank line to seperate different orbit with different distance
      
!read*

enddo ! this is for one orbit
  

 close(fx_fft); close(fy_fft);  close(fz_fft); close(fa_fft); close(ffqmx);
 
stop
end program main_pc



















  

program main_fftw  !--- not a good one , just to test fftww, do not use this one ...., the good one is main_fft.f90
! test the package fftww, which is the most efficient open library to do fftw

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
!   2. use fftw from numerical recipes to do the Fast Fourier Transformation
!   3. also start from H0, and then test the other energy

! 20160122 -- test fft 
!   1. test the function in Gerard's book, to see if we get the right magnitude and frequency

!    the function is f(t) = cos(2*pi*0.13*t) - 1/2*sin(2*pi*0.27*t) + 3/4 *sin(2*pi*0.41*t)

! since most of the time is spent on the integration of the orbit, 3s for 2^19 points, 6s for 2^20 points
! and the 

use lf_mod

! to use fftw3, first, we need to add these two lines
use, intrinsic :: iso_c_binding
include 'fftw3.f03'

implicit none

! number of sampling points on the orbit: 2**32, 2**24, 2**16 points should be good
! to obtain more precise result, we need: 
!  1. take longer sampling time :  lt 
!  2. take more sampling points (of power 2): nsm

integer, parameter ::   nsm = 2**20,  &  !number of sampling points 
                        smpl = 1, np = nsm*smpl,  & ! steps of integration  
                        
			neq = 42 ! number of dimensions for monodrony matrix computation 


real(kind=dp), parameter :: pi = 4*datan(1.d0), lt = 4.d2 * 2 * pi !! time length of sampling 
    
         ! the sampling frequency  fs = np / lt:  cycle/second 
         !                         fs = np/ lt * 2*pi : radian/second, in my case, i'm taking this one
         

! specify the case, equilibrium point and parameter beta that we are going to study
integer ::  cs0, ieq, & ! case of lf 
	    nthrd, ifw ! multi-thread fftw
	     
real(kind=dp) ::  beta0,  dlf(n,n)

! local variables
integer :: i, j, ix, iy, icj, iscj,  imf, imax, ncrs, feqmf, feqpc, feqcj, &
	   isign, f_fftw1, f_fftw2, f_fftw3, f_fftw4, f_fftw !, nsm2 !--fftw

real(kind=dp) ::  wr(n),wi(n), vr(n,n), y0(neq), &   ! stability of equilibrium points
                  x0(6), vr_usd(6,3), dx0(6), epsl, epsl_dx(8), cj_tar, & !manifold of eq 
                  xf(6), tf, xpc(6), tmax, tfpl, & !  poinc_n of eq
                  cj_eq, cj, cj2, aaa, & ! debug auxillary variables
                  xl(np,6), data1(nsm),data2(nsm),  h, t, & !, h fftw  
                  f, amp1, amp2,  amp3, amp4 ! write to file    
      
                  
external :: gr_lf ! the differential of the vector field and the variational matrix


! ********* specified to multithread fftw *********
!Second, before calling any FFTW routines, you should call the function:

ifw = fftw_init_threads() 

!This function, which need only be called once, performs any one-time initialization required to use threads on your system. It returns zero if there was some error (which should not happen under normal circumstances) and a non-zero value otherwise.

!Third, before creating a plan that you want to parallelize, you should call:

!     void fftw_plan_with_nthreads(int nthreads);
nthrd = 4
call fftw_plan_with_nthreads( nthrd) 
!The nthreads argument indicates the number of threads you want FFTW to use (or actually, the maximum number). All plans subsequently created with any planner routine will use that many threads. You can call fftw_plan_with_nthreads, create some plans, call fftw_plan_with_nthreads again with a different argument, and create some more plans for a new number of threads. Plans already created before a call to fftw_plan_with_nthreads are unaffected. If you pass an nthreads argument of 1 (the default), threads are disabled for subsequent plans.

!Given a plan, you then execute it as usual with fftw_execute(plan), and the execution will use the number of threads specified when the plan was created. When done, you destroy it as usual with fftw_destroy_plan. As described in Thread safety, plan execution is thread-safe, but plan creation and destruction are not: you should create/destroy plans only from a single thread, but can safely execute multiple plans in parallel.

!There is one additional routine: if you want to get rid of all memory and other resources allocated internally by FFTW, you can call:
!   void fftw_cleanup_threads(void);
!which is much like the fftw_cleanup() function except that it also gets rid of threads-related data. You must not execute any previously created plans after calling this function.


!For example, here is how we would allocate an L × M 2d real array:

!       real(C_DOUBLE), pointer :: arr(:,:)
!       type(C_PTR) :: p
!       p = fftw_alloc_real(int(L * M, C_SIZE_T))
!       call c_f_pointer(p, arr, [L,M])
!       ...use arr and arr(i,j) as usual...
!       call fftw_free(p)
!and here is an L × M × N 3d complex array:

!       complex(C_DOUBLE_COMPLEX), pointer :: arr(:,:,:)
!       type(C_PTR) :: p
!       p = fftw_alloc_complex(int(L * M * N, C_SIZE_T))
!       call c_f_pointer(p, arr, [L,M,N])
!       ...use arr and arr(i,j,k) as usual...
!       call fftw_free(p)

!   Alternatively, for an in-place r2c transform, as described in the C documentation we must pad the first dimension of the real input with an extra two entries (which are ignored by FFTW) so as to leave enough space for the complex output. The input is allocated as a 2[L/2+1] × M × N array, even though only L × M × N of it is actually used. In this example, we will allocate the array as a pointer type, using ‘fftw_alloc’ to ensure aligned memory for maximum performance (see Allocating aligned memory in Fortran); this also makes it easy to reference the same memory as both a real array and a complex array.
  
!       real(C_DOUBLE), pointer :: in(:,:,:)
!       complex(C_DOUBLE_COMPLEX), pointer :: out(:,:,:)
!       type(C_PTR) :: plan, data
!       data = fftw_alloc_complex(int((L/2+1) * M * N, C_SIZE_T))
!       call c_f_pointer(data, in, [2*(L/2+1),M,N])
!       call c_f_pointer(data, out, [L/2+1,M,N])
!       plan = fftw_plan_dft_r2c_3d(N,M,L, in,out, FFTW_ESTIMATE)
!       ...
!       call fftw_execute_dft_r2c(plan, in, out)
!       ...
!       call fftw_destroy_plan(plan)
!       call fftw_free(data)

!    void fftw_execute_dft_r2c(
!          const fftw_plan p,
!          double *in, fftw_complex *out);
!     

!*********************************************************
! use this to initialize all the equilibrium points and the respective sign of (q/m), 
! which are the private variable in module lf_mod
call init_lf_mod

! the fixed stepsize for the integration of the orbit
h = lt/np  ! time interval 

     
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

do i = 1, n
 write(*,'(6f8.4)') dlf(i,:) 
enddo

! check the energy 
!C1 = -3   
!C2 = 2.2894  
!C3 = 1.8899 

call gr_cjlf(eq, cj_eq)
print*,'check energy of equilibrium point!, cj, ieq,', cj_eq, ieq, eq
read(*,*) aaa

! set the target energy to be cj_eq
 cj_tar = cj_eq  


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
open(feqmf,file='./dat/mf_fftww.dat',access ='append',status='replace')  
open(feqpc,file='./dat/pc_fftww.dat',access ='append',status='replace')  
!open(feqcj,file='./dat/eqcj.dat',access ='append',status='replace')  

f_fftw1 = 30; f_fftw2 = 31;  f_fftw3 = 32;  f_fftw4 = 33; f_fftw = 35
open(f_fftw1, file='./dat/fftw1.dat',access ='append',status='replace')  
open(f_fftw2, file='./dat/fftw2.dat',access ='append',status='replace')  

open(f_fftw3, file='./dat/fftw3.dat',access ='append',status='replace')  
open(f_fftw4, file='./dat/fftw4.dat',access ='append',status='replace')  

! the magnitude which is greater than 1.d-3 and the respective frequency
open(f_fftw, file='./dat/fftw.dat',access ='append',status='replace')  

! test the same function, and compare with the result of matlab
! test f(t) = cos(2*pi*0.13*t) - 1/2*sin(2*pi*0.27*t) + 3/4 *sin(2*pi*0.41*t)

do i = 1, np ! tested!
  t =  (i-1) * h  !* 2 * pi ! we have 2 cycles
  
  data1(i)  = dcos(0.13*t) - 1./2*dsin(0.27*t) + 3./4. *dsin(0.41*t);
  data2(i)  = dcos(1.3*t) + 1./4*dsin(2.7*t) - 1./2. *dsin(4.1*t);


!  data1(i) = 3.d0*dsin(8.d0*t) 
!  data2(i) = dcos(13.d0*t) - 1.d0/2.d0*dsin(3.d0*t) + 3.d0/4.d0 *dsin(5.d0*t)
end do
 
print*, 'data sampling finished!'
call dtwofftw(data1,data2,fftw1,fftw2, nsm)


do i = 1, np/2
! if we use real time
  f  = i / lt  * pi ! ! for twofftw, amp = amp/(nsm/2), f = f/2, but if t is the real time, f = f*2*pi/2 = f*pi
  
! cycle per second
!  f  = i / lt / 2! ! for twofftw, amp = amp/(nsm/2), f = f/2, but if t is the real time, f = f*2*pi/2 = f*pi
    
  amp1 = dabs( fftw1(i) )/ nsm * 2
  amp2 = dabs( fftw2(i) )/ nsm * 2
 
  if(i == 1) then 
    amp1 = amp1 / 2 
    amp2 = amp2 / 2 
  endif
  
  write(f_fftw1, '(2f20.10)')  f, amp1 
  write(f_fftw2, '(2f20.10)')  f, amp2 
  
  if (dabs(amp1) > 1.d-2 ) then 
    print*,'f1:', f, amp1
    write(f_fftw, *) 'f1:', f, amp1
  endif
  
  if (dabs(amp2) > 1.d-2 ) then 
   print*, 'f2:', f, amp2
   write(f_fftw, *) 'f2:', f, amp2
  endif
enddo 
 

!!number of crossing needed to be complete
!imax =  1 ! 10! 800
imax = 800 ! 400
tmax = 1.d2  ! 100 as the maximum integration time for one orbit



! 5.d-2 is the allowable biggest value, 5.5d-2 begin to escape
!epsl_dx = (/5.d-3,  7.5d-3, 1.d-2, 2.5d-2, 5.d-2, 5.5d-2, 6.d-2, 1.d-1/)  ! the second is too close 

epsl_dx = (/5.d-3, 1.5d-2, 3.d-2, 4.d-2, 5.d-2, 5.5d-2, 6.d-2, 1.d-1/)   !
 
  
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
  do ix = 1, 1 ! 1,5  
      
    x0(1) =  eq(1) +  epsl_dx(ix) ! x 
    x0(2) =  eq(2) +  epsl_dx(ix) ! y 
    x0(3) =  eq(3) ! z
    
! choose xdot and ydot randomly, and modify the value of vz to obtain specified energy   
!    x0(4) = 0.d0
    x0(4) = 1.d-4 ! vx
    x0(5) = -1.d-5 ! vy
     
    x0(6) = 0.d0 ! vz
   
    call gr_cjlf(x0, cj)
   
    
 ! check the energy level before changing the value of vz 
    print*, 'before vz, cj, cj_tar, x0', cj, cj_tar, x0
!     read(*,*) aaa

    call lf2cj(x0, 6, cj_tar, iscj) ! 6: vz to be modified
      
    if(iscj ==0) then 
      print*, 'Vz*Vz < 0 in the energy level', cj_tar,'! icj=',icj,' ix=',ix,' iy=',iy
      read(*,*) aaa
      cycle
    endif
      
    call gr_cjlf(x0, cj)
   
    ! check the energy here, and then decide which energy level to take  
    print*, 'initial point x0 on energy level cj_tar', x0, cj, cj_tar
!    read(*,*) aaa
    
 ! integrate for enough long time? how long? fixed step size
 !  and save the data 
!    tf = 1.d2 
!    call plob(x0,0.d0, tf, feqmf, gr_lf, xf) 
!    print*,' fininsh plob'
!    read(*,*)
    
!subroutine plob_fxd(y0, t0,  h, np, ftag, deriv, yall) 
    call plob_fxd(x0, 0.d0, h/smpl, np, feqmf, gr_lf, xl) 
    
!! use the components of the state vector from the integration with fixed stepsize
!! 
    data1 = xl(1:np:smpl, 1) ! x
    data2 = xl(1:np:smpl, 2) ! y 
!    
    call dtwofftw(data1,data2,fftw1,fftw2, nsm)
!    

do i = 1, np/2
! if we use real time
  f  = i / lt  * pi ! ! for twofftw, amp = amp/(nsm/2), f = f/2, but if t is the real time, f = f*2*pi/2 = f*pi
  
! cycle per second
!  f  = i / lt / 2! ! for twofftw, amp = amp/(nsm/2), f = f/2, but if t is the real time, f = f*2*pi/2 = f*pi
    
  amp1 = dabs( fftw1(i) )/ nsm * 2
  amp2 = dabs( fftw2(i) )/ nsm * 2
  
  if(i == 1) then 
    amp1 = amp1 / 2 
    amp2 = amp2 / 2 
  endif
    
  write(f_fftw3, '(2f20.10)')  f, amp1 
  write(f_fftw4, '(2f20.10)')  f, amp2 
  
  if (dabs(amp1) > 1.d-3 ) then 
    print*,'f3:', f, amp1
    write(f_fftw, *) 'f3:', f, amp1
  endif
  
  if (dabs(amp2) > 1.d-3 ) then 
   print*, 'f4:', f, amp2
  write(f_fftw, *) 'f4:', f, amp1
  endif
  
enddo 

!    write(*,*) 'Fourier transform of first function:'
!    call prntft(fftw1, nsm2, f_fftw1)
!    write(f_fftw1,*) ; write(f_fftw1,*)
!    
!    write(*,*) 'Fourier transform of second function:'
!    call prntft(fftw2, nsm2, f_fftw2)

!!     invert transform
!    isign = -1
!    call dfour1(fftw1, nsm, isign)
!    
!    write(*,*) 'Inverted transform = first function:'
!    call prntft(fftw1, nsm2, f_fftw1m)
!   
!    call dfour1(fftw2, nsm, isign)
!    write(*,*) 'Inverted transform = second function:'
!    call prntft(fftw2, nsm2, f_fftw2m)

!    read(*,*)  

!      
      
      
  
! check the final state and stop time? 
    print*, 'check the final state'
    print*, np, xl(np,:)
  
    call gr_cjlf(xl(np,:), cj2)
  
    write(*, *)  cj, cj2 , cj2-cj
    read(*,*) aaa
  enddo
  
!  write(feqmf,*); write(feqmf,*) !add two blank lines
  
enddo 


 close(feqcj)
stop
end program main_fftw



















  

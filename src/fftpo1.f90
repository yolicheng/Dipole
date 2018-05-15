program fftpo1

! for a given p.o., with X0 + TP0, use FFT to get the analytical Fourier expression 
! TODO : not necessary to test the package fftww, which is the most efficient open library to do fftw


use lf_mod

implicit none

! number of sampling points on the orbit: 2**32, 2**24, 2**16 points should be good
! to obtain more precise result, we need: 
!  1. take longer sampling time :  lt 
!  2. take more sampling points (of power 2): nsm

! for 1 period: 2^12=4096 is good enough 


integer, parameter ::   nsm  = 2**12,  &  !number of sampling points 
                        nfq  = 6 ! 4 +1 : keep the first 4 frequencies with the biggest amplitude + the constant part                        

!real(kind=dp), parameter :: pi = 4*datan(1.d0) 

 ! lt = 10.d1 * 2 * pi !! time length of sampling 
 ! the sampling frequency  fs = nsm / lt:     cycle/second 
 !                         fs = nsm / lt * 2*pi :   radian/second, in my case, i'm taking this one
         
! module-based variables
integer :: cs, ieq 
real(kind=dp) :: beta

! local variables
integer::  i, j, isft, fob, fpo_fft, fx_fft, fy_fft, fz_fft, fvx_fft, fvy_fft, fvz_fft, ffqmx !  --fft

! module-based variables
!real(kind=dp) ::  runit, tunit, vuint, & ! the unit of dimensionless scaling

real(kind=dp) ::  lt, tp0, x0(6), xf(6), xl(nsm,6), x(nsm),y(nsm), z(nsm),  & 
                  vx(nsm), vy(nsm), vz(nsm), h,  xmax, & ! fft  
!                   a(nsm), f(nsm/2), amp1, amp2,  amp3, amp4, iamp, ramp, &! write to file  
                  tst, tend, &  ! evaluate time  
                  
                  !the domimant frequency which has the first nfq maximum amplitude
                  fqmx_x(nfq, 4), fqmx_y(nfq, 4), fqmx_z(nfq, 4), fqmx_vx(nfq, 4), fqmx_vy(nfq, 4), fqmx_vz(nfq,4) , &
                  fourx, foury, fourz,  t
      
DOUBLE COMPLEX   fft1(nsm),fft2(nsm), fft3(nsm),fft4(nsm), fft5(nsm),fft6(nsm) 
 character(len=70) ::  fnpo, fnx, fny, fnz, fnvx, fnvy, fnvz, fna, fnmx  

 
! specify the case we are studying 
! ***************** case = 3 , N =[0 0 1] 
! -----  ieq = 2,   for +t direction, x-y plane symmetry, z-axis, origin 
! (x,y, z, vx,vy, vz ) ---> (x,y, -z, vx,vy, -vz ) 
!                           (-x,-y, z, -vx,-vy, vz) 
 cs = 3
! --- ieq =2 , beta = 1.d0 
ieq = 2 ! 0, y,z, 0, 0,0 
beta  = 1.d0 
!  
! Provided the case, beta and ieq,  beta, cs, sgn and eq are initialized by subroutine init_lf from module lf_mod
! subroutine init_lf(beta0, cs0, ieq)
call init_lf(beta, cs, ieq) 

print*, 'check the assignment with module' !--ckd
print*, 'beta, cs, ieq, sgn, eq', beta, cs, ieq, sgn, eq
read*     

fpo_fft = 20; fx_fft = 30; fy_fft = 31;  fz_fft = 32;   ! save x-y-z and  vx-vy-vx
fvx_fft = 33; fvy_fft = 34;  fvz_fft = 35;
 
ffqmx = 36; fob = 37 ! save the dominant magnitude and frequency

write(fnpo, fmt='(a)') './dat/poformfly/fft/po_fft.dat' 
write(fnx,  fmt='(a)') './dat/poformfly/fft/x_fft.dat' 
write(fny,  fmt='(a)') './dat/poformfly/fft/y_fft.dat' 
write(fnz,  fmt='(a)') './dat/poformfly/fft/z_fft.dat' 
write(fnvx,  fmt='(a)') './dat/poformfly/fft/vx_fft.dat' 
write(fnvy,  fmt='(a)') './dat/poformfly/fft/vy_fft.dat' 
write(fnvz,  fmt='(a)') './dat/poformfly/fft/vz_fft.dat' 
write(fnmx, fmt='(a)') './dat/poformfly/fft/fqmx.dat' 

! open an empty file using the fname
open(fob, file = './dat/poformfly/fft/ob.dat', access='append', status='replace')
open(fpo_fft, file= fnpo, access ='append',status='replace')  
open(fx_fft, file= fnx, access ='append',status='replace')  
open(fy_fft, file= fny, access ='append',status='replace')
open(fz_fft, file= fnz, access ='append',status='replace')

open(fvx_fft, file= fnvx, access ='append',status='replace')  
open(fvy_fft, file= fnvy, access ='append',status='replace')
open(fvz_fft, file= fnvz, access ='append',status='replace')

open(ffqmx,  file= fnmx, access ='append',status='replace')  
write(ffqmx,*) '# fx	ax 	fy 	ay 	fz 	az'

! ---------------- test function sampling-------------- 
! test f(t) = cos(2*pi*0.13*t) - 1/2*sin(2*pi*0.27*t) + 3/4 *sin(2*pi*0.41*t)
!call cpu_time(tst)    
!do i = 1, nsm ! --tested! result is fine!
!  t =  (i-1) * h * smpl !* 2 * pi ! we have 2 cycles
!  a(i)  = dcos(0.13*t) - 1./2*dsin(0.27*t) + 3./4. *dsin(0.41*t);
!end do

!print*, 'Reference data sampling finished!' 
!read*
  
  
 ! 40-th  -- the one that is parallel to x-y plane ! --- not good, in fact are two orbits.... almost overlap for 2 
!    0.18432164858035E+01  0.00000000000000E+00  0.36664748068655E+00  0.10117245290528E+01  0.10317592130296E+01  0.00000000000000E+00  0.00000000000000E+00 -0.23038657110273E+01  0.44408920985006E-15  0.25600000    1
  tp0 = 0.18432164858035d1
  x0 = (/ 0.d0,  0.36664748068655d0,  0.10117245290528d1,  0.10317592130296d1 , 0.d0,  0.d0 /)
  
  
  ! 46-th  i=45, the one before almost perpendicular to x-y plane-- turns out it is not almost perpendicular, but the size is small
!    0.35864311809759E+01  0.00000000000000E+00 -0.55680992271445E+00  0.11417028016985E+01  0.44424859433677E+00  0.00000000000000E+00  0.00000000000000E+00 -0.18033810944970E+01 -0.44408920985006E-15  0.25600000   -1

!  tp0 = 0.35864311809759d1
!  x0 = (/ 0.d0,  -0.55680992271445d0,  0.11417028016985d1,  0.44424859433677d0,  0.d0,  0.d0 /) 

! 47th, i = 46, almost perpendicular to x-y plane -- no big different with 46-th....
!  0.38029658821541E+01  0.00000000000000E+00 -0.70640514618501E+00  0.12077964137697E+01  0.26294530968888E+00  0.00000000000000E+00  0.00000000000000E+00 -0.18922397316259E+01  0.00000000000000E+00  0.25600000   -1
  tp0 = 0.38029658821541d1
  x0 = (/ 0.d0, -0.70640514618501d0,  0.12077964137697d1,  0.26294530968888d0, 0.d0, 0.d0 /)
  
! 32-th, i=31, a big one in size, just for testing   
!    0.34858334376319E+01  0.00000000000000E+00  0.95416926423625E+00  0.54917304800717E+00  0.34667273261016E+00  0.00000000000000E+00  0.00000000000000E+00 -0.17863900149851E+01  0.11102230246252E-14  0.12800000    1
  tp0 = 0.34858334376319d1 
  x0 = (/ 0.d0,  0.95416926423625d0,  0.54917304800717d0,  0.34667273261016d0, 0.d0, 0.d0 /) 
 


  lt = tp0  * 1.d0 

! the fixed stepsize for the integration of the orbit
  h = lt / (nsm-1)  ! time interval 

 
! automatic step control for orbit integration --- use as reference for comparision 
  call plob(x0, 0.d0, lt, 6, 1, fob, 1, gr_lf, gr_cjlf,  xf) 


 ! ****************************** orbit sampling ******************************
  xl = 0.d0 ! to avoid mistake, for every orbit, reinitialize the state vector
  xmax = 2.d0 ! the maximum value to detect escape 
      
  call cpu_time(tst)    
 
!  subroutine plob_fxd(y0,t0, h0, nsm, xmax, ftag, yall,  deriv, gr_cj)
  call plob_fxd(x0, 0.d0, h, nsm, xmax, fpo_fft, xl,  gr_lf, gr_cjlf) 

  call cpu_time(tend)
  
  print*, 'Elapsed time for orbit sampling is ', tend-tst
  print * ,'Integration of orbit finished!'

! check the final state and stop time? 
  print*, 'check the final state'
  print*, nsm, xl(nsm,:)
  read*
    
    
!! use the components of the state vector from the integration with fixed stepsize    
  x = xl(:, 1) ! x
  y = xl(:, 2) ! y 
  z = xl(:, 3) ! z 
  
  vx = xl(:, 4) ! vx
  vy = xl(:, 5) ! vy 
  vz = xl(:, 6) ! vz 
     
  call cpu_time(tst)    
  call dtwofft(x, y,fft1,fft2, nsm) ! for x and y
  call dtwofft(z, vz,fft3,fft6, nsm) ! for z and vz
  call dtwofft(vx, vy,fft4,fft5, nsm) ! for vx and  vy
 
 


! compare the result by integration and the result recoveried from  FFT   
!! print  -- discard, I have no idea how to deal with the output of dtwofft
!  print*, 'For x:'
!  call prntft(fft1, nsm/2) 
!  read* 
!   
!   print*, 'For y:'
!  call prntft(fft2, nsm) 
!  read* 
!       
!  call dtwofft(z, a,fft3,fft4, nsm) ! for z and a
!  call cpu_time(tend)

!  print*, 'For z:'
!  call prntft(fft3, nsm) 
!  read* 
!   
 
  print*, 'Elapsed time for  2 * twofft is ', tend-tst
  
!   ! read*
! Extract real magnitude and frequency from fft1 returned by twofft, save them in seperate files   
  print*,'#********* for x **********'  
  write(*, *) '# ********* x *********'
 
 isft = 0 ! do not save all the frequencies and corresponding amplitude

!subroutine fqext(fft, nsm, lt, f_fft, nfq, fqmx)
  call fqext(fft1, nsm, lt, isft, fx_fft, nfq, fqmx_x) 
  
  write(*, *) '# ********* y *********'
  call fqext(fft2, nsm, lt, isft,  fy_fft, nfq, fqmx_y) 
  
  
  write(*, *) '# ********* z *********'
  call fqext(fft3, nsm, lt, isft, fz_fft, nfq, fqmx_z) 
  
  
  write(*, *) '# ********* vx *********'
  call fqext(fft4, nsm, lt, isft, fvx_fft, nfq, fqmx_vx)  
  
  write(*, *) '# ********* vy *********'
  call fqext(fft5, nsm, lt, isft, fvy_fft, nfq, fqmx_vy) 
  
  write(*, *) '# ********* vz *********'
  call fqext(fft6, nsm, lt, isft, fvz_fft, nfq, fqmx_vz) 

! use the first nfq frequencies that has the largest magnitude  to recover the Fourier series  
  do i =  1, nsm, 1
    t = (i-1) * h 
    call fourseri(t,  nfq,  fqmx_x, fourx )   !x
    call fourseri(t,  nfq,  fqmx_y, foury )   !y
    call fourseri(t,  nfq,  fqmx_z, fourz )   !z
    
    print*, 't, x, y, z',  t, x(i), y(i), z(i)
    print*, 'Fourier: ',   t, fourx, foury, fourz
    print*, 'dx=fourx-x',  fourx-x(i),  foury-y(i),  fourz-z(i) 
    read* 
  end do
  

  write(ffqmx, *) '#  x:  f - amp- cos - sin'
  do j = 1, nfq
    write(ffqmx,'(4f14.6)')   fqmx_x(j, :) 
    write(*, '(4f14.6)')      fqmx_x(j, :) 
  enddo  
  write(ffqmx, *);    print* 
  
  write(ffqmx, *) '#  y:  f - amp- cos - sin'
  do j = 1, nfq
    write(ffqmx,'(4f14.6)')   fqmx_y(j, :) 
    write(*, '(4f14.6)')      fqmx_y(j, :) 
  enddo  
  write(ffqmx, *);    print* 

  write(ffqmx, *) '#  z:  f - amp - cos - sin'
  do j = 1, nfq
    write(ffqmx,'(4f14.6)')   fqmx_z(j, :) 
    write(*, '(4f14.6)')      fqmx_z(j, :) 
  enddo  
  write(ffqmx, *);    print* 
  
  write(ffqmx, *) '#  vx:   f - amp - cos - sin'
  do j = 1, nfq
    write(ffqmx,'(4f14.6)')   fqmx_vx(j, :) 
    write(*, '(4f14.6)')      fqmx_vx(j, :) 
  enddo  
  write(ffqmx, *);    print* 
  
  write(ffqmx, *) '#  vy:   f - amp- cos - sin'
  do j = 1, nfq
    write(ffqmx,'(4f14.6)')   fqmx_vy(j, :) 
    write(*, '(4f14.6)')      fqmx_vy(j, :) 
  enddo  
  write(ffqmx, *);    print* 
  
    write(ffqmx, *) '#  vz:   f - amp - cos - sin'
  do j = 1, nfq
    write(ffqmx,'(4f14.6)')   fqmx_vz(j, :) 
    write(*, '(4f14.6)')      fqmx_vz(j, :) 
  enddo  
  write(ffqmx, *);    print* 
  
  
! Pick the first nfq=4 dominant frequencies and then save them in the form 
!! fx- ax - fy - ay - fz - az, keep 4 digits after the decimer, for convenient use in latex tabular

!  print*, 'fqmx:  fx- ax - fy - ay - fz - az ' 

! 
!  do j = 1, nfq
!    write(ffqmx,'(6f14.6)')   fqmx_x(j, 1:), fqmx_y(j, :), fqmx_z(j, :)
!    
!    write(*, '(6f14.6)')     fqmx_x(j, :), fqmx_y(j, :), fqmx_z(j, :) 
!    print* 
!  enddo    
  
  write(ffqmx,*) ! add a blank line to seperate different orbit with different distance
 read*


 close(fx_fft); close(fy_fft);  close(fz_fft);   close(ffqmx);
 
stop
end program fftpo1



















  

program main_fft 
!  this is to test fft, and explore which values of the time length and  number of samples are good for the extract of frequecies and magnitude

!  Time length:         LT 
!  Number of samples:   NP 
!  Time interval between 2 consecutive samples:  Delta = LT/ (NP-1)? or LT/NP
!  Sampling rate:       fs = 1 / Delta  
!  Nyquist frequency:   fc = 1/2 * fs  (half of the sampling rate)

!  Nyquist rate:   the lower bound for the sample rate for alias-free signal sampling

!    Assume the highest frequency component in x(t) is B. 
!    As a result, the Nyquist rate is R_N = 2*B

! Nyquist criteria:  fs > R_N, in other words, fc > Max(f) = B 

! So we could estimate the number of samples for a given time length LT, if we know roughly the value of B  
! fs = (NP-1) / LT > R_N = 2 * B      ===>      NP > 2 * LT * B + 1

! for the tested function, B = 27; and we need 1 digit for frequency, take LT = 10 
!    NP > 2 * 10 * 27 + 1 =  541,  with 2**10 = 1052,  2**9 = 512,  2**8 = 256

! For now, a rough conclusion: 
!  1.  the precision of the frequency is 1/ LT ( cycle per unit time), 
!      so for one digit, we could take LT =10,  and 2--100, 3--1000, etc. 

!  2.  If we want to recover the frequency in unit of  radian per unit time, we should multiply the above LT by 2*pi 
!
!  3.  when LT is taken an appropriate value to recover the frequency,  and NP satisfies Nyquist criteria
!      enlarge the value of LT doesn't make sense, 
!      we need to decrease the time interval, that means increase the number of sample points to better recover the magnitude. 
!      but there is an upper bound for the value of NP, the accuracy of the magnitude is not improved with vaule bigger than this upper bound

!  Que:  1- what is the value of the upper bound of NP ?

use fft_mod
use lf_mod

implicit none

! np = 2^21, lt = 2.d3  is the best match, we have x-z clearly extracted.


! number of sampling points on the orbit: 2**32, 2**24, 2**16 points should be good
integer, parameter ::   nsm = 2**20,  &  !number of sampling points 
                        nfq = 10   ! keep the first 4 frequencies with the biggest amplitude

real(kind=dp), parameter :: pi = 4.d0*datan(1.d0), pi2 = 2.d0*pi ! & !! time length of sampling 

! local variables
integer :: i, j, iha, kk, np,  fa_fft, fb_fft, ffqmx, fob_fft !, &
!	   isign, fob_fft, fx_fft, fy_fft, fz_fft  !  --fft

real(kind=dp) ::  lt,  h, t,  a1, a2, a3,  a(nsm), b(nsm),  &
                  tst, tend, &  ! evaluate time 
                  fqmx_a(nfq+1,2), fqmx_b(nfq+1,2)  , & !the domimant frequency which has the first nfq maximum amplitude
                  feqm1(nsm/2+1, 4 ), feqm2(nsm/2+1, 4 ), feqm3(nsm/2+1, 4 )

! character(len=70) ::  fnmf, fna, fnb, fnmx  !  fnx, fny, fnz,

real(kind=dp), external :: dnrm2 ! from BLAS-LEVEL1

! fx_fft = 30; fy_fft = 31;  fz_fft = 32;  ! save x-y-z and
fa_fft = 33;   fa_fft = 34 ! test function a and b 
fob_fft = 35;  ffqmx = 36  ! save the dominant magnitude and frequency

! open an empty file using the fname
open(fob_fft, file= './dat/fft/ob.dat' , access ='append',status='replace')  
open(fa_fft, file= './dat/fft/a.dat' , access ='append',status='replace')
open(fb_fft, file= './dat/fft/b.dat' , access ='append',status='replace')
open(ffqmx,  file= './dat/fft/fqmx.dat', access ='append',status='replace')  
!write(ffqmx,*) '# fx	ax 	fy 	ay 	fz 	az'

! ---------------- test function sampling-------------- 
! test f(t) = cos(pi2*0.13*t) - 1/2*sin(pi2*0.27*t) + 3/4 *sin(pi2*0.41*t)

do kk  =  1, 100 

  print*, 'Time Length,  Number of samples (log_2):  '
  read(*,*) lt, np, iha 

  
  np = 2 ** np
  h = lt / 

  print*,'lt, h, np: ',  lt, h, np  
  print*; !read*

  open(fob_fft, file= './dat/ob_fft.dat' ) !, access ='append',status='replace')  

  a1 = 1.8d0 ; ! 10
  a2 = 0.d0 ! -2.d0;  ! 27
  a3 = 11.25d0 ! 4.1

! --------- sampling ---------
  call cpu_time(tst) 
  do i = 1, np ! --tested! result is fine!
    t =  (i-1) * h     !* 2 * pi ! we have 2 cycles
  
    a(i)  = a1 * dcos(pi2*10d0*t) + a2*dsin(pi2*27d0*t) + a3 *dsin(pi2*4.1d0*t) 
    b(i)  = 3.78d0+ a1 * dcos(pi2*10d0*t) + a2*dsin(pi2*27d0*t) + a3 *dsin(pi2*4.1d0*t) 
 
    write(fob_fft, *)  t, a(i), b(i)
  
!  a(i)  = 1.8d0 * dcos(10d0*t) - 2.d0/1.d0*dsin(27d0*t) + 45.d0/4.d0 *dsin(4.1d0*t);
!  b(i)  = 1.8d0 * dcos(10d0*t) - 2.d0/1.d0*dsin(27d0*t) + 45.d0/4.d0 *dsin(4.1d0*t);
!  b(i)  = dcos(pi2*1.3d0*t) - 1.d0/2.d0*dsin(pi2*4.1d0*t) + 5.d0/4.d0 *dsin(pi2*6.7*t);
  end do
!----------------------------------

  close(fob_fft )

  call  hannin(b,np, iha)

  call cpu_time(tend)
  print*, 'Elapsed time for sampling is ', tend-tst!;  print*
 

  ! -------- Fourier analysis -----------  
  call dtwofft(a, b, fft1, fft2, np) ! for a and b
  call cpu_time(tst)
  print*, 'Elapsed time for dtwofft ',  tend-tst; print* !; read*

  
  write(*, *) '# ********* a: f -- amp -- cos -- sin *********'
  call mod_twofft( fft1, np,  lt, feqm1)
  call fqext(feqm1, np/2+1, lt, 0, fa_fft, nfq, fqmx_a) 
  
  write(*, *) '# ********* b: f -- amp -- cos -- sin *********'
  call mod_twofft( fft2, np,  lt, feqm2)
  call fqext(feqm2, np/2+1, lt, 0,  fb_fft, nfq, fqmx_b) 
  
  print*, 'finish dtwofft!'; print*
!enddo 


! ********************** test fftsc provide by Gerard ********************
! the result is the same with the one obtained by dtwofft,  except the data process for modulus is a little bit different

print*, '******** test fftsc*************';

call cpu_time(tst)
! SUBROUTINE FFTSC  (A,N,ST,CT,IWK,WK,CWK)
call FFTSC(b, np, st, ct, iwk, wk, cwk )
call cpu_time(tend)
print*, 'Time elapsed for fftsc is: ', tend-tst; print*; read*! 

call mod_fftsc(st, ct, np, lt, feqm3)
call fqext(feqm3, np/2+1, lt, 0,  fb_fft, nfq, fqmx_b) 

enddo 
 
stop
end program main_fft





  

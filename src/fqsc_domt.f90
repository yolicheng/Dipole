subroutine fqsc_domt(st, ct, nsm, lt, isft, f_fft, nfq,  fqmx)
!  this subroutine is extract the dominant frequency and amplitude from the data returned by fftsc 
!  and save them to file ftag if isft == 1, and return the first nfq dominant frequencies, to be written to file later

!
! REMARKS  1.  FFTSC COMPUTES THE SINE TRANSFORM, ST, ACCORDING
!              TO THE FOLLOWING FORMULA;
!
!                ST(K+1) = 2.0 * SUM FROM J = 0 TO N-1 OF
!                          A(J+1)*SIN(2.0*PI*J*K/N)
!                FOR K=0,1,...,N/2 AND PI=3.1415...
!
!              FFTSC COMPUTES THE COSINE TRANSFORM, CT, ACCORDING
!              TO THE FOLLOWING FORMULA;
!
!                CT(K+1) = 2.0 * SUM FROM J = 0 TO N-1 OF
!                          A(J+1)*COS(2.0*PI*J*K/N)
!                FOR K=0,1,...,N/2 AND PI=3.1415...
!          2.  THE FOLLOWING RELATIONSHIP EXISTS BETWEEN THE DATA
!              AND THE COEFFICIENTS OF THE SINE AND COSINE TRANSFORM
!
!                A(J+1) = CT(1)/(2*N) + CT(N/2+1)/(2*N)*(-1)**J +
!                         SUM FROM K = 1 TO N/2-1 OF
!                          (CT(K+1)/N*COS((2.0*PI*J*K)/N) +
!                           ST(K+1)/N*SIN((2.0*PI*J*K)/N))
!                FOR J=0,1,...,N-1 AND PI=3.1415...
!


! the structure of fqmx 
!  1st row:             0, a0, a0, 0 
!  the other rows:      f, amp, real, imag
! 	Input
! nsm 	:	number of sampling points
! fft 	:	complex variable from twofft
! lt	: 	length of signal, used to compute frequency
! isft  :       flag to indict to save all the results from dtwofft or not ? 
! f_fft :	file tag to save the frequency and amplitude for fft1.dat
! nfq	:	number of  frequencies to keep, we only consider the number less than 10(included) 
!                       normally we take nfq = 4

!      Onput
! fqmx(nfq, 4):  the selected frequency(the first column), and amplitude(the second column) - real -imag
! 

implicit none 
integer, parameter:: dp = kind(1.d0) 
real(kind=dp), parameter :: pi = 4.d0*datan(1.d0)


integer, intent(in) ::  nsm, isft, f_fft, nfq    
real(kind=dp), intent(inout)  :: st(nsm/2+1), ct(nsm/2+1)

real(kind=dp), intent(in)   :: lt 
real(kind=dp), intent(out)  ::  fqmx(nfq+1, 4) 

! local variables
integer :: i, j
real(kind=dp) :: f, amp, iamp, ramp, x_curr(4), fqmxm1(nfq, 4) 


 CHARACTER(LEN=*), PARAMETER  :: fmt  = "(4f22.12)" ! the format for magnitude output
real(kind=dp), external :: dnrm2 ! from BLAS-LEVEL1

! the initial value should be zero
fqmxm1 = 0.d0  

do i = 1, nsm/2+1

! radian per unit time ! without 2*pi, it is cycle per unit time
!   f = (i-1) / lt  * pi * 2 ! ! for twofft, amp = amp/(nsm/2), f = f/2, but if t is the real time, f = f*2*pi/2 = f*pi
  
! I think it is better to take the general f, without multiplied by 2*pi, meaning: cycle per unit time  
   f =  (i-1) / lt 
   
! is this the right way to get the frequency?! Yes! frequency is the number of sampling per unit time 
 
! the amplitude, the real and imaginary part, which are the coefficient of cos and sin respectively.     
  ct(i) = ct(i) / nsm 
  st(i) = st(i) / nsm  
  
  amp =  dsqrt( ct(i)**2 +  st(i) **2 )
    
  if(i == 1 .or. i == nsm/2+1)  then  ! for f=0, the constant shift should be devided by 2
!    print*, 'the imaginary should be zero, i, iamp'!ck
!    print*, i-1, iamp ; print*; read*
    
    amp   = amp  / 2.d0
    ct(i) = ct(i) / 2.d0 
    
    ! save the costant part as the first row
    if(i == 1) fqmx(1, :) = (/f, amp, ct(i), st(i)/)
  endif
    
  if(isft == 1) write(f_fft, fmt)  f, amp,  ct(i), st(i)
!  write(*, fmt)  f, amp,  ramp, iamp  ! frequency domain + phase domain(polar radius + phase angle)
   
  if (f .lt. 1.d-5)  cycle ! the constant term, corresponding to f = 0 
  
  x_curr = (/f, amp, ct(i), st(i)/)
  call fqmax(x_curr, nfq, fqmxm1)
  
enddo 

! the second and latter rows are the frequencies with largest amplitude  
fqmx(2:nfq+1, :) = fqmxm1

 print*, nfq, 'dominant frequencies';  print*; 
  do j = 1, nfq 
    write(*,'(4f20.10)') fqmxm1(j, :)
  enddo 
  print* 

! add a blank line to seperate orbit as block
if(isft == 1)  write(f_fft, *) 
 
return
end    
    
    

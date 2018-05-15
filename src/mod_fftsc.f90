subroutine mod_fftsc(st, ct, nsm, lt, feqm)

!subroutine fqsc_domt(st, ct, nsm, lt, isft, f_fft, nfq,  fqmx)

!  this subroutine transform the output of fftsc into real frequencies and magnitude, and the coefficients of cosine and sine 

! the structure of fqmx 
! 	Input
! ct, st :      coefficients returned by fftsc 
! nsm 	:	number of sampling points
! fft 	:	complex variable from twofft
! lt	: 	length of signal, used to compute frequency

!      Onput
! feqm(nsm/2, 4):  frequency-amplitude- real -imag


implicit none 
integer, parameter:: dp = kind(1.d0) 
real(kind=dp), parameter :: pi = 4.d0*datan(1.d0)

integer, intent(in) ::  nsm   
real(kind=dp), intent(inout) :: st(nsm/2+1), ct(nsm/2+1)


real(kind=dp), intent(in)    ::  lt 
real(kind=dp), intent(out)   ::  feqm(nsm/2+1, 4) 

! local variables
integer :: i 
real(kind=dp) :: f, amp, iamp, ramp 


! question: what is the point of clain the nsm/2+1 -th frequency
! by test, it seems not necessary for  the nsm/2+1 -th frequency... 

do i = 1, nsm/2+1 

! radian per unit time ! without 2*pi, it is cycle per unit time
!   f = (i-1) / lt  * pi * 2 ! ! for twofft, amp = amp/(nsm/2), f = f/2, but if t is the real time, f = f*2*pi/2 = f*pi
  
! I think it is better to take the general f, without multiplied by 2*pi, meaning: cycle per unit time  
   f =  (i-1) / lt 
   
  ct(i) = ct(i) / nsm 
  st(i) = st(i) / nsm  
  
  amp =  dsqrt( ct(i)**2 +  st(i) **2 )
   
! the amplitude, the real and imaginary part, which are the coefficient of cos and sin respectively. 
  if(i == 1 .or. i == nsm/2+1)  then  ! for f=0, the constant shift should be devided by 2
    print*, 'the imaginary should be zero, i, ramp, iamp'!ck
    print*, i, ct(i), st(i) ; print*; read*
    
    amp   = amp   / 2.d0
    ct(i) = ct(i) / 2.d0 
  endif 
        
    
  ! save the costant part as the first row
  feqm(i, :) = (/f, amp, ct(i), st(i)/)
enddo 
 
return
end    
    
    

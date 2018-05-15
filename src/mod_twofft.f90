subroutine mod_twofft(fft, nsm, lt, feqm)
!  this subroutine transform the output of dtwofft into real frequencies and magnitude, and the coefficients of cosine and sine 

! the structure of fqmx 
! 	Input
! nsm 	:	number of sampling points
! fft 	:	complex variable from twofft
! lt	: 	length of signal, used to compute frequency

!      Onput
! feqm(nsm/2+1, 4):  frequency-amplitude- real -imag


implicit none 
integer, parameter:: dp = kind(1.d0) 
real(kind=dp), parameter :: pi = 4.d0*datan(1.d0)

integer, intent(in) ::  nsm   
double complex, intent(in)  :: fft(nsm)

real(kind=dp), intent(in)   ::  lt 
real(kind=dp), intent(out)  ::  feqm(nsm/2+1, 4) 

! local variables
integer :: i 
real(kind=dp) :: f, amp, iamp, ramp 


do i = 1, nsm/2+1 

! radian per unit time ! without 2*pi, it is cycle per unit time
!   f = (i-1) / lt  * pi * 2 ! ! for twofft, amp = amp/(nsm/2), f = f/2, but if t is the real time, f = f*2*pi/2 = f*pi
  
! I think it is better to take the general f, without multiplied by 2*pi, meaning: cycle per unit time  
   f =  (i-1) / lt 
   
! is this the right way to get the frequency?! Yes! frequency is the number of sampling per unit time 
 
! the amplitude, the real and imaginary part, which are the coefficient of cos and sin respectively.     
  amp =  cdabs( fft(i) )/ nsm * 2.d0
  iamp = dimag( fft(i))  / nsm * 2.d0
  ramp = dble( fft(i))   / nsm * 2.d0 
      
  if(i == 1 .or. i == nsm/2+1) then  ! for f=0, the constant shift should be devided by 2
    amp  = amp  / 2
    ramp = ramp / 2 
    print*, 'i, ramp, iamp', i, ramp, iamp
  endif
    
  ! save the costant part as the first row
  feqm(i, :) = (/f, amp, ramp, iamp/)
enddo 
 
return
end    
    
    

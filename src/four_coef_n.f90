subroutine four_coef_n( nf, fft, nsm, lt,  fcoef ) 
!discard... it is of no use to use the first few frequencies without sorting...
 
! This routine is to recover the first nf terms of Fourier series using the complex output of dtwofft 

! F = a0 + Sigma _from 1 to n { real * cos(2*pi*fi * n * t) + imag * sin(2*pi*fi * n * t) }

!       Input Variables 
!  fft          the output of dtwofft 
!  nsm          number of rows in array fft
!  nf           number of terms in Fourier we are going to keep 


!       Output Variables 
!  fourt        the value of Fourier series at time t with the first nf terms           

!  Routine Used:
!     

! Finally Revised by Yu -- 20160530
!----------------------------------------

implicit none
integer, parameter :: dp = kind(1.d0)
  
! Input  and Output Declaration   
integer, intent(in)        ::  nsm, nf 
double complex, intent(in)        ::  fft(nsm) 

real(kind=dp), intent(in)       ::  lt   
real(kind=dp), intent(out)        ::  fcoef(nf, 5) 
 
! Local Variable
real(kind=dp)  :: i, amp, iamp, ramp, f  
  
 
do i = 1, nf
  
! I think it is better to take the general f, without multiplied by 2*pi, meaning: cycle per unit time  
  f =  (i-1) !/ lt 
   
! is this the right way to get the frequency?! Yes! frequency is the number of sampling per unit time 
 
! the amplitude, the real and imaginary part, which are the coefficient of cos and sin respectively.     
  amp =  cdabs(fft(i))/ nsm * 2.d0
  ramp = dble(fft(i))/ nsm * 2.d0 
  iamp = dimag(fft(i))/ nsm * 2.d0 
  
    
  if(i == 1) then  ! for f=0, the constant shift should be devided by 2
    amp  = amp / 2
    ramp = ramp/ 2  
  endif
    
  fcoef(i, :) = (/i*1.d0, f, amp, ramp, iamp/)  
  print*, 'The i-th row in fcoef', fcoef(i,:) 
  read* 
enddo 

return  
end subroutine four_coef_n


subroutine fourseri(t,  nfq, fcoef, fourt )

! This routine is to recover the Fourier series using the first nfq frequencies 

! f(t) = a0 + Sigma _ i from 1 to nfq { real * cos(2*pi*f(i)  * t) + imag * sin(2*pi*f(i) * t) }

!       Input Variables 
!  t           epoch 
!  tp          the period 
!  nfq         number of terms in Fourier we are going to keep 
!  fcoef       the coefficients of the corresponding frequencies  n-f-amp-real-imag

!       Output Variables 
!  fourt        the value of Fourier series at time t with the first nfq terms           

!  Routine Used:
!    None   

! Finally Revised by Yu -- 20160530
!----------------------------------------

implicit none
integer, parameter :: dp = kind(1.d0)
real(kind=dp), parameter :: pi2 = 4.d0*datan(1.d0)

! Input  and Output Declaration   
integer, intent(in)             ::  nfq 
real(kind=dp), intent(in)       ::  t,  fcoef(nfq, 4) 
real(kind=dp), intent(out)      ::  fourt 
 
! Local Variable
integer :: i 
real(kind=dp)  ::  iamp, ramp, f, n   
  
! initialization 
fourt = 0.d0 
 
do i = 1, nfq

  f = fcoef(i, 1)
  ramp = fcoef(i, 3) 
  iamp = fcoef(i, 4)
      
  if(i == 1) then  ! for f = 0, the constant shift should be devided by 2
    
    fourt = fourt + ramp 
    
  else 
  
    fourt = fourt + ramp * dcos( 2.d0 * pi * f * t  ) + iamp * dsin( 2.d0 * pi * f  * t )
      
  endif
enddo 

!print*, 't, fourt', t, fourt  !ckd
!read* 
  
return  
end subroutine fourseri


!***********************************************************************
!  Evaluate a periodic function or signal using a truncated Fourier series with the first m terms, 
!  the Fourier coefficients csf(k) and sif(k), k = 0, ..., m are computed by fourier.f 
!    the period is 1 in time, 2pi in radian 

! -- 1. as a function of time t,  ---  flag = 2
!       f(t) = csf(0) + Sigma for k=1, ..., m of { csf(k) * cos( 2pi* k * t) +  ss(k+1) * sin( 2pi* k * t) }

! -- 2. in terms of the argument of the evaluated point, of which the value range is [0,2pi] -- the normally used one, flag = 1
!       since at time t the argument theta = 2*pi*t, we have

!       f(theta) = csf(0) + Sigma for k=1, ..., m of { csf(k)) * cos( theta * k ) +  sif(k) * sin( theta * k ) }   

!       Input Variables 
!  phi      the parameter according to which the point is evalutated using Fourier coefficients
!            time (flag = 1) or  angle (flag = 2) 
!  flag     specify the type of the input  phi--- discard, since we can process the input before we call four_seri 
!  m        number of Fourier modes
!  csf      cosinus coefficients for cosine(k*  phi ), k =0,...,m
!  sif      sinus  coefficients for  sine(k* phi), k =0,...,m 

!  Routine Used:
!    trigrec_c, trigrec_s from grlib.a 

!  Finally Revised by Yu -- 20160816
!***********************************************************************
function four_seri( phi, m, csf, sif)  result(fun)

use dp_mod
use pi_mod
implicit none

! Input  and Output Declaration
integer, intent(in)               :: m
real(kind=dp), intent(in)         :: phi 
real(kind=dp), dimension(0: m), intent(in)   ::  csf, sif 
real(kind=dp)      ::  fun
 
! Local Variable
real(kind=dp) :: theta, cn(m), sn(m), &
                 cts(m+1), sts(m+1), fun2 ! for test when theta > pi  
integer :: i
  
!  if( flag == 1 ) then
!    theta = phi
!  elseif(flag == 2) then
!    theta = pi2 * phi
!  end if 
   
  fun = csf(0) ! the constant term 

  ! cos(k*phi) and sin(k*sin) computed recurrently
  call trigrec_c( phi, m, cn)  
  call trigrec_s( phi, m, sn) 
  
  do i = 1, m, 1
    fun = fun + csf(i) * cn(i) + sif(i) * sn(i)
  end do
  
  return  
  
end function four_seri 



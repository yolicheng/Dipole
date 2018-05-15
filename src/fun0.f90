!***********************************************************************
!  Evaluate a periodic function or signal using a truncated Fourier series with the first m terms, 
!  the Fourier coefficients cof(k) and sif(k), k = 0, ..., m are computed by fourier.f 
!    the period is 1 in time, 2pi in radian 

! -- 1. as a function of time t,  --- the normally used one, flag = 1
!       f(t) = cof(0) + Sigma for k=1, ..., m of { cof(k) * cos( 2pi* k * t) +  ss(k+1) * sin( 2pi* k * t) }

! -- 2. in terms of the argument of the evaluated point, of which the value range is [0,2pi] -- flag = 2
!       since at time t the argument theta = 2*pi*t, we have

!       f(theta) = cof(0) + Sigma for k=1, ..., m of { cof(k)) * cos( theta * k ) +  sif(k) * sin( theta * k ) }   

!       Input Variables 
!  phi       the time (flag = 1) or the angle (flag = 2)of the point to be evalutated using Fourier coefficients
!  flag     specify the type of the input t 
!  m        number of Fourier modes
!  cof      cosinus coefficients for cosine(k*  phi ), k =0,...,m
!  sif      sinus  coefficients for  sine(k* phi), k =0,...,m 

!  Routine Used:
!    trigrec_c, trigrec_s from grlib.a 

! TODO: test if this is the right expression of Fourier series 

!  Finally Revised by Yu -- 20160815  --hope everything is fine now...          

!  Routine Used:
!     

!  Finally Revised by Yu -- 20160816
!***********************************************************************
double precision function four_seri( phi, flag, m, cof, sif), result fun

implicit none
integer, parameter        :: dp = kind(1.d0)   
real(kind=dp), parameter  :: pi2 = 8.d0*datan(1.d0)

! Input  and Output Declaration
integer, intent(in)               :: flag, m
real(kind=dp), intent(in)         :: phi 
real(kind=dp), dimension(0: m), intent(in)   ::  cof, sif 
 
! Local Variable

real(kind=dp) :: theta, cn(m), sn(m)  
integer :: i
  
  if( flag == 1 ) then
    theta = pi2 * phi
  elseif(flag == 2) 
    theta = phi
  end if 
  
   
  fun = cof(0) ! the constant term 
  
  ! cos(k*phi) and sin(k*sin) computed recurrently
  call trigrec_c( theta, m, cn)  
  call trigrec_s( theta, m, sn) 

  do i = 1, m, 1
    fun = fun + cof(i) * cn(i) + sif(i) * sn(i)
  end do

  return  
  
end function four_seri 

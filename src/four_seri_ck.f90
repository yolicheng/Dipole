!***********************************************************************
!  Evaluate a periodic function or signal using a truncated Fourier series with the first m terms, 
!  the Fourier coefficients cof(k) and sif(k), k = 0, ..., m are computed by fourier.f 
!    the period is 1 in time, 2pi in radian 

! -- 1. as a function of time t,  ---  flag = 2
!       f(t) = cof(0) + Sigma for k=1, ..., m of { cof(k) * cos( 2pi* k * t) +  ss(k+1) * sin( 2pi* k * t) }

! -- 2. in terms of the argument of the evaluated point, of which the value range is [0,2pi] -- the normally used one, flag = 1
!       since at time t the argument theta = 2*pi*t, we have

!       f(theta) = cof(0) + Sigma for k=1, ..., m of { cof(k)) * cos( theta * k ) +  sif(k) * sin( theta * k ) }   

!       Input Variables 
!  phi      the parameter according to which the point is evalutated using Fourier coefficients
!            time (flag = 1) or  angle (flag = 2) 
!  flag     specify the type of the input  phi
!  m        number of Fourier modes
!  csf      cosinus coefficients for cosine(k*  phi ), k =0,...,m
!  sif      sinus  coefficients for  sine(k* phi), k =0,...,m 

!  Routine Used:
!    trigrec_c, trigrec_s from grlib.a 

! TODO: test if this is the right expression of Fourier series 

!  Finally Revised by Yu -- 20160815  --hope everything is fine now...          

!  Routine Used:
!     

!  Finally Revised by Yu -- 20160816
!***********************************************************************
function four_seri( phi, flag, m, csf, sif)  result(fun)

use dp_mod
use pi_mod
implicit none

! Input  and Output Declaration
integer, intent(in)               :: flag,  m
real(kind=dp), intent(in)         :: phi 
real(kind=dp), dimension(0: m), intent(in)   ::  csf, sif 
real(kind=dp)      ::  fun
 
! Local Variable
real(kind=dp) :: theta, cn(m), sn(m), &
                 cts(m+1), sts(m+1), fun2 ! for test when theta > pi  
integer :: i
  
  if( flag == 1 ) then
    theta = phi
  elseif(flag == 2) then
    theta = pi2 * phi
  end if 
   
  fun = csf(0) ! the constant term 
  
  !check the constant term, csf, sif, m and flag, phi 
!  print*, 'constant term: ', fun 
!  print*, 'm=', m, flag 

  ! cos(k*phi) and sin(k*sin) computed recurrently
  call trigrec_c( theta, m, cn)  
  call trigrec_s( theta, m, sn) 
  
  
  do i = 1, m, 1
    if(theta > 3.d0)then 
    print*,  'i, theta', i, theta
    print*, 'check cos(i*theta), sin(i*theta) by intrinsic function'
    print*,  dcos(i*theta), dsin(i*theta)
    print*, 'check by my own routine trigrec_c and trigrec_s:'
    print*, cn(i), sn(i)
    read*
    endif
    fun = fun + csf(i) * cn(i) + sif(i) * sn(i)
  end do
  
 ! *** check here why for theta > pi, there is something wrong ***
 !      1.0000123759    0.4998386448    0.0000140138    0.0000039269    0.0000019001    0.0000011377    0.0000007637    0.0000005509    0.0000004177    0.0000003286    0.0000002659
!    0.0000000000    0.2499162099    0.0000056190    0.0000011803    0.0000004471    0.0000002182    0.0000001231    0.0000000763    0.0000000506    0.0000000353    0.0000000256
  cts = (/1.0000123759d0,  0.4998386448d0,  0.0000140138d0,  0.0000039269d0,  0.0000019001d0,  0.0000011377d0, &
          0.0000007637d0, 0.0000005509d0,  0.0000004177d0,  0.0000003286d0,  0.0000002659d0/)
  sts = (/0.0000000000d0,  0.2499162099d0,   0.0000056190d0,  0.0000011803d0,  0.0000004471d0,  0.0000002182d0, &
          0.0000001231d0,  0.0000000763d0,   0.0000000506d0,  0.0000000353d0,  0.0000000256d0/)

  if(theta > pi) then 
     ! check the coefficients
     print*, 'csf, sif'
     print*, csf 
     print*, sif 
     
     print*, 'cts, sts'
     print*, cts
     print*, sts
     read*
     fun2 = cts(1)
     do i = 1, m, 1
       fun2 = fun2 + cts(i+1)*dcos(i*theta) + sts(i+1)*dsin(i*theta)
     end do
     print*, 'fun, fun2', fun, fun2 
     read*
  endif 
  
  return  
  
end function four_seri 

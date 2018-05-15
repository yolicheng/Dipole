!***********************************************************************
! This function  evaluate the values of coordinates by Fourier series for a given frequency w  
!   y  = c(0) + Sigma for k=1,...,m of {c(k)*cos(k* t*2pi/w) + ss(k)* sin(k* t*2pi/w) }

!       Input Variables 
!  t        
!  m        number of Fourier modes
!  cs       coefficients of cosine(k*darg), k =1,...,m 
!  ss       coefficients of sine(k*darg), k =1,...,m             
!  w        the frequency 

!  Routine Used:
!    trigrec_c, trigrec_s from grlib.a 

! TODO: what is the parameter t??? 

!  Finally Revised by Yu -- 20160808
!***********************************************************************
double precision function fun( t, m, cs, ss, w)

!use fft_mod

implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration
integer, intent(in)          :: m 
real(kind=dp), dimension(0:11), intent(in)      ::  cs, ss 
real(kind=dp), intent(in)        ::  t, w
 
! Local Variable
real(kind=dp), parameter  :: pi2 = 8.d0*datan(1.d0)
real(kind=dp) ::  darg, cn(m), sn(m)  
integer :: i
    
  fun = cs(0) ! the constant term 
  
  darg = t*pi2/w
   
  call trigrec_c( darg, m, cn)  
  call trigrec_s( darg, m, sn) 

  do i = 1, m, 1
    fun = fun + cs(i) * cn(i) + ss(i) * sn(i)
  end do

  return  
  
end function fun


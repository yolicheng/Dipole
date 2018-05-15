!***********************************************************************
! This function  evaluate the values of coordinates by Fourier series for 
! give coefficients  csf(k) and sif(k)
!  
! the terms  related to the frequecies and argument phi should be: 
!    (sin cos) ( k * 2pi /prd * phi  ), where prd  = 2pi by default 
!   
! So finally we have the Fourier series
!   y(phi)  = c(0) + Sigma for k=1,...,m of {csf(k)* cos( k * phi) + sif(k)* sin( k * phi) }

!       Input Variables 
!   phi      angle of the point to be evalutated using Fourier coefficients
!!  rho      rotation number of the periodic curve, unit: cycle per return 
!!           so the period of the invariant curve is what
!            so we need 1/rho returns to finish 1 cycle of the curve   

!  ** NOTE **   prd = 2*pi 
!   answer: depend on the with respect to which parameter we define the  period, keep coherent w.r.t. the unit of the rotation number 
!           if we take angle in radians,  it is 2pi  --- at this moment, we take this one, it is more straightforward. 
!           if we take rotation number in cycles, it is 1  

!  m        number of Fourier modes
!  csf      cosinus coefficients for cosine(k*  phi ), k =1,...,m
!  sif      sinus  coefficients for  sine(k* phi), k =1,...,m 

!  Routine Used:
!    trigrec_c, trigrec_s from grlib.a 


!  Finally Revised by Yu -- 20160815  --hope everything is fine now...
!***********************************************************************
double precision function funcs( phi, m, csf0, csf, sif)

!use fft_mod

implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration
integer, intent(in)               :: m 
real(kind=dp), intent(in)         :: phi, csf0, csf(m), sif(m)
 
! Local Variable
real(kind=dp), parameter  :: pi2 = 8.d0*datan(1.d0)
real(kind=dp) ::   cn(m), sn(m)  
integer :: i
    
  funcs = csf0 ! the constant term 
  
  ! cos(k*phi) and sin(k*sin) computed recurrently
  call trigrec_c( phi, m, cn)  
  call trigrec_s( phi, m, sn) 

  do i = 1, m, 1
    funcs = funcs + csf(i) * cn(i) + sif(i) * sn(i)
  end do

  return  
  
end function funcs


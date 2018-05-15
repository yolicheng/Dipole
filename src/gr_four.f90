!****************************************************************
!  Computes the Fourier series of a periodic function F evaluated at N
!  equally spaced points in a period [0, 2*pi]

!  -- the other routine by Gerard:   gr_foun in fourier.f  
!  -- index start from 0, to keep coherent with gr_foun in fourier.f   

!  It gives:
!   csf0, CSF(i) i=1, M containig the cosinus coefficients.
!   SIF(i) i=1, M containig the sinus coefficients.

!     both of the term (SIN COS)( i * x ),  and M is the number of harmonics

! The minimum value of N is N = 2*m +1
!    since we have 2*m+1 unknowns: cof(0_m) and sif(0_m), so we need at least 2*m+1 equations(points) 
!    that is why in fourier.f, the explanation says N has to be odd 


! Assume we have n points ( theta(n) and f(n) ), then 
!   theta_k = (k-1) * 2pi /n, k = 1, ..., n 

!   csf0 =  1 / N * Sigma for k=1_n of { f(theta_k) }

! And for j = 1,...,m, we have

!   csf(i) =  2 / N * Sigma for k = 1_n of { f(theta_k) * cos(j* theta_k) }
!   sif(i) =  2 / N * Sigma for k = 1_n of { f(theta_k) * sin(j* theta_k) }

! So it is convenient to use j* 2pi /n as the basic theta for each mode,
! and call trigrec_c(s) to compute k*theta, k = 1,...,n-1


! -- Finally, the evaluated functions can be given by the inverse of Fourier transformation:

!   F(x) = csf0 + Sigma  i = 1_M { CSF(i)* cos(i * x) + SIF(i)* sin(i * x) } 


!       Input Variables 
!  f        dimension n, the n sample points equally spaced in [0, 2pi]          
!  n        number of sample points 
!  m        number of harmonics to compute
!  
!       Output Variables 
!  csf0, csf     cosinus coefficients,  to keep coherent with the ones in curve_mod, save c0, ck separately. 
!  sif           sinum coefficients,  

! 
!  Routine Used:
!     trigrec_c, trigrec_s  

!  Finally Revised by Yu -- 20160816
!***********************************************************************
subroutine gr_four(f, n, m, csf0, csf, sif)

use dp_mod
use pi_mod  ! for pi and pi2 
implicit none 

! Input  and Output Declaration   
integer, intent(in)            ::  n, m
real(kind=dp), intent(in)      ::  f(n)  
real(kind=dp), intent(out)     ::  csf0, csf(m), sif(m)   
 
! Local Variable
integer ::   i, k, k_modn
real(kind=dp)  ::  dtheta, cn(n), sn(n),  fpc, fps, &
                   cnm1(n-1), snm1(n-1)  ! for k * theta, k=1,...,n-1

!real(kind=dp)  ::  theta, fsum ! check
                  
  dtheta = pi2 / n 
  
  ! compute the basic angle for different Fourier modes j, that is j*theta
!  call trigrec_c(theta, m, cm_bas) 
!  call trigrec_s(theta, m, sm_bas)
!  
  cn(1) = 1.d0;  sn(1) = 0.d0    ! j=0, cos(0) -- sin(0)
  
  ! -- compute csf(0) =  1 / N * Sigma for k=1_n of { f(theta_k) }
  csf0  = sum(f) / n ! the costant term
  
  print*, 'the costant term', csf0
    
!  fsum = 0.d0  ! -- ckd, the same as csf(0), so keep only the above array operation 
!  do i = 1, n, 1
!    fsum = fsum + f(i)
!  end do
!  print*, 'csf(0) computed component-by-component', fsum/n   !ckd 
!  read*

  ! the angles are of j*k * dtheta, so we only need to compute i*dtheta, i = 0,...,n-1 
  ! if we do mod(j*k, n)
   
  ! k = 1,...,n-1
  call trigrec_c(dtheta, n-1, cnm1) 
  call trigrec_s(dtheta, n-1, snm1) 
  
  ! k = 0,...,n, the first element are the same for all the values of theta 
  cn(2:n) = cnm1
  sn(2:n) = snm1
    
  ! --- for m Fourier modes ---  since 
  do i = 1, m, 1
    ! comptute the trigometric functions for  k * theta
    ! theta  = i * dtheta 
    
    ! Initialize f-by-cos(k*theta) to 0 
    fpc = 0.d0;     fps = 0.d0 
    
    do k = 1, n, 1
      k_modn = mod(i*(k-1), n) + 1 
      
!      print*, 'i,k, i*k, mod(i*k,n)', i,k, i*k, mod(i*k,n)!ck
      
      fpc = fpc + f(k) * cn(k_modn)
      fps = fps + f(k) * sn(k_modn)
    end do
    
!    read*
    
    csf(i) = fpc * 2.d0 / n 
    sif(i) = fps * 2.d0 / n
     
  end do  

  return  
end subroutine gr_four 



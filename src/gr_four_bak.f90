!****************************************************************
!  Computes the Fourier series of a periodic function F evaluated at N
!  equally spaced points in a period [0, 2*pi]

!  X_k = k* (PRD / N),  k=0..N-1 (N must be odd). TODO: why N must be odd? 

!  It gives:
!   CSF(i) i=0, M containig the cosinus coefficients.
!   SIF(i) i=0, M containig the sinus coefficients.

!     both of the term (SIN COS)( i * x ),  M is the number of harmonics


! Assume we have n points ( theta(n) and f(n) ), then 
!   theta_k = (k-1) * 2pi /n

!   csf(0) =  1 / N * Sigma for k=1_n of { f(theta_k) }

! And for j = 1,...,m, we have

!   csf(i) =  2 / N * Sigma for k = 1_n of { f(theta_k) * cos(j* theta_k) }
!   sif(i) =  2 / N * Sigma for k = 1_n of { f(theta_k) * sin(j* theta_k) }

! So it is convenient to use j* 2pi /n as the basic theta for each mode,
! and call trigrec_c(s) to compute k*theta, k = 1,...,n-1


! -- Finally, the evaluated functions can be given by the inverse of Fourier transformation:

!   F(x) = csf(0) + Sigma  i = 1_M { CSF(i)* cos(i * x) + SIF(i)* sin(i * x) } 


!       Input Variables 
!  f        dimension n, the n sample points equally spaced in [0, 2pi]          
!  n        number of sample points 
!  m        number of harmonics to compute
!  
!       Output Variables 
!  csf      cosinus coefficients,  dimension(0:m), first element is the constant term 
!  sif      sinum coefficients, dimension(0:m), first element is 0, no meansing 
!           index start from 0, to keep coherent with gr_foun in fourier.f   
! 
!  Routine Used:
!     trigrec_c, trigrec_s  

!  Finally Revised by Yu -- 20160816
!***********************************************************************
subroutine gr_four(f, n, m, csf, sif)

implicit none 
integer, parameter :: dp = kind(1.d0)
real(kind=dp), parameter ::   pi2 = 8.d0*datan(1.d0)

! Input  and Output Declaration   
integer, intent(in)     ::  n, m
real(kind=dp), intent(in)      ::  f(n)  
real(kind=dp), dimension(0:m), intent(out)     ::  csf, sif  ! csf(m+1), sif(m+1)  
 
! Local Variable
integer :: i, k
real(kind=dp)  :: dtheta, cn(n), sn(n), theta, fpc, fps, &
                  cnm1(n-1), snm1(n-1) ! for k * theta, k=1,...,n-1
  
  dtheta = pi2 / n 
  
  ! compute the basic angle for different Fourier modes j, that is j*theta
!  call trigrec_c(theta, m, cm_bas) 
!  call trigrec_s(theta, m, sm_bas)
!  
  cn(1) = 1.d0;  sn(1) = 0.d0    ! j=0, cos(0) -- sin(0)
  
  ! -- compute csf(0) =  1 / N * Sigma for k=1_n of { f(theta_k) }
  csf(0) = sum(f) / n ! the costant term
  
  print*, 'the costant term', csf(1)
    
  ! --- for m Fourier modes --- 
  do i = 1, m, 1
    ! comptute the trigometric functions for  k * theta
    theta  = i * dtheta 
    
    ! k = 1,...,n-1
    call trigrec_c(theta, n-1, cnm1) 
    call trigrec_s(theta, n-1, snm1)
    
    ! k = 0,...,n, the first element are the same for all the values of theta 
    cn(2:n) = cnm1
    sn(2:n) = snm1
    
    ! Initialize f-by-cos(k*theta) to 0 
    fpc = 0.d0;     fps = 0.d0 
    
    do k = 1, n, 1
      fpc = fpc + f(k) * cn(k)
      fps = fps + f(k) * sn(k)
    end do
    
    csf(i) = fpc * 2.d0 / n 
    sif(i) = fps * 2.d0 / n
     
  end do  

  return  
end subroutine gr_four 



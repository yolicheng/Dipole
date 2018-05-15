!***********************************************************************
!     ****   num_nf   ****
! Modify the number of Fourier modes according to the maximum norm of the last half of the CK, SK
! which are ckm and skm, if |ckm| < tol &&  |skm| < tol, we decrease the size by half

! repeat the same process until we have the appropriate value for minimum number of nf that keeps
! the accuracy of the Fourier representation

! if |ckm| < tol || |skm| < tol, we double the size: nf = 2*nf? 
! TODO: this constraint is too strict??? 
! I think we only need to check the last 2.... 



!       Input Variables 
!  nf     the initial value for nf 
!         which will be updated later 

! ck, sk  dimension n*nf, the Foureier coefficients         

!       Output Variables 
!                   

!  Routine Used:
!     

!  Finally Revised by Yu -- 20160
!***********************************************************************
subroutine num_nf( nf, out  optional)

use dp_mod
implicit none
 
! Input  and Output Declaration   
integer, intent(in)     ::
real(kind=dp), intent(in)      ::    
real(kind=dp), intent(out)     ::    
 
! Local Variable
real(kind=dp)  ::
  
    ! -- check the Maximum norm of csf and sif, and choose an appropriate value for nf 
! FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
inf = nf_init / 2
 ckm = dlange('m', n, inf, csfn(:, inf+1 : nf_init), n, 0.d0)
 skm = dlange('m', n, inf, sifn(:, inf+1 : nf_init), n, 0.d0)
 
print*, 'maximum norm of CK and SK for the last ', nf_init/2, 'Fourier harmonics'
print*,  ckm, skm 

print*; read*

if(nf_init < tol / 10.d0) then 
  nf = nf_init / 2
else 
  nf = 2 * nf_init
endif 
  source

  return  
end subroutine num_nf


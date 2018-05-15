!***********************************************************************
!     ****   tau_fc   ****
! Compute the Fourier coefficients of tau(xi) from the small-divisor equation
! \tau(xi) - \tau(xi+rho)  =  T(\xi) - \hat T   (1)
!   both the left hand side is are zero-average, so is the right-hand side
!   so we have 
!      \hat T =   1/2pi * Sigma for xi = 0, ..., 2pi of { T(xi) * dxi }


!  T(xi) can be represented by a truncated Fourier series, with coefficients as CT_0, and CT_K, ST_K
!  where rho is the fixed rotation number, k = 0_Nf, and 

!  T(xi) =   CT_0 + Sigma for k = 1, ..., nf of { CT(k)) * cos( xi * k ) +  ST(k) * sin( xi * k ) }   
! 
! Similarily, we have 
! \tau(xi) =   C0_tau + Sigma for k = 1, ..., nf of { Ctau(k) * cos( xi * k ) +  Stau(k) * sin( xi * k ) }   
! 
! By imposing (1), we could obtain  CK_tau, Sk_tau as a function of CK_T and Sk_T, for k= 1_Nf
!  C_tau = CT(k) * cos(k*rho) - ST(K) * sin(k*rho) 
!  S_tau = CT(k) * sin(k*rho) + ST(K) * cos(k*rho) 

! TODO: Que: C_tau0 = 0 Alex suggests to be 0 
 
!       Input Variables 
! CT, SK    the coefficients for T(xi), which is the return time to the Poincare section
!           dimension 0:m, just to keep coherent with all the other routines           
! rho       the rotation number, unit 2*pi 
! nf        number of Fourier modes 


!       Output Variables 
! Ctau , Stau      dimension (0:nf)the coefficients for \tau(xi), k=1_Nf         

!  Routine Used:
!     trigrec_cot from libgr.a

!  Finally Revised by Yu -- 20160825
!***********************************************************************
subroutine tau_fc( ct, st, rho, nf, ctau, stau)

use dp_mod
implicit none
 
! Input  and Output Declaration   
integer, intent(in)     ::  nf
real(kind=dp), intent(in)    ::  rho, ct(nf), st(nf)
real(kind=dp), intent(out)   ::  ctau(nf), stau(nf)   
 
! Local Variable
integer :: i
real(kind=dp)  :: cotn(nf), rhohf
                  
  rhohf = rho / 2.d0
  
  ! cos(k* rho/2) and sin(k* rho/2) computed recurrently  by trigrec_cot, k = 1_nf
  call trigrec_cot( rhohf, nf, cotn)  
  
  do i = 1, nf, 1
    ctau(i) =   .5d0*( ct(i) + st(i) * cotn(i) )
    stau(i) =   .5d0*( st(i) - ct(i) * cotn(i) )
  end do
  
  return  
end subroutine tau_fc


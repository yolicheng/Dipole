!***********************************************************************
!     ****   rot_fcs   ****

!  This routine is to avoid the indetermination of the phase shift in  the 
!  computation of invariant curve. 

!  The phase shift indetermination is: given \varphi(xi) as a solution that satisfies
!  the invariance equation \varphi(xi+rho) - P( \varphi(xi) ) = 0, then for any xi0 lies in R
!  \varphi(xi + xi0) also does. 

!  We choose to compute xi0 such that we can make c1(or S1)_J = 0 

!  The Foureier coefficients CK and SK  are transformed to 

!  [ck_hat, sk_hat] ^ T = |  cos(k*xi0)   sin(k*xi0) |  * [ck, sk] ^ T
!                         | -sin(k*xi0)   cos(k*xi0) |

!  For k = 1, we have 
!  c1_hat =  cos(xi0) * c1 + sin(xi0) * s1 
!  s1_hat = -sin(xi0) * c1 + cos(xi0) * s1 

! which is actually a rotation around the normal through an angle xi0 

! if is_ck = 0, we need to s1_hat = 0,  so xi0 = -theta 
! where theta is the phase anlge of [c1, s1]  in polar coordinate 

! else if  is_ck = 1, we need to C1_hat = 0,  so xi0 = pi/2 - theta
 
 
! Now the problem becomes the computation of xi0 by c1 and s1, be careful with the 
! quadrant 

!  where r1 = sqrt(c1^2 + s1^2), (c1, s1) is the first element in ck and sk respectively

! 
!       Input  Variables 
!  nf        dimension of ck and sk 
 
!       Input and Output Variables 
!  ck, sk    Foureier coefficients CK and SK and the updated ones such that sk(1) = 0          

!  Routine Used:
!     

!  Finally Revised by Yu -- 20160903
!***********************************************************************
subroutine rot_fcs( ck, sk, n, nf)

use dp_mod
use pi_mod
implicit none
 
! Input  and Output Declaration   
integer, intent(in)      :: n,  nf  
real(kind=dp), intent(inout)   ::  ck(n, nf), sk(n, nf)  
 
! Local Variable
integer :: i, k 
real(kind=dp)  :: ckcopy(n, nf), skcopy(n, nf), c1, s1, r1, theta, xi0, cn(nf), sn(nf)

  ! make a copy 
  ckcopy = ck
  skcopy = sk
    
  ! extract cos(alpha) and sin(alpha) from ck(1) and sk(1)
  ! assume we fix s1_1, from the first components 
  c1 = ck(1, 1)
  s1 = sk(1, 1)
  r1 = dsqrt(c1**2 + s1**2)
  
  !  compute the phase angle of (c1, s1), and be careful with the quadrant
  ! the range of value for dacos is [0, pi]
  theta = dacos(c1 / r1  ) ! 
  
  ! if (c1, s1) lies in the third or the fourth quadrant, add pi to the result
  if( s1 < 0.d0) theta = theta + pi  
  
  ! compute xi0, according to the  value of is_ck
  if(is_ck == 1) then 
    xi0 = pi/2 - theta
  else 
    xi0 = -theta 
  endif  
    
  call trigrec_c(xi0, nf, cn)  
  call trigrec_s(xi0, nf, sn)  
  
  ! test with element-by-element computation  --- keep this for the moment!
  ! ck_hat =  cos( k*xi0 ) * ck + sin(k*xi0 ) * sk 
  ! sk_hat = -sin( k*xi0 ) * ck + cos(k*xi0 ) * sk 
  
  do i = 1, n 
    do k= 1, nf, 1
      ck(i, k) =  cn(k) * ckcopy(i, k) + sn(k) * skcopy(i, k)
      sk(i, k) = -sn(k) * ckcopy(i, k) + cn(k) * skcopy(i, k)
    enddo    
  end do
  
  print*, 'check carefully if ', ind_cs1, '-th coefficient = 0? ';
  print*, 'is_ck = ', is_ck
  print*, 'if is_ck = 0, S(', ind_cs1, ',1) = 0' 
  print*, 'if is_ck = 1, C(', ind_cs1, ',1) = 0'
  print*; read* 
  
  print*, ' ck \\ sk'
  do i = 1, n, 1
    write(*, '(10e20.10)') ck(i, :); print* 
    write(*, '(10e20.10)') sk(i, :)
    print*; print*; read*
  end do
 
  return  
end subroutine rot_fcs






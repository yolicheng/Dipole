subroutine detmat( a, n, det)
! This routine is to compute the determinant of a square matrix A,  based on the Lapack routines DGETRF for real matrices 
! 	Input
!  a		The input matrix to compute the determinant 
!  n 		number of column in A

! 	Output
!  det 		the determinant of A


implicit none 
integer, parameter :: dp=kind(1.d0)

! Input and Output declaration
integer, intent(in)  ::  n
!integer, intent(out) ::  info
real(kind=dp), intent(in) :: a(n,n)
real(kind=dp), intent(out) :: det 
external :: dgetrf

!  Local 
integer  :: i, info, ipiv(n)
real(kind=dp) :: sgn, lu(n,n)


ipiv = 0
lu = a 

! do the LU decomposition 
call dgetrf(n, n, lu, n, ipiv, info )
 
 
! TODO, --- no good... but any way it works   
if(info /= 0) then 
  print*, 'Fail to do LU decomposition!'
  det = 0.d0
  return
endif 

! check the output of dgetrf
!print*, 'check if a is diagonal'
!do i = 1, n
!  print*, lu(i,:) 
!enddo 
!read*


! Determinant of U is the product of the diagonal elements
det = 1.d0
do i = 1, n
  det = det * lu(i, i)
end do


!-- Adjust sign based on number of row exchanges, odd: -1; even: 1 
sgn = 1.d0

do i = 1, n

  if(ipiv(i) /= i) then
    sgn = -sgn
  end if

end do

det = sgn * det   

end subroutine detmat


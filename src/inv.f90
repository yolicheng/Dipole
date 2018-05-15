! Returns the inverse of a square matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
! From the internet, http://fortranwiki.org/fortran/show/inv

subroutine inv(A, nrow, ncol, Ainv)
  implicit none 
  integer, parameter :: dp = kind(1.d0)
 
  integer, intent(in) ::  nrow, ncol 
  real(kind=dp), intent(in) :: A(nrow, ncol)
  real(kind=dp), intent(out):: Ainv(nrow, ncol)


  ! Local Varaibles
  real(kind=dp)  :: work(nrow)  ! work array for LAPACK
  integer :: ipiv(nrow)   ! pivot indices
  integer ::  info !, i

  ! External procedures defined in LAPACK
  external :: DGETRF, DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
!  n = size(A,1)
!  
!  print*, 'check A' ! ckd 
!  print*, A
!  do i = 1, nrow
!    print*, A(i,:)
!  enddo
!  read*
  
  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(nrow, nrow, Ainv, nrow, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(nrow, Ainv, nrow, ipiv, work, nrow, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
  
return
end subroutine inv


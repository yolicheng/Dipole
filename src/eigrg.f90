subroutine eigrg(a,n,isv, wr,wi,vr)

!  The routine computes for a general, real,  n-by-n nonsymmetric square matrix A, the
!  eigenvalues and, optionally, the right eigenvectors. 
!  

!The right
!  eigenvector v(j) of A satisfies
!
!  A*v(j)= lambda(j)*v(j)
!
!  where lambda(j) is its eigenvalue. The left eigenvector u(j) of A satisfies
!
!  u(j)H*A = lambda(j)*u(j)H
!
!  where u(j)H denotes the conjugate transpose of u(j). The computed
!  eigenvectors are normalized to have Euclidean norm equal to 1 and
!  largest component real.

! TODO, non-square matrix, compute the singular values of 
 
!Non square matrices can indeed have eigenvalues. Well... not exactly. The definition of an eigenvalue is for square matrices. For non square matrices, we can define singular values:

!Definition: The singular values of a m×nm×n matrix AA are the positive square roots of the nonzero eigenvalues of the corresponding matrix ATAATA. The corresponding eigenvectors are called the singular vectors.

!Of course, these have certain properties, that may or may not be useful for what you are trying to study.


!    Input
!  a	      n-by-n real general matrix
!  n	      dimension of matrix a
!  isv	1:  compute also the right eigenvector, 0:  only eigenvalue

!	Output
!  wr,wi the real and imaginary part the eigenvalues
!  vr    right eigenvector, for real 1-by-1 map, 
!        for complex(only keep the eigenvalue with positive imaginary part), 1-real part, 2-imaginary part
!  
!     .. Parameters ..
implicit none
integer, parameter :: dp = kind(1.d0), lwmax = 1000

integer, intent(in) :: n, isv
real(kind=dp), intent(in) :: a(n,n)
real(kind=dp), intent(out) :: wr(n), wi(n), vr(n,n)

character:: jobvl, jobvr

integer  ::  lda, ldvl, ldvr, info, lwork, i, j, shp(2)
!
!     .. Local Arrays ..
real(kind=dp) :: acopy(n,n), vl(n, n), work(lwmax) 

!     .. External Subroutines ..
EXTERNAL         DGEEV
EXTERNAL         PRT_EIGVAL, PRT_EIGVEC
!
!     .. Intrinsic Functions ..
INTRINSIC     INT, MIN

! just to keep a unmodified 
acopy = a 

lda = n
ldvl = n
ldvr = n

if(isv == 1) then 
  jobvr = 'V'
else
  jobvr = 'N'
endif

jobvl = 'N'

!     Query the optimal workspace.
!
LWORK = -1
CALL DGEEV( jobvl, jobvr, N, A, LDA, WR, WI, VL, LDVL, &
               VR, LDVR, WORK, LWORK, INFO )
LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!
!     Solve eigenproblem.
CALL DGEEV( jobvl, jobvr, N, A, LDA, WR, WI, VL, LDVL, &
                 VR, LDVR, WORK, LWORK, INFO )
!
! Check for convergence.
IF( INFO.GT.0 ) THEN
  WRITE(*,*)'The algorithm failed to compute eigenvalues.'
  STOP
END IF

! Print eigenvalues.
! cancel the numerical residual 

do i = 1, n
  if (dabs( wr(i) ) < 1.d-15) wr(i) = 0.d0
  if (dabs( wi(i) ) < 1.d-15) wi(i) = 0.d0
enddo    

shp = shape(vr) 

do i = 1, shp(1)
  do j = 1, shp(2)
    if (dabs( vr(i,j) ) < 1.d-15) vr(i,j) = 0.d0
  enddo
enddo


print*, 'Eigenvalues'
CALL PRT_EIGVAL(  N, 6, WR, WI )
!
! Print left eigenvectors.
! CALL PRINT_EIGENVECTORS( 'Left eigenvectors', N, WI, VL, LDVL )
!


if (isv == 1)  then 
  print*, 'Eigenvectors'
  CALL PRT_EIGVEC(N, 6, WI, VR )
endif
      
return
END

 

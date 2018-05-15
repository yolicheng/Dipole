! Compute the null space of a general rectangal matrix A, which maybe rank-deficient
! leading dimension LDA:  is used to define the distance in memory between elements of two consecutive columns which have the same row index

! --1. use DGEQP3 for QR fracterization

! --2. use  DORMQR to get the orthogonal basis and the complement for A, 
!      as the first rank(A) columns and rank(A)+1 to N columns, respectively. 


! ** NOTE **  Be careful with the dimension
!
! -- for DORMQR:  SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
!                   SIDE = 'L'     SIDE = 'R'
!       TRANS = 'N':      Q * C          C * Q
!       TRANS = 'T':      Q**T * C       C * Q**T
!    C is DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N matrix C. 
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. LDC >= max(1,M)
!          while Q is of dimension m*n, if we want to recover the whole Q, set C identity 

 !   C has to to identity matrix of size N-by-max(M,N),  if we want to recover Q
 !   the second dimension of C has to be at least N, ldc = max(m,n) 
 

! --  for DGEQP3:  subroutine dgeqp3	(M,N,A,LDA,JPVT,TAU,WORK,LWORK,INFO )
!  DGEQP3 computes a QR factorization with column pivoting of a
!  matrix A:  A*P = Q*R  using Level 3 BLAS.
!
!   A is DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit, the upper triangle of the array contains the
!          min(M,N)-by-N upper trapezoidal matrix R; the elements below
!          the diagonal, together with the array TAU, represent the
!          orthogonal matrix Q as a product of min(M,N) elementary
!          reflectors.


! -- ref3: http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=5&t=4506
!          http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=5&t=4542
! The best way is to do a FULL QR factorization of A^T and get the Null space from the last column of Q 
! To compute its null-space, the matrix to factorize has to be transposed 


!-- ref2: http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=3866
! Look for the orthogonal basis for {u: Au=0}, where u is in R^n, A is a matrix m*n (m<<n), with rank m.
! One way is to do QR factorization for A^T (transpose of A), the last n-m columns of Q will be the 
! orthogonal basis for the null space. That is, if written Q=[Q1,Q2]. Q2 is what I need.

! QR factorization can be done with LAPCK subroutine DGEQRF, then the orthogonal Q matrix can be found with DORMQR.

! However, after done this, I tested that when apply A to one column of Q2, it's not zeros (as large as ~1.e-2). 
! So, Q2 is not exactly the basis for the null space. 

 
! -- ref1: http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=227
! One way to solve the problem would be to compute a QR factorization with
! Householder transformation of the initial matrix. 
! To do this you call DGEQRF on A where A is M-by-N where M>=N. 
!     CALL DGEQRF( M, N, A, M, TAU, WORK, LWORK, INFO ) obsoleted and replaced by  DGEQP3
! Then you call DORMQR 
!      CALL DORMQR( 'L', 'N', M, M, N, A, M, TAU, C, M, WORK, LWORK, INFO )
!  Where C has previously been initialized to the identity matrix of size M-by-M. 
!  The first rank columns of C represent an orthogonal basis for B. 
!  The next rank+1 to N columns of C represent a basis of the orthogonal complement of B.

!	Input 
!  nr		number of rows in A 
!  nc		number of columns in A
!  A		left hand side coefficient matrix   
!  b		right-hand-side value function
!   
!	Output 
!   a_null 	the null space of A, dimension: nc - by - nc-rank 

!  Finaly revised by Yu 20170611
!----------------------------------------------------------------------                                                                             

subroutine null_cal(nr, nc, a, rank, a_base)
  use dp_mod 
  implicit none 
  
  integer, intent(in)        :: nr, nc 
  integer, intent(out)       :: rank
  real(kind=dp), intent(in)  :: a(nr, nc) 
  real(kind=dp), intent(out) :: a_base(nc, nc)
  
    
! Local Varaibles
  integer, parameter :: lwork = 2**18 !  3n+1
  integer :: jpvt(nr), info, i, nrhs, ldaT, mn, mx, ldc, debug  
!  the dimension of the array JPVT must be at least max(1,n). but here we work with nc-nr transpose 

!  TAU is DOUBLE PRECISION array, dimension (min(M,N))
  real(kind=dp)  :: acopy(nr,nc), aT(nc,nr), work(lwork),  tol,   & 
                    tau(min0(nr,nc)), c_work(nc,max0(nr,nc)),  amula_null(nr,nc)
                    
                    
! real(kind=dp)  :: aT2(nc, max0(nr,nc))m amula_null2(nr, max0(nr,nc)) ! debug using DORGQR, discard finally
                    
  
  EXTERNAL  ::  DGEQP3, DORMQR, DORGQR ! Lapack LLS solver
  
  debug = 0
!  debug = 2 ! only print the Singular values of A,check the dimension of the kernel 
  
  
  if(debug == 1) then 
    print*,'debug null_cal'
    print*, 'shape(A)', shape(a); print*;read*
  
    print*,'last 4 columns origin A:'
    do i = 1, nr,  1
      print*,  a(i,nr-3:nc)
    end do
    print*;read*
  endif    

  mn = min0(nr,nc)
  mx = max0(nr,nc)
  
  ! copy a and b, because after the call of dgelsy, a and b will be modified
  acopy = a 
  aT = transpose(acopy)
  
  nrhs = 1
  ldaT = nc  ! number of rows in aT 
  ldc  = nc  ! number of rows in c_work, leading dimension 
  
! Initialize JPVT to be zero so that all columns are free
  jpvt = 0
 
 
! Choose tol to reflect the relative accuracy of the input data, for double precision
! use a small one, a good option is one smaller than the cpu precision. 
  TOL = 1.d-9; !5.d-5
  if(debug == 2) TOL = 1.d-11
  
  
  TOL = 1.d-11 ! for tori = 22
  
! SUBROUTINE DGEQP3 (M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO)
  CALL DGEQP3(nc, nr, aT, ldaT, JPVT, TAU, WORK, LWORK, INFO)

!  On exit: if m >= n, the elements below the diagonal are overwritten by details of the orthogonal
!  matrix Q and the upper triangle is overwritten by the corresponding elements of the n by n upper
!  triangular matrix R.
!  If m < n, the strictly lower triangular part is overwritten by details of the orthogonal matrix Q and
!  the remaining elements are overwritten by the corresponding elements of the m by n upper
!  trapezoidal matrix R.

 
!  Determine and print the rank, K, of R relative to TOL
  do i = 1, mn, 1  
    IF (dabs(aT(i,i)) .LE. TOL*dabs(aT(1,1)) ) exit 
    if(debug == 2)      print*, aT(i,i)
  enddo 
  rank = i - 1
  
  if(debug == 1) then 
    print*, 'Tolerance used to estimate the rank of A'
    print*,  TOL
    print*, 'Estimated rank of A^T,  shape of A^T'
    print*,  rank,  shape(aT)
    print*; read*;  read* 
  endif 
    
  
  ! check which are the basis and kernel of A 
!  aT2 = 0.d0
!  aT2(:, 1:nr) = aT 
  
!  --- ! discard DORGQR, since the orthogonal matrix Q cannot be recovered completely if n>m (at most n-columns )
  !     and use DORMQR instead, which is far better 
  !  use this one to check  ! N>=K>=0 
!  SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )

!      A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!              On entry, the i-th column must contain the vector
!              which defines the elementary reflector H(i), for i =
!              1,2,...,k, as returned by DGEQRF in the first k
!              columns of its array argument A.  On exit, the M-
!              by-N matrix Q.

!  print*, 'dimension of aT2 before dorgqr:', shape(aT2); print*
!  CALL DORGQR (nc, nr, mx,    aT2, ldaT, TAU, WORK, LWORK, INFO)
!  but note that the second dimension of the array A must be at least M, which may be larger than was
!  required by DGEQP3
!  print*, 'dimension of aT2 after dorgqr:', shape(aT2); print*; read*
! ---------------------------------------------------------------------



 ! C has to to identity matrix of size N-by-N, if we want to recover Q
  c_work = 0.d0
  do i = 1, nc, 1
    c_work(i,i) = 1.d0
  end do
  
!  print*, 'dimension of aT and c_work before dormqr:', shape(aT), shape(c_work)
!  print*; read*
  

!  subroutine DORMQR(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO)
!  dimension: m X n for C, lda X k for A, k <=M for 'L' and K<=N for 'R'
  CALL DORMQR('L','N', nc, mx, nc, aT, ldaT, TAU, c_work, ldc, WORK, LWORK, INFO)
    
! the first 1-rank columns of C_work will be the orthogonal base of A, while 
! the rank+1-N columns correspond to the null space 
  a_base = c_work(1:nc, 1:nc) 
  
  if(debug == 1) then 
    print*,' The null space of A'
    do i = 1, nc,  1
      print*,  a_base(i, rank+1:nc)  
    end do
    print*; read*
  endif 
    
!  
!  print*,'origin A:'
!  do i = 1, nr,  1
!    print*,  a(i,:)
!  end do
!  
!  print*,'Base of  A:'
!  do i = 1, nc,  1
!    print*,  a_base(i,:)
!  end do
!  print*; read*
  
  ! check if A is orthogonal to the null space, the last nc-rank columns
  amula_null = 0.d0  
  amula_null  = matmul(a, a_base) 
 
  ! check if the null space is orthogonal with A 
  print*, 'check if a and a_null are orthogonal: maxval(a*a_null)', maxval(dabs(amula_null(:, rank+1:nc))); print*;
!  do i = 1, nc,  1
!    print*,  amula_null(i,rank+1:nc)  
!  end do
  
   if(debug == 2)   read*
  
  return
end subroutine null_cal
   

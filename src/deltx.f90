!***********************************************************************
!     ****   deltx   ****
!  Linear equation solver from Lapack: deglsd
!        A* X = B, min || A * X - B||

!  since so far we only need to solve the linear system with nrhs == 1, so we keep b as a column vector 

!-- TODO: for the future possible case in which nrhs >  1, process the dimension of corresponding variables 
!         in the callee by the 'reshape' function. 
  
!       Input Variables 
! nr      number of row in A               
! nc      number of columns in A, also the dimension of dx
! a       dimension nr-by-nc, the coefficient matrix, for Newton Method, it is usually the Jacobi matix 
! b       dimesnion nr*nrhs, the right hand array (always the error to be cancelleds) 
!         in fact, lapack can deal with more columns, but in our application, 1 column is enough

!       Output Variables 
! x       dimension nc-by-nrhs, the solution of A*X = B        

!  Routine Used:
!    deglsd  

!  Finally Revised by Yu -- 20160803 
!***********************************************************************
subroutine deltx( nr, nc, nrhs, a, b, x, info)

implicit none 
integer, parameter :: dp = kind(1.d0)

! Input  and Output Declaration   
integer, intent(in)         ::  nr, nc, nrhs 
integer, intent(out)        ::  info
real(kind=dp), intent(in)   ::  a(nr, nc), b(nr*nrhs) 
real(kind=dp), intent(out)  ::  x(nc*nrhs)
 
 
!  --- explation of the variables usde in dgelsd -------------
! dgelsd (M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, IWORK, INFO)
!  M       number of rows in A
!  N       number of columns in A

! NRHS     number of right hand sizes, i.e., number of column of X and B 
! A        LDA-by-N matrix, the coefficient matrix, use as input, so it's more flexible to call deltx 
! LDA      leading dimension of A,  >=max(1,M)

! B        LDB-by-NRHS matrix, the right hand side matrix is DOUBLE PRECISION array,
!          dimension (min(M,N)). On entry, the M-by-NRHS right hand side matrix B.
!          On exit, B is overwritten by the N-by-NRHS solution matrix X.
! LDB      leading dimension of B,  >=max( 1, max(M,N) )

! S        The singular values of A in decreasing order. dimension (min(M,N))
!          The condition number of A in the 2-norm = S(1)/S(min(m,n)). 

!RCOND     RCOND is DOUBLE PRECISION
!          RCOND is used to determine the effective rank of A.
!          Singular values S(i) <= RCOND*S(1) are treated as zero.
!          If RCOND < 0, machine precision is used instead.

! RANK     The effective rank of A, i.e., the number of singular values
!          which are greater than RCOND*S(1).

! WORK     WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
! LWORK    For good performance, LWORK should generally be larger.
!           12*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2
!            NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )
!            SNLSIZ ~= 25

! IWORK     dimension (MAX(1,LIWORK))
!           LIWORK >= max(1, 3 * MINMN * NLVL + 11 * MINMN),MINMN = MIN( M,N ).
!           On exit, if INFO = 0, IWORK(1) returns the minimum LIWORK.           

! INFO     = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  the algorithm for computing the SVD failed to converge;
!                if INFO = i, i off-diagonal elements of an intermediate
!                bidiagonal form did not converge to zero

! -------- for deglsd  ---------- 
!external  ::  dgelsd  ! lapack lls solver
      
! length of work,  3*n + 64*(m+1), use a big value to allocate enough memory
!integer            ::  lwork = 3*nc + 64*(nr+1) !2**20 

integer, parameter     ::  lwork =  2**22 
integer            ::  rank, iwork(lwork)    
real(kind=dp)      ::  s(nr), rcond, work(lwork)  
! -----------------------------------------------------------
      
! Local Variable
integer :: ldb, debug, allocatestatus!, lwork
real(kind=dp)  :: acopy(nr,nc), awork(nr,nc), bcopy(nr*nrhs)
real(kind=dp), allocatable  ::  bwork(:, :) 
                   
! ---  debug  --- 
integer :: i 
real(kind=dp)  :: df(nr), dbm0, dfm,  dxm 
!real(kind=dp)  :: ainv(nc, nr), x2(nr), dx2m, df2(nr), df2m


real(kind=dp)  :: dnrm2
  
! checked everything, to keep a and b untouched, use the copy to do any computation 

  debug = 0 ! if we want to debug, put 2
  
  ! B is of dimension ldb-by-nrhs, where ldb = max(1, max(m,n))
  ldb   = max0(nr,nc) 
  allocate ( bwork(ldb, nrhs), stat = allocatestatus)
  if (allocatestatus /= 0)  stop "*** not enough memory ***"
 
 
 ! first of all, to avoid the value of a and b being modified by any chance, 
 ! we put a copy, and keep the input as input, not being touched at all. 
  acopy = a
  bcopy = b
  
!  if(debug == 1) then 
!    print*, 'nr, nc, nrhs =', nr, nc,  nrhs; read*
!!    print*, ' A = '
!!    do i = 1, nr, 1
!!      write(*, '(10e20.10)') A(i, :)
!!    end do

!  print*, ' b = '
!  write(*, '(10e20.10)') b 
!  print*; read*
!  endif 
   
 
  dbm0 = dnrm2(nr*nrhs, bcopy, 1) / dsqrt( dble(nr*nrhs) )
  print*, 'error to cancel:   ', dbm0;  !read*
  
  
! LLS  slover initialization 
!  nrhs = 1  ! only one column in B  -- discard, and declare as input 
 
! Choose RCOND to reflect the relative accuracy of the input data, for double precision, 
! use a small one, a good option is one smaller than the cpu precision. 
  RCOND = 1.d-14 ! double precision 
  
  awork = acopy 
  bwork = 0.d0 
  bwork(1:nr, :) = reshape(bcopy, (/nr, nrhs/) )    ! working array, not to change b... 


! dgelsd (M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, IWORK, INFO)
  CALL DGELSD(nr, nc, nrhs, awork, nr, bwork, LDB, S, RCOND, RANK, WORK, LWORK, IWORK, INFO)
  
  ! use reshape to get a 1d array X, but make sure the second input of 'reshape' function is an array 
  ! NOTE: be careful !!! the size of bwork on exit of dgelsd is  nc-by-nrhs, the solution x 
  x = reshape( bwork(1:nc, :), (/nc*nrhs/) )
  
  if(debug == 1) then 
    WRITE (*,*) 'DGELSD-- Least norm solution'
    WRITE (*, '(1X,6e24.14)')  (bwork(I,1),I=1, 12)
    print*, 'rank(A) = ', rank, '         shape(A):', shape(a) ;
  
    df  = matmul(a, x ) - b
    dfm = dnrm2(nr, df, 1) / dsqrt( dble(nr) )
    print*, 'dfm = ', dfm;  
    read* 
    
    print*, 'nc, nrhs', nc, nrhs ;     read*
    print*, 'singular values of A :'
    write(*, '(8e24.14)') s 
    
    print*, 'rcond = ', rcond ;     read* 
    print*, 'rank(A) = ', rank ;
    read* 
  
    df  = matmul(a, x ) - b
    print*, 'df =  A * X - b : ' !, df 
    print* 
  endif 
   
   
  ! ---- debug, use direct matrix multiplication to check, failed for non-square matrix   
!  dfm = dnrm2(nr, df, 1) / dsqrt( dble(nr) )
!  print*, 'dfm = ', dfm; read* 
!  
!  dxm = dnrm2(nc*nrhs, x, 1) / dsqrt( dble(nc) )
!  
!  print*, 'x by lapack, dxm = ', dxm
! 
!  write(*, '(10e20.10)')  x ; print*;    read* 
!  
!!  try the general approach to slove A*X=B by least norm,  X = inv(A)*B
!  acopy = a 
!  bcopy = b
!  call inv(acopy, nr, nc, ainv)
!  x2 = matmul(ainv, bcopy)
!  
!  dx2m = dnrm2(nc, x2, 1) / dsqrt( dble(nc) ) 
!  print*; print*, 'x = A^(-1)*b, dxm = ', dx2m
!  write(*, '(7e24.14)')  x2 
!  
!  df2  = matmul(A,  x2) - bcopy 
!  df2m = dnrm2(nr, df2, 1) / dsqrt( dble(nr) )
!  print*, 'df2m = ', df2m; read* 
!  
!   by check,  df2m =    1.7023887882498004E-017, dx2m is of order 1.d-3
!   while  dfm =    2.6770179825759835E-003, dxm is of order 1.d-4
!   but we tried also replacing x by x2, the result is even worse...
  
!  read*
  
  return  
end subroutine deltx


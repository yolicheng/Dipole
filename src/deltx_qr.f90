!******************************************************************************
! Use Lapack driver dgelsy to compute the  solve the least norm problem  Ax = b, 
! where A may be rank-deficient, based on DGEQP3 for the QR fracterization with column pivoting  

! while dgelsd used SVD fracterization 

! SUBROUTINE DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,  WORK, LWORK, INFO ) ! used in deltx.f90  

!	Input 
!  nr		number of rows in A 
!  nc		number of columns in A
!  A		left hand side coefficient matrix   
!  b		right-hand-side value function
!   
!	Output 
!  dx 	the least-norm solution for null space 

!  Finaly revised by Yu 20170611
!*****************************************************************************

subroutine deltx_qr(nr, nc, a, b,  dx, info)
  use dp_mod 
  implicit none 
  
  integer, intent(in)      :: nr, nc 
  integer, intent(out)     ::  info 
  real(kind=dp), intent(in)  :: a(nr, nc), b(nr)
  real(kind=dp), intent(out) :: dx(nc)
  
    
! Local Varaibles
  integer, parameter :: lwork = 2**18 !  3*nc + 64*(nr+1)
  integer :: jpvt(nc), rank, i, nrhs, ldb , allocatestatus, debug  
  
  ! -- dgelsy
  real(kind=dp)  :: acopy(nr,nc), bcopy(nr), a_work(nr, nc), rcond, work(lwork) 
  real(kind=dp), allocatable  ::  b_work(:, :) 
  
  ! -- dgelsd
!  integer :: iwork(lwork)
!  real(kind=dp) :: s(nr), a2_work(nr, nc) 
!   real(kind=dp), allocatable  ::   b2_work(:, :)
  
  EXTERNAL  ::  DGELSY, DGElSD  ! Lapack LLS solver
  
  acopy = a 
  bcopy = b 
  
  nrhs = 1
 
  debug = 0
  
  if(debug == 1) then 
    print*, 'the last three columns in A'
    do i = 1, nr 
      print*, a(i, nc-3:nc)
    enddo
    print*; read*
  endif 
    
   ! B is of dimension ldb-by-nrhs, where ldb = max(1, max(m,n))
  ldb   = max0(nr,nc) 
  allocate ( b_work(ldb, nrhs), stat = allocatestatus)
  if (allocatestatus /= 0)  stop "*** not enough memory ***"
  
  ! copy a and b, because after the call of dgelsy, a and b will be modified
  a_work   = acopy 
  
  b_work = 0.d0 
  b_work(1:nr, :) = reshape(bcopy, (/nr, nrhs/) )    ! working array, not to change b... 
  
  
  
  
! Initialize JPVT to be zero so that all columns are free
  jpvt = 0
 
! Choose RCOND to reflect the relative accuracy of the input data, for double precision
! use a small one, a good option is one smaller than the cpu precision. 
  RCOND = 1.d-10
     
!  SUBROUTINE DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,  WORK, LWORK, INFO )
!   LDB >= max(1,M,N)     
! Solve the least squares problem min( norm2(b - Ax) ) for the x of minimum norm.
  CALL DGELSY(nr, nc, nrhs, a_work, nr, b_work, ldb, JPVT,RCOND, RANK, WORK, LWORK, INFO)
 
  if(debug == 1) then 
    WRITE (*,*) 'Estimated rank of A'
    WRITE (*,'(1X,I6)') RANK
  endif   
    
   ! checked! both solver return the same results, but we keep deglsy for this routine, and dgelsd for deltx.f90  
!  allocate ( b2_work(ldb, nrhs), stat = allocatestatus)
!  if (allocatestatus /= 0)  stop "*** not enough memory ***"
!  a2_work  = acopy 
!  b2_work  = b_work
!  CALL DGELSD(nr, nc,nrhs, a2_work, nr, b2_work, ldb, S, RCOND, RANK, WORK, LWORK, IWORK, INFO)

! Print solution
  if( info == 0 ) then 
  
!    WRITE (*,*) 'DGELSD-- Least norm solution'
!    WRITE (*, '(1X,6e24.14)')  (b2_work(I,1),I=1,12)
  
    ! assign the solution to dx  
    dx = reshape( b_work(1:nc, :), (/nc*nrhs/) )
    
    if(debug == 1) then 
      WRITE (*,*) 'DGELSY-- Least norm solution'
      WRITE (*, '(1X,6e24.14)')  (b_work(I,1),I=1, 12)
    
     ! Print the effective rank of A
      WRITE (*,*)
      WRITE (*,*) 'Tolerance used to estimate the rank of A'
    
      WRITE (*,'(3X,1P,E11.2)')  RCOND
    
      WRITE (*,*) 'Estimated rank of A'
      WRITE (*,'(1X,I6)') RANK
    endif   

  else 
  
    WRITE (*,*) 'The algorithm failed to converge! deltx_qr'; read*
    
  endif
 
  if(debug == 1) then 
    print*, 'the correction: dx:'
    do i = 1, nc 
      print*, dx(i)
    enddo
    print*; read*
    
    print*, 'the last two components in dx:'
    print*, dx(nc-1:nc); print*; read*
  endif 
   
  return
end subroutine deltx_qr
   

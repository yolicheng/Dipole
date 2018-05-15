!******************************************************************************
! Use Lapack driver dgelsy (QR fracterization)to solve the least norm problem  Ax = b, 
! where A may be rank-deficient. 
! set b=0 to compute the null space 

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

subroutine null_cal(nr, nc, a, b,  dx)
  use dp_mod 
  implicit none 
  
  integer, intent(in)        :: nr, nc 
  real(kind=dp), intent(in)  :: a(nr, nc), b(nr)
  real(kind=dp), intent(out) :: dx(nc)
  
    
! Local Varaibles
  integer, parameter :: lwork = 2**14 !  3*nc + 64*(nr+1)
  integer :: jpvt(nc), rank, info, i, nrhs 
  
  real(kind=dp)  ::  a_work(nr, nc), b_work(nr+nc), rcond, work(lwork) !dgelsy
  
  ! dgelsd
  integer :: iwork(2**14) 
  real(kind=dp) :: s(nr), a2_work(nr, nc), b2_work(nr+nc)
  
  EXTERNAL  ::  DGELSY, DGElSD  ! Lapack LLS solver
  
  
  nrhs = 1
 
 ! copy g and f, because after the call of dgelsy, a and b will be modified
  a_work       = a 
  b_work(1:nr) = b
 
  
! Initialize JPVT to be zero so that all columns are free
  jpvt = 0
 
! Choose RCOND to reflect the relative accuracy of the input data, for double precision, 
! use a small one... 
  RCOND = 1.d-6
     
! Solve the least squares problem min( norm2(b - Ax) ) for the x of minimum norm.
  CALL DGELSY(nr, nc, nrhs, a_work, nr, b_work,  nc, JPVT,RCOND, RANK, WORK, LWORK, INFO)
 
  a2_work       = a 
  b2_work(1:nr) = b
  CALL DGELSD(nr, nc,nrhs, a2_work, nr, b2_work, nc, S, RCOND, RANK, WORK, LWORK, IWORK, INFO)

! Print solution
  if( info == 0 ) then 
    WRITE (*,*) 'DGELSY-- Least norm solution'
    WRITE (*, '(1X,6e24.14)')  (b_work(I),I=1, 12)
  
    WRITE (*,*) 'DGELSD-- Least norm solution'
    WRITE (*, '(1X,6e24.14)')  (b2_work(I),I=1,12)
  
  
    ! assign the solution to dx  
    dx = b_work(1:nc)
    read*
    
! Print the effective rank of A
    WRITE (*,*)
    WRITE (*,*) 'Tolerance used to estimate the rank of A'
    
    WRITE (*,'(3X,1P,E11.2)')  RCOND
    
    WRITE (*,*) 'Estimated rank of A'
    WRITE (*,'(1X,I6)') RANK

  else 
  
    WRITE (*,*) 'The algorithm failed to converge'
    
  endif
  
!    
  end subroutine null_cal
   

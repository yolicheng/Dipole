!******************************************************************************
! Use Lapack driver dgelsy (QR fracterization)to solve the least norm problem  Ax = b, 
! where A may be rank-deficient. 
! set b=0 to compute the null space 


! while dgelsd used SVD fracterization 

! which is better than computing the inverse of the matrix followed by the matrix multiplication 
!   here it is G * X = -F   (A=G : 6 X 7, B = -F)

!	Input 
!  nr		number of rows in G 
!  nc		number of columns in G
!  g		Jacobi Matrix of F w.r.t (X0 + T) 
!  f		Target function
!   
!	Output 
!  dx 		the least-norm correction by sloving G * DX = -F 

!  Finaly revised by Yu 20160407
!*****************************************************************************

  subroutine deltx_lns(f, g, nr,  nc,  dx)
  
  implicit none 
  integer, parameter :: dp = kind(1.d0)
  
  integer, intent(in) :: nr, nc 
  real(kind=dp), intent(in) :: g(nr, nc), f(nr)
  real(kind=dp), intent(out) :: dx(nc)
  
    
! Local Varaibles
! for linear equation solver
  integer, parameter :: lwork = 2**14 !  3*nc + 64*(nr+1)
  integer :: jpvt(nc), rank, info, i, nrhs 
  
  real(kind=dp)  ::  a(nr, nc), b(nr+nc), rcond, work(lwork) !dgelsy
  
  ! dgelsd
  integer :: iwork(2**14) 
  real(kind=dp) :: s(nr), a2(nr, nc), b2(nr+nc)
  
  EXTERNAL  ::  DGELSY, DGElSD  ! Lapack LLS solver
  
  
  nrhs = 1
 
 ! copy g and f, because after the call of dgelsy, a and b will be modified
  a = g 
  b(1:nr) =  -f
 
 ! check this least-norm solver with the symmetric approach, 2 equations and 3 unknowns
  
! Initialize JPVT to be zero so that all columns are free
  jpvt = 0
 
! Choose RCOND to reflect the relative accuracy of the input data, for double precision, 
! use a small one... 
  RCOND = 1.d-6
     
! Solve the least squares problem min( norm2(b - Ax) ) for the x of minimum norm.
!  SUBROUTINE DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,  WORK, LWORK, INFO )  
  CALL DGELSY(nr, nc, nrhs, A, nr, B, nc, JPVT,RCOND, RANK, WORK, LWORK, INFO)
 
  a2 = g 
  b2(1:nr) =  -f
  CALL DGELSD(nr, nc,nrhs, A2, nr, B2, nc, S, RCOND, RANK, WORK, LWORK, IWORK, INFO)

! Print solution
  if( info == 0 ) then 
    WRITE (*,*) 'DGELSY-- Least norm solution'
    WRITE (*, '(1X,7e24.14)')  (B(I),I=1, nc)
  
    WRITE (*,*) 'DGELSD-- Least norm solution'
    WRITE (*, '(1X,7e24.14)')  (B2(I),I=1, nc)
  
  
    ! assign the solution to dx  
    dx = b(1:nc)
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
  end subroutine deltx_lns
   

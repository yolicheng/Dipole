! ***************************************************************************
  subroutine deltx(f,g, dx)
! 20160407 
!  discard this approach by computing the inverse of the matrix followed by matrix multiplication
!  replace with the least-norm solution by Lapack  Linear least square (LLS) solve deglsy  

!  Nov.26,2014 -- general solution for the Modified Newton method with 3 free variables(control), 2 target equations
!  we have Xk = X(k-1) + deltX(k-1)
!  where we require G( deltX(k-1) ) = - F( X(k-1) )
!  It is better to just do the computation here, and check if we reach to the tolerance in the upper subroutine

!  minimize the norm
!     (deltX)T * Q * deltX
!      Q is a given diagonal positive definite weight matrix,  here,just set it to identity

!  so, the solution of this problem of minima is given by
!     (deltXk-1)= -Q-1 (G)T * (G Q-1 (G)T)-1 * F( X(k-1) ) 
!     (deltXk-1)= -(G)T* (G*(G)T)-1 * F( X(k-1) ) 

!    Input Variables
!  f(*)		target variables to be set 0, (here vx,vz )
!  g(*,*) 	Jacobian matrix of f(*) 

!    Output Variables
!  dx   	the difference for the initial state to meet the final target

!    Module-based Varaible
!  ntar, nctr 

!  subroutine used: none
! --------------------------------------------------------------------------
  
  implicit none 
  
  real(kind=dp), intent(in)  :: f(ntar), g(ntar, nctr) 
  real(kind=dp), intent(out) :: dx(nctr)

! Local Variables  
  integer :: i, j, k, &
  	     n, info ! for the norm computed by lapack
  real(kind=dp) :: g2(ntar, ntar), g2inv(ntar,ntar), gnew(nctr, ntar),  & 
                   gnm(nctr),  gnormi, g2m, g2inv2(ntar, ntar), gnew2(nctr, ntar), dx2(nctr)! compute the matrix norm to see if the equation is in well condition

! checked minv, the inv using lapack is not accurate enough, use direct computation subroutine minv instead

! g2 = g*g^T   
  g2 = matmul( g, transpose(g) )
 
  call minv(g2, g2inv, ntar)
!  read*
  
  gnew  = matmul( transpose(g), g2inv )


! Dx = - Gnew * F :   nctr-by-ntar  X  ntar-by-1 ===> nctr-by-1
  dx = - matmul(gnew, f)
  
  if(debug == 1) then ! check the matrix norm of g^T(g*g^T)^(-1)
!   DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )

    gnormi = dlange('I', nctr,  ntar, gnew,  nctr, gnm)
    print*, 'g^T(g*g^T)^(-1), : norm_inf,  norm_row(3)'
    print*,  gnormi, gnm
    write(fgnew,*) gnormi, gnm
    print*
    
!    print*, 'g^T(g*g^T)^(-1)'
!    do i = 1, nctr
!      print*, gnew(i,:)
!    enddo

    read*
    print*, 'dx', dx  ! ckd, exactly the same
  endif 
     
  return      
  end subroutine deltx  

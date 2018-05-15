!***********************************************************************
!     ****   diffpoinc   ****
! The general routine to numercially compute the differential of Poincare map   
! See Carles' note:  On the Analytical and Numerical Approximation of Invariant Manifolds, P14

! We consider a simple situation: the standard first(imax-th) return map to a surface of section  Sigma

! Starting at x0 \in Sigma (the Poincare section), one defines P:  Sigma --> Sigma by 
! P( x0 ) = xf = \varphi( t(x0), x0 ) such that \psi(t) = g( \varphi( t(x0), x0 ) ) = 0
! the arrival time depends, of course, on the initial point. The refinement of the
! intersection is done by routine 'dt_poinc'
 
! let xf \in Sigma be the arrival point. Then, by differentiation w.r.t. x0 one has 

!  d P / d x0 = f(xf) d t / d x0 + d \varphi( t(x0), x0 )  / d x0 

!  where d t / d x0  is computed by the differentiation of g( \varphi( t(x0), x0 ) ) = 0
!  and we have 
!     D g(xf) * ( f(xf) dt / d x0 + d \varphi / d x0 ) = 0

!     where f is the vector field at current point (\varphi(t, x0),  
!     D g(xf) is the differential of g w.r.t. the state of the arrival point 

! d P / d x0 = - f(xf) / ( \nabla g(xf), f(xf) ) D g(xf) * d \varphi / d x0 + d \varphi / d x0

! if we take g(xf) = xf(ind) = 0, then the donominator ( \nabla g(xf), f(xf) ) = f(xf)_ind
! We assume the matrix of the first variational equations is given by d \varphi / d x0  = a_ij

! then the (n-1)-by-(n-1) matrix d P / d x0 | (x0)  restricted to the n-1 components expect ind-th one is given by 

!  d P / d x0 |k,j =  a_k,j - f(xf)_k / f(xf)_ind * a_ind,j ,  1<=j,k<=n && j,k /= ind 

! In order to make it more general, we only compute the necessary components of d P / d x0 


!       Input Variables 
! phi       variational matrix at time t, (State transition matrix)   
! f         vector filed at the arrival point  xf, dimension 6, velocity+acceleration
! ndim      the dimension of the phrase, in order to make general routines available for 4 (prtbp) and 6(RTBP)
! ind       the index of the component of xf set to be zero to define 
!           the Poincare section:  xf(ind) =0 , ie. y=0, we have ind=2
! nr, nc    dimension of dpdx, nr-by-nc 
! para      index of the free parameter, dimension nc
! fun       index of target fcuntions,  dimension nr
!            
! ***Comment***: for asymmetric p.o., para = fun = [1,2,3,4,5] if ind=6
!                for symmetric p.o., it depends, the values may vary a lot
 
!       Output Variables 
!  dpdx     differential of Poincare map associated with target components specified by fun,
!           w.r.t. the parameters specified by para.  dimension nr-by-nc
 
 
!  Routine Used:
!    None  

!  TODO: unchecked.... 
!  Finally Revised by Yu -- 20160808 
!***********************************************************************
subroutine diffpoinc( phi, f,  ndim, ind, nr, nc, para, fun, dpdx)
implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration  
integer, intent(in)            ::  ndim, ind, nr, nc, para(nc), fun(nr) 
real(kind=dp), intent(in)      ::  phi(ndim, ndim), f(ndim)  
real(kind=dp), intent(out)     ::  dpdx(nr, nc)
 
! Local Variable 
integer :: i, j
  
  dpdx = 0.d0
  
!   TODO  problem sloved!!!!
!  print*, 'STM' 
!  do i = 1, ndim, 1
!    write(*, '(10f18.8)')  phi(i, :)
!  enddo 
!  read* 
! 
!  print*, 'vector field' 
!  write(*, '(10f14.8)')  f 
!  read* 
!  
!  print*, 'dpdx'
  ! compute the component one-by-one, the component of dpdx with index i,j corresponds to 
  ! the one specified by fun(i) and para(j) in phi 
  
  do i = 1, nr 
  
      if(fun(i) == ind) cycle 
     
    do j = 1, nc 
    
       if(para(j)  == ind) cycle 
      
      dpdx(i,j) = phi( fun(i), para(j) ) -  f( fun(i) ) / f(ind) * phi( ind, para(j) ) 
    enddo 
    
!    write(*, '(10f14.8)')  dpdx(i, :) 
  enddo 
   
!  read*
  return  
end subroutine diffpoinc  

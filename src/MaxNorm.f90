!***********************************************************************
!     ****   MaxNorm   ****
! Compute the maximum norm of a vector of dimenison n
! For matrix case, we could simply call Lapack function 'dlange'

!       Input Variables 
!  A      1d array 
!  n      dimension of A        

!       Function return Variable 
!  maxA            

!  Routine Used:
!   dabs, maxval  

!  Finally Revised by Yu -- 20160907
!***********************************************************************
function MaxNorm(A, n) result(maxA) 

use dp_mod
implicit none
 
! Input  and Output Declaration   
integer, intent(in)            ::  n
real(kind=dp), intent(in)      ::  A(n)   

real(kind=dp)   ::    maxA
 
! Local Variable
real(kind=dp)  :: absA(n)

  absA = dabs(A)
  maxA = maxval(absA)
    
  return  
end function MaxNorm


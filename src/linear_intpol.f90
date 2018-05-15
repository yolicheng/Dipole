!***********************************************************************
!     ****   linear_intpol   ****
! The linear interpolation, given 2 points with (x1, y1) and (x2, y2)
! compute the value of y for a value of x satisfies is x1 < x < x2 

! Note: we assme x2 is not equal to x1 in this case.... 

!       Input Variables 
!  x0,y0    dimension 2, the input data             
!  x        the value of x to be interpolated 

!       Output Variables 
!  y        the linear interpolation for y at x           

!  Routine Used: None
!     
!  Finally Revised by Yu -- 20160808
!***********************************************************************
subroutine linear_intpol(x0, y0, x, y)
implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration   
real(kind=dp), intent(in)       :: x0(2), y0(2), x   
real(kind=dp), intent(out)      :: y 
 
 
! Local Variable
real(kind=dp)  ::  dx 
  ! we need to check if the denominator is close to zero 
  dx =   x0(2) - x0(1)
  
  if(dabs(dx) < 1.d-4)  then 
    print*, 'Denominator too small! dx=', dx 
    read*
  endif  
  
  y = y1 + ( y0(2) - y0(1) ) / dx * (x - x0(1))
  
  return  
end subroutine linear_intpol


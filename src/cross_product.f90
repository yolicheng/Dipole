subroutine cross_product(a, b, cross )
!! compute the cross production of two array, only work with array of size 3 

implicit none
integer, parameter :: dp = kind(1.d0)

real(kind=dp), intent(in)   :: a(3), b(3)  
real(kind=dp), intent(out)  :: cross(3)

 
  cross(1) = a(2)*b(3) - a(3)*b(2)
  cross(2) = a(3)*b(1) - a(1)*b(3)
  cross(3) = a(1)*b(2) - b(1)*a(2)

  return
end subroutine cross_product


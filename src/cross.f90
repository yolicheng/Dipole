!subroutine cross_product(a, b, cross )
!!! compute the cross production of two array, only work with array of size 3 

!implicit none
!integer, parameter :: dp = kind(1.d0)

!real(kind=dp), intent(in) :: a(3), b(3)  
!real(kind=dp), intent(out)  :: cross(3)

! 
!  cross(1) = a(2)*b(3) - a(3)*b(2)
!  cross(2) = a(3)*b(1) - a(1)*b(3)
!  cross(3) = a(1)*b(2) - b(1)*a(2)

!  return
!end subroutine cross_product


!---discard the function approach, tried several ways to declare the size of cross but failed. .... 
!double precision, dimension(3)  function cross(a, b) 
!! compute the cross production of two array, only work with array of size 3 
!! TODO: untested
!implicit none
!integer, parameter :: dp = kind(1.d0)
!  
!!real(kind=dp), dimension(3)  :: cross
!real(kind=dp), intent(in)    :: a(3), b(3)
! 
!  cross(1) = a(2)*b(3) - a(3)*b(2)
!  cross(2) = a(3)*b(1) - a(1)*b(3)
!  cross(3) = a(1)*b(2) - b(1)*a(2)
!  
!end function cross


 
!function s3_product(a, b, c)
!integer, parameter :: dp = kind(1.d0)
!  real(kind=dp) :: s3_product
!  real(kind=dp), dimension(3), intent(in) :: a, b, c
! 
!  s3_product = dot_product(a, cross_product(b, c))
!end function s3_product


! 
!function v3_product(a, b, c)
!integer, parameter :: dp = kind(1.d0)
!  real(kind=dp), dimension(3) :: v3_product
!  real(kind=dp), dimension(3), intent(in) :: a, b, c
! 
!  v3_product = cross_product(a, cross_product(b, c))
!end function v3_product
  
  

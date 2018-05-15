! This subroutine computes the angle between 3 consecutive points 
!  We check the angle formed by the three last tori 
!  (to compute this angle we take into account h, δ, ρ and the Fourier coefficients of order zero). 
!  If it is greater than a given tolerance (we have typically used 15◦), 
!  we divide by 2 the continuation step and restart
!  In some situations, it is necessary to restart the process from the first of the three last tori.
  
! 	Input Parameters
!  curv(3,*)	the coordinates of three points 
!  ncol		number of columns in curv 
  
! 	Output Parameters
!  csangle	the cosine of the angle beween the two vectors that connects the 3 points 

!  Last modified  2017-06-17 21:14:30 
! ---------------------------------------------------------------------- 
  
subroutine curv_angle( curv, ncol,  csangle)  
use dp_mod 
implicit none 

integer, intent(in) :: ncol 
real(kind=dp), intent(in)  :: curv(3, ncol)
real(kind=dp), intent(out) :: csangle

real(kind=dp) :: dnrm2 


!	Local Parameters
integer :: i, debug 
real(kind=dp) :: v1(ncol), v2(ncol )


debug = 0 ! --ckd

! check the three point 
if( debug == 1) then  !ckd 
  print*, 'check curv for angle compuation'
  do i = 1 ,3 
    print*, curv(i, :)
  enddo 
  print* ;read*
endif 


v1 = curv(2, : ) - curv(1, : )
v2 = curv(3, : ) - curv(2, : )

! normalize 
v1 = v1 / dnrm2(ncol, v1, 1)
v2 = v2 / dnrm2(ncol, v2, 1)

 csangle = dot_product( v1, v2 ) ! v1, v2 are both unit vector, so do not need to devide the norm 

if(debug == 1) then 
  print*, 'cos<v1,v2>=', csangle 
  read* 
endif   

print*, 'curv_angle:  ', csangle; print*


return
end subroutine curv_angle

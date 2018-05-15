! This subroutine computes the angle between 3 consecutive curves 

! 	Input Parameters
!  curv(3,*)	the coordinates of three points 
!  ncol		number of column in curv 
  
! 	Output Parameters
!  csangle	the cosine of the angle beween the two vectors that connects the 3 points 
  
subroutine curv_angle( curv, ncol,  csangle)  

implicit none 
integer, parameter :: dp=kind(1.d0)

integer, intent(in) :: ncol 
real(kind=dp), intent(in)  :: curv(3, ncol )
real(kind=dp), intent(out) :: csangle

!	Local Parameters
integer :: i
real(kind=dp) :: v1(ncol), v2(ncol )

! check the three point 
if( debug == 1) then  !ckd 
  print*, 'check curv for angle compuation'
  do i = 1 ,3 
    print*, curv(i,:)
  enddo 
  print* ;read*
endif 


v1 = curv(2, : ) - curv(1, : )
v2 = curv(3, : ) - curv(2, : )

! normalize 
v1 = v1 / dnrm2(nctr, v1, 1)
v2 = v2 / dnrm2(nctr, v2, 1)

 csangle = dot_product( v1, v2 ) ! v1, v2 are both unit vector, so do not need to devide the norm 

if(debug == 1) then 
  print*, 'cos<v1,v2>=', csangle 
  read* 
endif   


return
end subroutine curv_angle

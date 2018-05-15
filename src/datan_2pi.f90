!***********************************************************************
!     ****   datan_2pi   ****
! Compute the angle between a vector (x,y) and the x-axis, the value range is [0,2pi]
! using atan(y/x), and combined with the fact of the signs of coordinates

! ** NOTE ** 
!  the notation of the quadrants depends on the sense of the curve 
!  if it moves clockwise, it should be     IV-III-II-I,      angle = datan(-y/x)
!  if it moves counterclockwise, then we follow the classical notation: I-II-III-IV
!      angle = datan(y/x)

!   since the intrinsic function datan has the value range in [-pi/2, pi/2], 
!   acos: [0, pi];   asin: [-pi/2: pi/2]
!   to obtain values in the intereval [0, 2*pi], we need to do for points in the 
!   -- first quadrant (x>0,y>0), do nothing
!   -- second and third quadrant (x<0),  datan(y/x) + pi 
!   -- fourth quadrant, datan(y/x) + 2*pi 

!       Input Variables 
!  x,y      the coordinates of point that is the end of a vector
!  dir      the sense of the curve, 
!                 -1: clockwise, with fourth quadrant as the actual first one 
!                  1: conterclockwise, the classical quadrants notation       

!       Output Variables 
!  arg      the angle, value range is [0, 2pi], based on which quadrant the point is in          

!  Routine Used:
!    atan  

!  Finally Revised by Yu -- 20160811
!***********************************************************************
real(kind=8) function datan_2pi( x, y, dir) ! result arg

implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration   
integer, intent(in)            ::  dir
real(kind=dp), intent(in)      ::  x, y    
 
! Local Variable
real(kind=dp), parameter ::  pi = 4.d0*datan(1.d0) 
  

  datan_2pi = datan(dir*y/x) 
  
!  print*, 'x, dir*y', x, dir*y, 'arg = ', datan_2pi 
!  read*
  
  ! -- 1st quadrant,  arg = arg
  
  ! for the 2nd and 3rd quadrant, arg = arg + pi 
  if( x <  0.d0 )  then 
     datan_2pi = datan_2pi + pi
          
  elseif ( dble(dir*y) < 0.d0  ) then ! this is x<0, y<0, the fourth quadrant
  
    ! 4-th quadrant, datan(y/x) is in [-pi/2, 0], so we need to add 2*pi
    datan_2pi = datan_2pi + 2.d0*pi 
     
  end if  
  
  if(dsqrt(x*x+y*y) < 1.d-10) then 
    datan_2pi = 0; return 
  endif 
  
  if(dabs(x) < 1.d-10) then  ! x = 0 
     datan_2pi = pi/2.d0 
  endif 
  
  return  
end function datan_2pi

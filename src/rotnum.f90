!***********************************************************************
!     ****   rotnum   ****
! Compute the rotation number of an invariant curve, which is a set of points
! See Carles' note, P63. 

! And return the value of the ordered index and the coorsponding angles, in order to do interpolation later 

! here, we only need two coordinates, and we introduce an angle, with the 
! mean value of both coordinates as the origin

! 1. for each point with coordinate (x,y), we compute the angle
!   alpha = atan( (y - ymean) / (x-xmean) )

! 2. assume we now have alpha0, alpha1, ..., alpha_n, the number of full revolution 
!    n = 0, n1, n2, ..., nj,  n_j+1 = n_j, if alpha_j+1 > a_j, otherwise n_j+1 = n_j + 1

! 3. Improvement: 
!    3.1 we have alpha0, alpha_n1, ..., alpha_nn \in [0, 2pi], corresponding to 
!               n0 (=0),       n1, ...,     nn 
!                 
!    3.2 Let us order the values of alpha_j increasingly 
!        alpha_i0 <= alpha_i1  <= ... <= alpha_i_j <= alpha_i_j+1 <= ... <= alpha_inn          
!           n_i0       n_i1                 n_i_j      n_i_j+1                n_inn 
!           
!    3.3 Assume  that looking at alpha_i_j <= alpha_i_j+1, one has n_i_j < n_i_j+1 in i_j+1 - i_j iterates, 
!        the angle has increased by alpha_i_j+1 - alpha_i_j+1 + sth.

!    3.4 Assume alpha_i_j <= alpha_i_j+1, but n_i_j >= n_i_j+1  (=> i_j >= i+j+1)
!        that means going from iterate i_j+1 to i_j.
 
!        In this case, we start with rho_sup = 1, rho_inf = 0 
!        Look at alpha_i0 <= alpha_i1 
!        if( n_i0 < n_i1),  rho_inf = max{ rho_inf, (n_i_j+1 - n_i_j) / (i_j+1 - i_j)}
!       if( n_i0 > n_i1),   rho_sup = min{ rho_sup, (n_i_j+1 - n_i_j) / (i_j+1 -
! 
!     3.5 At the end, we have a value of rho_inf and rho_sup 
!         we know rho_inf < rho < rho_sup 
!         if rho_sup - rho_inf < 1.d-8, we know rho with an accuracy of 1.d-8 

! **Note** 
!   1.  what happens if rho_sup < rho_inf, initial point NOT IN A CURVE!!!

!   2.  Expected accuracy:
!     First method without improvement: rho = n_m / m, after m iterates, errors = O(1/m)    
!     second method: errors = O ( 1 / m^2)
!   3.  arg is the angle between the position vector (origin at the center of all the points) and the x-axis,
!       arg is not equal to rho*2pi, arg is just an auxillary angle to sort all the points to compute the rotation number, 
!      since rh0=0 for the first point, which is a new estimate to evalution the invariant curve 
               

!       Input Variables 
!  pt       dimension   np-by-2, in principle, are the x-y coordinates of the points            
!  np       number of points

!  TODO : try to assign the value of dir automatically, 
!         compute the angles of the first two points by default quadrants, and compare the value 
!  dir      the sense how the points move along the curve, 1 : counterclockwise, -1: clockwise

!       Output Variables 
!  rho      the rotaion number, unit:  cycle per revolution 
!           **NOTE** if we want the unit to be radian, multiply by pi2             

!  Routine Used:
!     sort_1d (from sort_mod module)

!  Finally Revised by Yu -- 20160
!***********************************************************************
subroutine rotnum( pt, np, dir, rho)

use sort_mod

implicit none
integer, parameter :: dp = kind(1.d0)   
real(kind=dp), parameter :: pi2 = 8.d0*datan(1.d0)

! Input  and Output Declaration   
integer, intent(in)        ::  np, dir

real(kind=dp), intent(inout)    ::  pt(np,2)
real(kind=dp), intent(out)      ::  rho
 
! Local Variable
integer       :: i, nrev(np), nj, njp1, ind_sorted(np)
real(kind=dp) :: xmean, ymean,  x, y,  rho_inf, rho_sup, rho_temp, arg(np)

real(kind=dp)  :: arg2pi  ! declare the type of the function 
  
  
 ! the mean of both coordinates, seen as the origin to compute the argument   
  xmean = sum( pt(:,1) ) / np
  ymean = sum( pt(:,2) ) / np
  
  write(*,*)' average =', xmean,  ymean; read*; !ck
  
  ! -- compute the angle for each point, remember that we want the value range to be [0, 2pi]
  ! so be careful with the sign of x and y components and decide the quadrant respectively.
  
  ! and also the number of full revolutions of each point  based on the values of the angle at current and previous points
  nrev(1) = 0
  
  ! by default, the curve goes counterclockwise, dir = 1, the quadrants are the standard one  --- dicard, and set dir manually. 
  ! TODO: this one is hard to do..... keep the previous one, observe the motion of the curve and set
  !       the value of dir by hand ....
!  dir = 1
!  
!  do i = 1, 2, 1
!    x = pt(i,1) - xmean 
!    y = pt(i,2) - ymean
!    
!    arg(i) = arg2pi(x, y, dir) 
!  endif 
  
  ! check the direction of the curve, since value range is [0, 2pi]
  ! if arg(2)> arg(1), it goes counterclockwise, dir = 1
  
  
  do i = 1, np, 1
    x = pt(i,1) - xmean 
    y = pt(i,2) - ymean
    
    arg(i) = arg2pi(x, y, dir) 
    
!    print*, '(x0,y0) = ', pt(i,:), '(x,y)=', x,y, 'arg=', arg(i)
!    read*
!  
    ! -- compute the number of full revolutions of each point by compare the angle with the one of the previous point  
    !    start from i=2, the first point has nrev(1) = 0
    
    if (i == 1 )  cycle  ! skip the first point( already initialized as 0)
    
    ! if the angle increases, the current point is at the same revolution as the previous point 
    ! if it decreases, we go to the next revolution 
      
    if( arg(i) > arg(i-1) ) then 
      nrev(i) =  nrev(i-1)
    else 
      nrev(i) = nrev(i-1) + 1
    endif 
    
  end do    
 
  ! -- sort the angles in increasing order --  
  !  subroutine sort_1d( a, na, order, ind_sorted)
  call sort_1d( arg, np, 1, ind_sorted)
  
  ! update the revolutions, to keep coherent with the values of the angles
  nrev = nrev(ind_sorted)
  
  ! the order of point start from 0, so decrease the order by 1 for all the points 
  ! in order to compute rotation number later
  ind_sorted = ind_sorted - 1
  
  ! check the angles and the order   !-- ckd
!  print*, 'angle -- order  -- revolution'
!  write(*, '(5e24.14)')  arg         ; print*
!  write(*, '(5I5)')      ind_sorted  ; print*
!  write(*, '(5I5)')      nrev 
!  read*
   
  ! Initialize the value of rho_inf, rho_sup
  rho_inf = 0.d0
  rho_sup = 1.d0
  
  ! Use the method in Carles' class note, compute the estimate of the rotation number 
  ! be careful with the order, start from 0, the first point doesn't count 
  do i = 1, np-1, 1
  
    ! check each step all the relevent values  -- ckd 
!    print*, i, '-th point!'
!    print*, ind_sorted(i:i+1), 'iterate!' 
!    print*, 'angles = ', arg(i:i+1), 'revolutions = ', nrev(i:i+1)
!    read*
    
    nj   = nrev(i)
    njp1 = nrev(i+1) 
    
    ! we only update rho when the revolutions of the two adjacent points are different 
    ! when they are equal, do nothing
    if(nj == njp1) cycle 
    
    ! convert the integer value to double precision 
    rho_temp = dble( ( njp1 - nj ) ) / dble( ( ind_sorted(i+1) - ind_sorted(i) ) )
    

    if(nj < njp1) then 
      rho_inf = dmax1(rho_inf, rho_temp)
      
    else
      rho_sup = dmin1(rho_sup, rho_temp)
    endif  
    
    ! check the rotation number for ij-i_j+1 iterate  ! --  ckd
!    print*, i, '-th point!'
!    print*, 'njp, njp1, i_j, i_j+1: ', nj, njp1, ind_sorted(i:i+1)
!    print*, 'Rotation number: rho_inf, rho_temp, rho_sup', rho_inf, rho_temp, rho_sup
!    read*
    
  end do
  
  ! The final guess of the estimate of the rotation number 
  rho = ( rho_inf + rho_sup ) / 2.d0
  print*, 'rho = ', rho, 'rho_inf=', rho_inf, 'rho_sup =',  rho_sup
  
  ! check if arg ind_sorted rho, if they satisfy mod( rho * ind_sorted*2pi, 2pi) = arg? 
  ! NO!  --- ckd
  ! arg /= rho * (ind_sorted) * 2pi, they are not of the same estimate, 
  ! see more details in the comment 
!  print*, 'check mod( rho * ind_sorted) = arg for the first 6 components?'
!  
!  print*, 'rho = ', rho 
!  print*, 'ind_sorted(1:6) = ', ind_sorted(1:6) 
!  print*, 'rho * ind_sorted(1:6) = ', rho*ind_sorted(1:6) 
!  print*; read*
!  
!  print*, 'rho * ind_sorted * 2pi ',  rho * ind_sorted(1:6) * pi2 
!  print*, 'mod( rho * ind_sorted*2pi, 2pi)', dmod( rho * ind_sorted(1:6)*pi2, pi2)
!  print*, 'arg', arg(1:6)
!  print*; read*
  
  return  
end subroutine rotnum

!***********************************************************************
!     ****   arg2pi   ****
! Compute the angle between a vector (x,y) and the x-axis, the value range is [0,2pi]
! using atan(y/x), and combined with the fact of the signs of coordinates

! ** NOTE ** 
!  the notation of the quadrants depends on the sense of the curve 
!  if it moves clockwise, it should be     IV-III-II-I,      angle = datan(-y/x)
!  if it moves counterclockwise, then we follow the classical notation: I-II-III-IV
!      angle = datan(-y/x)

!   since the intrinsic function datan has the value range in [-pi/2, pi/2], 
!   to obtain values in the intereval [0, 2*pi], we need to do for points in the 
!   -- first quadrant (x>0,y>0), do nothing
!   -- second and third quadrant (x<0),  datan(y/x) + pi 
!   -- fourth quadrant, datan(y/x) + 2*pi 

!       Input Variables 
!  x,y      the coordinates of point that is the end of a vector
!  dir      the sense of the curve, 1: counterclockwise, with fourth quadrant as the actual first one 
!                                  -1: clockwise, the classical quadrants notation       

!       Output Variables 
!  arg      the angle, value range is [0, 2pi], based on which quadrant the point is in          

!  Routine Used:
!    atan  

!  Finally Revised by Yu -- 20160811
!***********************************************************************
real(kind=8) function arg2pi( x, y, dir ) ! result arg

implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration   
integer, intent(in)            ::  dir
real(kind=dp), intent(in)      ::  x, y    
 
! Local Variable
real(kind=dp), parameter ::  pi = 4.d0*datan(1.d0) 
  

  arg2pi = datan(dir*y/x) 
  
!  print*, 'x, dir*y', x, dir*y, 'arg = ', arg2pi 
!  read*
  
  ! -- 1st quadrant,  arg = arg
  
  ! for the 2nd and 3rd quadrant, arg = arg + pi 
  if( x <  0.d0 )  then 
     arg2pi = arg2pi + pi
          
  elseif ( dble(dir*y) < 0.d0  ) then ! this is x<0, y<0, the fourth quadrant
  
    ! 4-th quadrant, datan(y/x) is in [-pi/2, 0], so we need to add 2*pi
    arg2pi = arg2pi + 2.d0*pi 
     
  end if  
  
  return  
end function arg2pi



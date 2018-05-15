subroutine gr_plf_n(t, y, nvar, f)
!  Compute the vector field of n points simultaneously in  planar lorentz force problem 
!  for the purpose to integrate several points simultaneously or obtain the time-t map for a curve 

! ** NOT ** in order to do integration using rk78, we have to pass y as an  array with only 1 row.... 

!  so we have 
!  1         - n(=4),    the first point 
!  n+1       - 2n,       the second point 
! ......
!  n*(i-1)+1 - 6*i       for i-th point
! where n=4 for plf, is the dimension of the vector field 

!  include this subroutine into the respective module, so we can directly use the system-based paramters
!  without adding to the input parameter list,  keep the form deriv(t, y, n, f) to be called in gr_rk78

!     Input parameters:
!  t              independent variable (typically time)
!  y(*,*)         state of point, of dimension nrow X 4
!  nvar             number of variables (np/6 points) 

!  Output parameters:
!     f(*)       vector field
!                in y(*) and f(*) the first 6 components correspond to the position and velocity
!
! Subroutine use:  gr_lf 

! Finaly revised by Yu -- 20160426
! -----------------------------------------------
use dp_mod
implicit none 
 
!     Input and Output
integer,  parameter ::  n = 4
integer, intent(in) ::  nvar  
real(kind=dp), intent(in)  :: t, y(nvar)  
real(kind=dp), intent(out) :: f(nvar)

!     Local Variable
integer       :: i, nrow, debug  
real(kind=dp) :: yi(n), fi(n) 

  debug = 0
!     Resource 
  nrow = nvar / n

  if(debug == 1) then 
    print*, 'nvar=', nvar, 'nrow=', nrow; read*
  endif   

 ! make sure the reshape  of the array is done in the same way, which are are column-wise in fortran 
  do i = 1, nrow 

  ! the index of i-th point saved in y
    yi = y( n*(i-1)+1 : n*i )
  
    call  gr_plf(t, yi, n, fi)
  
    f( n*(i-1)+1 : n*i ) = fi
    
    if(debug == 1 ) then  !ckd
      print*, i,'-th point:'
      print*,  yi
      print*,  fi
    endif
     
  enddo 

  if(debug == 1) then  
    print*; read*
  endif   


  return
end subroutine gr_plf_n


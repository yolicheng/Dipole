subroutine gr_lf_n(t, y, nvar, f)
!  computation of the vector field to be integrated for lorentz force problem 
!  designed we can integrate several points simultaneously 
!  but in order to do integration using rk78, we have to pass y as an  array with only 1 row.... 

!  so we have 
!  1 - 6,  the first point 
!  7 - 12, the second point 
! ......
!   6*(i-1)+1 - 6*i    	for i-th point

!  include this model subroutine into the respective module, so we can directly use the system-based paramters
!  without adding to the input parameter list,  keep the form deriv(t, y, n, f) to be called in gr_rk78

!  	Input parameters:
!  t  		independent variable (typically time)
!  y(*,*)	state of point, of dimension nrow X 6
!  np 		number of variables (np/6 points) 

!  Output parameters:
!     f(*)       vector field
!                in y(*) and f(*) the first 6 components correspond to the position and velocity
!
! Subroutine use:  gr_lf 

! Finaly revised by Yu -- 20160426
! -----------------------------------------------

implicit none 
 
! 	Input and Output
integer, intent(in) ::  nvar  
real(kind=dp), intent(in)  :: t, y(nvar)  
real(kind=dp), intent(out) :: f(nvar)

! 	Local Variable
integer :: i, nrow 
real(kind=dp) :: yi(6), fi(6) 

nrow = nvar / 6

!print*, 'nvar=', nvar, 'nrow=', nrow
!read* 

do i = 1, nrow 

! the index of i-th point saved in y
  yi = y( 6*(i-1)+1 : 6*i )
  
  call  gr_lf(t, yi, 6, fi)
  
  f( 6*(i-1)+1 : 6*i ) = fi
  
!  print*, i,'-th point:'
!  print*,  yi
!  print*,  fi
enddo 

!print*; read*

return
end subroutine gr_lf_n


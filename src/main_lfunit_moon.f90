program main_lfunit_moon

implicit none
integer, parameter :: dp = kind(1.d0)

! Variables
integer    ::  i 
real(kind=dp) ::  beta_array(4), beta, runit, tunit, vunit



!subroutine lfunit(beta, runit, tunit, vunit)

beta_array = [1.d-1, 1.d0, 1.d1, 1.d2] !%, 1.d2, 1.d3] 

do i = 1, 3
  beta = beta_array(i)
  call lfunit(beta, runit, tunit, vunit)
enddo 


end program main_lfunit_moon

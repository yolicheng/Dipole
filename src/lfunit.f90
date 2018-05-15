 
!*******************************************************************************  
! Compute the unit of distance for given beta, the parameters B0 and q2m are specified already
! runit = (abs( B0/n * q2m /beta ) ) **(1./3.) 

! TODO: at this moment, we are not sure if the unit is right. check really carefully later ... 

! this file saves all the rescaling unit, including distance, time, velocity for Lorentz force problem
! Note: The unit of distance for HCW equations is 1 

!  	Input variables  
!   beta   :    n / wc

! 	Output Variables
!   runit : unit of distance  : km
!   tunit : unit of time      : s
!   vunit : unit of velocity  : km/s
  
!  ROUTINE USED:  none
!  Finally revised by Yu  -20160428
!*************************************************************************
subroutine lfunit(beta, runit, tunit, vunit)

implicit none
integer, parameter:: dp = kind(1.d0) 

real(kind=dp), intent(in):: beta  
real(kind=dp), intent(out):: runit, tunit, vunit 
 
! 	Local variables
real(kind=dp) ::  muE, rE, pi, B0, alt, q2m, rc, n  

! a = (B0/n * q/m * 1/beta)^(1/3)  ---  the unit of distance
!  where n is the mean orbital motion of the chief around the Earth, which is in a keplerian,circular orbit 

muE  = 3.986d5   ! u_earth -- 		km^3/s^2
rE   = 6371d0    ! radius of Earth --	km

pi =  3.141592653589793 ! common used parameter

!  magnetic moment, take the value of Geomagnetic moment at the surface of the Earth's equator 
!  Peck's paper is 8e15 Wb-m, is the same value in different unit
!  T(tesla) = 1Wb/m^2
! TODO: we have to figure out what is an appropriate value for this parameter.... 
B0 = 8.d6; ! T km^3

! -----  variables needed to be specified ------------------------
!alt = 20200d0            ! altitude -- km  --- Medium Earth Orbit
alt =  35786d0 - rE       ! GEO 35786 km, geostationary orbit, circular orbit around the Earth with the same period 

! q2m = [1e-6 1e-5 1e-4 1e-3] -- take the first one to get the minimum  
q2m = 1.d-6

rc  = rE + alt           ! the radius of the chief from the Earth's center
n   = dsqrt(muE/rc**3)   ! rad/s

! -- time 
tunit = 1.d0 / n 

!tunit = 6860.30148727185      ! s    -- time  !ckd
!print*, 'tunit=', tunit; ! read*
 

! v = n*rc   ! test with r_geo=42,164(alt=35,786), proves right 

! r_star  = (B0/n * q/m * 1/beta)^(1/3)  --- test with the function assignment, fine!
runit = (dabs( B0/n * q2m / beta ) ) ** (1.d0/3.d0) 

! -- distance  
! runit = 38.0024033332631    ! km -- beta = 1
! runit = 30.1625275142708    ! km -- beta = 2
! runit = 17.6391530962123    ! km -- beta = 10  

vunit = runit / tunit         ! km/s -- velocity - pay attention here, because there 

print*, 'runit(km), tunit(s), vunit(km/s)',  runit, tunit, vunit  !ckd
print* !read*

return
end subroutine lfunit

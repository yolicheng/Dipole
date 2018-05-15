!*******************************************************************************  
! Compute the unit of length for the thesis, where the Leader is in a lunar orbit 
!  Finally revised by Yu  -20180307  for the feasibility check of the new model 


! Debye length :   Stubbs, order of 1 m 
!                 Properties of plasma near the moon in the magnetotail
!                 however, could be 100 m 
 
!the parameters B0 and q2m are specified already
! runit = (abs( B0/n * q2m /beta ) ) **(1./3.) 
! 

! this file saves all the rescaling unit, including distance, time, velocity for Lorentz force problem
! Note: The unit of distance for HCW equations is 1 

!  	Input variables  
!   beta   :    n / wc

! 	Output Variables
!   runit : unit of distance  : km
!   tunit : unit of time      : s
!   vunit : unit of velocity  : km/s
  
!  ROUTINE USED:  none

!*************************************************************************
subroutine lfunit(beta, runit, tunit, vunit)

implicit none
integer, parameter:: dp = kind(1.d0) 

real(kind=dp), intent(in):: beta  
real(kind=dp), intent(out):: runit, tunit, vunit 
 
! 	Local variables
!B0 = mu_0/4/pi*n_c*i_c*pi*rc^2 
integer    ::  n_c
real(kind=dp) ::  muM, rM, pi, B0, alt, q2m, rc, n, &
                  mu_0, i_c, r_c   

!  a = (B0/n * q/m * 1/beta)^(1/3)  ---  the unit of distance
!  where n is the mean orbital motion of the chief around the Earth, which is in a keplerian,circular orbit 

!6.674×10−11 m3⋅kg−1⋅s−2. the gravitional constant 

! the gravitational constant of Earth and Moon: m^3/s^2
!Earth	3.986004418(9)×10^14
!Moon	      4.9048695(9)×10^12

!muE  = 3.986d5   ! u_earth -- 		km^3/s^2
!rE   = 6371d0    ! radius of Earth --	km


!M_moon = 7.342×10^22 kg ---- wikipedia 
!muM  = 4.905d3   
! u_Moon = G*M_moon = 6.674d-11 * 7.342d22  = 4.9001d12 m^3/s^2 
muM = 4.9001d12 !  m^3/s^2
rM  = 1737.1d3  ! radius of Moon 1737.1km  	m


pi =  3.141592653589793 ! common used parameter

!  magnetic moment, take the value of Geomagnetic moment at the surface of the Earth's equator 
!  Peck's paper is 8e15 Wb-m, is the same value in different unit
!  T(tesla) = 1Wb/m^2
! TODO: we have to figure out what is an appropriate value for this parameter.... 


!B0 = 8.d6; ! T km^3 % this is for GEO--- not appropriate here 

! B0 is produced by the artificial magnetic field. 
!B0 = mu_0/4/pi*n_c*i_c*pi*r_c^2 

mu_0 = 4.d0*pi*1.d-7;  ! N/A^2
n_c  = 1000 

! critical current density of more than 1 MA/cm^2 https://en.wikipedia.org/wiki/Superconducting_wire
i_c =  1.d6*1.d4 ! A  
r_c  = 1.d0  ! m

B0 = mu_0/4.d0/pi*n_c*i_c*pi*r_c*r_c 
print*, 'B0=', B0; 
print*; read*


! -----  variables needed to be specified ------------------------
!alt = 20200d0            ! altitude -- km  --- Medium Earth Orbit
!alt =  35786d0 - rE       ! GEO 35786 km, geostationary orbit, circular orbit around the Earth with the same period 

! q2m = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1] -- take the first one to get the minimum  
!q2m = 1.d0
!q2m = 1.d-2 
!q2m = 1.d-4
q2m = 1.d-6
  

!rc  = rE + alt           ! the radius of the chief from the Moon's center

! take a 200km lunar orbit as an example 
alt = 200.d3 ! m
rc  = rM + alt           ! the radius of the chief from the Moon's center

!n   = dsqrt(muE/rc**3)   ! rad/s for earth orbit 
n   = dsqrt(muM/rc**3)    ! rad/s

print*, 'n=',n 
print*; read*


! -- time --- I don't think this one is right
tunit =   1.d0 / n 


!tunit = 6860.30148727185      ! s    -- time  !ckd
!print*, 'tunit=', tunit; ! read*
 


! v = n*rc   ! test with r_geo=42,164(alt=35,786), proves right 

! r_star  = (B0/n * q/m * 1/beta)^(1/3)  --- test with the function assignment, fine!
runit = (dabs( B0/n * q2m / beta ) ) ** (1.d0/3.d0) 

!we need to increase B0 and q2m 

! -- distance  -- for geostationary orbit with B_0 
! runit = 38.0024033332631    ! km -- beta = 1
! runit = 30.1625275142708    ! km -- beta = 2
! runit = 17.6391530962123    ! km -- beta = 10  

vunit = runit / tunit         ! km/s -- velocity - pay attention here, because there 

print*, 'runit(m), tunit(s), vunit(m/s)',  runit, tunit, vunit  !ckd
print* !read*

return
end subroutine lfunit

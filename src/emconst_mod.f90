! Put all the constant parameters of EM system into this module for polar problem
! avoid the case that we have to pass a lot of shared arguments down to all
!  subroutines that use them

! ** Note **
! in module, we cannot assign the variables directly other than the parameter type 
!  we have to use a subroutine call init_emconst for the initialization of the parameters


! for the safe side, declare the variables 'private' attribute that are only used in this module

!   NOTE: 
!    ---to keep coherent with the model we are using. Take a reasonable value for mean distance from the 
!       Moon to the Earth, we take the one for Andreu's thesis:   dsem = 384403.6782590720 km
!    I think if we use the CRTBP model, then the orbit of the Moon around the Earth is a circle
!    if we know the mean distance of the Moon to the Earth dsem, we could compute the period as 2*pi 
!       prd = sqrt( dsem**3 / gmm )

!  So we cannot use the real values for both dsem and prd, they are related.... 
!  while the real values take into account all the perturbations. 
! 
! Question:   which values to choose??????   ephemeris?  or prd = sqrt(dsem**3 / gmm )????
 

!    Varaibles 
!  mu	      mass ratio 
!  rmoon	radius of the Moon 
!  dsem	unit of distance, defined as the distance between the Earth and the Moon
!  tuem	unit of time,  so the period of the Moon orbiting around the Earth is 2*pi ! about  1 month/2*pi 
!           so the angular velocity of the Moon around the Earth is n=1
!  vuem	unit of velocity, km/s

!  tmat     the rotation matrix from the initial epoch 2020 Jan.1st. Recovied from the reveiw.tex 

!  0.99338553935461127 & -6.3603242008194415E-003  & -0.11465040984317662   \\ 
!  2.5088190255703746E-003 & 0.99942862516935727 & -3.3706513008396610E-002 \\
!  0.11479928583508195 & 3.3195925475105303E-002 & 0.99283390076266331

! --- Gravitational costant:   -- by data from wikipedia
!  G        = 6.67408(31)  10^{-11}  m^3 kg^{-1} s^{-2}
!  M_{moon}  = 7.34767309 × 10^22  kg 
!  mu_{moon} = 4903.89807942938    km^3 s^{-2}
!  mu_{moon} = 0.0121862743144954 (non-dimensional)

! --- For DE421:  
!  AU:      0.149597870699626200D+09         -- km
!  EMRAT:   0.813005690699153000D+02         -- M_E / M_M
!  GMB  :   0.899701140826804900D-09         -- G * (M_E + M_M)  for Earth-Moon Barycenter [au**3/day**2]

! From which we derive:  
!  GMM  :   4902.80007622774                 -- GMB/(EMRAT+1) * au**3 / day**2   km^3/s^{–2} 

! merat = 1 / emrat  = 0.0123000369055232   -- M_M / M_E

!  --------- The mean distance from the Moon to the Earth   ------------  

! While from wikipedia: the Moon's orbit
!    Astronomers O'Keefe and Anderson calculated the lunar distance by observing 4 occultations from 9 locations in 1952.[31] They calculated a mean distance of 384407.6±4.7 km, however the value was refined by in 1962 by Irene Fischer, who incorporated updated geodetic data to produce a value of 384403.7±2 km.[7]

!   'The Constitution and Structure of the Lunar Interior'.  P305, semi-major axis:  384399 km
!    Orbital period (days)                 --- 27.321582 (2)
!    Mean planetary radius (km)            --- 1737.103 ± 0.015(3)
!    GM_Moon (km3 s−2)                     --- 4902.801076 (4,5)
 
 
! 'DE430 Lunar Orbit, Physical Librations, and Surface Coordinates',  P10   
!    For an elliptical orbit, the time-averaged mean distance is a(1+e2/2). 
!    For the Moon a = 384,399 km and the mean distance is 385,000 km including solar perturbations

! dsem = 385,000 km   ! discard this one, because it takes into account of sun's pertubation

module emconst_mod
 
implicit none
integer, parameter, private :: dp = kind(1.d0) 

real(kind=dp), parameter, private ::  pi = 4.d0*dtan(1.d0)
real(kind=dp) ::  emrat, rmoon, dsem, tuem,vuem , tmat(3,3), dtmat(3,3), &
                  mu_m, gmm  
                  
! The Earth-Moon System, from Elisa Maria's paper: two manoeuver 
!the unit of distance equals 384,400 km, 
!the unit of velocity equals 1.02316 km/s and 
!the dimensionless mass of the Moon is l ¼ 0:012150582

!  ----------  mass ratio of Moon to Earth  ---------------
! --TODO:  which to use?
! mu = 0.1215058609624d-1       ! mu = 0.012150582  from Gerard's book,   0.012150586D from DE405

! while in DE421, emrat = 0.813005690699153000D+02, and mu = 1/emrat = 0.0123000369055232
 
contains 


  subroutine init_emconst
!  implicit none 
                                
  emrat = 0.0123000369055232d0   ! use this one for the moment ...                           

  !  ----------  mean radius of the Moon  ---------------              
  ! TODO: which value to use for the distance?                               
  rmoon = 1.7371d3                  ! km  1737.4 \pm 1,  from DE421;   DE430: 1738 km

 !  ----------  mean distance from the Moon to the Earth  ---------------              
  dsem = 3.8440d5       ! km    wikipedia:  3.8440d5

  !377496.086808103


  !   ---------  orbital period ---------------------
  ! which to use????  -- orbital period, or synodic period??? here I take the latter
  ! Orbital period  --- in fact we should take this one ...
  ! 27.321661 d
  ! (27 d 7 h 43.19 min 11.5 s[1])
  ! Synodic period
  ! 29.530589 d
  ! (29 d 12 h 44 min 2.9 s)

  ! tuem = 29.530589*24*3600 ! second  - Period: 1 month  --- discard this one 

  ! TODO:  mean period of the Moon around the Earth from wikipedia, orbital period. 
  tuem = 27.321661*24*3600/2.d0/pi   ! 375699.807500923 second  - 1 orbital Period is 2*pi

! but if we use the dsem as the radius of the circular orbit of the Moon around the Earth 
! prd = sqrt( dsem**3 / gmm )  -- unit : s 

! use the value of dsem as 384402.1 km, and the above gmm =  from DE421, we get  
!  prd   =  377493.729937798 sec  
!  vuem  = dsem/prd   =  1.01830062200858

!  DE405:  tuem = 375698.7248459793 sec,  from paper1
!  dsem = 384400 km,    vuem = 1.02316024670456

!  From Horizon:  http://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
!  mean orbital velocity:   1.022 km/s 

!   --- Unit of the velocity --- -
  vuem = dsem/tuem  !  1.02316 km/s  


  print *, 'Unit of distance(km) :     ', dsem 
  print *, 'Unit of time (s/2/pi):     ', tuem
  print *, 'Unit of distance(km/s):    ', vuem
  read*

  ! gravitational parameter of the Moon:    km^3/s^{–2} 
  !  GMM  :   4902.80007622774       -- DE421:  GMB/(EMRAT+1) * au**3 / day**2  
 !  GMM  :   4902.801076            -- ref: Lunar Interior Constitution & Structure: 2006 
 
 
  gmm =   4902.80007622774  ! DE421, keep this one....               
  mu_m =  gmm/vuem/dsem     ! rescale to non-dimensional unit

  !  0.0121862743144954   -- gmm computed by values from wikipedia
  !  0.0121836101005017   -- gmm computed with DE421

  !	alt = .3d3 ! the altitude of target orbit 
 ! radius of the Moon, 1737.1km 

  tmat(1,:) = (/0.99338553935461127d0,  -6.3603242008194415d-3, -0.11465040984317662d0/)
  tmat(2,:) = (/2.5088190255703746d-3, 0.99942862516935727d0,   -3.3706513008396610d-2/)
  tmat(3,:) = (/0.11479928583508195d0,   3.3195925475105303d-2,  0.99283390076266331d0/)
  
  dtmat = 0.d0
  
  end subroutine init_emconst    
 
end module emconst_mod 
	
	

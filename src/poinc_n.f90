!  A general routine to compute npoinc prescribed returns to the Poincare section, used to generate Poincare representation!
!  For the moment, we only focus on the state vector, take ndim = 6. Consider not the variational matrix. 

!  For each poincare map, the maximal time and domain is tmax and xmax, respectively.

!  the poincare section is defined by y(ind) = z0,  with the componnet of velocity y(ind+3)  in the direction specified by dir 


! 20160202 -- Add domain constraint, since for the Poincare map, in principle, we are only interested in specific domain. 
!             if x,y and z component are beyong a certain value xmax, in principle, 1 in adimensional unit. The orbit escapes, stop the integration  


!     Input Varaibles
!  xi(ndim)       dimension, ndim,  the initial state 
!  ndim           the dimension of the state: 6 for spatial problem; 4 for planar ones
!  nvar           the dimension of the vector to integrate, ndim*(ndim+1) if we want to variational matrix
!                 = ndim, if we only need to integration the state vector 
!  npoinc         number of Poincare maps to compute
!  tdir           time sense for integration, 1: forward; -1: backward
!  fpc            file tag to save tpc--ypc of each map  

!        Output Variables
!  tf             time spent by the orbit to go from xi to xf
!  xf(*)         the state at ncrs-th crossing
!  ncrs           real number of crossing through the poincare section (for escape name, the value is smaller than n)


!  ---------- use poinc_mod to pass these parameters ----------
!  imax           the number of time of the intersecion with the Poinare section to be considered 
!  tmax           the maximum integration time for one Poincare map (allowed to be updated for the continuation of p.o.s )
!  xmax           the maximum domain for the orbit, beyond that the orbit is treated as escape 
!  ind,p0         parameters that define the poincare section:  st(ind) = p0
!  dir            the direction of the velocity considered when cross Poincare section 
!                     1: v(ind)>0; -1: v(ind)<0; 0: keep both 
! ***********************************************************************
subroutine poinc_n(xi, ndim, nvar, tdir, fpc, npoinc, tf,xf, ncrs, deriv, gr_cj) 

use dp_mod
use poinc_mod, only :   tmax  ! to pass parameters
implicit none


! Input and Output declarations
integer, intent(in)  :: ndim, nvar, tdir, npoinc, fpc 
integer, intent(out) :: ncrs

real(kind=dp), intent(in)        :: xi(ndim)  
real(kind=dp), intent(out)       :: xf(ndim), tf

external :: deriv, gr_cj  ! the vector field and the differential of variational equations


! local variables
integer ::  ispc   

real(kind=dp) ::  x(nvar), tpc, tpci, xpc(nvar), cj, hminim 
                  
! Initialization for pv (or variational matrix if nvar > ndim)
x(1:ndim) = xi
if(nvar > ndim) then 
  x(ndim+1:nvar)            = 0.d0
  x(ndim+1 : nvar : ndim+1) = 1.d0
endif 
 

print*,'poinc_n, x0', xi;  read* !ckd

! Initialization of the Poincare maps 
tpci = 0.d0  ! time 
xpc  = 0.d0  ! state 
ncrs = 0     ! number of obtained maps

!write(fpc, *) ! add one blank line 
 
! Idea: use poinc to obtain one crossing, and then use this point as the new initial point to get the next poincare map
! save the time and state of each intersection point to file fpc(eqpoinc.dat)

do  
!subroutine poinc(sti, ndim, nvar, tdir, tf, stf, hminim, ispc, deriv, gr_cj) 
  call poinc(x, ndim, nvar, tdir, tpc, xpc, hminim, ispc, deriv, gr_cj) 
  
  ! TODO: hminim is necessary to pass as an output parameter? 
  
  if(ispc == 0 ) exit ! tpc, xpc not updated if there is no intersecion

  ncrs  = ncrs + 1    ! record the number of maps we obtained 
  
  ! TODO: the return time or the elapsed time??  Probably the latter one is better
!  tpci = tpc 
  
  tpci  = tpci + tpc  ! time at the intersection: previous time + time spent for the current intersection  

  x     = xpc         ! as the new starting point for the next Poincare map
  
  ! the energy 
  call gr_cj(xpc, cj)
  write(fpc, *)  tpci, xpc(1:ndim), cj, ncrs  
  
!  write(*, *)    tpci, xpc(1:ndim), cj, ncrs   !ck
!  read*
    
!  write(fpc, '(8e24.14, 1I5)')  tpci, xpc, cj, ncrs 
  
!  write(*, '(8f16.8, 1I5)')  tpci, xpc, cj, ncrs ! ck 
!  read*  
  
  ! update maximum integration time, which is used as the stop condition
  if(5*tpc > tmax) tmax = 5*tpc  
  
  ! use the maximal intersecion as the termination condition
  if(ncrs .ge. npoinc) exit
enddo 

! return the state at the last crossing
tf = tpci
xf = xpc

write(fpc, *) ! add one blank line 

!  print*, 'npoinc =', npoinc, 'ncrs =', ncrs; read*
  
return
end subroutine poinc_n



















  

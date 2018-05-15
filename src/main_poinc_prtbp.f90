program poinc_prtbp

! This is to compute the poincare representation of a bunch of initial points, results saved in torusall.dat and curveall.dat 

!  Study the planar Restricted Three Body Problem, and choose L4, which is of type 
!  center X center, so around L4 we have a plenty 2D tori. 
!    1. start from a point that is close to L4, and integrate the orbit 
!    2. compute the Poincare maps, npoinc = 10000, with Poincare section taken as x = p0 (or y= p0)
!    3. compute the rotation number of numerically computed invariant curve 
!    4. use the linear interpolation to obtain equally spaced points w.r.t. the angle between the relative position vector from L4 and the x-axis!
!    5. do a simple Fourier analysis to compute the Fourier coefficients as the initial guess to refine the curve 
!    6. Impose the invariance equations, and specify the energy level, to compute the invariant curve 
!    7. Globalize the curve to an invariant torus 


!  --- data saved in ./dat/curve_prtbp/  subfolder 
!  1. torus.dat         -- the original orbit of the 2D torus 
!  2. curve.dat         -- the invariant curve obtained by Poincare map 
!  3. ob_interpol. dat  -- linear interpolated points for Fourier analysis
!  4. fun.dat           -- approximated orbits by the Fourier representation 
!  
use emconst_mod 
use rtbpconst_mod

implicit none
integer, parameter  ::  dp = kind(1.d0)
integer, parameter  ::  npoinc =  200 , & !20 !1000  ! number of Poincare maps, to start test, 11 is enough?
                        ndim0   = 4, & ! for PRTBP
                        npvar   = 4, & ! 20 for variational matrix, 4 for state vector 
                        np_interpol = npoinc+1, &  ! in order to use fourier, we need odd number of points 
                        nf   =  100 ! how many Fourier modes ???
                         
real(kind=dp)       ::  pi2 = 8.d0*datan(1.d0),  day = 24.d0 * 60.d0 * 60.d0  


! Local Variables
real(kind=dp) :: p0, tmax, xmax, pvi(npvar), tf, pvf(npvar), tic, toc, cj0,   &  ! torus
                 xl4, yl4, & ! L4 equilibirum point
                 tpc(npoinc), xpc(npoinc, ndim0), cj, t, pv(ndim0), & !curve 
                 pt(npoinc, 2), rho, arg(npoinc), & ! rotnum
                 ptnew(np_interpol, ndim0), pv_in(np_interpol), dt, fun6(ndim0)  ! interpol + fourier
  
real(kind=dp), dimension(0:nf, ndim0) ::  CSF6, SIF6 ! fourier coefs of  all six compoments
real(kind=dp), dimension(0:nf)        ::  CSF, SIF   ! fourier coefs of 1 component
 
                  
integer ::  i, j, k, l,  ncrs, ind_sorted(npoinc), &
            ftorus, ispl, fcurve, fcurve_interpol, ffc, ffun, ffcbas, & 
            ind, dir, imax, pt_col, isvy, np ! poinc 
            
character(len=70) :: fntorus, fncurve, fncurve_interpol, fnfc, fnfun,  fnfcbas   


external :: gr_cjprtbp, gr_prtbp ! Jacobi constant + vector field 

real(kind=dp) ::  funcs

! -- assign the value of mass ratio for EM system and put the values into the common module --- 
call init_emconst
print*, 'emrat = ', emrat 
call init_rtbpconst(emrat, ndim0)

print*, 'check the mass ration, mu =', mu, 'ndim=', ndim 
read*


! take a point that is close the L4 point,  with the big primary at (mu, 0, 0), 
! the triangular libration point L4 form a equilateral triangle with the two primaries. 
! so the location of L4 is ( -1/2+mu, sqrt(3)/2  ) and  L5 located at ( -1/2+mu, -sqrt(3)/2  )
xl4 = -.5d0 + mu 
yl4 = dsqrt(3.d0) / 2.d0

! the initial point --- instead of take random points, we fix the energy and poincare section 
pvi = (/xl4 + 1.d-3, yl4 + 2.d-3, 1.d-4, -1.d-4/) ! add a small deviation on L4 to obtain an initial point
! cj0 = 3.0000046704015602  

!this is a test, to fix the energy level, and x0 = xl4,  and change the value of y, to have a look of a lot tori. 


! files to save the original torus + invariant curve + interpolated one + approximated one 
ftorus = 21;  fntorus    =  './dat/curve_prtbp/torusall.dat'
open(ftorus  ,file = fntorus, access ='append', status = 'replace' )
write(ftorus ,*)  ' # Original 2D torus: t      (x y vx vy)   cj'

! first integrate the orbit for a long time, and  observe the geometry to select an appropriate Poincare section
!subroutine plob(y0,t0,tf, n, tdir, ftag, ispl, deriv, gr_cj,  y) 

!tf = 2.d4 !2.d5 !1.d4 ! test a long time interval  -- already computed, do not overlap it....
!call plob(pvi, 0.d0, tf, 4, 1, ftorus, 1, gr_prtbp, gr_cjprtbp,  pvf) 

!stop 


!-- curve by computating Poincare map ----
fcurve   = 25;  fncurve  =  './dat/curve_prtbp/curveall.dat'
open(fcurve  ,file = fncurve  , access ='append',status = 'replace' )
write(fcurve  ,*)  ' # Invariant curve by Poincare map:  t      (x y vx vy)   cj  npc'


! ---------- Computation -----------------------------
! outer bound for escape, nondimensinal unit 1 is big enough....
xmax =  2.d0 
tmax =  5.d1 ! a rough guess?? or not...

! Poincare section: x=-0.49 by observing the plot of the torus 
ind  = 1           ! x component 
p0   = xl4   !yl4 !-0.49d0     ! rough guess, by observing directly at the plot of the tori  
dir  = 1           ! take the positive velocity 
imax = 1           ! consider every intersection, and choose by the direction of velocity


! -- compute the Poincare map   
!subroutine poinc_n(xi, ndim, npvar, ind, p0,  dir, ispl, fob fpc, npoinc, imax, tmax, xmax,  tf,xf, ncrs, deriv, gr_cj) 
!call poinc_n(pvi,ndim, ndim, ind, p0,  dir, 1, ftorus, fcurve, npoinc, imax, tmax, xmax, tf, pvf, ncrs, gr_prtbp, gr_cjprtbp) 

!this is a test, to fix the energy level, and x0 = xl4,  and change the value of y, to have a look of a lot tori. 
 cj0 = 3.0000046704015602  

! fix also the Poincare section, and take points from this section
pvi(ind) = p0 
 
 
ispl = 1

! first we try 10-by-10 points with different values of y and vy 
np = 10 
do i = 1, np, 1
  ! the initial point --- instead of take random points, we fix the energy and poincare section 
  !    
  if(ind == 1) pvi(2) = (i-5)*5.d-3 +  yl4 
  if(ind == 2) pvi(1) = (i-5)*5.d-3 +  xl4 
  
  ! the value of vx, step size 2.d-3
  do j = 1, np, 1
    pvi(3)  = (j-5)*5.d-3
    
    ! --  the value of vy, positive, index of vy is 4
    call  cj2v_prtbp(pvi, cj0,  4, isvy)
    
    if( isvy == 0 ) cycle
    
    
  !  pvi = (/xl4 + 1.d-3, yl4 + 2.d-3, 1.d-4, -1.d-4/) ! add a small deviation on L4 to obtain an initial point

    call poinc_n( pvi,ndim, ndim, ind, p0,  dir, fcurve, ispl, ftorus,  npoinc, & 
                  imax, tmax, xmax, tf, pvf, ncrs, gr_prtbp, gr_cjprtbp)   
   
  ! add two blank lines to seperate different initial points
    write(fcurve, *);   write(fcurve, *);  
    write(ftorus, *);   write(ftorus, *);
  
  end do  
end do 

 close(fcurve) 
stop  

! read to Poincare maps from the file  
open(fcurve, file = fncurve) 
read(fcurve, *) ! skip the first line 
!stop  
! read the Poincare maps from the data file 
do i = 1, npoinc, 1  
!  write(fpc, '(8e24.14, 1I5)')  tpci, xpc, cj, ncrs !-- from poinc_n.f90
  read(fcurve, *) t, pv, cj, ncrs
!  read(fcurve, '(6e24.14, 1i5)') t, pv, cj, ncrs
!  print*, t, pv, cj, ncrs
!  read*
  
  tpc(i) = t
  xpc(i, :) = pv
enddo   
  
! -- compute the rotation number -- 
! be careful, which two columns to use to compute the rotation number...
! x-vx or y-vy 
if(ind == 1) pt_col = 2
if(ind == 2) pt_col = 1
pt = xpc(1:npoinc,  (/pt_col, pt_col+ndim/2/) ) 

call  rotnum( pt, npoinc, rho, arg, ind_sorted)

print*, 'rho= ', rho 
! update the points by an increasing order of arg 
xpc = xpc(ind_sorted+1, :) 

write(*,'(10f14.10)')  arg
read* 



! -- interpolated curve (equally spaced in rho) --
fcurve_interpol = 26;   fncurve_interpol  =  './dat/curve_prtbp/curve_interpol.dat'
open(fcurve_interpol  ,file = fncurve_interpol)
write(fcurve_interpol, *)  ' # Interpolated curve: argument      (x,y,z,vx,vy,vz)'  ! the explanation line 

! -- Fourier analysis : Fourier coefficients + approximated curve --    
! the coefficients computed by general Fourier analysis:  fourier.f + fun.f 
ffc = 111;  fnfc = './dat/curve_prtbp/fcs.dat'
open(ffc, file = fnfc) 
write(ffc,*) '# For each column: (x,y, vx,vy),  for each line: the Four Coef ck(i), sk(i), i=0,...,', np_interpol 

! the approximated Fourier function, to check if they match the original data 
ffun = 222;  fnfun = './dat/curve_prtbp/fun.dat'
open(ffun,  file = fnfun) 
write(ffun, *) ' # argument    (x,y, vx,vy)' ! TODO -- check the energy? cj'

ffcbas = 223;  fnfcbas = './dat/curve_prtbp/fcbas.dat'
open(ffcbas  , file=fnfcbas,  access ='append',status='replace')
! TODO: need to explain this file... not used yet

print*, 'Files names for data read and write'
print*, fntorus, fncurve,  fncurve_interpol, fnfc,   fnfcbas;  print*   !ck


! ----- linear interpolation to get equally spaced points in angle -----------
! do the linear interpolation by the value of arg, to obtain npoinc points equally spaced in arg within the interval [0, 2pi]

call  interpol( xpc, npoinc, ndim, arg, np_interpol, fcurve_interpol, ptnew)
 
! -- save all the points to file ob_interpol.dat, better to do this in the routine interpol  


! -- general Fourier analysis for nf Fourier modes, with np_interpol points  

! deal with the six components one by one 
do k = 1, ndim, 1
  pv_in =  ptnew(:, k) ! one component, the k-th column  
 
  call gr_foun( pv_in, np_interpol, nf, csf, sif)
  
  csf6(:, k) = csf
  sif6(:, k) = sif
  
end do
  
  
! write the coefficients C's and S's into file fcs.dat (x,y,z,vx,vy,vz)  cf // sf !x  // cfy-sfy -- ......
do k = 0, nf, 1
  
   write(ffc, *)  csf6(k, :) 
   write(ffc, *)  sif6(k, :) 
   write(ffc,*)
!     sif6(1,k), csf6(2,k), sif6(2,k), csf6(3,k), sif6(3,k), &
!                     csf6(4, k), sif6(4,k), csf6(5,k), sif6(5,k), csf6(6,k), sif6(6,k)
enddo


! Period = 2*pi/rho, the relation between the period and the rotation number 
dt = pi2/rho/np_interpol ! take this value to keep coherent with the original points  

! Evaluate the truncated Fourier series with the above coefficients
do l = 1, np_interpol*2 
 
    t = (l-1) * dt
    
    ! the 6 components
    do k = 1, ndim, 1
      fun6(k) = funcs( t, rho, nf, csf6(:, k), sif6(:, k) )
    end do
    
!    call gr_cjlf(fun6, cj)
    write(ffun, *) t, fun6  !, cj 
    
  enddo 
  
  write(ffc, *)
  write(ffun, *) 



stop

end 


 




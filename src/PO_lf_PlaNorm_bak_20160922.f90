program po_plf

! just save as a backup

! This is to compute the periodic orbits in the planar LF problem of the normal case 

! Take a grid of points 20*20 for  discretisized in (x, vx)
! and the other two components are zeros: y=0; z=0 
! check the first and the second intersections to see if there is periodic orbits
 

!  --- data saved in ./dat/curve_plf/  subfolder 

!  1. torus.dat         -- the original orbit of the 2D torus 
!  2. curve.dat         -- the invariant curve obtained by Poincare map 
! 
use dp_mod
use pi_mod 
use poinc_mod
use plf_mod 
use po_plf_mod
implicit none

integer, parameter  ::  ndim   = 4,  & ! for lf
                        npvar  = 6, &  !42 for variational matrix, 6 for state vector 
                        npoinc = 100 !0
! Local Variables
integer       ::  debug, i
real(kind=dp) ::  beta0 

! eigenvalues and eigenvectors around the equilibria                 
real(kind=dp) ::  dlf(ndim, ndim), wr(ndim), wi(ndim), vr(ndim, ndim) 


! Poincare map 
integer       :: ind0, dir0, imax0
real(kind=dp) :: p00,  tmax0, xmax0, cj0 

! Po 
integer ::  ispo 
real(kind=dp) ::  y0(ndim-2),   TP, yf(ndim), dg(1,1), cj  

! curve
integer :: ftorus, fcurve , ftorus_poinc, ind_curv, isv, ispc  
real(kind=dp) :: pvi(ndim), pvf(ndim), cj_i, tf, hminim, tf_total
character(len=70) :: fntorus, fncurve 

! for fft_ob 
integer           ::  istol, nfmax, nsm,  nitmax, isfour 
real(kind=dp)     ::  tolres,  tolfc 

! take eq1  in planar normal case as an example, type center X saddle  equilibria
! with beta = 1.d0, there is resonance when beta is around 1, avoid this for the moment 
debug = 1

call init_plf     ! to assign the value of  eq 
 
beta0 = 2.d0   ! use 2 to check 
call init_beta(beta0)   
     
print*, 'check,  beta',  beta0 

call gr_cjplf(eq, cj)
print*, 'check energy! ', cj, eq

! -- the energy level we r taking, smaller than the one at the equilibria, so we have a 
! closed UFO 

 cj0 = cj + 0.05d0
 
read*
    
! -- the real units-----
print*, 'runit(km), tunit(s), vuint(km/s): ', runit, tunit, vunit

! Jacobi matrix of the lorentz force with respective to the state 
call dflrtz2(eq, dlf)

print*, 'Check DG in the planar case: ' 
do i = 1, ndim 
  write(*,'(4e18.8)') dlf(i,:) 
enddo
print*; read* 

call eigrg(dlf, ndim,1, wr,wi, vr)
print*, 'Check if we have a center X center X saddle equilibria:' 


! Take initial point on the x-axis, with velocity perpendicular to x-axis, and compute the P.O. 
! by Poincare map with the section as y=0, take vy as a function of f(cj, x, y=0, vx)
! so we have 2 dimension of freedom (x, vx )

y0 = 0.d0 
y0(1) =  0.61d0 

!! ---------- Computation -----------------------------
! For poincare map limitation 
xmax0 =  0.7d0  !
tmax0 =  3.d1 ! 

! Poincare section:  y=0  
ind0  = 2           ! x component 
p00   = 0.d0        ! rough guess, by observing directly at the plot of the tori  
dir0  = 1           ! no requirement in the direction, just compute the first turn 
imax0 = 1           ! consider every intersection, and choose by the direction of velocity

!subroutine init_poinc( ind0, p00, dir0, imax0, tmax0, xmax0, ndim0, h00) 
call init_poinc(ind0, p00, dir0, imax0, tmax0, xmax0, ndim, cj0)


! instead of refining for po orbit, we just take a lot of points on x-axis, and do the poincare section 
fcurve   = 25;  fncurve  =  './dat/curve_plf/curve_poinc.dat'
open(fcurve  ,file = fncurve  , access ='append',status = 'replace' )
write(fcurve  ,*)  ' # Invariant curve by Poincare map:  t      (x y  vx vy )   cj'

ftorus = 21;  fntorus    =  './dat/curve_plf/torusall.dat'
open(ftorus  ,file = fntorus, access ='append',  status = 'replace' )
write(ftorus ,*)  ' # Original  tori: t      (x y  vx vy )   cj'

ftorus_poinc = 24
open(ftorus_poinc, file = './dat/curve_plf/torus_poinc.dat', access ='append',  status = 'replace' )
write(ftorus_poinc ,*)  ' # the orbit by poinc: t      (x y  vx vy )   cj'


ind_curv = 0
do 
  
  ! for the moment, only compute 10 points
  if(ind_curv >= 3 ) exit
  
  pvi = 0.d0 
  pvi(1) = y0(1) - 1.d-3 * ind_curv ! x
  
!  pvi(1) = -0.60740924112296302d0
!  pvi = (/ 0.56784837604607341 ,      -1.1119702540638688E-022,   8.7175931844212136E-003 , 0.33555563012250067/)

  call cj2v_plf(pvi, cj0, ind+2, isv)
   
!  -0.60740924112296302        0.0000000000000000       -7.2296819316566494E-002  0.13242131288321832
!  the routine stops with this point, test this one  
  
!  print*, 'State of initial point: ', pvi; read*
  
!   check the torus 
  if(ind_curv > 0) open(ftorus, file=fntorus, access='append', status='old')
!  
  call plob(pvi, 0.d0, 1.d3, ndim, 1, ftorus, 1, gr_plf, gr_cjplf,  pvf) 
  write(ftorus,*); write(ftorus,*)
  close(ftorus)
  print*, 'Finish plob, check the orbit!'; read* 

  ind_curv = ind_curv +1
  
  open(fcurve, file=fncurve, access='append', status='old')
  tf_total = 0.d0
  do i = 1, npoinc
!   subroutine poinc(sti, ndim, nvar, ispl, fob, tf,stf, hminim, ispc, deriv, gr_cj) 

    call poinc(pvi, ndim, ndim, 0, ftorus_poinc, tf, pvf, hminim, ispc, gr_plf, gr_cjplf)
    if(ispc == 0) exit 
    
    tf_total = tf_total + tf
    
    if(tf_total > 2.d3) exit 
    
    ! we only care about the velocity vy is small enough... that will be close to a periodic orbit.
    
    call gr_cjplf(pvf, cj_i)
    write(fcurve, *) tf, pvf, cj_i
    
!   if( mod(i, 100) == 0)  
   write(*, *) tf_total, tf, pvf, cj_i; print*; read*
    pvi = pvf 
  enddo 
  
  print*, 'finish one initial point!'
  read* 
  
  write(fcurve, *); write(fcurve, *)
  write(ftorus_poinc, *); write(ftorus_poinc, *)
  close(fcurve)
enddo 




!subroutine init_po_plf( ind_zero )
! we want to take vx0 = 0, which is the second element in (x, vx), so ind_zero=2
call init_po_plf( (/2/) )


!subroutine refn_po(y0, tf, yf, dg, ispo)  
call refn_po(y0, TP, yf, dg, ispo) 

print*, 'After the refinment of P.O.'
print*, 'y0: ', y0
print*, 'tp, yf', TP, yf
print*; read* 

if(ispo == 1) then 
  open(100, file = './dat/curve_plf/po_refn.dat', status='replace')
!  TP = ntp * TP
  pvi = 0.d0
  pvi(1) = y0(1) 
  call cj2v_plf(pvi, cj0, ind+2, isv)
  
  call plob(pvi, 0.d0, 2*TP, ndim, 1, 100, 1, gr_plf, gr_cjplf,  pvf) 
endif 

! -- do fourier analysis on this orbit ---
write(*,*) 'Maximum amplitude of the residual //  Tolerance in correction for refinement? '
print*, '( 1 : Default: 1.d-4, 1.d-6;  Otherwisze input new values)'
read(*,*) istol

if(istol == 1) then 
  tolres = 1.d-4; tolfc = 1.d-6
else 
 read(*,*)  tolres,  tolfc
endif 


tf = 1.d3 
nitmax = 10  
nsm = 2**18
nfmax = 10 

open(60  ,file = 'dat/curve_plf/ob.dat', access ='append',status = 'replace')
write(60, *) ' #   t      (x y z vx vy vz)    cj '

open(61  ,file = 'dat/curve_plf/fcs.dat', access ='append',status = 'replace')
open(62  ,file = 'dat/curve_plf/fcsbas.dat', access ='append',status = 'replace')

!subroutine fft_ob(x0, n, lt, np, tol, tolres, nitmax, nfmax, xmax, fob, ffc, ffbas, deriv, gr_cj) 
call fft_ob(pvi, ndim, tf, nsm, tolfc, tolres, nitmax, nfmax, xmax, 60, 61, 62, gr_plf,  gr_cjplf)  

read*
 
 

stop 

 
 
end 
 



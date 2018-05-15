program poinc_lf

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


!  --- data saved in ./dat/curve_lf/  subfolder 
!  1. torus.dat         -- the original orbit of the 2D torus 
!  2. curve.dat         -- the invariant curve obtained by Poincare map 
! 
use dp_mod
use pi_mod 
use lf_mod 
!use poinc_mod
implicit none

integer, parameter  ::  npoinc =  5000 , & !20 !1000  ! number of Poincare maps, to start test, 11 is enough?
                        ndim   = 6, & ! for lf
                        npvar  = 6  !42 for variational matrix, 6 for state vector 

! Local Variables
integer       ::  debug, isv , ispc, istol
real(kind=dp) ::  beta0

real(kind=dp) :: x0, p0, tmax, xmax, pvi(npvar), tf, pvf(npvar), tic, toc, cj0,   &  ! torus
                 tpc(npoinc), xpc(npoinc, ndim), cj, t, pv(ndim), hminim, cj_i

! for fft_ob 
integer           ::  nfmax, nsm,  nitmax, isfour, ind_curv !, fob, ffc, ffcbas
real(kind=dp)     ::  tolres,  tol 

! for the curve 
real(kind=dp)     ::   w1, w2, dh, pvall(npoinc, ndim), pv0(ndim)
  
! eigenvalues and eigenvectors around the equilibria                 
real(kind=dp) ::  dlf(ndim, ndim), wr(ndim), wi(ndim), vr(ndim, ndim), &
                  dlf_p(ndim-2, ndim-2), wr_p(ndim-2), wi_p(ndim-2), vr_p(ndim-2, ndim-2)

                  
integer ::  i, j, k, l,  ncrs,  &
            ftorus, ispl, fcurve,  & 
            ind, dir, imax,  np ! poinc 
            
character(len=70) :: fntorus, fncurve   


!external :: gr_cjlf, gr_lf ! Jacobi constant + vector field 

! take eq3 in normal case as an example, which is of type center X center X saddle 
! with beta = 2.d0, there is resonance when beta is around 1, avoid this for the moment 
debug = 1
isfour = 0

! to assign the value of cs, ieq, eq 
call init_lf 

!  read from the screen the value of beta  
if(cs == 2) then 
   ! if(ieq==1)  beta = 10.d0 ! TODO 
   if(ieq == 2)  beta0 = 1.d0
   if(ieq == 3)  beta0 = 10.d0 
  
 elseif(cs == 1) then 
   if(ieq == 1) beta0 = 2.d0 
   if(ieq == 2) beta0 = 1.d0
 endif  

!beta0 = 2.d0  
call init_beta(beta0)   
     
print*, 'check, cs, ieq,  beta', cs, ieq, beta0 

call gr_cjlf(eq, cj)
print*, 'check energy! ', cj, eq
read*
    
! Compute the unit of distance, time and velocity, which are all function of beta 
 call lfunit(beta, runit, tunit, vunit)
    
! -- the real units-----
print*, 'runit(km), tunit(s), vuint(km/s): ', runit, tunit, vunit

! Jacobi matrix of the lorentz force with respective to the state 
call dflrtz(eq, dlf)
dlf_p = dlf( (/1,2, 4,5/), (/1,2, 4,5/) )

print*, 'DG : ' 
do i = 1, ndim
  write(*,'(6f18.8)') dlf(i,:) 
enddo
read* 

! compute the eigenvalues and eigenvectors  of dlf
!subroutine eigrg(a,n,isv, wr,wi,vr)
call eigrg(dlf, ndim,1, wr,wi, vr)
print*, 'Check if we have a center X center X saddle equilibria:' 


print*, 'Check the planar case: '
do i = 1, ndim-2
  write(*,'(4f18.8)') dlf_p(i,:) 
enddo
read*

call eigrg(dlf_p, 4, 1, wr_p , wi_p, vr_p)

read*  

! add a small deviation on equilibria to obtain an initial point for a torus 
!pvi =  eq + (/0.d0, 0.d0, 0.d0, 1.d-5, 0.d0, 2.d-5 /) 
!pvi = (/0.6707d0, -0.02d0, 0.d0, 1.d-4,  0.d0, 2.d-5/)

pv0 = (/0.507d0, -0.02d0, 0.d0, 2.d-1,  0.d0, 1.d-1/)
 
 
!pvi = eq + (/1.d-4, 1.d-4, 1.d-4, -1.d-3,  0.d0, 2.d-3/)
! this is the one where the zero velocity surface has a close region, we start from this case 

 cj = cj + 0.05d0

 
! take 10 initial points on the z=0 plane, with velocity perpendicular to z=0 plane 
! and compute the Poincare map, check if there is fixed point, which means there will be periodic orbit 
! from which we can go to the torus 

!! ---------- Computation -----------------------------
! outer bound for escape, nondimensinal unit 1 is big enough....
xmax =  2.d0 
tmax =  8.d2 ! 

! Poincare section:  z=0 by observing the plot of the torus
!  
ind  = 3           ! x component 
p0   = 0.d0     ! rough guess, by observing directly at the plot of the tori  
dir  = 1           ! take the positive velocity 
imax = 1           ! consider every intersection, and choose by the direction of velocity


! -- compute the Poincare map   

! fix also the Poincare section, and take points from this section
pvi = 0.d0  
pvi(ind) = p0 
ispl = 0

fcurve   = 25;  fncurve  =  './dat/curve_lf/curve_poinc.dat'
open(fcurve  ,file = fncurve  , access ='append',status = 'replace' )
write(fcurve  ,*)  ' # Invariant curve by Poincare map:  t      (x y z vx vy vz)   cj  npc'

ftorus = 21;  fntorus    =  './dat/curve_lf/torusall.dat'
open(ftorus  ,file = fntorus, access ='append',  status = 'replace' )
write(ftorus ,*)  ' # Original  tori: t      (x y z vx vy vz)   cj'

! --- we are going to study the planar normal case, where z-direction has invariance

x0       = 0.61d0 
ind_curv = 0
do 
  
  ! for the moment, only compute 10 points
  if(ind_curv >= 10 ) exit
  
  ! we have z=0, vx=0, vz=0, 
  pvi = 0.d0 
  pvi(1) = x0 - 3.d-2 * ind_curv ! x
  pvi(2) = 0.d0          ! y
  pvi(ind) = p0          ! z
!  pvi(ndim+1:npvar:ndim+1) = 1.d0 
  
  print*,' pvi before cj2v_lf:', pvi; read* 
    
  call cj2v_lf(pvi, cj, cs, 5, isv)
 
  if(isv == 0) cycle  
  
  print*, 'isv, ind_curv, pvi', isv, ind_curv, pvi; read* 
  
!  ! check the torus 
!  if(ind_curv > 0) open(ftorus, file=fntorus, access='append', status='old')
!  call plob(pvi, 0.d0, 1.d3, 42, 1, ftorus, 1, gr_lf, gr_cjlf,  pvf) 
!  write(ftorus,*); write(ftorus,*)
!  close(ftorus)
!  read* 

  ind_curv = ind_curv +1
  
  
  open(fcurve, file=fncurve, access='append', status='old')
  do i = 1, npoinc

    call poinc(pvi, ndim, 42, 1, tf, pvf, hminim, ispc, gr_lf, gr_cjlf)
    if(ispc == 0) then 
      print*, 'Fail to reach to Poincare section!'
      exit 
    endif 
  
    if(i == 1) cycle 
    call gr_cjlf(pvi, cj_i)
    write(fcurve, *) tf, pvi(1:ndim), cj_i
  
    pvi = pvf 
  enddo 
  
  read* 
  
  write(fcurve, *); write(fcurve, *)
  close(fcurve)
enddo 

stop 

! files to save the original torus + invariant curve + interpolated one + approximated one 
!ftorus = 21;  fntorus    =  './dat/curve_lf/torusall.dat'
!open(ftorus  ,file = fntorus, access ='append',  status = 'replace' )
!write(ftorus ,*)  ' # Original 2D torus: t      (x y z vx vy vz)   cj'

!! first integrate the orbit for a long time, and  observe the geometry to select an appropriate Poincare section
!!subroutine plob(y0,t0,tf, n, tdir, ftag, ispl, deriv, gr_cj,  y) 

!tf = 1.d3 !2.d5 !1.d4 ! test a long time interval  -- already computed, do not overlap it....
!call plob(pvi, 0.d0, tf, 6, 1, ftorus, 1, gr_lf, gr_cjlf,  pvf) 
!write(ftorus,*); write(ftorus,*)
!stop 

!-- curve by computating Poincare map ----
!!do a Fourier analysis for the supposed torus

if (isfour == 1) then 
write(*,*) 'Maximum amplitude of the residual //  Tolerance in correction for refinement? '
print*, '( 1 : Default: 1.d-4, 1.d-6;  Otherwisze input new values)'

read(*,*) istol
if(istol == 1) then 
  tolres = 1.d-4; tol = 1.d-6
else 
 read(*,*)  tolres,  tol 
endif 


tf = 1.d3 
xmax =  2.d0 
nitmax = 15  
nsm = 2**18
nfmax = 25 

open(60  ,file = 'dat/curve_lf/ob.dat', access ='append',status = 'replace')
write(60, *) ' #   t      (x y z vx vy vz)    cj '

open(61  ,file = 'dat/curve_lf/fcs.dat', access ='append',status = 'replace')
open(62  ,file = 'dat/curve_lf/fcsbas.dat', access ='append',status = 'replace')

!subroutine fft_ob(x0, n, lt, np, tol, tolres, nitmax, nfmax, xmax, fob, ffc, ffbas, deriv, gr_cj) 
call fft_ob(pvi, 6, tf, nsm, tol, tolres, nitmax, nfmax, xmax, 60, 61, 62, gr_lf,  gr_cjlf)  

read*
 
! now I have two basic frequecies, with Maximum order (sum(abs(coefs)) = 51 
!    and    Tolerace for the freq identification = 1.d-5

!4.2050014392525877E-002
!  2.5727192937821428
w1 = 4.2050014392525877E-002
w2 = 2.5727192937821428
!dh = 2*pi/w2 
dh = 1.d0 / w2 


! take the one with short peroid, t2 = pi2/w1, and compute the Time T map  
open(222,  file = './dat/curve_lf/curve.dat') 
write(222, *) ' # t    (x,y,z, vx,vy,vz)  cj'

! ****************************** orbit sampling ******************************
! Time interval for discretisize the invariant curve with in one period along frequency w2 
! NOTE: be careful if we should multiply 2pi or not 
print*, 'start orbit sampling'; read*

!subroutine plob_fxd(y0, n, t0, h0, np, xmax, ftag, yall,  deriv, gr_cj) 
call plob_fxd(pv0, ndim, 0.d0, dh, npoinc, xmax, 222, pvall, gr_lf, gr_cjlf)   
print*, 'Finish the Time T map : curve.dat'  
read*  

endif 
! ------------------ fourier analysis for basic frequecies of the torus ------------------

              
!! ---------- Computation -----------------------------
! outer bound for escape, nondimensinal unit 1 is big enough....
xmax =  2.d0 
tmax =  5.d1 ! 

! Poincare section:  z=0 by observing the plot of the torus
!  
ind  = 3           ! x component 
p0   = 0.d0     ! rough guess, by observing directly at the plot of the tori  
dir  = 1           ! take the positive velocity 
imax = 1           ! consider every intersection, and choose by the direction of velocity


! -- compute the Poincare map   

! fix also the Poincare section, and take points from this section
pvi(ind) = p0 
ispl = 0

fcurve   = 25;  fncurve  =  './dat/curve_lf/curve_poinc.dat'
open(fcurve  ,file = fncurve  , access ='append',status = 'replace' )
write(fcurve  ,*)  ' # Invariant curve by Poincare map:  t      (x y z vx vy vz)   cj  npc'

do i = 1, npoinc

  call poinc(pvi, ndim, ndim, 1, tf, pvf, hminim, ispc, gr_lf, gr_cjlf)
  if(ispc == 0) then 
    print*, 'Fail to reach to Poincare section!'
    exit 
  endif 
  
  if(i == 1) cycle 
  call gr_cjlf(pvi, cj_i)
  write(fcurve, *) tf, pvi(1:ndim), cj_i
  
!  write(*, *) tf, pvi(1:ndim), cj_i; read*  
  pvi = pvf 
enddo 


write(fcurve, *);  write(fcurve, *)
!call poinc_n( pvi, ndim, ndim, ind, p0,  dir, fcurve, ispl, ftorus, npoinc, & 
!                   imax, tmax, xmax, tf, pvf, ncrs, gr_lf, gr_cjlf)   
print*, 'tf = ', tf 




stop 
   
!  ! add two blank lines to seperate different initial points
!    write(fcurve, *);   write(fcurve, *);   
!ispl = 1

!! first we try 10-by-10 points with different values of y and vy 
!np = 10 
!do i = 1, np, 1
!  ! the initial point --- instead of take random points, we fix the energy and poincare section 
!  !    
!  if(ind == 1) pvi(2) = (i-5)*5.d-3 +  yl4 
!  if(ind == 2) pvi(1) = (i-5)*5.d-3 +  xl4 
!  
!  ! the value of vx, step size 2.d-3
!  do j = 1, np, 1
!    pvi(3)  = (j-5)*5.d-3
!    
!    ! --  the value of vy, positive, index of vy is 4
!    call  cj2v_lf(pvi, cj0,  cs, 5, isv)
!    
!    if( isvy == 0 ) cycle
!    
!    
!  !  pvi = (/xl4 + 1.d-3, yl4 + 2.d-3, 1.d-4, -1.d-4/) ! add a small deviation on L4 to obtain an initial point

!    call poinc_n( pvi,ndim, ndim, 1, npoinc, tf, pvf, ncrs, gr_lf, gr_cjlf)   
!   
!  ! add two blank lines to seperate different initial points
!    write(fcurve, *);   write(fcurve, *);  
!    write(ftorus, *);   write(ftorus, *);
!  
!  end do  
!end do 

! close(fcurve) 
!stop  


!stop

end 


! 




program pofam_plf

! This is to compute the periodic orbits and the continuation along the energy level  in the planar LF problem of the normal case 

! NOTE: while po_plf deals with one family instead of the family, and their stability(eigenvalues of the monodromy matrix)

! using the initial guess obtained by poinc_plf, pick the second map where vx=0
! which gives us v0 =0.319770, retro(vy  > 0)

! **NOTE***
!  y0 is only 2d, free parameters on the Poincare Map, so becareful 

! TODO:
! --1. refine the periodic orbit
! --2. do the continuation by changing energy 
! --3. update the continuation by using the parameters of the previous two p.o.s
!      introduce a 2 rows array to save the data, when the guess fails, we still have the data of the previous two 
!      since if we only use x0_pre, and cj_pre, only one previous p.o. remains when the continuation fails.
! --4. study also the stability of the po family

!  --- data saved in ./dat/map_plf/  subfolder 

!  1. po.dat         -- the refined periodic orbits
!  2. 
!  2. poit.dat       -- the intermediate orbits during refinement 
  
use dp_mod
use pi_mod 
use poinc_mod
use plf_mod 
use po_plf_mod
implicit none

integer, parameter  ::  ndim   = 4,  &  ! for lf
                        npvar  = 20, &  ! 20 for variational matrix, 4 for state vector 
                        npp    = 100  
! Local Variables
integer       ::  debug, i, k
real(kind=dp) ::  beta0, cj_eq 

! eigenvalues and eigenvectors around the equilibria                 
real(kind=dp) ::  dlf(ndim, ndim), wr(ndim), wi(ndim), vr(ndim, ndim), phi(ndim, ndim)


! Poincare map 
integer       :: ind0, dir0, imax0, ind_free
real(kind=dp) :: p00,  tmax0, xmax0, cj0 

! Po 
integer ::  ispo, isv, ispl, niter, ipo, fob, finit, feig
real(kind=dp) ::  y0(ndim-2), TP, yf(ndim), dg(1,1), cj, kern_dg , x0,  &
                  x0_pre(2), cj0_pre(2), dcj, pvi(npvar), pvf(npvar), tf, dhmax, dhmin
                  
character(len=70) :: fnob, fninit, fneig

! for fft_ob 
integer           ::  istol, nfmax, nsm,  nitmax, isfour 
real(kind=dp)     ::  tolres,  tolfc 

! take eq1  in planar normal case as an example, type center X saddle  equilibria
! with beta = 1.d0, there is resonance when beta is around 1, avoid this for the moment 
debug = 0

! limits of step size for continuation
dhmax =  2.d1 !1.d-1  ! for y
!dhmax = 1.d-1 
dhmin = 1.d-7 
 
! the variation in energy, we need to decrease to go from the closed region to open one
dcj  =  -1.d-3  
    
! to assign the value of  eq     
call init_plf     
 
beta0 = 2.d0   ! use 2 to check 
call init_beta(beta0)   
     

call gr_cjplf(eq, cj_eq)
print*, 'check energy! ', cj_eq, eq

 
! Take initial point on the x-axis, with velocity perpendicular to x-axis, and compute the P.O. 
! by Poincare map with the section as y=0, take vy as a function of f(cj, x, y=0, vx)
! so we have 2 dimension of freedom (x, vx )

!! ---------- Computation -----------------------------
! For poincare map limitation 
!xmax0 =  2.d2  !
!tmax0 =  4.d2 ! 
xmax0 =  6.d4  !
tmax0 =  4.d4 ! 

! Poincare section:  y=0  
!ind0     = 1           ! x(y) component 
!ind_free = 1        ! doubly-symmetric p.o. 

ind0  = 1            ! x=0,  symmetric w.r.t. y-axis 
!ind0  = 2           ! y=0,  symmetric w.r.t. x-axis
  
ind_free = 3-ind0   ! y(x) component, the free parameter
p00   = 0.d0        ! rough guess, by observing directly at the plot of the tori  
dir0  = 0           ! no requirement in the direction, just compute the first turn 
imax0 = 1           ! consider every intersection, and choose by the direction of velocity

!subroutine init_poinc( ind0, p00, dir0, imax0, tmax0, xmax0, ndim0, h00) 
call init_poinc(ind0, p00, dir0, imax0, tmax0, xmax0, ndim, cj0)

ispl = 1

! ---------- compute periodic orbit --------------------------------
! by observing the first two maps, we find the second return with vy>0, retrograde orbits

!for x0, p.o. symmetric w.r.t. x-axis
!x0 = 0.319776d0 
!  0.12425502665701038       0.48572777020884461       -6.4623485355705287E-025  -4.3730346033019958E-015 -0.96200794605542128        3.8998678590602029     

! another family, the bifuraction detected 
!x0  = 0.45143999792d0 
! cj0 = 3.8270713160552936d0
! 
!! check if there r two families, at energy level  cj0 = 3.8270713160552936 + 1.d-3 
!! x0_1 = 0.4465
!!  x0_2 = 0.45863995
!  
!  cj0 = 3.8270713160552936 + 1.d-3
!  x0  = 0.4465
!  x0  = 0.45863995
!   the energy level is within the interval [4.320714 : 3.827071361]
 ! periodic orbit in the choatic region

! for y0, p.o. symmetric w.r.t. y-axis
!x0 = 0.36714425247  
! ! where we stop and cannot continue further  
!!   0.54362985216868198        0.0000000000000000       0.70000004456192400       0.92987365800236521        0.0000000000000000        1.9924776554107650     
! x0 = 0.70000004456192400
! cj0 =  1.9924776554107650d0

!subroutine init_po_plf( ind_zero )
! we want to take vx0 = 0, which is the second element in (x, vx), so ind_zero=2
call init_po_plf( (/2/) )

fob = 100; finit = 200; feig = 88
! better to use the general name, and copy to the specified folder later
write(fnob,    fmt='(a)')  './dat/map_plf/po/poy_ob.dat'   ! general name for   p.o. orbits
write (fninit, fmt='(a)')  './dat/map_plf/po/poy_init.dat' ! I.C.
write (fneig,  fmt='(a)')  './dat/map_plf/po/poy_eig.dat'  ! eigenvalues values

!write(fnob, fmt='(a)')      './dat/map_plf/poxy_ob.dat' !general name for 1 p.o.
!write (fninit,  fmt='(a)')  './dat/map_plf/poxy_init.dat'

!write(fnob, fmt='(a)')      './dat/map_plf/h3/poy_ob.dat'
!write (fninit,  fmt='(a)')  './dat/map_plf/h3/poy_init.dat'

open(fob,   file = fnob,   status='replace', access='append')
open(finit, file = fninit, status='replace', access='append')
open(feig,  file = fneig,  status='replace', access='append')


! --- Initialization for the first p.o. ------------ 
ipo  = 0 ! p.o. counter
ispo = 0 ! flag of success of the refinement 

! y0 is the initial point on the Poincare section, (x0, vx0 = 0.d0)
y0 = 0.d0;      y0(1) = x0 
x0_pre = 0.d0;  cj0_pre = 0.d0 

! we take one that is symmetric w.r.t. to both x-axis and y-axis 

! we failed using only 1 symmetry, and end up with the previously detected p.o. x0=0.32, although we want to start with 0.24 
!x0 = 0.24d0     ! npoinc = 35  !

x0 = 0.3194d0   ! -- for y0 in the same energy level. 

! Type 3 at H3 energy level for pox symmetric w.r.t. x-axis
!x0 = 0.6182491304E+00   ! good one, outermost one
!x0 = 0.602d0        ! 2 with loops
!x0 = 0.584d0        ! 3 with loops
!x0 = 0.5705d0      ! 4 with loops
!x0 = 0.5604d0       ! 5 with loops
!x0 = 0.5516d0       ! 6 with loops
!x0 = 0.544d0       ! 7 with loops
!x0 = 0.537d0       ! 8 with loops
!x0 = 0.531d0       ! 9 with loops
!x0 = 0.5255d0       ! 10 with loops
!x0 = 0.5204d0       ! 11 with loops


 cj0 =  cj_eq + 0.05d0
print*, 'cj0 = ', cj0, '= 4.3767487109222252 ' ; read*
dcj = -dabs(dcj) ! we need increase the energy, go backwards

call init_h0(cj0)
print*, 'check the initial energy, cj_eq: ', cj0, cj_eq; read*

!4.3767487109222252 
! -------------- continuation---------------------------
do  
!  if(ipo > 500) exit

!    if(ipo > 0) exit

  ! update the energy and the free parameter   
  ! y0 = (x0, vx0=0), the initial point on the Poincare section, with zero velocity(in the direction of the free parameter)
  
  call init_h0(cj0) 
  y0 = 0.d0  
  y0(1) = x0    !x0
!  y0(3) = 1.d-1 ! vx
   
!   print*, 'y0 =', y0; read*
  ! refinement using Newton method 
  call refn_po(y0, TP, yf, dg, niter, ispo) 
  
  if(ispo == 0)  then 
    dcj = dcj / 2
    print*, 'Too big step in continuation, fail! Decrease by half', dcj
    
    if (dabs(dcj) < 1.d-8) then 
      print*, 'Too small step of H for continuation! <1.d-8, stop the continuation', dcj;
      read* 
      exit
    endif  
    
    ! use the two previous good p.o.s for predicition with smaller stepsize
    call pred_h(x0_pre, cj0_pre, dcj, ipo, x0, cj0)
    cycle 
    
  else 
    
    ! we obtain a new periodic orbit
    ipo = ipo + 1
    x0 = y0(1) 
    tp = 2.d0*tp
    
    ! -- update x0_pre and cj0_pre 
    if(ipo > 2) then 
      x0_pre(1)  = x0_pre(2)
      cj0_pre(1) = cj0_pre(2)
      
      x0_pre(2)  = x0
      cj0_pre(2) = cj0 
    else 
      x0_pre(ipo)  = x0
      cj0_pre(ipo) = cj0
    endif 
    
    ! -- modify the stepsize based on the iterations of Newton method 
    if(niter < 3)  dcj = dmin1( dabs(2*dcj), dhmax) * dcj / dabs(dcj)
    if(niter > 5)  dcj = dcj / 2 
    
    ! --- plot po for 1 period -----
    pvi = 0.d0
    pvi(ind_free)= y0(1) 
    call cj2v_plf(pvi(1:ndim), cj0, ind+2, isv)
    
    print*, 'check the energy:', cj0 
    print*, 'After the refinment of P.O., tp =', tp
    print*, 'y0', pvi(1:ndim)
    print*, 'yf', yf(1:ndim)
!    print*; read*
    
    open(finit, file=fninit, status='old', access='append')
    write(finit, *) tp, pvi(1:ndim), cj0; close(200)
    
!    write(*, *) tp, pvi(1:ndim), cj0 ; print*; read*
    
    open(fob, file=fnob, status='old', access='append')
    
    
    ! --- plot the orbit for one revolution and compute the Eigenvalues of the Monodromy Matrix to check the stability
    pvi(5:20) = 0.d0
    pvi(5:20:5) = 1.d0
    call plob(pvi, 0.d0, 1*TP, ndim, npvar, 1, 1, fob, gr_plf, gr_cjplf,  pvf) 
    
    phi = reshape(pvf(5:20), (/4,4/))
    
!    do i = 1, ndim, 1
!      print*, phi(i,:)
!    end do
!    print*; read*
    
    call eigrg(phi, ndim,1, wr, wi, vr)
 
   ! Save the eigenvalues to files + x0 as the first parameter  
    write(feig, '(2e24.14)', advance='no')  x0, cj
    do i = 1, ndim, 1
      write(feig, '(2e24.14)', advance='no')  wr( i ), wi( i ) 
    end do  
    write(feig, *)
    
!    call prt_eigval( ndim, feig, wr, wi )
    
    
    ! --- we only integrate the orbit 
!    call plob(pvi, 0.d0, TP, ndim, 1, 100, 1, gr_plf, gr_cjplf,  pvf) 
!    plob(y0, t0, tf, n, tdir, ftag, ispl, deriv, gr_cj,  y) 
    
     close(fob); !read*
    
     ! ---- continuation ------------
    ! since we have only 1 free parameter to continue along: x 
    ! we take a step forward and update the energy level  
  
    ! ---- prediction for the guess of a new p.o -------
    !  using linear extrapolation based on the previous 2 p.o.s
    call pred_h(x0_pre, cj0_pre, dcj, ipo, x0, cj0)
    
    if( debug == 1 ) then
      print*, 'Initial guess for', ipo, '-th new p.o.:'
      print*, x0, cj0 ; read*
    end if
    
  endif  ! ispo == 1
 
enddo 

 close(feig)
stop 




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

open(60  ,file = 'dat/map_plf/ob.dat', access ='append',status = 'replace')
write(60, *) ' #   t      (x y z vx vy vz)    cj '

open(61  ,file = 'dat/map_plf/fcs.dat', access ='append',status = 'replace')
open(62  ,file = 'dat/map_plf/fcsbas.dat', access ='append',status = 'replace')

!subroutine fft_ob(x0, n, lt, np, tol, tolres, nitmax, nfmax, xmax, fob, ffc, ffbas, deriv, gr_cj) 
call fft_ob(pvi, ndim, tf, nsm, tolfc, tolres, nitmax, nfmax, xmax, 60, 61, 62, gr_plf,  gr_cjplf)  

read*
 
 

stop 

 
 
end 
 



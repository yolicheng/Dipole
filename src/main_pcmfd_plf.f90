program pcmfd_plf

! TODO 2016-12-01 09:16:49 
! check the first poincare map, it should be a loop, For the idea of the intersections of the 
! manifold, we have to take several points from the initial p.o. as the initial point of the mfd
! the idea of taking one point and do a lot of intersections is wrong! 

! This is to compute the first few returns of the invariant manfold to Poincare section y=0 (x-axis)
! just to check the manifolds as the seperatrice of different motions

! Take 100 points as the initial manifold, try the first 10 returns

!  --- data saved in ./dat/pcmfd_plf/h  subfolder  ! for general use, and after the check, we can copy the files to the specified energy level  

!  1. mfd_init**.dat   --  I.C. of the manifold  
!  2. pc**.dat         --  Poincare map of the manifold cut by Poincare section
!  2. ob**.dat         --  orbits on the manifold 
 
use dp_mod
use pi_mod 
use poinc_mod
use plf_mod 
implicit none

integer, parameter  ::  ndim   = 4,    &    ! for lf
                        npvar  = 20,   &    ! 20 to compute MM and mfd 
                        npoinc = 20,   &    ! number of returns to Poincare section  
                        nmf    = 200       ! number of orbits on the manifold for the poincare map 

! Local Variables
integer       ::  debug, i, k, imfd_debug
real(kind=dp) ::  beta0, dlf(ndim, ndim)

! Poincare map 
integer       :: ind0, dir0, imax0, ind_free   
! 
real(kind=dp)     :: p00, tmax0, xmax0, cj0 
character(len=70) :: fnmap 

! Poinc_n 
integer       ::  ncrs  
real(kind=dp) ::  yi(ndim), tf, hminim , cj, ylim

! P.O. 
integer           ::  fob, fpo_init, istop
real(kind=dp)     ::  y0(ndim), TP,  yf(ndim) 
character(len=70) ::  fnob, fnpo_init 

! mfd -- Eigenvalues and eigenvectors of MM       
integer         ::  col, finit, fmfd_pc, fmfd, ispl, tdir     
real(kind=dp)   ::  pvi(npvar), pvf(npvar), phi(ndim, ndim), wr(ndim), wi(ndim), vr(ndim, ndim), &
                    vepu(ndim), veps(ndim), ymfd(nmf, ndim),  epsl, yf2(ndim), &
                    yfmfd(nmf, ndim), yn(nmf, ndim), dt, yf3(ndim), t0 ! ck plob


! take eq1  in planar normal case as an example, type center X saddle  equilibria
! with beta = 1.d0, there is resonance when beta is around 1, avoid this for the moment 
debug = 0

epsl  = 1.d-6 ! the commonly used one

!epsl = .5d-6 ! use smaller ones to check 

call init_plf     ! to assign the value of  eq 

! TODO: the value of beta matters! Be careful to use the same value when we do exploration 
!       At this moment, we take beta=2 

beta0 = 2.d0   ! use 2 to check, but 1 is better?  

print*, 'Check the value of beta (= 2?)', beta0; read*; print*

!print*, 'Input the value of beta (n/w_c):'
!read*, beta0 

call init_beta(beta0)   
print*, 'check,  beta', beta0 

call gr_cjplf(eq, cj)
print*, 'check energy! ', cj, eq

! -- the energy level we r taking, smaller than the one at the equilibria, so we have a 
! closed UFO 
! for x ----
 cj0 = cj + 0.05d0
 print*, 'check the energy! H0 = ', cj, ' H3 =', cj0
 print*, 'cj0 =', cj0; read*
 
 ! 4.3267487109222253  for equilibirum in plf
 ! 4.3767487109222252  for h3 !--- the current studied one --- 
  
 ! ** NOTE ** 
 ! ----- for p.o. symmetric to x-axis, we have ....
 !  the energy level is within the interval [4.320714 : 3.827071361], we have  two families of periodic orbits 
 
 !  for energy above 4.320714,  one family exists and goes to the origin, the other family diminish, we have choatic region
 !  this case is what we are studing now. H3 = 4.3767487109222252 
 
 ! for energy level less than  3.827071361, there will be no periodic orbit.

 ! and we also have one family that is symmetric w.r.t. y-axis, which goes to infinity....
 
! cj0 = 4.1d0 
  
! cj0 = 3.8270713160552936          ! the one the continutation stops 
 
! cj0 = 3.8270713160552936 - 1.d-5  ! the one there is no intersection? 
 
! cj0 = 3.8270713160552936 + 1.d-3  ! the one with 2 periodic orbits 
   
!! try another one, where the continuation stops 
!  cj0 =  3.8270713125978135    ! 3.8270713125975111
!! 
!!! try a even smaller one 
!! cj0 =  3.8270713125978135 - 1.d-5
! 
! ! for continuation of Type 2 p.o. 
!!  cj0 =  1.9970714282989488 -  1.d0 !5.d-1 !- 1.d-6    
!  cj0 = 1.9924776554107650
read*

    
! -- the real units-----
print*, 'runit(km), tunit(s), vuint(km/s): ', runit, tunit, vunit

! Jacobi matrix of the lorentz force with respective to the state 
call dflrtz2(eq, dlf)

!print*, 'Check DG in the planar case: ' 
!do i = 1, ndim 
!  write(*,'(4f18.8)') dlf(i,:) 
!enddo
!print*; read* 

call eigrg(dlf, ndim,1, wr,wi, vr)
print*, 'Check the type of equilibria!' 
print*; read* 

! Take initial point on the x-axis, with velocity perpendicular to x-axis, and compute the P.O. 
! by Poincare map with the section as y=0, take vy as a function of f(cj, x, y=0, vx)
! so we have 2 dimension of freedom (x, vx )


!! ---------- Computation -----------------------------
! For poincare map limitation 
xmax0 =  2.d0  !
tmax0 =  4.d2 ! 

! TODO: Poincare section:  x=0
!print*, 'check the value of ind, 1 for p.o. symmetric w.r.t. x-axis' 
!ind0  = 1       !  y-axis, x = 0
ind0  = 2       !  x-axis, y = 0

p00   = 0.d0        ! rough guess, by observing directly at the plot of the tori  

! TODO: NOTE: the first crossing is 
!dir0  = 1          ! For this purpose, we only take the plus crossing, use dir = 1
print*, 'Input the direction of the intersection: 1 (v>0) -1: (v<0)'
read*, dir0 

imax0 = 1           ! consider every intersection, and choose by the direction of velocity

!subroutine init_poinc( ind0, p00, dir0, imax0, tmax0, xmax0, ndim0, h00) 
call init_poinc(ind0, p00, dir0, imax0, tmax0, xmax0, ndim, cj0)

! TODO: remember to modify the name for specific initial condition
fpo_init = 200  
write(fnpo_init,  fmt='(a)')  './dat/pcmfd_plf/h3/poy_init.dat' ! 

! do the transformation: (x,y,vx,vy) --> (x, -y, -vx, vy)


! Read from the refined periodic orbits.
open(fpo_init, file = fnpo_init, status='old')
print*; print*, 'Energy level: ', cj0 

read(fpo_init, *)  tp, y0(1:ndim), cj0 

! do the transformation: (x,y,vx,vy) --> (x, -y, -vx, vy)
print*, 'Input 1: Take the top p.o., others: the bottom one'
read(*,*) istop 
if(istop /= 1)  y0(2:3) = -y0(2:3)

write(*, *)     'I.C. of p.o.: tp, y0, cj0'
print*, tp, y0(1:ndim), cj0;   read* ! ck


! -- the initial condition for the manifolds, use different names  for different manifolds
! --- compute the manifold, using the eigenvector of Monodromy Matrix as linear approximation
! --- plot the orbit for one revolution and compute the Eigenvalues of the Monodromy Matrix to check the stability
pvi(1:ndim) = y0  
pvi(5:20)   = 0.d0
pvi(5:20:5) = 1.d0

! do not plot, only use this for Monodromy matrix computation
fob = 532
open(fob, file = './dat/pcmfd_plf/h/po2.dat', status = 'replace')
call plob(pvi, 0.d0, TP, ndim, npvar, 1, 1, fob, gr_plf, gr_cjplf,  pvf) 
 close(fob) 
 
! print*, 'check the periodic orbit!'  !--ckd!
!stop 
 
! check the final state 
print*, 'I.C. :', pvi(1:4) !ck
print*, 'F.C. :', pvf(1:4) 
print*; read*

! Monodromy Matrix  
phi   = reshape(pvf(5:20), (/4,4/))
print*, 'Monodromy matrix: '
do i = 1, ndim, 1
   print*, phi(i,:)
end do
print*; read*
    
call eigrg(phi, ndim, 1, wr,wi, vr)
print*, 'Pick the dominant unstable eigenvector ( Re(lambda_max) > 1)! Input the column'
read(*, *)  col 
vepu =  vr(:, col)
print*, 'eigenvector :', vepu; print*; read*

print*, 'Pick the stable eigenvector ( 1 / Re(lambda_max) < 1)! Input the column'
read(*, *)  col 
veps = vr(:, col)
print*, 'eigenvector :', veps; print*; read*


! Compute the  initial conditions of manifold  and the first npoinc returns together with the orbit 
!subroutine gr_mfdinst(ypo, n, tp, mmat, nmf, epsl, ftag, ymfd, deriv, gr_cj)
!epsl = 1.d-5 or 1.5-6 ! doesn't matter too much... 

finit = 100; fmfd_pc = 110; fmfd = 120  

! --- right half of stable manifold of the downstairs p.o. ---- 
tdir = -1 

! -- I.C. 
open(finit, file = './dat/pcmfd_plf/h/mfd_initsp.dat', status='replace', access='append')
write(finit,*) '# I.C. of manifold  WS+: dt    (x, y, vx, vy)   cj'

epsl = -dabs(epsl)
call gr_mfdinst(y0,  ndim, tp, veps, nmf, epsl, finit, ymfd, gr_plf, gr_cjplf) 
 close(finit)
 print*, 'Finish WS+'; print*; read*; !ck
 
! -- maps for manifold 
open(fmfd_pc, file = './dat/pcmfd_plf/h/pcsp.dat', status='replace', access='append')
write(fmfd_pc,*) '# Map of manifold WS+:  tpc (x, y, vx, vy)   ncrs'

! --do the poincare map for the first npoinc returns for all the nmf points
do k = 1, nmf, 1 
!subroutine poinc_n(xi, ndim, np, fpc, npoinc, tf,yf, ncrs, deriv, gr_cj) 
  yi = ymfd(k, 1:ndim)
  call  poinc_n(yi, ndim, ndim, tdir, fmfd_pc, npoinc, tf, yf, ncrs, gr_plf, gr_cjplf) 
enddo  
 close(fmfd_pc)
stop  



! --- right half of unstable manifold ----
tdir = 1 
open(finit, file = './dat/pcmfd_plf/h/mfd_initup.dat', status='replace', access='append')
write(finit,*) '# I.C. of manifold  WU+:  dt    (x, y, vx, vy)   cj'

epsl = dabs(epsl)
call gr_mfdinst(y0, ndim, tp, vepu, nmf, epsl, finit, ymfd, gr_plf, gr_cjplf) 
 close(finit)

!  ** NOTE** 
!  we only plot |y| < 0.1, when the orbit is close to x-axis, to avoid huge data file 
ylim = 2.d-1    
ispl = 0
if( ispl == 1 ) then
  open(fmfd, file = './dat/pcmfd_plf/h/obup.dat', status='replace', access='append')
  write(fmfd,*) '# Orbits on the manifold  WU+:  t (x, y, vx, vy)'
end if

open(fmfd_pc, file = './dat/pcmfd_plf/h/pcup.dat', status='replace', access='append')
write(fmfd_pc,*) '# Map of manifold  WU+:  tpc (x, y, vx, vy)   ncrs'

!--------------ckd!! debuged k=403 for iter>10 ------------------------
!  403-- 1.2250807938099942E-002  0.28434843467464005       -1.2673055018281816       -1.0221269553536949     
! remove the appendix 2 after the debug
!imfd_debug = 92 !1009 
!do i = imfd_debug, 95, 1
!  yi = ymfd(i, 1:ndim); tf = 105.d0
!  call plob(yi, 0.d0, tf, ndim, ndim, 1, ispl,  fmfd, gr_plf, gr_cjplf, yf2) 
!end do
!stop
!--------------

! do the poincare map for the first 20 returns for all the nmf points with vy>0

! do k = imfd_debug, nmf, 1   ! debug for iter>10, k=403
!   if(k > 96) stop

ispl = 0
do k = 1, nmf, 1 

  yi = ymfd(k, 1:ndim)
  call  poinc_n(yi, ndim, ndim, tdir, fmfd_pc, npoinc, tf, yf, ncrs, gr_plf, gr_cjplf) 
  
  if( ncrs == 0 ) then 
    print*, 'Fail to obtain the crossing!'
    print*, k, '-th  orbit on the manifold'
    print*, yi(1:ndim) 
    print*; read*
  endif 
    
  if(ispl == 1) then 
    call plob_y0(yi, 0.d0, tf, ndim, ndim, 1, ispl, ind0, ylim, fmfd, gr_plf, gr_cjplf, yf2) 
  endif 
enddo  

 close(fmfd_pc)
print*, 'finish WU+'; read*

stop 

! --- right half of stable manifold of the downstairs p.o. ---- 
tdir = -1 

! -- I.C. 
open(finit, file = './dat/pcmfd_plf/h/mfd_initsp.dat', status='replace', access='append')
write(finit,*) '# I.C. of manifold  WS+: dt    (x, y, vx, vy)   cj'

epsl = -dabs(epsl)
call gr_mfdinst(y0,  ndim, tp, veps, nmf, epsl, finit, ymfd, gr_plf, gr_cjplf) 
 close(finit)
 print*, 'Finish WS+'; print*; read*; !ck
 
! -- maps for manifold 
open(fmfd_pc, file = './dat/pcmfd_plf/h/pcsp.dat', status='replace', access='append')
write(fmfd_pc,*) '# Map of manifold WS+:  tpc (x, y, vx, vy)   ncrs'

! --do the poincare map for the first npoinc returns for all the nmf points
do k = 1, nmf, 1 
!subroutine poinc_n(xi, ndim, np, fpc, npoinc, tf,yf, ncrs, deriv, gr_cj) 
  yi = ymfd(k, 1:ndim)
  call  poinc_n(yi, ndim, ndim, tdir, fmfd_pc, npoinc, tf, yf, ncrs, gr_plf, gr_cjplf) 
enddo  
 close(fmfd_pc)
!stop  




! --- left half of unstable manifold ---- 
tdir = 1
! -- I.C. 
open(finit, file = './dat/pcmfd_plf/h/mfd_initum.dat', status='replace', access='append')
write(finit,*) '# I.C. of manifold WU-:  dt   (x, y, vx, vy)   cj'

! -- compute nmf = 100 initial points on the manifold
epsl = -dabs(epsl)
call gr_mfdinst(y0,  ndim, tp, vepu, nmf, epsl, finit, ymfd, gr_plf, gr_cjplf) 
 close(finit)
print*, 'finish WU-'; read*

! -- manifold 
open(fmfd_pc, file = './dat/pcmfd_plf/h/pcum.dat', status='replace', access='append')
write(fmfd_pc,*) '# Map of manifold WU-:  tpc (x, y, vx, vy)   ncrs'

! do the poincare map for the first 10 returns for all the nmf points
do k = 1, nmf, 1 
!subroutine poinc_n(xi, ndim, np, fpc, npoinc, tf,xf, ncrs, deriv, gr_cj) 
  yi = ymfd(k, 1:ndim)
  call  poinc_n(yi, ndim, ndim, tdir, fmfd_pc, npoinc, tf, yf, ncrs, gr_plf, gr_cjplf) 
enddo  
 close(fmfd_pc)



! --- left half of  stable manifold ---- 
tdir = -1
! -- I.C.
open(finit, file = './dat/pcmfd_plf/h/mfd_initsm.dat', status='replace', access='append')
write(finit,*) '# I.C. of manifold WS-: dt   (x, y, vx, vy)   cj'
! -- compute nmf = 100 initial points on the manifold

epsl = dabs(epsl)
call gr_mfdinst(y0,  ndim, tp, veps, nmf, epsl, finit, ymfd, gr_plf, gr_cjplf) 
 close(finit)
print*, 'Finish WS-'; print*; read*; !ck

! -- maps for manifold 
open(fmfd_pc, file = './dat/pcmfd_plf/h/pcsm.dat', status='replace', access='append')
write(fmfd_pc,*) '# Map of manifold WS-:  tpc (x, y, vx, vy)   ncrs'

! do the poincare map for the first 10 returns for all the nmf points
do k = 1, nmf, 1 
!subroutine poinc_n(xi, ndim, np, fpc, npoinc, tf,xf, ncrs, deriv, gr_cj) 
  yi = ymfd(k, 1:ndim)
  call  poinc_n(yi, ndim, ndim, tdir, fmfd_pc, npoinc, tf, yf, ncrs, gr_plf, gr_cjplf) 
enddo  

 close(fmfd_pc) 
 
end 
 



program mfd_plf

! 2016-12-15 16:15:59 
! In order not to mix the computation of the manifold with the poincare map, 
! we do a standalone routine only for the  computation of the manifold, 
! 16 orbits will be enough 

! This is only for the manifold, to avoid huge data file, we only take 16 initial points

!  --- data saved in ./dat/mfd_plf/h  subfolder  
! ! for general use, and after the check, we can copy the files to the specified energy level  

!  1. mfd_init**.dat   --  I.C. of the manifold  
!  2. ob**.dat         --  orbits on the manifold 
 
use dp_mod
use pi_mod 
use plf_mod 
implicit none

integer, parameter  ::  ndim   = 4,    &    ! for lf
                        npvar  = 20,   &    ! 20 to compute MM and mfd 
                        nmf    = 8       ! number of orbits on the manifold  for the shape 
                      
! Local Variables
integer       ::  debug, i, k, imfd_debug
real(kind=dp) ::  beta0, dlf(ndim, ndim), cj_eq, cj0

! P.O. 
integer           ::  fob, fpo_init, istop
real(kind=dp)     ::  y0(ndim), TP,  yf(ndim) 
character(len=70) ::  fnob, fnpo_init 

! mfd -- Eigenvalues and eigenvectors of MM       
integer         ::  col, finit, fmfd_pc, fmfd, ispl, tdir     
real(kind=dp)   ::  cj, yi(ndim), t0, tf, pvi(npvar), pvf(npvar), phi(ndim, ndim), wr(ndim), wi(ndim), vr(ndim, ndim), &
                    vepu(ndim), veps(ndim), ymfd(nmf, ndim),  epsl, yf2(ndim) 


! take eq1 in planar normal case as an example, type center X saddle  equilibria
! with beta = 1.d0, there is resonance when beta is around 1, avoid this for the moment 
debug = 0

epsl  = 1.d-6 ! the commonly used one

!epsl = .5d-6 ! use smaller ones to check, obtain the same result as with 1.d-6

call init_plf     ! to assign the value of  eq 

! TODO: the value of beta matters! Be careful to use the same value when we do exploration 
!       At this moment, we take beta=2 

beta0 = 2.d0   ! use 2 to check, but 1 is better?  

print*, 'Check the value of beta (= 2?)', beta0; read*; print*

!print*, 'Input the value of beta (n/w_c):'
!read*, beta0 

call init_beta(beta0)   
print*, 'check,  beta', beta0 

call gr_cjplf(eq, cj_eq)
print*, 'check energy! ', cj_eq, eq

! -- the energy level we r taking, smaller than the one at the equilibria, so we have a 
! closed UFO 
! for x ----
 cj0 = cj_eq + 0.05d0
 print*, 'check the energy! H0 = ', cj_eq, ' H3 =', cj0
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

! TODO: remember to modify the name for specific initial condition
fpo_init = 200  
write(fnpo_init,  fmt='(a)')  './dat/map_plf/h3/poy_init.dat' ! 

! do the transformation: (x,y,vx,vy) --> (x, -y, -vx, vy)


! Read from the refined periodic orbits.
open(fpo_init, file = fnpo_init, status='old')
read(fpo_init, *)  tp, y0, cj0 

! do the transformation: (x,y,vx,vy) --> (x, -y, -vx, vy)
print*, 'Input 1: Take the top p.o., others: the bottom one'
read(*,*) istop 
if(istop /= 1) y0(2:3) = -y0(2:3)


write(*, *)         'I.C. of p.o.: tp, y0, cj0'
print*, tp, y0, cj0;   read* ! ck


! -- the initial condition for the manifolds, use different names  for different manifolds
! --- compute the manifold, using the eigenvector of Monodromy Matrix as linear approximation
! --- plot the orbit for one revolution and compute the Eigenvalues of the Monodromy Matrix to check the stability
pvi(1:ndim) = y0  
pvi(5:20)   = 0.d0
pvi(5:20:5) = 1.d0

! do not plot, only use this for Monodromy matrix computation
fob = 532
open(fob, file = './dat/mfd_plf/h/po2.dat', status = 'replace')
call plob(pvi, 0.d0, TP, ndim, npvar, 1, 1, fob, gr_plf, gr_cjplf,  pvf) 
 close(fob) 
 
! print*, 'check the periodic orbit!'  !--ckd!
!stop 
 
 
! check the final state 
print*, 'I.C. :', pvi(1:4) !ck
print*, 'F.C. :', pvf(1:4) 
print*; read*

! Monodromy Matrix  
phi = reshape(pvf(5:20), (/4,4/))
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


finit = 100;   fmfd = 120  

t0 = 0.d0 !75.d0
tf = 150.d0
  
! --- right half of stable manifold of the downstairs p.o. ---- 
tdir = -1 
! -- I.C. 
open(finit, file = './dat/mfd_plf/h/mfd_initsp8.dat', status='replace', access='append')
write(finit,*) '# I.C. of manifold  WS+: dt    (x, y, vx, vy)   cj'

! -- compute nmf = 100 initial points on the manifold
epsl = -dabs(epsl)
call gr_mfdinst(y0,  ndim, tp, veps, nmf, epsl, finit, ymfd, gr_plf, gr_cjplf) 
 close(finit)
 print*, 'Finish WS+'; print*; read*; !ck
 
! --- plot the manifold of 16 orbits, take tf to cover  
ispl = 1
print*, 'Plot the orbit of WS+'
if( ispl == 1 ) then
  open(fmfd, file = './dat/mfd_plf/h/obsp8.dat', status='replace', access='append')
  write(fmfd,*) '# Orbits on the manifold  WS+:  t (x, y, vx, vy)'
  
 ! check the first 12 Poincare maps which complete a curve, tf is more or less 105
 ! we do not need to compute the whole curve, just the ones close to the x-axis, we take t \in [70-130]
  do k = 1, nmf, 1 
    yi = ymfd(k, 1:ndim)
!    call plob(yi,  0.d0, t0, ndim, ndim, tdir, 0, 6, gr_plf, gr_cjplf,   yf2) ! do not plot 
    call plob(yi,  t0, tf, ndim, ndim, tdir, 1, fmfd, gr_plf, gr_cjplf, yf2) 
  enddo 
end if

!stop
 
 
! --- right half of unstable manifold ----
tdir = 1 
open(finit, file = './dat/mfd_plf/h/mfd_initup8.dat', status='replace', access='append')
write(finit,*) '# I.C. of manifold  WU+:  dt    (x, y, vx, vy)   cj'

! -- compute nmf = 100 initial points on the manifold
epsl = dabs(epsl)
call gr_mfdinst(y0, ndim, tp, vepu, nmf, epsl, finit, ymfd, gr_plf, gr_cjplf) 
 close(finit)

! --- plot the manifold of 16 orbits, take tf to cover  
ispl = 1
print*, 'Plot the orbit of WU+'
if( ispl == 1 ) then
  open(fmfd, file = './dat/mfd_plf/h/obup8.dat', status='replace', access='append')
  write(fmfd,*) '# Orbits on the manifold  WU+:  t (x, y, vx, vy)'
  do k = 1, nmf, 1 
    yi = ymfd(k, 1:ndim)
    call plob(yi, t0, tf, ndim, ndim, tdir, 1, fmfd, gr_plf, gr_cjplf, yf2) 
  enddo 
end if

print*, 'finish WU+'; read*
!stop 

! --- left half of unstable manifold ---- 
tdir = 1
! -- I.C. 
open(finit, file = './dat/mfd_plf/h/mfd_initum8.dat', status='replace', access='append')
write(finit,*) '# I.C. of manifold WU-:  dt   (x, y, vx, vy)   cj'

! -- compute nmf = 100 initial points on the manifold
epsl = -dabs(epsl)
call gr_mfdinst(y0,  ndim, tp, vepu, nmf, epsl, finit, ymfd, gr_plf, gr_cjplf) 
 close(finit)
print*, 'finish WU-'; read*

ispl = 1
print*, 'Plot the orbit of WU-'
if( ispl == 1 ) then
  open(fmfd, file = './dat/mfd_plf/h/obum8.dat', status='replace', access='append')
  write(fmfd,*) '# Orbits on the manifold  WU-:  t (x, y, vx, vy)'
  do k = 1, nmf, 1 
    yi = ymfd(k, 1:ndim)
    call plob(yi, t0, tf, ndim, ndim, tdir, ispl, fmfd, gr_plf, gr_cjplf, yf2) 
  enddo 
end if

print*, 'finish WU-'; read*
!stop 


! --- left half of  stable manifold ---- 
tdir = -1
! -- I.C.
open(finit, file = './dat/mfd_plf/h/mfd_initsm8.dat', status='replace', access='append')
write(finit,*) '# I.C. of manifold WS-: dt   (x, y, vx, vy)   cj'
! -- compute nmf = 100 initial points on the manifold

epsl = dabs(epsl)
call gr_mfdinst(y0,  ndim, tp, veps, nmf, epsl, finit, ymfd, gr_plf, gr_cjplf) 
 close(finit)
print*, 'Finish WS-'; print*; read*; !ck

ispl = 1
print*, 'Plot the orbit of WS-'
if( ispl == 1 ) then
  open(fmfd, file = './dat/mfd_plf/h/obsm8.dat', status='replace', access='append')
  write(fmfd,*) '# Orbits on the manifold  WS-:  t (x, y, vx, vy)'
  do k = 1, nmf, 1 
    yi = ymfd(k, 1:ndim)
    call plob(yi, t0,   tf, ndim, ndim, tdir, ispl, fmfd, gr_plf, gr_cjplf, yf2) 
  enddo 
end if

 close(fmfd_pc) 
 
end 
 



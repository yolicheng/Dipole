program tmapmfd_plf

! This is to compute the tmap of the invariant manfold, as a comparation with the Poincare map pcmfd_plf
! just to check how the manifolds evolves, an intuitive guess is that they keep as a closed curve. 

! ** NOTE ** 
! -- 1. We are using linear approaximation of the manifold, so we cannot trust too much the result after a long time interval 

! -- 2. Focus on the time around first crossing of x-axis, take shorter time interval to observe the evolution of Time-T map 

! -- 3. data saved in ./dat/map_plf/h  subfolder  with 'tmap' to classify

    !  for general use, and after the check, we can copy the files to the specified energy level  
    !  --3.1. mfd_init.dat        -- the I.C. of the manifold  
    !  --3.2. map_mfd.dat         -- the Poincare map of the manifold cut by Poincare section
 
 
use dp_mod
use pi_mod 
use plf_mod 
implicit none

integer, parameter  ::  ndim   = 4,  &    ! for lf
                        npvar  = 20, &    ! 20 for variational matrix, to compute Monodromy matrix 
                        nmf    = 2000 !2000     ! numbeer of orbits on the manifold 

! we have computed nmf = 2000 and npoinc = 20, saved in h_h3 subfolder   

! take npoinc = 1000, nmf=1 to get the full poince map, and                        
! Local Variables
integer       ::  debug, i 
real(kind=dp) ::  beta0, dlf(ndim, ndim)


! P.O. 
integer           ::  fob, fpo_init
real(kind=dp)     ::  y0(ndim), TP, cj0, cj 
character(len=70) ::  fnpo_init 

! mfd -- Eigenvalues and eigenvectors of MM       
integer         ::  col, finit, ftmap, fmfd, ispl, tdir    
real(kind=dp)   ::  pvi(npvar), pvf(npvar), phi(ndim, ndim), wr(ndim), wi(ndim), vr(ndim, ndim), & !mm
                    vepu(ndim), veps(ndim),  epsl,  t0, tf,  dt, &  !tmap 
                    eqmfd_init(ndim) ! eq_mfd

real(kind=dp) ::  dnrm2 

integer ::  isup, isum, issp, issm 

isup = 0; isum = 0 ! already well computed, do not over write 

! TODO: debug, for t = 130, there is some problem
issp = 1
issm = 1


! take eq1  in planar normal case as an example, type center X saddle  equilibria
! with beta = 1.d0, there is resonance when beta is around 1, avoid this for the moment 
debug = 0

epsl = 1.d-6

call init_plf     ! to assign the value of  eq 

! TODO: the value of beta matters! Be careful to use the same value when we do exploration 
!       At this moment, we take beta=2 

beta0 = 2.d0   ! use 2 to check, but 1 is better?  

print*, 'Check the value of beta :', beta0; read*; print*

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
 ! 4.3767487109222252  for h3 !--- the current studied closed one --- 
  
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

! Compute the stable and unstable manfold of the equilibria. check the energy level 
print*, 'Pick the dominant unstable eigenvector ( lambda > 0)! Input the column'
read(*, *)  col 
vepu =  vr(:, col)
print*, 'eigenvector :', vepu; print*; read*

print*, 'Pick the stable eigenvector (   lambda_max < 0)! Input the column'
read(*, *)  col 
veps = vr(:, col)
print*, 'eigenvector :', veps; print*; read*


print*, 'Compute the stable and unstable manfold of the equilibria'; read*; !ck

epsl = 1.d-6; tf = 1.d2

! unstable, right branch
open(211, file='./dat/map_plf/eqmfd_up.dat', status='replace', access='append')
eqmfd_init = eq + epsl * vepu/ dnrm2(ndim, vepu, 1)
call plob(eqmfd_init, 0.d0, tf, ndim, ndim, 1, 1, 211, gr_plf, gr_cjplf,  pvf) 
 close(211)

! unstable, left branch
open(211, file='./dat/map_plf/eqmfd_um.dat', status='replace', access='append')
eqmfd_init = eq - epsl * vepu / dnrm2(ndim, vepu, 1)
call plob(eqmfd_init, 0.d0, tf, ndim, ndim, 1, 1, 211, gr_plf, gr_cjplf,  pvf) 
 close(211)
 
open(211, file='./dat/map_plf/eqmfd_sp.dat', status='replace', access='append')
! unstable, right branch
eqmfd_init = eq + epsl * veps / dnrm2(ndim, veps, 1)
call plob(eqmfd_init, 0.d0, tf, ndim, ndim, -1, 1, 211, gr_plf, gr_cjplf,  pvf) 
 close(211)
 
! unstable, left branch
open(211, file='./dat/map_plf/eqmfd_sm.dat', status='replace', access='append')
eqmfd_init = eq - epsl * veps/ dnrm2(ndim, veps, 1)
call plob(eqmfd_init, 0.d0, tf, ndim, ndim, -1, 1, 211, gr_plf, gr_cjplf,  pvf) 
 close(211)

stop 


! TODO: remember to modify the name for specific initial condition
fpo_init = 200  
write(fnpo_init,  fmt='(a)')  './dat/map_plf/h3/poy_init.dat' ! 

! Read from the refined periodic orbits.
open(fpo_init, file = fnpo_init, status='old')
read(fpo_init, *)  tp, y0, cj0 
write(*, *)        'I.C. of p.o.: tp, y0, cj0'
print*, tp, y0, cj0;   read* ! ck

! -- the initial condition for the manifolds, use different names  for different manifolds

! --- compute the manifold, using the eigenvector of Monodromy Matrix as linear approximation
! --- plot the orbit for one revolution and compute the Eigenvalues of the Monodromy Matrix to check the stability
pvi(1:ndim) = y0  
pvi(5:20) = 0.d0
pvi(5:20:5) = 1.d0
! do not plot, only use this for Monodromy matrix computation
fob = 532
open(fob, file = './dat/map_plf/h/po.dat', status = 'replace')
call plob(pvi, 0.d0, TP, ndim, npvar, 1, 1, fob, gr_plf, gr_cjplf,  pvf) 
 close(fob) 
 
 
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
print*, 'Pick the dominant unstable eigenvector ( lambda_max > 1)! Input the column'
read(*, *)  col 
vepu =  vr(:, col)
print*, 'eigenvector :', vepu; print*; read*

print*, 'Pick the stable eigenvector ( 1 / lambda_max < 1)! Input the column'
read(*, *)  col 
veps = vr(:, col)
print*, 'eigenvector :', veps; print*; read*


! Compute the  initial conditions of manifold  and the first npoinc returns together with the orbit 
!subroutine gr_mfdinst(ypo, n, tp, mmat, nmf, epsl, ftag, ymfd, deriv, gr_cj)
!epsl = 1.d-5 or 1.5-6 ! doesn't matter too much... 

finit = 100; ftmap = 110; fmfd = 120  
epsl = 1.d-6

! ** NOTE **
! by checking the first return time  from mfd_mapup.dat and  mfd_mapum.dat
! we have tpc1 = 105 unit time  approximately; 
! take  t0 = 80, tf=130,   dt = 2.d0  to check

! From mfd_mapsp.dat and  mfd_mapsm.dat
! we have tpc1 = 155 unit time  approximately
! take  t0 = 120, tf=190 , dt = 2.d0  to check

! they are symmetric w.r.t. the x-axis

! ************************  Unstable Manifold  *************************
! check tf = 80:  120, make sure the middle point is more or less the time across the x-axis
dt = 2.d0
t0 = 80.d0
tf = 130.d0 

! ****************** Right half of unstable manifold *******************
tdir = 1
epsl = dabs(epsl)

ispl = 0

if(isup == 1) then 
open(finit, file = './dat/tmapmfd_plf/h/mfd_initup.dat', status='replace', access='append')
write(finit,*) '# I.C. of manifold  WU+:   (x, y, vx, vy)   cj'

open(ftmap, file='./dat/tmapmfd_plf/h/tmapup.dat', access='append', status='replace')
write(ftmap, *) '# Time-T map of WU+:  T       (x, y, vx, vy)    cj'

open(fmfd, file='./dat/tmapmfd_plf/h/obup.dat', access='append', status='replace')
write(fmfd,*) '# Orbits on the manifold WU+:  t       (x, y, vx, vy)    cj'

! subroutine tmap_mfd(y0, tp, ndim, nmf, vep, epsl, tdir, t0, tf, dt, finit, ftmap, ispl, fob)
call tmap_mfd(y0, tp, ndim, nmf, vepu, epsl, tdir, t0, tf, dt, finit, ftmap, ispl, fmfd, &
               gr_plf, gr_plf_n, gr_cjplf)


print*, 'Finish Time-T map of WU+'; print*; read*; !ck
endif 
stop 

! **********************************************************************



! ****************** Left half of unstable manifold *******************
tdir = 1 ! ustable
epsl = -dabs(epsl)

ispl = 0
if(isum == 1) then 
open(finit, file = './dat/map_plf/h/mfd_initum.dat', status='replace', access='append')
write(finit,*) '# I.C. of manifold  WU-:   (x, y, vx, vy)   cj'

open(ftmap, file='./dat/map_plf/h/mfd_tmapum.dat', access='append', status='replace')
write(ftmap,*) '# Time-T map of WU-:  T       (x, y, vx, vy)    cj'

open(fmfd, file='./dat/map_plf/h/mfd_tmap_obum.dat', access='append', status='replace')
write(fmfd,*) '# Orbits on the manifold WU-:  t       (x, y, vx, vy)    cj'

! -- I.C. of manifold +  Time-T map 
!subroutine tmap_mfd(y0, tp, ndim, nmf, vep, epsl, tdir, t0, tf, dt, finit, ftmap, ispl, fob, deriv, gr_cj)
call tmap_mfd(y0, tp, ndim, nmf, vepu, epsl, tdir, t0, tf, dt, finit, ftmap, ispl, fmfd, &
              gr_plf, gr_plf_n, gr_cjplf)

print*, 'Finish WU-'; print*; read*; !ck
endif
! **********************************************************************



! *************** Stable Manifold *********************

! From mfd_mapsp.dat and  mfd_mapsm.dat
! we have tpc1 = 155 unit time  approximately
! take  t0 = 120, tf=190 , dt = 2.d0  to check

! Check tf = 120: 190 to see how the loops appear
dt = 3.d0
t0 = 120.d0 
tf = 190.d0 

! ****************** Right branch of stable manifold *******************
tdir = -1;   ! stable manifold, integrate backward in time 
epsl = dabs(epsl)
ispl = 0

!! ** TODO ** for debug 
!ispl = 1
!t0 = 130.d0
!tf = 140.d0

!open(finit, file = './dat/map_plf/h/mfd_initsp1.dat', status='replace', access='append')
!write(finit,*) '# I.C. of manifold  WS+ :   (x, y, vx, vy)   cj'
! 
!open(ftmap, file='./dat/map_plf/h/mfd_tmapsp1.dat', access='append', status='replace')
!write(ftmap,*) '# Time-T map WS+ :  T       (x, y, vx, vy)    cj'

!open(fmfd, file='./dat/map_plf/h/mfd_tmap_obsp1.dat', access='append', status='replace')
!write(fmfd,*) '# Orbits on the manifold WS+ :  t       (x, y, vx, vy)    cj'

!!subroutine tmap_mfd(y0, tp, ndim, nmf, vep, epsl, tdir, t0, tf, dt, finit, ftmap, ispl, fob)
!call tmap_mfd(y0, tp, ndim, nmf, veps, epsl, tdir, t0, tf, dt, finit, ftmap, ispl, fmfd, &
!              gr_plf, gr_plf_n, gr_cjplf)

!print*, 'Finish WS+'; print*; read*; !ck

!! ***************************


if(issp == 1) then 
open(finit, file = './dat/map_plf/h/mfd_initsp.dat', status='replace', access='append')
write(finit,*) '# I.C. of manifold  WS+ :   (x, y, vx, vy)   cj'
 
open(ftmap, file='./dat/map_plf/h/mfd_tmapsp.dat', access='append', status='replace')
write(ftmap,*) '# Time-T map WS+ :  T       (x, y, vx, vy)    cj'

open(fmfd, file='./dat/map_plf/h/mfd_tmap_obsp.dat', access='append', status='replace')
write(fmfd,*) '# Orbits on the manifold WS+ :  t       (x, y, vx, vy)    cj'

!subroutine tmap_mfd(y0, tp, ndim, nmf, vep, epsl, tdir, t0, tf, dt, finit, ftmap, ispl, fob)
call tmap_mfd(y0, tp, ndim, nmf, veps, epsl, tdir, t0, tf, dt, finit, ftmap, ispl, fmfd, &
              gr_plf, gr_plf_n, gr_cjplf)

print*, 'Finish WS+'; print*; read*; !ck
endif
! **********************************************************************


! ****************** Left branch of stable manifold *******************
tdir = -1
epsl = -dabs(epsl)
ispl = 0 

if(issm == 1) then 
open(finit, file = './dat/map_plf/h/mfd_initsm.dat', status='replace', access='append')
write(finit,*) '# I.C. of manifold  WS-:   (x, y, vx, vy)   cj'

open(ftmap, file='./dat/map_plf/h/mfd_tmapsm.dat', access='append', status='replace')
write(ftmap,*) '# Time-T map WS-:  T       (x, y, vx, vy)    cj'

open(fmfd, file='./dat/map_plf/h/mfd_tmap_obsm.dat', access='append', status='replace')
write(fmfd,*) '# Orbits on the manifold WS-:  t       (x, y, vx, vy)    cj'

!subroutine tmap_mfd(y0, tp, ndim, nmf, vep, epsl, tdir, t0, tf, dt, finit, ftmap, ispl, fob)
call tmap_mfd(y0, tp, ndim, nmf, veps, epsl, tdir, t0, tf, dt, finit, ftmap, ispl, fmfd, &
              gr_plf,  gr_plf_n, gr_cjplf)

print*, 'Finish WS-'; print*; read*; !ck
endif
! **********************************************************************

end




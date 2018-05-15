program poinc_plf

! This is to compute the Poincare map representation of the phase space,  
! we do for both right half x-axis and y-axis 
! 
! whereas map_plf is for the first two returns for symmetric p.o.s, do not mix the two routines

! we fix the energy level and the Poincare section y=0, so we only have two degree of freedom 
! to make it simple, we take vx0 = 0 
! and treat vy as a function of (h, x, y=0, vx=0)

!  we explore the interval  x \in [0, 0.72],  and take only the crossing with positive velocity vy 

! that means for both x and y axis
!  --- data saved in ./dat/poinc_plf2/  subfolder 

! **NOTE** 
!  be careful of the retrograde and prograde definition here,  instead of the one used for search of symmetric periodic orbit
!  we consider as cirteria the direction of vx when they pass the Poincare section.

!  1. map_pro_x.dat         -- the  poincare map with positive vy
!  2. map_retro_x.dat       -- the  poincare map with negative vy

use dp_mod
use pi_mod 
use poinc_mod
use plf_mod 
implicit none

integer, parameter  ::  ndim   = 4,  & ! for lf
                        npvar  = 6  !42 for variational matrix, 6 for state vector 
integer ::              npoinc = 1000 
! Local Variables
integer       ::  debug, i, j, k
real(kind=dp) ::  beta0 

! eigenvalues and eigenvectors around the equilibria                 
real(kind=dp) ::  dlf(ndim, ndim), wr(ndim), wi(ndim), vr(ndim, ndim) 


! Poincare map 
integer       :: ind0, dir0, imax0, ispro, ind_free, fmap_pro, fmap_retro, fmap 
!                ind will be the Poincare section, ind_free will be the free parameter left
!                first we take ind=2, ind_free = 1 to detect p.o. symmetric w.r.t. x-axis
!                then we take ind=1, ind_free = 2 to detect p.o. symmetric w.r.t. y-axis
! 
real(kind=dp) :: p00,  tmax0, xmax0, cj0, vstep, v0 
character(len=70) :: fnpro, fnretro, fnmap

! Po 
integer ::  ispo 
real(kind=dp) ::  y0(ndim-2),   TP, yf(ndim), dg(1,1), cj  

! curve
integer :: ftorus,foball, ind_curv, isv, ispc, ispl, np
real(kind=dp) :: pvi(ndim), pvf(ndim), cj_i, tf, hminim, xstep, x0!, tf_total
character(len=70) :: fntorus, fncurve 

! for fft_ob 
integer           ::  istol, nfmax, nsm,  nitmax, isfour 
real(kind=dp)     ::  tolres,  tolfc 

! take eq1  in planar normal case as an example, type center X saddle  equilibria
! with beta = 1.d0, there is resonance when beta is around 1, avoid this for the moment 
debug = 0

call init_plf     ! to assign the value of  eq 

! TODO: the value of beta matters! at this moment, we take beta=2 
beta0 = 2.d0   ! use 2 to check, but 1 is better?  
!print*, 'Input the value of beta (n/w_c):'
!read*, beta0 
call init_beta(beta0)   
     
print*, 'check,  beta',  beta0 

call gr_cjplf(eq, cj) 
print*, 'check energy! ', cj, eq 
! 4.3267487109222253  for equilibirum in plf


 ! ** NOTE ** 
 ! ----- for p.o. symmetric to x-axis, we have ....
 !  the energy level is within the interval [4.320714 : 3.827071361], we have  two families of periodic orbits 
 
 !  for energy above 4.320714,  one family exists and goes to the origin, the other family diminish, we have choatic region
 !  this case is what we are studing now. H3 = 4.3767487109222252 
 
 ! for energy level less than  3.827071361, there will be no periodic orbit.

 ! and we also have one family that is symmetric w.r.t. y-axis, which goes to infinity....
 
 
! -- the energy level we r taking, smaller than the one at the equilibria, so we have a 
! closed UFO 
! for x ----
 cj0 = cj + 0.05d0
 print*, 'cj0 =', cj0; read*  
 ! cj0 =   4.3767487109222252  for h3 ! a little bit larger than the one of equilibirum point (cj)
 
 ! try 4.15, when we have no chaotic region
 cj0 = 4.15d0
 
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

! compute the zero velocity curve 

!call zvc_plf(cj0, 2.d0)

! -- the real units-----
print*, 'runit(km), tunit(s), vuint(km/s): ', runit, tunit, vunit

! Jacobi matrix of the lorentz force with respective to the state 
call dflrtz2(eq, dlf)

print*, 'Check DG in the planar case: ' 
do i = 1, ndim 
  write(*,'(4f18.8)') dlf(i,:) 
enddo
print*; read* 

call eigrg(dlf, ndim,1, wr,wi, vr)
print*, 'Check if we have a center X center X saddle equilibria:' 
print*; read* 

! Take initial point on the x-axis, with velocity perpendicular to x-axis, and compute the P.O. 
! by Poincare map with the section as y=0, take vy as a function of f(cj, x, y=0, vx)
! so we have 2 dimension of freedom (x, vx )


!! ---------- Computation -----------------------------
! For poincare map limitation 
xmax0 =  8.d0  !
tmax0 =  4.d1 ! 

! Poincare section:  y=0  
ind_free  = 1        ! 1: move along x - axis
!ind_free  = 2       ! 2: move along  y - axis
ind0  = 3-ind_free  ! 2: Poincare section: y (or 1: x) = 0 

p00   = 0.d0        ! rough guess, by observing directly at the plot of the tori  
dir0  = 0           ! no requirement in the direction, just compute the first turn 
imax0 = 1           ! consider every intersection, and choose by the direction of velocity

!subroutine init_poinc( ind0, p00, dir0, imax0, tmax0, xmax0, ndim0, h00) 
call init_poinc(ind0, p00, dir0, imax0, tmax0, xmax0, ndim, cj0)
  
! instead of refining for po orbit, we just take a lot of points on x-axis, and do the poincare section 
! observe the first and second return to the Poincare section to see if there exist p.o. 
! ** NOTE *** 
! Definition of prograde and retrograde orbit 

! Prograde: The normal to orbit is the same with the normal to the chief's orbit around the Earth 
! retrograde: .... the contrary 

! ** Criteria **
! Since we have the symmetry about the origin in '+t' time sense, once we have (x, y, vx, vy)
! as a solution, we will have also (-x, -y, -vx, -vy), so we only need to study the first qudrant (x>0, y>0)
 
! we take initial point on x-axis with an initial velocity perpendicular to x-axis, (x0, 0, 0, vy0), with x0>0
! and check the first return to x-axis(y=0 section), take vy0>0. 

! We could have two kinds of p.o.,   
! -- 1. one is prograde,  so xf < x0 after half revolution. 

! -- 2. one is retrograde, so xf > x0 after half revolution 

! Now we have to make clear if it is necessary to study also the case when vy0 < 0  
! for the above two kinds of p.o.s
! -- 1. if xf > 0 for all the grids, there is no need to take vy0 < 0, since for every point (x0>0, vy0<0), its first return to y=0 section 
!       corresponds to a point (xf > 0, vyf > 0) 
! -- 2. if xf < 0, we obtain a p.o. symmetric w.r.t. to both x-axis and y-axis, which ends up with a circle. 

fmap_pro = 25;    fmap_retro   = 52;  fmap = 50

if(ind_free == 1) then  ! for x axis
!  if(dir == 0)  write(fnmap, fmt='(a)')   './dat/poinc_plf2/map_x.dat'
!  if(dir == 1)  write (fnpro,  fmt='(a)')  './dat/poinc_plf2/map_pro_x.dat'
!  if(dir == -1)  write(fnretro,fmt='(a)')  './dat/poinc_plf2/map_retro_x.dat'
!  write(fnmap,  fmt='(a)')   './dat/poinc_plf2/map_x.dat'

  write (fnpro, fmt='(a)')  './dat/poinc_plf2/map_pro_x.dat'
  write(fnretro,fmt='(a)')  './dat/poinc_plf2/map_retro_x.dat'
!  
elseif(ind_free == 2) then ! for y axis
!  if(dir == 0)  write(fnmap, fmt='(a)')   './dat/poinc_plf2/map_y.dat'
!  if(dir == 1)  write(fnpro,  fmt='(a)')  './dat/poinc_plf2/map_pro_y.dat'
!  if(dir == -1) write(fnretro,fmt='(a)')  './dat/poinc_plf2/map_retro_y.dat'
  
!  write(fnmap, fmt='(a)')   './dat/poinc_plf2/map_y.dat'
  write(fnpro,  fmt='(a)')  './dat/poinc_plf2/map_pro_y.dat'
  write(fnretro,fmt='(a)')  './dat/poinc_plf2/map_retro_y.dat'
endif  

!open(fmap, file = fnmap, access ='append',status = 'replace')
!write(fmap, *)  ' # Poincare map +-:  t      (x y  vx vy )   h'

open(fmap_pro, file = fnpro, access ='append',status = 'replace')
write(fmap_pro, *)  ' # Pro Poincare map v+:  t      (x y  vx vy )   h'

open(fmap_retro, file = fnretro  , access ='append',status = 'replace' )
write(fmap_retro, *)  ' # Retro Poincare map v-:  t      (x y  vx vy )   h'

! plot the orbit, just to check what happened with the poincare map close to the origin
!open(888, file = './dat/poinc_plf2/pcob.dat', access ='append',status = 'replace')
!write(888, *)  ' # Orbits to check:  t      (x y  vx vy )   h'

ispl = 0
np   = 40 ! number of initial points for the poincare map of the phase space ...

! -- only consider the positvie half x-axis, AND vary x and vx in grid 
! for x = [-1, 1], because we also want to see what happens outside the zvc 
! for vx, instead of taking 0, which gives us perpendicular orbit, we look inside [-8,8]

! update finally -- by applying the symmetry, we can only explore the positive half axis, with 
x0 = 1.d0
v0 = 4.d0

! -- only the first quadrant, and take advantage of the symmetry
xstep = x0 / np
vstep = v0 / np 

! this is for general plots, and we do local with smaller stepsize 
xstep = 1.d-2  ! for the region of [0.32:0.38], we do more plots 
x0 = 0.d0 
do  i =  0, 100000, 1
  x0  = x0 + xstep   ! x
!  if( (x0 > 0.30d0 .and. x0 < 0.40d0)  .or. (x0 > 0.50d0 .and. x0<0.70d0) ) then 

! for x axis, 0.2410 -- 0.4656 is the region for tori around the elliptic p.o.
! and 0.32 is more or less a fixed point in the Poincare map 

! so, we can skip 0.32-0.46
!  if(ind_free == 1 .and. (x0>0.34d0 .and. x0 < 0.4d0) ) cycle
  
  !  if(ind_free == 1 .and. (x0>0.350 .and. x0 < 0.45d0) ) cycle ! for energy level 0.41
  
  if(  x0 > 0.24d0 .and. x0 < 0.70d0 ) then 
    xstep =  5.d-3 
    npoinc = 1000
  else 
    xstep  = 2.d-2
    npoinc = 400
  endif  
   
  

 ! for the region of [0.32:0.38], we do more plots 
!  x0  = -0.38d0 + i * xstep   ! x
!  x0  = 0.32d0 + i * xstep 
!  if(dabs(x0) < 0.32d0 .or. dabs(x0) >0.38d0) cycle
 
! this is for a particular zone...  close to the chaotic edges
!  x0  = 0.42d0 + i * xstep 
!  xstep = 4.d-2
!  npoinc = 20000
!  if(dabs(x0) < 0.4d0 .or. dabs(x0) >0.65d0) cycle
  
! --- start from -0.72  --- discard--- we only need to consider the x>0, vx>0 (and y>0, vy>0) due to the symmetry
! xstep = x0 * 2 / np /4
!!vstep = v0 * 2 / np/10
! vstep = v0   / np/10
!do  i = 0, 2*np, 1 ! vary x
! !!   the whole axis
!  x0  =  -1.d0 +  i * xstep   ! x \in  [-1, 1]

  print*, 'x0 = ', x0 
  if(dabs(x0) > 1.d0) exit
 
  ! skip the region that is close to the origin
  if( dabs(x0) <1.d-1 ) cycle
  

  !--  vary vx 
  !  np = 0
  do j = 0,  np , 1  
  
    pvi = 0.d0           ! y  = 0  ! vx = 0 
    
!    v0  = -2.d0 + j * vstep   ! vx \in [-8, 8] -- discard 

    v0  = j * vstep   ! only positive velocity vx0

!    v0 = 0.d0 ! this is for velocity perpendicular to the starting axis
    
    pvi(ind_free)    = x0 
    pvi(ind_free +2) = v0 
    
    ! the other  velocity: vy
    call cj2v_plf(pvi, cj0, ind+2, isv)
   
    ! TODO: this is with negative vy0 
!    pvi(ind+2) = -pvi(ind+2)
  
    if(isv == 0 ) then 
      print*, pvi 
      cycle
    endif  
  
    ! check the energy
    call gr_cjplf(pvi, cj_i)
    if (debug == 1) print*, 'State of initial point: ', pvi, 'cj = ', cj_i, 'dcj', dabs(cj_i-cj0) ; ! read*
  
    ! before doing Poincare map, we check the orbits
!    open(888, file = './dat/poinc_plf2/pcob.dat', access ='append',status = 'replace')
!    call plob(pvi, 0.d0, 2.d1, 4, 1, 888, 1, gr_plf, gr_cjplf,  pvf) 
!    close(888)
!    print*, 'Finish one orbit!'; 
!    read*; read*
    
    do k = 1, npoinc, 1
    !   subroutine poinc(sti, ndim, nvar,  tf,stf, hminim, ispc, deriv, gr_cj) 
      call poinc(pvi, ndim, ndim,   tf, pvf, hminim, ispc, gr_plf, gr_cjplf)

      if(ispc == 0) exit 
      
      ! check the energy 
      call gr_cjplf(pvf, cj_i)
      
      if(dir == 0) then ! keep intersetions in both directions 
!        write(fmap, *)    tf, pvf, cj_i ! we'd better seperate them... 
        if(pvf(ind+2) > 0.d0 ) then 
          write(fmap_pro, *)    tf, pvf, cj_i
        else 
          write(fmap_retro, *)  tf, pvf, cj_i
        endif
      endif       
      
      if(dir ==  1)  write(fmap_pro, *)    tf, pvf, cj_i
      if(dir == -1)  write(fmap_retro, *)  tf, pvf, cj_i
    
      if(debug == 1) then 
         print*, 'i,k, dir', i, k, dir;  read* 
      endif    
    
    
      if (debug == 1) then 
        print*, k, '-th return to Poincare section!'
        write(*, *) k, tf, pvf, cj_i;  print*;  read*
      endif 
    
      pvi = pvf  ! as the new start point 
    enddo 
  
    if(dir == 0) then
!      write(fmap, *);
      write(fmap_pro, *);  write(fmap_retro, *)  
    endif 
     
    if(dir ==  1)   write(fmap_pro, *)  
    if(dir == -1)   write(fmap_retro, *)  
 
  enddo 
enddo   

! close(fmap_pro); 
 close(fmap_retro); close(fmap)
 stop 
 
 
end 
 



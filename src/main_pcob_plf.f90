program pcob_plf

! This is to compute the interesting orbits from the Poincare map representation of the phase space by poinc_plf 
! For simplicity, we take vx0 = 0 for the initial points. 

! By choosing carefully the initial points, we can generate the Poincare map of the phase space 

! we fix the energy level and the Poincare section y=0, so we only have two degree of freedom 
!  to make it simple, we take vx0 = 0 
!  and treat vy as a function of (h, x, y=0, vx=0)

!  we explore the interval  x \in [0, 0.72],  and take only the crossing with positive velocity vy 

!  --- data saved in ./dat/pcob_plf/  subfolder 

! **NOTE** 
!  be careful of the retrograde and prograde definition here,  instead of the one used for search of symmetric periodic orbit
!  we consider as cirteria the direction of v(ind_free) when the orbit passes the Poincare section.

!  1. pc_pro_x.dat         -- the  poincare map with positive vy
!  2. pc_retro_x.dat       -- the  poincare map with negative vy

! instead of save the return time for each map, we save the time elapsed from the initial point, good to pruduce 

use dp_mod
use pi_mod 
use poinc_mod
use plf_mod 
implicit none

integer, parameter  ::  ndim   = 4,  &  ! for lf
                        npvar  = 4      ! 20 for variational matrix, 4 for state vector 

 
! the problem is the value of npoinc for one revolution of different tori might be different, 
! so we have to deal with the orbits one by one ...      
! and since we are clear about the symmetry, at this moment, we only need to keep the positive crossing  dir = 1
! and if we want to obtain the Poincare map representation, we take dir = -1 then 

! take the initial point to be npc0+1 iteration, but I think it is not necessary, discard at the moment  
! since any curve can be obtained with vx0=0 . Original, npc0=2000    
integer, parameter  ::  npc0 = 0! 2000
              
integer, parameter  ::  npoinc = 1000 !*2 ! 10000   ! for the full poincare map 

! the rest are only for one revolution 
!integer, parameter  ::  npoinc = 20 !22  !19! 26! 31! 35! 40! 48! 60 !82  

!integer, parameter  ::  npoinc = 131! 259 ! 131  
!integer, parameter  ::  npoinc = 259  

! see file pc_init.dat for the number of iterats
! x0 = 0.30,    npoinc = 60
! x0 = 0.305,   npoinc = 82
! x0 = 0.31,   npoinc = 131 
! x0 = 0.315   npoinc = 259

!0.20, 22
! if we want to produce the Poincare portrait, it would be nice to use more points here 
!integer, parameter  ::  npoinc = 1000                                                                                          

! Local Variables
integer       ::  debug, i, k ! j,
real(kind=dp) ::  beta0 

! eigenvalues and eigenvectors around the equilibria                 
real(kind=dp) ::  dlf(ndim, ndim), wr(ndim), wi(ndim), vr(ndim, ndim) 


! Poincare map 
integer       :: ind0, dir0, imax0, ind_free, fpc_pro, fpc_retro, fpc  ! ispro, 
!                ind will be the Poincare section, ind_free will be the free parameter left
!                first we take ind=2, ind_free = 1 to detect p.o. symmetric w.r.t. x-axis
!                then we take ind=1, ind_free = 2 to detect p.o. symmetric w.r.t. y-axis
! 
real(kind=dp) :: p00,  tmax0, xmax0, cj0,   xi(1000), pv0(4), tftot !, tfhf
character(len=70) :: fnpro, fnretro, fnpc


! Poinc
integer ::  foball, isv, ispc, ispl, nob !np, 
real(kind=dp) :: pvi(ndim), pvf(ndim), cj_i, tf, hminim, x0, dx, x0i, x0f,  cj, yf(ndim)

! for fft_ob 
!integer           ::  istol, nfmax, nsm,  nitmax, isfour 

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
print*, 'check energy! ', cj, eq  ! eq: 0.6933612744, 0

! -- the energy level we r taking, smaller than the one at the equilibria, so we have a 
! closed UFO 
! for x ----
 cj0 = cj + 0.05d0  ! -- h3 

! cj0 = cj  ! -- h0 
 print*, 'cj0, cj_eq', cj0, cj; read*
 
 ! 4.3267487109222253  for equilibirum in plf
 ! 4.3767487109222252  for h3 !--- the current studied one 
 
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



! Jacobi matrix of the lorentz force with respective to the state 
call dflrtz2(eq, dlf)

print*, 'Check DG in the planar case: ' 
do i = 1, ndim 
  write(*,'(4f18.8)') dlf(i,:) 
enddo
print*; read* 

call eigrg(dlf, ndim,1, wr,wi, vr)
print*, 'Check if we have a center X saddle equilibria:' 
print*; read* 

! Take initial point on the x-axis, with velocity perpendicular to x-axis, and compute the P.O. 
! by Poincare map with the section as y=0, take vy as a function of f(cj, x, y=0, vx)
! so we have 2 dimension of freedom ( x, vx )


!! ---------- Computation -----------------------------
! For poincare map limitation 
xmax0 =  2.d0  !
tmax0 =  4.d1 ! 

! Poincare section:  y=0  
ind_free  = 1           ! 1: x (or 2: y) component
ind0      = 3-ind_free  ! 2: y (or 1: x) component, x(ind0) = 0 as poincare section 

p00   = 0.d0        ! rough guess, by observing directly at the plot of the tori  
dir0  = 1           ! no requirement in the direction, just compute the first turn 
imax0 = 1           ! consider every intersection, and choose by the direction of velocity

!subroutine init_poinc( ind0, p00, dir0, imax0, tmax0, xmax0, ndim0, h00) 
call init_poinc(ind0, p00, dir0, imax0, tmax0, xmax0, ndim, cj0)

fpc_pro = 25;    fpc_retro   = 52;  fpc = 50

if(ind_free == 1) then  ! for x axis crossing
  write(fnpc,   fmt='(a)')   './dat/pcob_plf/h/pc_x.dat'
  write(fnpro,  fmt='(a)')  './dat/pcob_plf/h/pc_pro_x.dat'
  write(fnretro,fmt='(a)')  './dat/pcob_plf/h/pc_retro_x.dat'
  
elseif(ind_free == 2) then ! for y axis crossing
   write(fnpc, fmt='(a)')   './dat/pcob_plf/h/pc_y.dat'
   write(fnpro,  fmt='(a)')  './dat/pcob_plf/h/pc_pro_y.dat'
   write(fnretro,fmt='(a)')  './dat/pcob_plf/h/pc_retro_y.dat'
  
endif  

open(fpc, file = fnpc, access ='append',status = 'replace')
write(fpc, *)  ' # Poincare map +-:  t      (x y  vx vy )   h'

open(fpc_pro, file = fnpro, access ='append',status = 'replace')
write(fpc_pro, *)  ' # Pro Poincare map :  t      (x y  vx vy )   h'

open(fpc_retro, file = fnretro  , access ='append',status = 'replace' )
write(fpc_retro, *)  ' # Retro Poincare map:  t      (x y  vx vy )   h'

! For several initial points, plot the orbits 
foball = 700
open(foball, file = './dat/pcob_plf/h/pcob.dat', status='replace')
write(foball, *) '# orbit: :  t      (x y  vx vy )   h'

ispl = 1

!xi(1:11) = (/0.31975d0,  0.372d0, 0.371d0, 0.37d0, 0.365d0, 0.353018d0,  &
!             0.334841d0, 0.330610d0, 0.326723d0, 0.319541d0, 0.31d0/)

!xi(1:9) = (/(0.31d0+i*1.d-2, i=0,8)/)  

! to make it more clear, we need to add three more orbits. 0.315, 0.33, 0.35 

! pc_pro_x_1m2.dat:   x0 \in [0.62: 0.01 : 0.1] + (/0.3125d0, 0.315d0, 0.595d0, 0.613d0 /) 
!x0i = 0.62d0
!x0f = 1.d-1
!xi(1:4) = (/0.3125d0, 0.315d0, 0.595d0, 0.613d0 /)
!dx = -1.d-2   ! step size


! Finally taken value: 
! h3/pc_pro_x_2m2.dat , 2.d-2 is enough
!xi(1:2) = (/0.315d0, 0.613d0/)
!x0i = 0.6d0   
!x0f = 1.d-1
!dx = -2.d-2   ! we take more orbits
!nob = int(  (x0f - x0i) / dx  )
!do i = 3, nob+2, 1
!  xi(i) = x0i + (i-3) * dx
!end do

! For h0 and h1 --- we can try other values 
xi(1:2) = (/0.319d0, 0.455d0/)  ! for y, 0.465-0.315
x0i = 0.46d0   
x0f = 1.d-1
dx = -5.d-3   ! we take more orbits

nob = int(  (x0f - x0i) / dx  )

do i = 3, nob+2, 1
  xi(i) = x0i + (i-3) * dx
end do

nob = nob +2 

print*, 'nob=', nob; print*; read*

! we need to add three more points to make the picture more clear 
! 1- 0.618  the outermost P.O.    
! 2- 0.32   the elliptic p.o. at x0 = 0.32 
! 3- 0.315, just to make clear the elliptic p.o. 

! the seperatrice are at x0 = 0.241084 and x0 = 0.3994123 
! so within [0.32:0.40] we can skip 

! and for x < 0.24.... we need to take more iterations to complete the curve 


! do the orbit one-by-one to check the trend
! so better to start from the middle point x0=0.32, first explore the left-hand side space 

! for every initial points, we have to check the poincare map and see when they finish one evolution... 
! xi(1) = 0.32d0  !; npoinc = 259
! xi(1) = 0.315d0 !; npoinc = 259
! xi(1) = 0.31d0  !;  npoinc = 131
! xi(1) = 0.305d0 !;  npoinc = 82
! xi(1) = 0.30d0  !;  npoinc = 60
! xi(1) = 0.295d0 !;  npoinc = 48
! xi(1) = 0.29d0  !;  npoinc = 40
! xi(1) = 0.285d0 !;  npoinc = 35
! xi(1) = 0.28d0  !;  npoinc = 31 
! xi(1) = 0.27d0  !;  npoinc = 26 
! xi(1) = 0.25d0  !;  npoinc = 19  
! ! xi(1) = 0.2475d0 !;  npoinc = 19  ! --almost the same as 0.25, so skip this one 

!! the transit must be within 0.245 and 0.2475
!!  xi(1) = 0.24625d0  ! 
! xi(1) = 0.24535d0  ! 
! xi(1) = 0.245d0 !;  npoinc = 19 
!   
!   ---- clearly goes to a difference torus
! xi(1) = 0.24d0 !;  npoinc = 35  

! xi(1) = 0.2d0 !;  npoinc = 22 

!xi(1) = 0.325d0 ! discard, since the next crossing is more or less symmetric w.r.t. 0.32, so we only need to explore even more right 
! we start from the left in 0.31, so we could continue with 0.335 

!xi(1) = 0.2d0 !335d0 

! save the x0, H0 and number of iterations for Poincare maps 
open(444, file='./dat/pcob_plf/h/pc_init.dat',status='replace',access='append')
write(444,*) ' # H0 = ', cj0 
write(444,*) ' # npc0(start)     npoinc(N.O.)    (x0, y0, vx0, vy0) '


do  i =  1, nob, 1
  pvi = 0.d0          ! y  = 0,  vx = 0 
  x0  = xi(i)

  
  ! for another torus at x0 = 0.615 
  x0 = 0.615d0    ! torus2 -- around center X center p.o.
!  x0 = 0.2  ! ob11
!  x0 = 0.48 ! torus4 -- close to p.o.s of type III 
  if(i>1) stop 
  
  print*, 'i, x0: ', i, x0; print*; read*
!  if(x0 > 0.325 .and. x0 < 0.395) cycle 

  pvi(ind_free) = x0    !x0
!  pvi(3) = 1.d-1 ! vx    -- if we want to try another velocity other than the one perpendicular to x-axis

  ! the other  velocity: vy
  call cj2v_plf(pvi, cj0, ind+2, isv)
   
  ! TODO: this is with negative vy0 
!  pvi(ind+2) = -pvi(ind+2)
  
  if(isv == 0 ) then 
    print*, pvi 
    cycle
  endif  
  
  ! save the initial data in pc_init.dat 
!  write(444, *)
  write(444, *) npc0, npoinc, pvi 
  
  pv0 = pvi
  
  ! check the energy
  call gr_cjplf(pvi, cj_i)
  if (debug == 1) print*, 'State of initial point: ', pvi, 'cj = ', cj_i, 'dcj', dabs(cj_i-cj0) ; ! read*

  tftot = 0.d0 
  
 ! put a blank line before each block, also easier to check the data for index or ev:::i::i when we r not sure how we save the data
  if(dir == 0) then
    write(fpc_pro, *);  write(fpc_retro, *)  
  endif  
  
  if(dir == 1)   write(fpc_pro, *)  
  if(dir == -1)   write(fpc_retro, *)  
    
  
  ! save also the first point, since it is on the x-axis, a nice start for every curve 
  do k = 1, npoinc+npc0, 1
  
   !   subroutine poinc(sti, ndim, nvar, tdir,  tf,stf, hminim, ispc, deriv, gr_cj) 
    call poinc(pvi, ndim, ndim, 1, tf, pvf, hminim, ispc, gr_plf, gr_cjplf)

    if(ispc == 0) exit 
    
    ! we do more iterations after we complete one revolution... 
    if(k <= npc0)  cycle  ! skip the first npc0 iterations, in fact it doesn't affect the result
    
    if(k == npc0+1)  then 
      pv0 = pvi
    ! save also the first point, since it is on the x-axis, a nice start for every curve 
      write(fpc_pro, *)  0.d0, pvi, cj_i
      
      !print, check with npc0=0, we shoulde have pvi of form (x0, 0,0,vy+) !--ckd
!      write(*, *) 'pvi=', pvi, 'tftot=', tftot, 'cj=', cj_i, 'x0=', x0
!      print*; read*
    endif  
    
    tftot = tftot + tf
    
   ! check the energy 
    call gr_cjplf(pvf, cj_i)
      
    if(dir == 0) then 
!      write(fpc, *)    tf, pvf, cj_i
      if(pvf(ind+2) > 0.d0 ) then 
        write(fpc_pro, *)    tftot, pvf, cj_i
      else 
        write(fpc_retro, *)  tftot, pvf, cj_i
      endif   
    endif  
    
    if(dir ==  1)    write(fpc_pro, *)    tftot, pvf, cj_i
    if(dir == -1)    write(fpc_retro, *)  tftot, pvf, cj_i
   
    if(debug == 1) then 
       print*, 'i,k, dir', i, k, dir;  read* 
    endif    
    
    if (debug == 1) then 
      print*, k, '-th return to Poincare section!'
      write(*, *) k, tftot, pvf, cj_i;  print*;  read*
    endif 
   
    pvi = pvf  ! as the new start point 
  enddo 
  
  if(dir == 0) then
    write(fpc_pro, *);  write(fpc_retro, *)
  endif  
  if( dir == 1)     write(fpc_pro, *)  
  if(dir == -1)     write(fpc_retro, *)  
 
  print*, 'Integrated time for orbit: ', tftot; read*
  print*, 'pv0 : ', pv0
!  tftot = 1.d3
  
!subroutine plob(y0, t0, tf, ndim, nvar, tdir, ispl, ftag, deriv, gr_cj,  y) 
  if(ispl == 1)  call plob(pv0, 0.d0, tftot, 4, 4, 1, ispl, foball, gr_plf, gr_cjplf,  yf) 
  print*, 'Finish the orbit plot!'  
  
  
  close(fpc_retro); close(fpc)
 

enddo 


end 
 



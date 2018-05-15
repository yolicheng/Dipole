program fft_plf

! This is to do the Fourier analysis for several orbits to observe how the frequecies envolves
! we take x0 \in [0.32: -0.005 : 0.24] 

! we fix the energy level and the Poincare section y=0, so we only have two degree of freedom 
!  to make it simple, we take vx0 = 0 
!  and treat vy as a function of (h, x, y=0, vx=0)

!  we explore the interval  x \in [0, 0.72],  and take only the crossing with positive velocity vy 

!  --- data saved in ./dat/fft_plf/  subfolder 

! **NOTE** 
!  be careful of the retrograde and prograde definition here,  instead of the one used for search of symmetric periodic orbit
!  we consider as cirteria the direction of vx when they pass the Poincare section.

!  1. map_pro_x.dat         -- the  poincare map with positive vy
!  2. map_retro_x.dat       -- the  poincare map with negative vy

! instead of save the return time for each map, we save the time elapsed from the initial point, good to pruduce 


use dp_mod
use pi_mod 
use poinc_mod
use plf_mod 
!use po_plf_mod
implicit none

integer, parameter  ::  ndim   = 4,  & ! for lf
                        npvar  = 6   !42 for variational matrix, 6 for state vector 

 
! the problem is the value of npoinc for one revolution of different tori might be different, 
! so we have to deal with the orbits one by one ...      
! and since we are clear about the symmetry, at this moment, we only need to keep the positive crossing  dir = 1
! and if we want to obtain the Poincare map representation, we take dir = -1 then 

! take the initial point to be npc0+1 iteration, but I think it is not necessary, discard at the moment  
! since any curve can be obtained with vx0=0 . Original, npc0=2000    
integer, parameter  ::  npc0 = 0! 2000
              
!integer, parameter  ::  npoinc = 1000 !*2 ! 10000   ! for the initial try   

! the rest are only for one revolution 
integer, parameter  ::  npoinc = 100 !19! 26! 31! 35! 40! 48! 60 !82  
!it seems between 0.30 and 0.295, the topology of the orbits changes

!integer, parameter  ::  npoinc = 131! 259 ! 131  

!integer, parameter  ::  npoinc = 259  

! see file pc_init.dat for the number of iterats
! x0 = 0.30,   npoinc = 60
! x0 = 0.305,   npoinc = 82
! x0 = 0.31,   npoinc = 131 
! x0 = 0.315   npoinc = 259

! if we want to produce the Poincare portrait, it would be nice to use more points here 
!integer, parameter  ::  npoinc = 1000                                                                                          

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
real(kind=dp) :: p00,  tmax0, xmax0, cj0, vstep, v0, xi(100), pv0(4), tftot, tfhf
character(len=70) :: fnpro, fnretro, fnmap

! Po 
integer ::  ispo 
real(kind=dp) ::  y0(ndim-2),   TP, yf(ndim), dg(1,1), cj  

! Poinc
integer :: ftorus,foball, ind_curv, isv, ispc, ispl, np
real(kind=dp) :: pvi(ndim), pvf(ndim), cj_i, tpc, hminim, xstep, x0 !, tf_total
character(len=70) :: fntorus, fncurve 

! for fft_ob 
integer       ::  istol, nfmax, nsm,  nitmax,  nob, nf_bas, fx0_fabas, isrpt 
real(kind=dp) ::  dx,  tf,  fa_bas(5, ndim)

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

! -- the energy level we r taking, smaller than the one at the equilibria, so we have a 
! closed UFO 
! for x ----
 cj0 = cj + 0.05d0
 print*, 'cj0 =', cj0; read*
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

! compute the zero velocity curve 

!call zvc_plf(cj0, 2.d0)


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
xmax0 =  2.d0  !
tmax0 =  4.d1 ! 

! Poincare section:  y=0  
ind_free  = 1       ! 1: x (or 2: y) component
ind0  = 3-ind_free  ! 2: y (or 1: x) component 

p00   = 0.d0        ! rough guess, by observing directly at the plot of the tori  
dir0  = 1           ! no requirement in the direction, just compute the first turn 
imax0 = 1           ! consider every intersection, and choose by the direction of velocity

!subroutine init_poinc( ind0, p00, dir0, imax0, tmax0, xmax0, ndim0, h00) 
call init_poinc(ind0, p00, dir0, imax0, tmax0, xmax0, ndim, cj0)

fmap_pro = 25;    fmap_retro   = 52;  fmap = 50

if(ind_free == 1) then  ! for x axis
!  if(dir == 0)   write(fnmap, fmt='(a)')   './dat/fft_plf/pc_x.dat'
!  if(dir == 1)   write (fnpro,  fmt='(a)')  './dat/fft_plf/map_pro_x.dat'
!  if(dir == -1)  write(fnretro,fmt='(a)')  './dat/fft_plf/map_retro_x.dat'

  write(fnmap, fmt='(a)')   './dat/fft_plf/pc_x.dat'
  write(fnpro,  fmt='(a)')  './dat/fft_plf/pc_pro_x.dat'
  write(fnretro,fmt='(a)')  './dat/fft_plf/pc_retro_x.dat'
  
elseif(ind_free == 2) then ! for y axis
  if(dir == 0)  write(fnmap, fmt='(a)')   './dat/fft_plf/map_y.dat'
  if(dir == 1)  write(fnpro,  fmt='(a)')  './dat/fft_plf/map_pro_y.dat'
  if(dir == -1) write(fnretro,fmt='(a)')  './dat/fft_plf/map_retro_y.dat'
  
endif  

open(fmap, file = fnmap, access ='append',status = 'replace')
write(fmap, *)  ' # Poincare map +-:  t      (x y  vx vy )   h'

open(fmap_pro, file = fnpro, access ='append',status = 'replace')
write(fmap_pro, *)  ' # Pro Poincare map :  t      (x y  vx vy )   h'

open(fmap_retro, file = fnretro  , access ='append',status = 'replace' )
write(fmap_retro, *)  ' # Retro Poincare map:  t      (x y  vx vy )   h'

! For several initial points, plot the orbits 
foball = 100
open(foball, file = './dat/fft_plf/pcob.dat', status='replace')
write(foball, *) '# orbit: :  t      (x y  vx vy )   h'

! save the x0, H0 and number of iterations for Poincare maps 
open(444, file='./dat/fft_plf/pc_init.dat',status='replace',access='append')
write(444,*) ' # H0 = ', cj0 
write(444,*) ' # npc0(start)     npoinc(N.O.)    (x0, y0, vx0, vy0) '

! -- do fourier analysis on this orbit ---
tf     = 2.d2
nitmax = 20  
nsm    = 2**16
nfmax  = 5  ! we have 3 basic frequencies at this moment, so let's assume 5 will be enough for the detection of basic frequecies
xmax   = 8.d0
  
open(60  ,file = './dat/fft_plf/fob.dat', access    ='append',status = 'replace')
write(60, *) ' #   t      (x y z vx vy vz)    cj '

open(61  ,file = './dat/fft_plf/fcs.dat', access    ='append',status = 'replace')
open(62  ,file = './dat/fft_plf/fcsbas.dat', access ='append',status = 'replace')
fx0_fabas = 555
open(fx0_fabas,file = './dat/fft_plf/x0_fabas.dat', access ='append',status = 'replace')

ispl = 1

!xi(1:9) = (/(0.31d0+i*1.d-2, i=0,8)/)  !(/ (i, i=1,MAXDIM) /)

! -- Initialize the orbits to do Fourier analysis
! try x0 \in [0.32: 0.005 : 0.24], and check the behavior of the frequecies
dx = 1.d-2  !5.d-3
nob = (0.32d0 - 0.08d0) / dx + 1  ! number of orbits we want to study 
do i = 1, nob, 1
  xi(i) = 0.32d0 - (i-1) * dx
end do

! do the orbit one-by-one to check the trend
! so better to start from the middle point x0=0.32, first explore the left-hand side space 

! for every initial points, we have to check the poincare map and see when they finish one evolution... 
! xi(1) = 0.32d0  !;  npoinc = 259
! xi(1) = 0.315d0 !;  npoinc = 259
! xi(1) = 0.31d0  !;  npoinc = 131
! xi(1) = 0.305d0 !;  npoinc = 82   
!   3.1067156584394899E-002
!   7.1394201558354116 
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
! xi(1) = 0.22d0 !;  npoinc = 35  
! xi(1) = 0.20d0 !;  npoinc = 35 
! xi(1) = 0.14d0 !;  npoinc = 35 
! xi(1) = 0.08d0 !;  npoinc = 35
!xi(1) = 0.325d0 ! discard, since the next crossing is more or less symmetric w.r.t. 0.32, so we only need to explore even more right 
! we start from the left in 0.31, so we could continue with 0.335 



print*, 'Finish file open!'; read*
!debug = 1
i = 1
do  !  i =  1, nob, 1
  if(i > nob) exit
  pvi = 0.d0          ! y  = 0,  vx = 0 
  x0  = xi(i)
  
  if(x0 <= 0.12d0) then 
   tf  = 2.5d1
   nsm = 2**19
  elseif(x0 < 0.24d0) then 
   tf = 1.d2
   nsm = 2**18
  endif    
  
  pvi(ind_free)  = x0 
  
  pvi(1) = x0    !x0
  ! the other  velocity: vy
  call cj2v_plf(pvi, cj0, ind+2, isv)
   
  ! TODO: this is with negative vy0 
!  pvi(ind+2) = -pvi(ind+2)
  
  if(isv == 0 ) then 
    print*, pvi 
    cycle
  endif  
  
  write(444, *)
  write(444, *) npc0, npoinc, pvi 
  
  
  ! -- do fourier analysis on this orbit ---
  print*; print*, '************  DO fourier analysis! ************';  print*;! read*
  
  
  !subroutine fft_ob(x0, n, lt, np, nitmax, nfmax, xmax, fob, ffc, ffbas, deriv, gr_cj) 
  open(60  ,file = './dat/fft_plf/fob.dat', access    ='append',status = 'old')
  call fft_ob(pvi, ndim, tf, nsm,  nitmax, nfmax, xmax, 60, 61, 62, nf_bas, fa_bas,  isrpt, gr_plf, gr_cjplf)  
  
  ! if we are not satisfied with the current detection of results, we can update nfmax to be a bigger value 
!  print*, 'Not satisfied with', nfmax, 'frequecies, use more frequecies ? (1 = Yes)' 
!  read*, isrpt
  if(isrpt == 1) cycle 
  
  
  i = i+1
  ! save x0 and the basic frequecies+amplitudes into file 
  open(fx0_fabas,file = './dat/fft_plf/x0_fabas.dat', access ='append',status = 'old')
  print*, 'we have detected ', nf_bas, 'basic frequecies'
  do k = 1, nf_bas, 1
    write(fx0_fabas, '(9e24.14)') x0, fa_bas(k, :)
  end do
  write(fx0_fabas, *) ! one blank line to seperate different tori
  close(fx0_fabas)
  
 
  
  ! **NOTE** at this moment, skip the Poincare map representation
  cycle 
   
   
  pv0 = pvi
  
  ! check the energy
  call gr_cjplf(pvi, cj_i)
  if (debug == 1) print*, 'State of initial point: ', pvi, 'cj = ', cj_i, 'dcj', dabs(cj_i-cj0) ; ! read*
 

  tftot = 0.d0 
  
 ! put a blank line before each block, also easier to check the data for index or ev:::i::i when we r not sure how we save the data
  if(dir == 0) then
    write(fmap_pro, *);  write(fmap_retro, *)  
  elseif( dir == 1)  then
    write(fmap_pro, *)  
  elseif(dir == -1) then 
    write(fmap_retro, *)  
  endif  
  
  ! save also the first point, since it is on the x-axis, a nice start for every curve 
  
  do k = 1, npoinc+npc0, 1
  
  ! discard the plot within poinc
  ! TODO: debug, some problem in ispl==1
   !   subroutine poinc(sti, ndim, nvar,  tf,stf, hminim, ispc, deriv, gr_cj) 
    call poinc(pvi, ndim, ndim,   tpc, pvf, hminim, ispc, gr_plf, gr_cjplf)

    if(ispc == 0) exit 
    
    ! we do more iterations after we complete one revolution... 
    if(k <= npc0)  cycle  ! skip the first npc0 iterations, in fact it doesn't affect the result
    
    if(k == npc0+1)  then 
      pv0 = pvi
    ! save also the first point, since it is on the x-axis, a nice start for every curve 
      write(fmap_pro, *)  0.d0, pvi, cj_i
      
      !print, check with npc0=0, we shoulde have pvi of form (x0, 0,0,vy+) !--ckd
!      write(*, *) 'pvi=', pvi, 'tftot=', tftot, 'cj=', cj_i, 'x0=', x0
!      print*; read*
    endif  
    
    tftot = tftot + tpc
    
    
   ! check the energy 
    call gr_cjplf(pvf, cj_i)
      
    if(dir == 0) then 
      if(pvf(ind+2) > 0.d0 ) then 
        write(fmap_pro, *)    tftot, pvf, cj_i
      else 
        write(fmap_retro, *)  tftot, pvf, cj_i
      endif   
      
    elseif( dir == 1)  then
      write(fmap_pro, *)    tftot, pvf, cj_i
      
    elseif(dir == -1) then  
      write(fmap_retro, *)  tftot, pvf, cj_i
    endif 
   
    if(debug == 1) then 
       print*, 'i,k, dir', i, k, dir;  read* 
    endif    
    
    
    if (debug == 1) then 
      print*, k, '-th return to Poincare section!'
      write(*, *) k, tftot, pvf, cj_i;  print*;  read*
    endif 
   
    pvi = pvf  ! as the new start point 
  enddo 
  
!  if(dir == 0) then
!!    write(fmap, *);  
!    write(fmap_pro, *);  write(fmap_retro, *)  
!  elseif( dir == 1)  then
!    write(fmap_pro, *)  
!  elseif(dir == -1) then 
!    write(fmap_retro, *)  
!  endif   
 

! subroutine plob(y0, t0, tf, n, tdir, ftag, ispl, deriv, gr_cj,  yf) 
!  call plob(pv0, 0.d0, tftot, 4, 1, foball, 1, gr_plf, gr_cjplf,  yf) 
!  print*, 'Finish the orbit plot!'  
  close(fmap_pro); 
  close(fmap_retro); close(fmap)

  read* 

enddo 


end 
 



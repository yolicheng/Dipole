program fft_plf_1ob
! 2017-03-14 22:30:03  at this moment, we do not need to know the basic frequency 
! we can use poincare map method directly for the computation of invariant tori 

! This is to do the Fourier analysis  for only 1 orbit in plf problem  

!  we fix the energy level and the Poincare section y=0, so we only have two degree of freedom 
!  to make it simple, we take vx0 = 0 
!  and treat vy as a function of (h, x, y=0, vx=0)

!  --- data saved in ./dat/fft_plf/1ob/1ob  subfolder 

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
implicit none

integer, parameter  ::  ndim   = 4  ! for plf

! take the initial point to be npc0+1 iteration, but I think it is not necessary, discard at the moment  
! since any curve can be obtained with vx0=0 . Original, npc0=2000    
integer, parameter  ::  npc0 = 0! 2000

! the rest are only for one revolution 
integer, parameter  ::  npoinc =  260
! x0 = 0.315   npoinc = 259

! Local Variables
integer       ::  debug, i, j, k
real(kind=dp) ::  beta0 

! eigenvalues and eigenvectors around the equilibria                 
real(kind=dp) ::  dlf(ndim, ndim), wr(ndim), wi(ndim), vr(ndim, ndim) 


! Poincare map 
integer       :: ind0, dir0, imax0, ispro, ind_free, fmap_pro, fmap_retro, fmap 
!                ind will be the Poincare section, ind_free will be the free parameter left
!                first we take ind=2, ind_free = 1, we fix y0 
!                then we take  ind=1, ind_free = 2, we fix x0 = p0   
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
 print*, 'cj0 =', cj0;   read*
 
 ! 4.3267487109222253  for equilibirum in plf
 ! 4.3767487109222252  for h3 !--- the current studied one 
 
read*

! 2017-03-14 22:27:28 

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
ind0      = 3-ind_free  ! 2: y (or 1: x) component 

! rough guess of poincare section, by observing directly at the plot of the tori 
! if ind = 2, we take y0 = 0, if ind = 1, we take x0 = 0.315 
if(ind == 1)  p0 = 0.315d0 
if(ind == 2)  p00 = 0.d0 
 
dir0  = 1           ! 1: with v>0, -1: v<0, 0: no constraint  
imax0 = 1           ! consider every intersection, and choose by the direction of velocity

! TODO: check the init_poinc 
!subroutine init_poinc( ind0, p00, dir0, imax0, tmax0, xmax0, ndim0, h00) 
call init_poinc(ind0, p00, dir0, imax0, tmax0, xmax0, ndim, cj0)

fmap_pro = 25;    fmap_retro = 52;  fmap = 50

if(ind_free == 1) then  ! for x axis

  write(fnmap, fmt='(a)')   './dat/fft_plf/1ob/pc_x.dat'
  write(fnpro,  fmt='(a)')  './dat/fft_plf/1ob/pc_pro_x.dat'
  write(fnretro,fmt='(a)')  './dat/fft_plf/1ob/pc_retro_x.dat'
  
elseif(ind_free == 2) then ! for y axis

  write(fnmap,  fmt='(a)')  './dat/fft_plf/1ob/map_y.dat'
  write(fnpro,  fmt='(a)')  './dat/fft_plf/1ob/map_pro_y.dat'
  write(fnretro,fmt='(a)')  './dat/fft_plf/1ob/map_retro_y.dat'
  
endif  

open(fmap, file = fnmap, access ='append',status = 'replace')
write(fmap, *)  ' # Poincare map +-:  t      (x y  vx vy )   h'

open(fmap_pro, file = fnpro, access ='append',status = 'replace')
write(fmap_pro, *)  ' # Pro Poincare map :  t      (x y  vx vy )   h'

open(fmap_retro, file = fnretro  , access ='append',status = 'replace' )
write(fmap_retro, *)  ' # Retro Poincare map:  t      (x y  vx vy )   h'

! For several initial points, plot the orbits 
foball = 100
open(foball, file = './dat/fft_plf/1ob/pcob.dat', status='replace')
write(foball, *) '# orbit: :  t      (x y  vx vy )   h'

! save the x0, H0 and number of iterations for Poincare maps 
open(444, file='./dat/fft_plf/1ob/pc_init.dat',status='replace',access='append')
write(444,*) ' # H0 = ', cj0 
write(444,*) ' # npc0(start)     npoinc(N.O.)    (x0, y0, vx0, vy0) '

! -- do fourier analysis on this orbit ---
tf     = 2.d2
nitmax = 20  
nsm    = 2**16
nfmax  = 5  ! we have 3 basic frequencies at this moment, so let's assume 5 will be enough for the detection of basic frequecies
xmax   = 8.d0
  
open(60  ,file = './dat/fft_plf/1ob/fob.dat', access    ='append',status = 'replace')
write(60, *) ' #   t      (x y z vx vy vz)    cj '

open(61  ,file = './dat/fft_plf/1ob/fcs.dat', access    ='append',status = 'replace')
open(62  ,file = './dat/fft_plf/1ob/fcsbas.dat', access ='append',status = 'replace')
fx0_fabas = 555
open(fx0_fabas,file = './dat/fft_plf/1ob/x0_fabas.dat', access ='append',status = 'replace')

ispl = 1

! -- Initialize the orbits to do Fourier analysis
! I.C. of the torus, from pic/normal/tori/torus1/ subfolder 
!    0.6932587044-309    0.3150000000E+00    0.0000000000E+00    0.0000000000E+00    0.1506695934E+01    0.4376748711E+01


print*, 'Finish file open!'; read*

! 2017-03-14 22:36:27  -- stop here 

!debug = 1
i = 1
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
  open(60  ,file = './dat/fft_plf/1ob/fob.dat', access    ='append',status = 'old')
  call fft_ob(pvi, ndim, tf, nsm,  nitmax, nfmax, xmax, 60, 61, 62, nf_bas, fa_bas,  isrpt, gr_plf, gr_cjplf)  
  
  ! if we are not satisfied with the current detection of results, we can update nfmax to be a bigger value 
!  print*, 'Not satisfied with', nfmax, 'frequecies, use more frequecies ? (1 = Yes)' 
!  read*, isrpt
  if(isrpt == 1) cycle 
  
  
  i = i+1
  ! save x0 and the basic frequecies+amplitudes into file 
  open(fx0_fabas,file = './dat/fft_plf/1ob/x0_fabas.dat', access ='append',status = 'old')
  print*, 'we have detected ', nf_bas, 'basic frequecies'
  do k = 1, nf_bas, 1
    write(fx0_fabas, '(9e24.14)')  x0, fa_bas(k, :)
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
   !   subroutine poinc(sti, ndim, nvar, 1, tf,stf, hminim, ispc, deriv, gr_cj) 
    call poinc(pvi, ndim, ndim, 1, tpc, pvf, hminim, ispc, gr_plf, gr_cjplf)

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
 
  close(fmap_pro); 
  close(fmap_retro); close(fmap)

  read* 

enddo 


end 
 



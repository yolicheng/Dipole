program fft_lf
! 2017-03-14 22:30:03  at this moment, we do not need to know the basic frequency 
! we can use poincare map method directly for the computation of invariant tori 

! This is to do the Fourier analysis  for only 1 orbit in lf problem  

!  we fix the energy level and the Poincare section y=0, so we only have two degree of freedom 
!  to make it simple, we take vx0 = 0 
!  and treat vy as a function of (h, x, y=0, vx=0)

!  --- data saved in ./dat/fft_lf/1ob/1ob  subfolder 

! **NOTE** 
!  be careful of the retrograde and prograde definition here,  instead of the one used for search of symmetric periodic orbit
!  we consider as cirteria the direction of vx when they pass the Poincare section.

!  1. map_pro_x.dat         -- the  poincare map with positive vy
!  2. map_retro_x.dat       -- the  poincare map with negative vy

! instead of save the return time for each map, we save the time elapsed from the initial point, good to pruduce 


use dp_mod
use pi_mod 
use lf_mod 
implicit none

integer, parameter  ::  ndim   = 6  ! for lf

! Local Variables
integer       ::  debug, i 
real(kind=dp) ::  beta0 

! eigenvalues and eigenvectors around the equilibria                 
real(kind=dp) ::  dlf(ndim, ndim), wr(ndim), wi(ndim), vr(ndim, ndim) 

! for fft_ob 
integer       ::  nfmax, nsm,  nitmax,   nf_bas,  isrpt 
real(kind=dp) ::  tf,  fa_bas(5, ndim), cj, cj0, pvi(ndim), xmax

! take eq1  in planar normal case as an example, type center X saddle  equilibria
! with beta = 1.d0, there is resonance when beta is around 1, avoid this for the moment 
debug = 0

call init_lf

! TODO: the value of beta matters! at this moment, we take beta=2 
beta0 = 2.d0   ! use 2 to check, but 1 is better?  
!print*, 'Input the value of beta (n/w_c):'
!read*, beta0 
call init_beta(beta0)   
     
print*, 'check,  beta',  beta0 

call gr_cjlf(eq, cj)
print*, 'check energy! ', cj, eq

! -- the energy level we r taking, smaller than the one at the equilibria, so we have a 
! closed UFO 
! for x ----
 cj0 = cj + 0.05d0
 print*, 'cj0 =', cj0;   read*
 
 ! 4.3267487109222253  for equilibirum in lf
 ! 4.3767487109222252  for h3 !--- the current studied one 
 
read*

! 2017-03-14 22:27:28 

! Jacobi matrix of the lorentz force with respective to the state 
call dflrtz(eq, dlf)

print*, 'Check DG in the planar case: ' 
do i = 1, ndim 
  write(*,'(6f18.8)') dlf(i,:) 
enddo
print*; read* 

call eigrg(dlf, ndim,1, wr,wi, vr)
print*, 'Check if we have a center X center X saddle equilibria:' 
print*; read* 


! -- do fourier analysis on this orbit ---
tf     = 2.d2
nitmax = 20  
nsm    = 2**16
nfmax  = 5  ! we have 3 basic frequencies at this moment, so let's assume 5 will be enough for the detection of basic frequecies
xmax   = 8.d0
  
open(60  ,file = './dat/fft_lf/1ob/fob.dat', access    ='append',status = 'replace')
write(60, *) ' #   t      (x y z vx vy vz)    cj ' !; close(60)

open(61  ,file = './dat/fft_lf/1ob/fcs.dat', access    ='append',status = 'replace')
open(62  ,file = './dat/fft_lf/1ob/fcsbas.dat', access ='append',status = 'replace')


! -- Initialize the orbits to do Fourier analysis
!    0.494069644968E+01    0.429219426309E+00    0.713275462622E+00    0.422292367972E-01   -0.215974185230E+00   -0.674220592223E+00    0.379991680857E-03    0.181287323131E+01   -0.222044604925E-14  0.1250E-03    1

pvi = (/0.429219426309d0, 0.713275462622d0, 0.422292367972d-1,  & 
        -0.21597418523d0,-0.674220592223d0, 0.379991680857d-3/)
  
    ! -- do fourier analysis on this orbit ---
print*; print*, '************  DO fourier analysis! ************';  print*;! read*
  
  !subroutine fft_ob(x0, n, lt, np, nitmax, nfmax, xmax, fob, ffc, ffbas, deriv, gr_cj) 
  call fft_ob(pvi, ndim, tf, nsm,  nitmax, nfmax, xmax, 60, 61, 62, nf_bas, fa_bas,  isrpt, gr_lf, gr_cjlf)  

stop 
end 
 



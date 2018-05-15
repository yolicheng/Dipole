program curve_tmap_prtbp
! check the Time T map method, and compare with the Poincare map method 
! use the curve by Poincare map method to obtain the initial curve for Time T map method 
!  --1. do the Fourier analysis for the return time T, and use the constant term C0 for Time T map 
!  --2. compute also the small deviation time \tau(xi), and integrate for this time to obtain points on curve for Time T map method
! 

!  Study the planar Restricted Three Body Problem, and choose L4, which is of type 
!  center X center, so around L4 we have a plenty 2D tor
!    1. start from a point that is close to L4, and integrate the orbit 
!    2. compute the Poincare maps, nmap = 10000, with Poincare section taken as x = p0 (or y= p0)
!    3. compute the rotation number of numerically computed invariant curve 
!    4. use the linear interpolation to obtain equally spaced points w.r.t. the angle between the relative position vector from L4 and the x-axis!
!    5. do a simple Fourier analysis to compute the Fourier coefficients as the initial guess to refine the curve 
!    6. Impose the invariance equations, and specify the energy level, to compute the invariant curve 
!    7. Globalize the curve to an invariant torus 


!  --- data saved in ./dat/curve_tmap_prtbp/  subfolder 
!  1. torus.dat            -- the original orbit of the 2D torus 
!  2. curve.dat            -- the invariant curve obtained by Poincare map 
!  3. curve_interpol. dat  -- linear interpolated points for Fourier analysis
!  4. fun.dat              -- approxiamted Fourier representation \varphi(xi+rho), the check \varphi(xi+rho), P( \varphi(xi) )
!  5. fcs.dat              -- fourier coefficients, For each coordinate, three parts: c0 // ck // sk


!  
use dp_mod
use pi_mod
use emconst_mod 
use rtbpconst_mod
use sort_mod

use curve_mod
!use poinc_mod  ! still we  use poinc for the initial guess of the return time T


implicit none



integer, parameter  ::  nmap  = 2000 , &   ! number of Time T map, keep the same as for Poincare map method
                        ndim0   = 4, & ! for PRTBP
                        npvar   = 20, & ! 20 for variational matrix, 4 for state vector 
                        np_interpol = 2000   ! in order to use fourier, we need odd number of points 

! --  ---- declare the size for Fourier coefficients of the original  Fourier analysis   ----
! how many Fourier modes ???, Alex suggests of power 2 
! Gerard said 60 something is enough ....
! use 2 to check... 
! TODO: apply Alex's strategy: claim a big enough array at first, then check the maximum norm 
!       of the last half, if it is greater than a given threshold, which should be a little bit smaller than 
!       the precision we ask for the Newton method,  we increase the size by 2

!       if it is smaller, then we decrease the size by half 
!       the arrays that store ck, sk are claimed dynamical arrays.                  

!       since they are just initial seed computed from a linear interpolated curve, not accurate enough as 
!       the estimate of the precision of the refined curve 
!       so we check again the refined Fourier coefficients by looking at the Maximum norm of the last half.

! TODO: JM's idea, paper  Physic D. P12. Section 2.2.6. Error estimate
!       increase the value of nf, until the error in Invariance equation is under a given tolerance. 1.d-10
!       where he uses multiple shooting with intermediate node M=50 
! 
integer, parameter  ::  nf_init  = 128 !  

! -- fourier coefs for Fourier analysis -- 
real(kind=dp), dimension(ndim0, 0:nf_init)  ::  csfn, sifn 
real(kind=dp), dimension(0: nf_init)        ::  csf, sif   ! fourier coefs of 1 component
real(kind=dp), dimension(0: nf_init)        ::  ctau, stau   ! fourier coefs for \tau


! -- fourier coefs for Refinement  -- 
! the number of Fourier modes we are going to use will be specified 
! after we check the maximum norm of ck and sk, all the following allocatable arrays have nf-columns  

integer                                      ::  nf, inf, nfmn
real(kind=dp)                                ::  c0(ndim0), ckm, skm   ! fixed size 
real(kind=dp), dimension(:, :), allocatable  ::  ck, sk


!  ------ Local Variables ----------- 
real(kind=dp) :: yi(ndim0), yf(ndim0), yf2(ndim0), tau1,  cj0, rho, pvi(ndim0), &   ! torus
                 tpc(nmap), xpc(nmap, ndim0), cj, t, T2, tpre, tau(nmap), &    ! curve 
                 ptold(nmap, ndim0), rho_new, theta(nmap),          &    ! rotnum
                 
                 ptnew(np_interpol, ndim0),  pt1(ndim0),   &    ! interpolation  
                 
                 pv_in(np_interpol), fun(ndim0),  fun_rho(ndim0),  &  ! interpol + fourier
                 darg, arg,                     & ! \varphi(xi+rho)
                 tol, tol_err, pv_curv(ndim0), xi, &         ! Netwon method
                 pv_refn(np_interpol, ndim0), rho_refn
                  
                  
integer ::  i, k, l,  ncrs, ind_sorted(nmap), dircurve,  &
            ftorus, fcurve, fcurve_poinc, fcurve_interpol, ffcs, ffun, ffun_rho,  & 
            
            iscurve, isinterpol, isfourier, & !  control which subroutine to call
            
!            ispc, &  ! P(\varphi(xi))
            nitmax, isref, opt, tdir 
             
character(len=70) :: fntorus, fncurve, fncurve_poinc, fncurve_interpol, fnfun_rho, fnfc, fnfun   

external :: TimeMap_prtbp         ! general 4d map for prtbp  
external :: gr_cjprtbp, gr_prtbp  ! Jacobi constant + vector field 

real(kind=dp) :: four_seri, dlange  !funcs

! which approach we want to take 
print*, 'Choose approach: (I prefer 3)'  
print*, '1: no constraint;  '
print*, '2: Avoid phase shift indetermination by adding one more equation'
print*, '3: Avoid phase shift indetermination by removing one cofficient C1_J = 0 or S1_J = 0'
 
read(*,*) opt 

iscurve    = 1
isinterpol = 1
isfourier  = 1

! -- assign the value of mass ratio for EM system and put the values into the common module --- 
call init_emconst
print*, 'emrat = ', emrat 
call init_rtbpconst(emrat, ndim0)

print*, 'check the mass ration, mu =', mu, 'ndim =', ndim0
read*

! files to save the original torus + invariant curve + interpolated one + approximated one 
ftorus = 21;  fntorus  =  './dat/curve_prtbp/torus.dat'
fcurve = 25;  fncurve  =  './dat/curve_tmap_prtbp/curve.dat'

fcurve_poinc = 26;  fncurve_poinc  =  './dat/curve_prtbp/curve.dat'

if(iscurve == 1) then 
  open(fcurve  ,file = fncurve  , access ='append', status = 'replace')
  write(fcurve  ,*)  ' # Invariant curve by Time T map:  t   (x y vx vy)   cj  np'
endif 

! ---------- Computation -----------------------------

! use the same value as with the poincare map method 
 cj0 = 3.0000046704015602d0 
 rho = 1.1140282949439548d0  
 
! -- Initialization for the module  curve_poinc_mod 
nitmax  = 15 
tol     = 1.d-10
tol_err = 1.d-10 


! read from the curve on Poincare section (fcurve_poinc  ./dat/curve_prtbp/curve.dat)
!open(100, file='./dat/curve_prtbp/para_poinc.dat'); read(100,*)
!read(100, *) h0, rho, t2_0, nf

pvi = (/-0.48769996309447677d0,  0.87102540378443860d0,  5.0000000000000001d-3, 5.1381386919162088d-3/) 

! read to Poincare maps from the file a 
open(fcurve_poinc, file = fncurve_poinc);  read(fcurve_poinc, *) ! skip the first line 

! -- read the Poincare maps from the data file -- 
do i = 1, nmap, 1  
  
  read(fcurve_poinc, *)  t, pt1,  cj,  ncrs
  
  if(i == 1) then 
    tpc(1) =  t ; ! xpc(1, :) = pvi  
  else
    tpc(i)  =  t - tpre
  endif 
  
  xpc(i, :) = pt1 
!  if(i < nmap)  xpc(i+1, :) = pt1 
  
  tpre = t
  
!  print*, 'tf =', tpc(i),  pt1; read*
enddo 

! -- check the rotation number of the poincare map curve....
dircurve = 1
ptold = xpc 
call rotnum( ptold(:, 2:4), nmap, dircurve, rho_new)

print*, 'rho_new = (unit cycle)', rho_new 
rho_new = rho_new * pi2
print*, 'transfrom rho to radian: ', rho_new, 'targeted rho = ', rho  
print*; read* 

! --- initial guess of return time T2 for Time T map method T2 = <T_PC> --- 

! TODO: since we save the time tpc spent to arrive at each xpc, so the return time tpc(i) is for the point xpc(i-1, :)

T2 = sum(tpc)/nmap 
print*, 'T2 = ', t2;  print*; read*

!  T2 =    20.928875635274963     ! by average 

!T2 = 20.929626993648455d0

! so the time deviation \tau should be T2- T_PC for each point 
! then for each segment of P i -- P i+1, applying the homotopology equation 
! \tau( xi_i + rho ) - \tau ( xi_i ) = TPC ( xi ) - T2
! tau(2) - tau(1) = TPC(1) - T2
! here, we take tau(1) = 0, so 
!     tau(2) = TPC(1) - T2
!     tau(3) = tau(2) + TPC(2) - T2
!   ... ... 

! for each point \hat P_i on the time T map, we can obtain from the point P_i on the Poincare map 
!  \hat P_i = P_i + \phi_ { tau(i) } (P_i)

tau(1) = 0.d0 
ptold(1, :) = xpc(1, :)
write(fcurve, *) xpc(1, :)

if(iscurve == 1) then

  do i = 1, nmap-1, 1
  
    ! we move yi for time T2-tpc(i), then the new yf should be one point on the 
    
!    if( i == 1)  yi =  pvi !xpc(1, :)
    if( i == 1)  yi =   xpc(1, :)
    call plob(yi, 0.d0, T2, ndim, 1 , 0, 0, gr_prtbp, gr_cjprtbp,  yf) 
    yi = yf 
    
    ! --- for i >  1 ----
    
    tau(i+1) = tpc(i+1) - T2 + tau(i)
    tau1 = tau(i)
    if(tau1 > 0.d0)  then 
      tdir = 1
    else 
      tdir = -1
    endif  
    
    ! subroutine plob(y0, t0, tf, n, tdir, ftag, ispl, deriv, gr_cj,  y) 
    call plob(xpc(i,:), 0.d0, tau1, ndim, tdir , 0, 0, gr_prtbp, gr_cjprtbp,  yf2) 
    
    ! check yf == yf2?? 
    print*, 'by \tau yf2,  by T2 Map yf= ' 
    print*, yf2 
    print*, yf 
    print*, 'Error: ', dabs(yf-yf2)  
    ! print*; read*
     
    ptold(i+1, :) = yf
    write(fcurve, *) yf, yf2
  end do 

  close(fcurve)
endif   



! for each point in Poincare map curve, we integrate \tau(xi) 
! the value of parameter xi for each point is defined as theta(1) = 0, theta(2) = rho, ..., 
!  given rho in radian, theta = mod(rho*ind_sorted, 2pi) 
theta = (/ ( dmod(rho*i, pi2) , i = 0, nmap-1) /) ! implicit-do 

! the new estimate angle of all the points, defined as theta(1) = 0, theta(2) = rho, ..., 
call sort_1d( theta, nmap, 1, ind_sorted)

! update the points by an increasing order of theta 
ptold = ptold(ind_sorted, :) 
  
!  call system( 'gnuplot /home/yu/Dropbox/dipole/dat/curve_tmap_prtbp/curve_rotnum.pl' )
!  print*, 'Input the sense of the curve: 1: counterclockwise, -1: clockwise'; print* 
!  read*, dircurve
!else   
!  dircurve = -1  ! save time 
!endif 

! -- interpolated curve (equally spaced in rho) --
fcurve_interpol = 26;   fncurve_interpol  =  './dat/curve_tmap_prtbp/curve_interpol.dat'

! ----- linear interpolation to get equally spaced points in angle: rho  -----------
! rho in radian is the new estiamte of the curve

if ( isinterpol == 1) then 
  open(fcurve_interpol,  file = fncurve_interpol)
  write(fcurve_interpol, *)  ' #  Interpolated curve: arg, (x, y, vx, vy)'  
  
  call  interpol( ptold, nmap, ndim, theta, np_interpol, ptnew)
  
  ! Period = 2*pi in radian w.r.t. the angle of the points 
  darg = pi2 / np_interpol ! take this value to keep coherent with the original points 
  do i = 1, np_interpol, 1
    arg = (i-1) * darg
    pt1 = ptnew(i, : ) ! the sorted equispaced points
    write(fcurve_interpol, *)  arg,  pt1 ! save all the four components? not necessary... 
  end do
  
  close(fcurve_interpol) 

else  ! ( isinterpol == 0) then 

  open(fcurve_interpol, file = fncurve_interpol) 
  read(fcurve_interpol,*) ! skip the comment line 
  do i = 1, np_interpol, 1
    read(fcurve_interpol,*) arg, pt1
    ptnew(i,:) = pt1
  end do
endif 

print*, 'check the interpolated curve!'
stop

! -- check ptnew by gnuplot --- discard at the monent 
!call system( 'gnuplot /home/yu/Dropbox/dipole/dat/curve_tmap_prtbp/curve_interpol.pl' )

! ---------------- Fourier analysis  ---------------------------
! Fourier coefficients + approximated curve, skip  the  component pos(ind) and vel(ind)

! the coefficients computed by general Fourier analysis:  gr_four.f90 
ffcs = 111;  fnfc = './dat/curve_tmap_prtbp/fcs.dat'
open(ffcs, file = fnfc) 
write(ffcs,*) '# For columns the components; for each line: Fourier Coef ck(i), sk(i), i=0,...,', nf 

sifn = 0.d0;   csfn = 0.d0 

!  we have also checked the Fouriere analysis for vx, and compare with the one obtained from energy 
do k = 1, ndim, 1 
  pv_in =  ptnew(:, k ) ! one component, the k-th column  

  call gr_four(pv_in, np_interpol, nf_init, csf, sif)
 
  csfn(k, :) = csf
  sifn(k, :) = sif
end do
 
! -- check the Maximum norm of csf and sif, and choose an appropriate value for nf 
! FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
inf = nf_init / 2
 ckm = dlange('m', ndim, inf, csfn(:, inf+1 : nf_init), ndim, 0.d0)
 skm = dlange('m', ndim, inf, sifn(:, inf+1 : nf_init), ndim, 0.d0)
 
print*, 'maximum norm of CK and SK for the last ', nf_init/2, 'Fourier harmonics'
print*,  ckm, skm 

print*; read*

if(nf_init < tol / 10.d0) then 
  nf = nf_init / 2
else 
  nf = 2 * nf_init
endif 


! -- initialization for module curve_poinc_mod 
call init_curve(ndim, nf, nitmax, tol, tol_err, opt) 

! -- save all the relevent parameter to a file 
!open(222,  file = './dat/curve_tmap_prtbp/para_poinc.dat', status = 'replace') 
!write(222, *) ' # n, nf, ind, p0, dir, cj0,  imax, tmax, ndim0, nitmax, tol, tol_err'
!write(222, *)  n, nf,  cj0,  ndim0, nitmax, tol, tol_err 
! close(222)
 
! -- declare memory for  ck, sk 
print*, 'Number of fourier modes used, nf=  ', nf, 'Original nf_init = ', nf_init 
allocate(ck(ndim, nf))
allocate(sk(ndim, nf))


sk = 0.d0; ck = 0.d0  
nfmn = min0(nf, nf_init)

print*, 'nf=', nf, 'nf_init=', nf_init, 'nfmn=', nfmn ; print*; read*
 c0 = csfn(:, 0)                   ! the constant term 
 ck(:, 1:nfmn) = csfn(:, 1: nfmn)
 sk(:, 1:nfmn) = sifn(:, 1: nfmn)

! - write the initial guess of Fourier coefficients --
do i = 1, ndim, 1
  write(ffcs, '(10e24.14)')  c0(i)   ; write(ffcs,*)
  write(ffcs, '(10e24.14)')  ck(i, :); write(ffcs,*)
  write(ffcs, '(10e24.14)')  sk(i, :) 
  write(ffcs,*);  write(ffcs,*) 
end do
 close(ffcs )

 
! ------------ Approximated Fourier representation --------
! the approximated Fourier function, to check if they match the original data 
ffun = 222;  fnfun = './dat/curve_tmap_prtbp/fun.dat'
open(ffun,  file = fnfun) 
write(ffun, *) ' # Approximated Fourier representation for (x, y, vx, vy)' 

! can check if \varphi(xi+rho) and  \varphi(xi) are the same curves in front of rho(radian)
ffun_rho = 25;  fnfun_rho  =  './dat/curve_tmap_prtbp/fun_rho.dat' 
open(ffun_rho  ,file = fnfun_rho, access ='append',status = 'replace' )
write(ffun_rho  ,*)  '# rho  \varphi(xi+rho)  P( \varphi(xi) )'  
  
   
! Evaluate the truncated Fourier series with the above coefficients
! for one coordinate, we can call four_seri
! for more than one coordinates, we need to call vaphi (in module curve_poinc_mod)

! Period = 2*pi in radian w.r.t. the angle of the points 

darg = pi2/np_interpol ! take this value to keep coherent with the original points 
do l = 1,  np_interpol, 1  
 
  arg = (l-1) * darg
  
  ! compare four_seri and varphi(from module curve_poinc_mod) 
  call varphi( arg,     ndim, nf, c0, ck, sk, fun)
  call varphi( arg+rho, ndim, nf, c0, ck, sk, fun_rho)
  
  write(ffun, *)  arg,  fun, fun_rho 
 
 
 ! TODO 
  ! check the error in the invariance equation of the first approximation of all the fourier coefficients
  ! do a poincare map to get P(\varphi(xi)) = ( \varphi(xi) ) 
  
!  subroutine poinc_tf_prtbp(pvin, tf, pvf)
!  call poinc_tf_prtbp(fun, tf, pvf )
!  tpc_fun(l) = tf 
  
!  ferr((l-1)*n+1 : l*n) = pvf(ind_fun)
  
!  write(ffun_rho, *) arg, pv_fun_rho,  pcf_rho(1:ndim), pvf

enddo   

 close(ffun)  
 
 
!  ---------- Refinement of the curve -------------
print*, 'before refine_curve_poinc, rho=', rho
read* 

!subroutine refine_curve_poinc(  rho, c0, ck, sk, isref, gr_poinc)
call refine_curve_poinc( rho, c0, ck, sk, isref, TimeMap_prtbp)

! --- use this updated c0, ck, sk to compute refined curve ---- 
open(177, file = './dat/curve_tmap_prtbp/fcs_refn.dat')
do i = 1, ndim, 1
  write(177, '(10e24.14)')  c0(i);     write(177,*) 
  write(177, '(10e24.14)')  ck(i, :);  write(177,*) 
  write(177, '(10e24.14)')  sk(i, :) 
  write(177,*); write(177,*) ! a blank line to seperate the component 
end do
 close(177)

! check the maximum norm of the refined Foureier coefficients 
 ckm = dlange('m', ndim, nf/2, ck(:, nf/2+1 : nf), ndim, 0.d0)
 skm = dlange('m', ndim, nf/2, sk(:, nf/2+1 : nf), ndim, 0.d0)
 
print*, 'maximum norm of CK and SK for the last', nf/2, 'Fourier harmonics'
print*,  ckm,  skm 
print*; read*

 
open(166, file = './dat/curve_tmap_prtbp/curve_refn.dat')
darg = pi2/np_interpol 
! Evaluate the truncated Fourier series for \hat T(xi)
do l = 1,  np_interpol, 1  
  xi = (l-1) * darg
  call varphi( xi, ndim, nf, c0, ck, sk, pv_curv)
  pv_refn(l, :) = pv_curv
  
  write(166, '(10e24.14)')  xi,  pv_curv
enddo 
 close(166)
 
! ---- compute the rotation number of the refined curve, check if it is 
!      the original fixed one...
! Evaluate the truncated Fourier series for \hat T(xi)

do l = 1,  np_interpol, 1  
  xi = (l-1) * rho
  call varphi( xi, ndim, nf, c0, ck, sk, pv_curv)
  pv_refn(l, :) = pv_curv
enddo 

call rotnum( pv_refn, np_interpol,  dircurve, rho_refn)
rho_refn = rho_refn*pi2

print*, 'Original rho: ', rho 
print*, 'rho of the refined curve: ', rho_refn 
print*; read*
 
!! ---- the final step, goes to the  torus!  
!open(188, file = './dat/curve_tmap_prtbp/torus_refn.dat')


!! TODO: fourier analysis for the return time \tilta T(xi)  
!!  we have to record each tf, which are equispaced in rho(unit radian)
!!  and store the coefficients as the ind column, just to save memory
!call gr_four(tpc_fun, np_interpol, nf,  csf, sif) 

!open(155, file = './dat/curve_tmap_prtbp/tpc_fun.dat')
!! Evaluate the truncated Fourier series for \hat T(xi)
!do l = 1,  np_interpol, 1  
! 
!  arg = (l-1) * darg
!  
!  ! the fourier representation of the return time of the approximated curve    
!  tf  = four_seri( arg,  nf, csf, sif)
!   
!! call gr_cjlf(fun, cj)
!  write(155, *)  arg, tpc_fun(l), tf  
!enddo 

  
stop

end 


 




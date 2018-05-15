program curve_prtbp

!  Study the planar Restricted Three Body Problem, and choose L4, which is of type 
!  center X center, that means we have 4D centre manifold and two families of periodic orbits embeded in it 
!  so around L4 we have plenty of 2D tori. 

!    1. start from a point that is close to L4, and integrate the orbit 
!    2. compute the Poincare maps, npoinc = 10000, with Poincare section taken as x = p0 (or y= p0)
!    3. compute the rotation number of numerically computed invariant curve 
!    4. use the linear interpolation to obtain equally spaced points w.r.t. the angle between the relative position vector from L4 and the x-axis!
!    5. do a simple Fourier analysis to compute the Fourier coefficients as the initial guess to refine the curve 
!    6. Impose the invariance equations, and specify the energy level, to compute the invariant curve 
!    7. Globalize the curve to an invariant torus 


!  --- data saved in ./dat/curve_prtbp/  subfolder 
!  1. torus.dat            -- the original orbit of the 2D torus 
!  2. curve.dat            -- the invariant curve obtained by Poincare map 
!  3. curve_interpol. dat  -- linear interpolated points for Fourier analysis
!  4. fun.dat              -- approxiamted Fourier representation \varphi(xi+rho), the check \varphi(xi+rho), P( \varphi(xi) )
!  5. fcs.dat              -- fourier coefficients, For each coordinate, three parts: c0 // ck // sk
!  6. curve_refn.dat       -- the curve recovied by Fourier coefficients
!  7. tori_refn.dat        -- the global tori by the curve 


!  
use dp_mod
use pi_mod
use emconst_mod 
use rtbpconst_mod
use sort_mod
use curve_mod
use poinc_mod 


implicit none


integer, parameter  ::  npoinc  = 2000 , & !20 !1000  ! number of Poincare maps, to start test, 11 is enough?
                        ndim0   = 4, & ! for PRTBP
                        
                        npvar   = 20, & ! 20 for variational matrix, 4 for state vector 
                        np_interpol = 1000   ! in order to use fourier, we need odd number of points 


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
integer, parameter  ::  nf_init   =  128 ! 128 !  !128 ! 64 !  !32 ! 32   !64 ! 32  !8 

! -- fourier coefs for Fourier analysis -- 
real(kind=dp), dimension(ndim0-2, 0:nf_init)  ::  csfn, sifn 
real(kind=dp), dimension(0: nf_init)          ::  csf, sif   ! fourier coefs of 1 component
     

! -- fourier coefs for Refinement  -- 
! the number of Fourier modes we are going to use will be specified 
! after we check the maximum norm of ck and sk, all the following allocatable arrays have nf-columns  

integer                                      ::  nf, inf, nfmn
real(kind=dp)                                ::  c0(ndim0-2), ckm, skm   ! fixed size 
real(kind=dp), dimension(:, :), allocatable  ::  ck, sk


!  ------ Local Variables ----------- 
real(kind=dp) :: p00, tmax0, xmax0, pvi(ndim0), tf, pvf(ndim0),  cj0,   &  ! torus
                 xl4, yl4, &                                              ! L4 equilibirum point
                 tpc(npoinc), xpc(npoinc, ndim0), cj, t, pv(ndim0), &     ! curve 
                 pt(npoinc, 2),  rho,  theta(npoinc), &                   ! rotnum
                 
                 ptold(npoinc, ndim0-2), ptnew_aux(np_interpol,ndim0-2),& ! interpolation
                 ptnew(np_interpol, ndim0),  &                            ! interpolation  for full state 
                 
                 pv_in(np_interpol), fun(ndim0-2), pv_fun(ndim0), &  ! interpol + fourier
                 arg1, ptnew1(ndim0), pvin(ndim0-2),  darg, arg,  &  ! gamm 
                 
                 fun_rho(ndim0-2), & ! \varphi(xi+rho)
                 pv_fun_rho(ndim0),  tpc_fun(np_interpol), &              ! P(\varphi(xi))
                

                 tol, tol_err, pv_curv(ndim0-2), xi, &                    ! Netwon method
                 pv_refn(np_interpol, ndim0-2), rho_refn
                  
                  
integer ::  i, k, l,  ncrs, ind_sorted(npoinc), dircurve,  &
            ftorus, ispl, fcurve, fcurve_interpol, ffcs, ffun, ffun_rho,  & 
            
            ind0, dir0, imax0,  & ! poinc 
            
            iscurve, isinterpol, isfourier, & !  control which subroutine to call
            
!            ispc, &  ! P(\varphi(xi))
            nitmax, isref, opt
             
character(len=70) :: fntorus, fncurve, fncurve_interpol, fnfun_rho, fnfc, fnfun   

external :: PoincMap_prtbp ! general 2d map for prtbp  
external :: gr_cjprtbp, gr_prtbp, cj2v_prtbp, dvind_dx_prtbp  ! Jacobi constant + vector field 

real(kind=dp) :: four_seri, dlange  !funcs
!real(kind=dp) :: MaxNorm   


! which approach we want to take 
print*, 'Choose approach: (I prefer 3)'  
print*, '1: no constraint;  '
print*, '2: Avoid phase shift indetermination by adding one more equation'
print*, '3: Avoid phase shift indetermination by removing one cofficient C1_J = 0 or S1_J = 0'
 
read(*,*) opt 


iscurve    = 0
isinterpol = 0
isfourier  = 1

! -- assign the value of mass ratio for EM system and put the values into the common module --- 
call init_emconst
print*, 'emrat = ', emrat 
call init_rtbpconst(emrat, ndim0)

print*, 'check the mass ration, mu =', mu, 'ndim =', ndim0
read*

! take a point that is close the L4 point,  with the big primary at (mu, 0, 0), 
! the triangular libration point L4 form a equilateral triangle with the two primaries. 
! so the location of L4 is ( -1/2+mu, sqrt(3)/2  ) and  L5 located at ( -1/2+mu, -sqrt(3)/2  )
xl4 = -.5d0 + mu 
yl4 = dsqrt(3.d0) / 2.d0

! the initial point --- instead of take random points, we fix the energy and poincare section 
!pvi = (/xl4 + 1.d-3, yl4 + 2.d-3, 1.d-4, -1.d-4/) ! add a small deviation on L4 to obtain an initial point
! cj0 = 3.0000046704015602  

!this is a test, to fix the energy level, and x0 = xl4,  and change the value of y, to have a look of a lot tori. 

! files to save the original torus + invariant curve + interpolated one + approximated one 
ftorus = 21;  fntorus  =  './dat/curve_prtbp/torus.dat'
fcurve = 25;  fncurve  =  './dat/curve_prtbp/curve.dat'

if(iscurve == 1) then 
  open(fcurve  ,file = fncurve  , access ='append', status = 'replace')
  write(fcurve  ,*)  ' # Invariant curve by Poincare map:  t   (x y vx vy)   cj  npc'
  
  ispl = 1;
  open(ftorus  ,file = fntorus, access ='append', status = 'replace' )
  write(ftorus ,*)  ' # Original 2D torus: t      (x y vx vy)   cj'
endif 

! ---------- Computation -----------------------------
! outer bound for escape, nondimensinal unit 1 is big enough....
xmax0  =  2.d0 
tmax0 =  5.d1 ! a rough guess?? or not...

! Poincare section: x=-0.49 by observing the plot of the torus 
ind0  = 1           ! x component 
p00   = xl4   !yl4 !-0.49d0     ! rough guess, by observing directly at the plot of the tori  

dir0  = 1           ! take the positive velocity 
imax0 = 1           ! consider every intersection, and choose by the direction of velocity

!this is a test, to fix the energy level, and x0 = xl4,  and change the value of y, to have a look of a lot tori. 
 cj0 = 3.0000046704015602  
call init_poinc( ind0, p00, dir0, imax0, tmax0, xmax0, ndim0, cj0) 


! ---- compute the Poincare map----   
! -- Initialization for the module  curve_mod 
nitmax  = 15 
tol     = 1.d-10
tol_err = 1.d-10 

!subroutine init_curve(n0, nitmax0, tol0, tol_err0, opt0) 
!call init_curve(n, nf, nitmax, tol, tol_err, opt) 

! fix also the Poincare section, and take points from this section
pvi(ind) = p00 
 
!  pvi = (/xl4 + 1.d-3, yl4 + 2.d-3, 1.d-4, -1.d-4/) ! add a small deviation on L4 to obtain an initial point
!-0.48769996309447677       0.87102540378443860        5.0000000000000001E-003   5.1381386919162088E-003   3.0000047683715820 
! test with this point 
pvi = (/-0.48769996309447677d0,  0.87102540378443860d0,  5.0000000000000001d-3, 5.1381386919162088d-3/)  

print*, 'before poinc_n, ind =', ind, 'p0 =', p0
read*

if(iscurve == 1) then 
   call poinc_n( pvi,ndim, ndim, ind, p0,  dir, fcurve, ispl, ftorus, npoinc, & 
                 imax, tmax, xmax, tf, pvf, ncrs, gr_prtbp, gr_cjprtbp)   
   close(fcurve) ! alll the data is saved already inside poinc_n routine 
endif 
   
! read to Poincare maps from the file  
open(fcurve, file = fncurve) 
read(fcurve, *) ! skip the first line 

! -- read the Poincare maps from the data file 
do i = 1, npoinc, 1  
  read(fcurve, *) t, pv,  cj,  ncrs
!  print*, t, pv, cj, ncrs; read*
  tpc(i) = t
  xpc(i, :) = pv
enddo 
  
! ------ compute the rotation number ------ 
! be careful, which two columns to use to compute the rotation number...
! the indice are specified by init_curve for ind_fun in module curve_mod 

! we have to specify the sense of the curve 
! -- check ptnew by gnuplot --- discard at the monent 
!print*
!print*, 'Observe the plot and determine the sense of the curve: 1: counterclockwise, -1: clockwise'
!read*

if(iscurve == 1) then 
  call system( 'gnuplot /home/yu/Dropbox/dipole/dat/curve_prtbp/curve_rotnum.pl' )
  print*, 'Input the sense of the curve: 1: counterclockwise, -1: clockwise'; print* 
  read*, dircurve
  
else   
  dircurve = -1  ! save time 
  
endif 

pt = xpc(:,  ind_fun) 
call rotnum( pt, npoinc, dircurve, rho )

print*, 'rho= (unit cycle)', rho  
rho = rho * pi2
print*, 'transfrom rho to radian: ', rho
read* 

! the new estimate angle of all the points, defined as theta(1) = 0, theta(2) = rho, ..., 
! given rho in radian
! theta = mod(rho*ind_sorted, 2pi) 
theta = (/ ( dmod(rho*i, pi2) , i = 0, npoinc-1) /) ! implicit-do 
call sort_1d( theta, npoinc, 1, ind_sorted)


! update the points by an increasing order of theta 
xpc = xpc(ind_sorted, :) 

! -- interpolated curve (equally spaced in rho) --
fcurve_interpol = 26;   fncurve_interpol  =  './dat/curve_prtbp/curve_interpol.dat'

! ----- linear interpolation to get equally spaced points in angle: rho  -----------
! rho in radian is the new estiamte of the curve
! -- TODO, here we only do the interpolation for y and vy, skip x and vx, 
!    since x is constant and vx can be computed as a function of vx = f(H0, x, y, vy) 
!    and we can do also Foureier analysis for vx and check if the value match with the function 

if( isinterpol == 1) then 
  open(fcurve_interpol,  file = fncurve_interpol)
  write(fcurve_interpol, *)  ' #  Interpolated curve: arg, PV,  the interpolated components (', ind_fun,'), '  
  print*,'The index of components to do interpolation:', ind_fun
  read*
  
  ptold = xpc(:, ind_fun)
  call  interpol( ptold, npoinc,n, theta, np_interpol, ptnew_aux)
  
  ! Period = 2*pi in radian w.r.t. the angle of the points 
  darg = pi2/np_interpol ! take this value to keep coherent with the original points 

  do i = 1, np_interpol, 1
    arg = (i-1) * darg
    pvin = ptnew_aux(i, : ) ! the sorted equispaced points
    
    ! subroutine gamm( pvin, n, ind_fun, h0, ind, p0, pv, cj2v) !--ckd
    call gamm(pvin, n, ind_fun, cj0, ind, p0, ptnew1, cj2v_prtbp)
    write(fcurve_interpol, *)  arg,  ptnew1 ! save all the four components? not necessary... 
    ptnew(i, :) = ptnew1
    
    !  ---  check gamm by the full state and energy    ! -- ckd, cj and gamm! ok!
!    call gr_cjprtbp(ptnew1, cj) ! to check  
!    print*, 'pv(ind_fun) = ', pvin  
!    print*, 'full state by gamm: ', ptnew1 ;  read*
!    print*, 'energy =', cj, 'H0 =', cj0; print*; read*
!    
  end do
  close(fcurve_interpol) 
endif 

print*, 'Files names for data read and write'
print*, fntorus, fncurve,  fncurve_interpol, fnfc;   print*   !ck

! -- general Fourier analysis for nf Fourier modes, with np_interpol points  
! --  read from fcurve_interpol for the array ptnew --- 
if( isinterpol == 0) then 
  open(fcurve_interpol, file = fncurve_interpol) 
  read(fcurve_interpol,*) ! skip the comment line 
  do i = 1, np_interpol, 1
    read(fcurve_interpol,*) arg1, ptnew1
    ptnew(i,:) = ptnew1
  end do
endif 

! -- check ptnew by gnuplot --- discard at the monent 
!call system( 'gnuplot /home/yu/Dropbox/dipole/dat/curve_prtbp/curve_interpol.pl' )

! ---------------- Fourier analysis  ---------------------------
! Fourier coefficients + approximated curve, skip  the  component pos(ind) and vel(ind)

! the coefficients computed by general Fourier analysis:  gr_four.f90 
ffcs = 111;  fnfc = './dat/curve_prtbp/fcs.dat'
open(ffcs, file = fnfc) 
write(ffcs,*) '# For columns the components; for each line: Fourier Coef ck(i), sk(i), i=0,...,', nf 

sifn(:, ind) = 0.d0;   csfn(:, ind) = 0.d0 

!  we have also checked the Fouriere analysis for vx, and compare with the one obtained from energy 
do k = 1, n, 1 
  pv_in =  ptnew(:, ind_fun(k) ) ! one component, the k-th column  
  call gr_foun( pv_in, np_interpol, nf_init, csf, sif)
  
  ! check gr_foun and gr_four(my routine) !--ckd, the same!
!  print*, 'csf and sif by gr_foun' !ckd
!  write(*,'(10f20.14)') csf;   write(*,'(10f20.14)') sif;   read* 
  
!  subroutine gr_four(f, n, m, csf, sif)
  call gr_four(pv_in, np_interpol, nf_init, csf, sif)
 
!  print*, 'csf and sif by my own routine gr_four' ! --ckd, the same 
!  write(*,'(10f20.14)') csf;   write(*,'(10f20.14)') sif;   read* 
  
  csfn(k, :) = csf
  sifn(k, :) = sif
end do
 
! -- check the Maximum norm of csf and sif, and choose an appropriate value for nf 
! FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
inf = nf_init / 2
 ckm = dlange('m', n, inf, csfn(:, inf+1 : nf_init), n, 0.d0)
 skm = dlange('m', n, inf, sifn(:, inf+1 : nf_init), n, 0.d0)
 
print*, 'maximum norm of CK and SK for the last ', nf_init/2, 'Fourier harmonics'
print*,  ckm, skm 

print*; read*

if(nf_init < tol / 10.d0) then 
  nf = nf_init / 2
else 
  nf = 2 * nf_init
endif 


! -- initialization for module curve_mod 
call init_curve(n, nf, nitmax, tol, tol_err, opt) 

! -- save all the relevent parameter to a file 
open(222,  file = './dat/curve_prtbp/para_poinc.dat', status = 'replace') 
write(222, *) ' # n, nf, ind, p0, dir, cj0,  imax, tmax, ndim0, nitmax, tol, tol_err'
write(222, *)  n, nf, ind, p0, dir, cj0,  imax, tmax, ndim0, nitmax, tol, tol_err 
 close(222)
 
! -- declare memory for  ck, sk 
print*, 'Number of fourier modes used, nf=  ', nf, 'Original nf_init = ', nf_init 
allocate(ck(n, nf))
allocate(sk(n, nf))


sk = 0.d0; ck = 0.d0  
nfmn = min0(nf, nf_init)

print*, 'nf=', nf, 'nf_init=', nf_init, 'nfmn=', nfmn ; print*; read*
 c0 = csfn(:, 0)                   ! the constant term 
 ck(:, 1:nfmn) = csfn(:, 1: nfmn)
 sk(:, 1:nfmn) = sifn(:, 1: nfmn)

! - write the initial guess of Fourier coefficients --
do i = 1, n, 1
  write(ffcs, '(10e24.14)')  c0(i)   ; write(ffcs,*)
  write(ffcs, '(10e24.14)')  ck(i, :); write(ffcs,*)
  write(ffcs, '(10e24.14)')  sk(i, :) 
  write(ffcs,*);  write(ffcs,*) 
end do
 close(ffcs )

 
! ------------ Approximated Fourier representation --------
! the approximated Fourier function, to check if they match the original data 
ffun = 222;  fnfun = './dat/curve_prtbp/fun.dat'
open(ffun,  file = fnfun) 
write(ffun, *) ' # Approximated Fourier representation, x(',ind_fun, ')' 

! can check if \varphi(xi+rho) and  \varphi(xi) are the same curves in front of rho(radian)
ffun_rho = 25;  fnfun_rho  =  './dat/curve_prtbp/fun_rho.dat' 
open(ffun_rho  ,file = fnfun_rho, access ='append',status = 'replace' )
write(ffun_rho  ,*)  '# rho  \varphi(xi+rho)  P( \varphi(xi) )'  
  
   
! Evaluate the truncated Fourier series with the above coefficients
! for one coordinate, we can call four_seri
! for more than one coordinates, we need to call vaphi (in module curve_mod)

! Period = 2*pi in radian w.r.t. the angle of the points 

darg = pi2/np_interpol ! take this value to keep coherent with the original points 
do l = 1,  np_interpol, 1  
 
  arg = (l-1) * darg
  
  ! compare four_seri and varphi(from module curve_mod) 
  call varphi( arg,     n, nf, c0, ck, sk, fun)
  call varphi( arg+rho, n, nf, c0, ck, sk, fun_rho)
 
  ! use gamma to get the full state 
  call gamm(fun,     n, ind_fun, cj0, ind, p0, pv_fun,     cj2v_prtbp)
  call gamm(fun_rho, n, ind_fun, cj0, ind, p0, pv_fun_rho, cj2v_prtbp)
  
  ! check if pv_fun and pv_fun_rho are the same curve... also there is phase shift w.r.t. xi 
  ! --ckd!! the same curve, and the curve of y, y(xi+rho) in front of xi has a shift in xi 
  
  write(ffun, *)  arg,  pv_fun, pv_fun_rho 
 ! 
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
call refine_curve_poinc(n, nf, rho, c0, ck, sk, isref, PoincMap_prtbp)



! --- use this updated c0, ck, sk to compute refined curve ---- 
open(177, file = './dat/curve_prtbp/fcs_refn.dat', status='replace')
write(177, *) '# refined curve: c0 // ck // sk  format -- (10e24.14)'
do i = 1, n, 1
  write(177, '(10e24.14)')  c0(i);     write(177,*) 
  write(177, '(10e24.14)')  ck(i, :);  write(177,*) 
  write(177, '(10e24.14)')  sk(i, :) 
  write(177,*); write(177,*) ! a blank line to seperate the component 
end do
 close(177)

! check the maximum norm of the refined Foureier coefficients 
 ckm = dlange('m', n, nf/2, ck(:, nf/2+1 : nf), n, 0.d0)
 skm = dlange('m', n, nf/2, sk(:, nf/2+1 : nf), n, 0.d0)
 
print*, 'maximum norm of CK and SK for the last', nf/2, 'Fourier harmonics'
print*,  ckm,  skm 
print*; read*

 
open(166, file = './dat/curve_prtbp/curve_refn.dat')
darg = pi2/np_interpol 
! Evaluate the truncated Fourier series for \hat T(xi)
do l = 1,  np_interpol, 1  
  xi = (l-1) * darg
  call varphi( xi, n, nf, c0, ck, sk, pv_curv)
  pv_refn(l, :) = pv_curv
  
!  subroutine gamm( pvin, n, ind_fun, h0, ind, p0, pv, cj2v)
  call gamm(pv_curv, n, ind_fun, cj0, ind, p0, pv, cj2v_prtbp)
  write(166, '(10e24.14)')  xi,  pv
enddo 
 close(166)
 
! ---- compute the rotation number of the refined curve, check if it is 
!      the original fixed one...
! Evaluate the truncated Fourier series for \hat T(xi)

do l = 1,  np_interpol, 1  
  xi = (l-1) * rho
  call varphi( xi, n, nf, c0, ck, sk, pv_curv)
  pv_refn(l, :) = pv_curv
enddo 

call rotnum( pv_refn, np_interpol,  dircurve, rho_refn)
rho_refn = rho_refn*pi2

print*, 'Original rho: ', rho 
print*, 'rho of the refined curve: ', rho_refn 
print*; read*
 
! ---- the final step, goes to the  torus!  
open(188, file = './dat/curve_prtbp/torus_refn.dat')

  
! TODO: fourier analysis for the return time \tilta T(xi)  
!  we have to record each tf, which are equispaced in rho(unit radian)
!  and store the coefficients as the ind column, just to save memory
call gr_four(tpc_fun, np_interpol, nf,  csf, sif) 
sifn(ind, :) = sif;  csfn(ind, :) = csf

open(155, file = './dat/curve_prtbp/tpc_fun.dat')
! Evaluate the truncated Fourier series for \hat T(xi)
do l = 1,  np_interpol, 1  
 
  arg = (l-1) * darg
  
  ! the fourier representation of the return time of the approximated curve    
  tf  = four_seri( arg,  nf, csf, sif)
   
! call gr_cjlf(fun, cj)
  write(155, *)  arg, tpc_fun(l), tf  
enddo 


! -- compute the Fourier coefficients ctau, stau for the small time deviation tau(xi)
! but at this moment, I don't know why we compute ctau and stau... 

!sif(0) = 0.d0;           csf(0) = c0(ind) 
!sif(1:nf) = sk(ind, :);  csf(1:nf) = ck(ind, : )

!!subroutine tau_fc( ct, st, rho,  nf, ctau, stau)
!call tau_fc( csf, sif, rho,  nf, ctau, stau)
  
stop

end 


 




program curve_plf
! old version with an arbitraray torus is saved in main_curve_plf_pc.f90, and probably discard forever.

! that takes the linear flow around an elliptic periodic orbit, more specific, centre X saddle (6D phase space) and centre (4D phase space), 
! as the initial guess of an invariant curve
! for Poincare map approach, we need one more step to integrate this curve for the next crossing with the Poincare section. 

! ** NOTE ** 
! DO not forget the initialization of the model parameters:   beta, sgn, ndim  

!  Study the planar norm case of lorentz force problem, take an elliptic periodic orbit, of type center X center
!     around it, we have plenty 2D tori, the one around x0 =0.35.  
 
!    1. start from a known torus, use fft to check the basic frequencies. 

!    2. compute the Poincare maps, np_interpoloinc = 10000, with Poincare section taken as x = p0 (or y= p0)
!    3. compute the rotation number of numerically computed invariant curve 
!    4. use the linear interpolation to obtain equally spaced points w.r.t. the angle between the relative position vector from L4 and the x-axis!
!    5. do a simple Fourier analysis to compute the Fourier coefficients as the initial guess to refine the curve 
!    6. Impose the invariance equations, and specify the energy level, to compute the invariant curve 
!    7. Globalize the curve to an invariant torus 


!  --- data saved in ./dat/curve_plf/  subfolder 
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
use sort_mod
use curve_mod
use poinc_mod 
use plf_mod 


implicit none

! ** NOTE ** 
! remember to update np_interpoloinc and np_interpol for different tori 
! torus1:  np_interpoloinc = 4000,  np_interpol = 2000 
! torus2:  np_interpoloinc = 1200,  np_interpol = 1000  

integer, parameter  ::  np_interpoloinc    = 4000 !1200   ! 4000 !  20 !1000  ! number of Poincare maps
integer, parameter  ::  np_interpol     =  2000 !1000 ! 2000   ! in order to use fourier, we need odd number of points 
 
integer, parameter  ::  ndim0   = 4  ! for plf
                        

real(kind=dp) ::  beta0 

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

 
integer, parameter  ::  nf_init   =  512 ! 128 !  !128 ! 64 !  !32 ! 32   !64 ! 32  !8 

! -- fourier coefs for Fourier analysis -- 
real(kind=dp), dimension(ndim0-2, 0:nf_init)  ::  csfn, sifn 
real(kind=dp), dimension(0: nf_init)          ::  csf, sif, ctau, stau   ! fourier coefs of 1 component
     

! -- fourier coefs for Refinement  -- 
! the number of Fourier modes we are going to use will be specified 
! after we check the maximum norm of ck and sk, all the following allocatable arrays have nf-columns  

integer                                      ::  nf, inf, nfmn, niter
real(kind=dp)                                ::  c0(ndim0-2), ckm, skm   ! fixed size 
real(kind=dp), dimension(:, :), allocatable  ::  ck, sk


!  ------ Local Variables ----------- 
real(kind=dp) :: p00, tmax0, xmax0, pvi(ndim0), tf, pvf(ndim0),  cj0,   &  ! torus
                 tpc, xpc(np_interpoloinc, ndim0), cj, t, pv(ndim0),     &  ! curve 
                 pt(np_interpoloinc, 2),  rho,  theta(np_interpoloinc),                   &  ! rotnum
                 
                 ptold(np_interpoloinc, ndim0-2), ptnew_aux(np_interpol,ndim0-2),&  ! interpolation
                 ptnew(np_interpol, ndim0),  &                             ! interpolation  for full state 
                 
                 pv_in(np_interpol), fun(ndim0-2), pv_fun(ndim0), &  ! interpol + fourier
                 arg1, ptnew1(ndim0), pvin(ndim0-2),  darg, arg,  &  ! gamm 
                 
                 fun_rho(ndim0-2), & ! \varphi(xi+rho)
                 pv_fun_rho(ndim0),  tpc_fun(np_interpol), &              ! P(\varphi(xi))
                

                 tol, tol_err, pv_curv(ndim0-2), xi, &                    ! Netwon method
                 pv_refn(np_interpol, ndim0-2), rho_refn 
                  
real(kind=dp) ::  dp_dx(2), tpc_refn(np_interpol) ,  pvpc(4), tau, T2                                  
                                                      
integer ::  i, j, k, l,  ncrs, ind_sorted(np_interpoloinc), dircurve, isv2,  &
            ftorus, ispl, fcurve, fcurve_interpol, ffcs, ffun, ffun_rho,  & 
            
            ind0, dir0, imax0, tdir,  & ! poinc 
            
            iscurve, isinterpol, isfourier, & !  control which subroutine to call
            
!            ispc, &  ! P(\varphi(xi))
            nitmax, isref, opt, ispc 
             
character(len=70) :: fntorus, fncurve, fncurve_interpol, fnfun_rho, fnfc, fnfun   

external :: PoincMap_plf ! general 2d map for plf  
!external :: gr_cjplf, gr_plf, cj2v_plf, dvind_dx_plf  ! Jacobi constant + vector field 

real(kind=dp) :: four_seri,  dlange  !funcs
!real(kind=dp) :: MaxNorm   

call init_plf
beta0 = 2.d0 
call init_beta(beta0)

! which approach we want to take 
print*, 'Choose approach: (I prefer 3)'  
print*, '1: no constraint;  '
print*, '2: Avoid phase shift indetermination by adding one more equation'
print*, '3: Avoid phase shift indetermination by removing one cofficient C1_J = 0 or S1_J = 0'
 
read(*,*) opt 

iscurve    = 1   ! 1: Do poincare map of the torus to get the curve; other: curve already computed. 
isinterpol = 1   ! 1: Do linear interpolation to get a curve equi-spaced in rho 
isfourier  = 1   ! 1: Do fourier analysis


! files to save the original torus + invariant curve + interpolated one + approximated one 
ftorus = 21;  fntorus  =  './dat/curve_plf/torus.dat'
fcurve = 25;  fncurve  =  './dat/curve_plf/curve.dat'


if(iscurve == 1) then 
  open(fcurve  ,file = fncurve  , access ='append', status = 'replace')
  write(fcurve  ,*)  ' # Invariant curve by Poincare map:  t   (x y vx vy)   cj  np_interpolc'
  
  ispl = 1;
  open(ftorus  ,file = fntorus, access ='append', status = 'replace' )
  write(ftorus ,*)  ' # Original 2D torus: t      (x y vx vy)   cj'
endif 

! ---------- Computation -----------------------------
! outer bound for escape, nondimensinal unit 1 is big enough....
xmax0 = 2.d0 
tmax0 = 5.d1        ! a rough guess, problem-dependent  
tdir  = 1           ! integrate forward.  
dir0  = 1           ! 1: positive velocity -1: v<0 
imax0 = 1           ! consider every intersection, and choose by the direction of velocity

! Poincare section:  PV(ind0) = p00 
ind0  = 2           ! 1: x,   2: y 

! rough guess of poincare section, by observing directly at the plot of the tori 
! if ind = 2, we take y0 = 0, if ind = 1, we take x0 = 0.315 
if(ind0 == 1)  p00 = 0.315d0 
if(ind0 == 2)  p00 = 0.d0 

!  ****   I.C. of torus 1 : x-y-vx-vy-hn  normal case  *****
 cj0 = 4.3767487109222252d0  !--- the current studied one, h3+0.05

pvi = (/0.315d0, 0.d0, 0.d0, 0.1506695934d1/)          ! torus1 
!pvi = (/0.615d0, 0.d0,  0.d0, 9.9793834493809713d-2/)    ! torus2 
     
! torus 3, four quadrants  --- discard this one, at different energy level. 
!  0.30322998672270002        0.0000000000000000        5.5018786752431001        0.0000000000000000       -23.399169929385398     
!pvi = (/0.30322998672270002d0, 0.d0,  5.5018786752431001d0, 0.14447870484784d1 /)

!subroutine cj2v_plf(X, cj0, indv, isv2)
call cj2v_plf(pvi, cj0,  4, isv2)
call gr_cjplf(pvi,  cj0)
print*, 'isv2, pvi, cj0', isv2, pvi, cj0; print*; read* 
 
! --- check if this the tori we want to compute? 
!subroutine plob(y0, t0, tf, ndim, nvar, tdir, ispl, ftag, deriv, gr_cj,  y) 
open(38, file = './dat/curve_plf/ob_test.dat', access ='append', status = 'replace' )
call plob(pvi, 0.d0, 35.d0, 4, 4, 1, 1, 38, gr_plf, gr_cjplf,  pvf) 
print*, 'Finish the orbit, check file ob_test.dat'; close(38)
print*; read*


!TODO: figure out yhy fix the energy level? for the compuation of initial v_y? 
!subroutine init_poinc( ind0, p00, dir0, imax0, tmax0, xmax0, ndim0, h00) 
call init_poinc( ind0, p00, dir0, imax0, tmax0, xmax0, ndim0, cj0) 


! ---- compute the Poincare map----   
! for torus1, np_interpoloinc = 260 to complete a closed poincare map curve. 

! -- Initialization for the module  curve_mod 
nitmax  = 15 
tol     = 1.d-10
tol_err = 1.d-10 

! fix also the Poincare section, and take points from this section
pvi(ind) = p0  

print*, 'before poinc_n, ind =', ind, 'p0 =', p0
read*


if(iscurve == 1) then 
!subroutine poinc_n(xi, ndim, np_interpol, tdir, fpc, np_interpoloinc, tf,xf, ncrs, deriv, gr_cj) 
   call poinc_n( pvi, ndim0, ndim0, dir, fcurve, np_interpoloinc, tf, pvf, ncrs, gr_plf, gr_cjplf)   
   close(fcurve) ! alll the data is saved already inside poinc_n routine 
endif 
   
! read to computed Poincare maps from the file  
open(fcurve, file = fncurve, status='old') 
read(fcurve, *) ! skip the first line 

! -- read the Poincare maps from the data file 
do i = 1, np_interpoloinc, 1  
  read(fcurve, *) t, pv,  cj,  ncrs
!  print*, t, pv, cj, ncrs; read*
!  tpc(i) = t
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
  call system( 'gnuplot /home/yu/Dropbox/dipole/dat/curve_plf/curve_rotnum.pl' )
  print*, 'Inp_interpolut the sense of the curve: 1: counterclockwise, -1: clockwise'; print* 
  read*, dircurve
  
else   
  dircurve = -1  ! to save time 
endif 


pt = xpc(:,  ind_fun) 
call rotnum( pt, np_interpoloinc, dircurve, rho )

print*, 'rho= (unit cycle)', rho  
rho = rho * pi2
print*, 'transfrom rho to radian: ', rho
read* 

! the new estimate angle of all the points, defined as theta(1) = 0, theta(2) = rho, ..., 
! given rho in radian
! theta = mod(rho*ind_sorted, 2pi) 
theta = (/ ( dmod(rho*i, pi2) , i = 0, np_interpoloinc-1) /) ! implicit-do 
call sort_1d( theta, np_interpoloinc, 1, ind_sorted)


! update the points by an increasing order of theta 
xpc = xpc(ind_sorted, :) 

! -- interpolated curve (equally spaced in rho) --
fcurve_interpol = 26;   fncurve_interpol  =  './dat/curve_plf/curve_interpol.dat'

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
  call  interpol( ptold, np_interpoloinc,n, theta, np_interpol, ptnew_aux)
  
  ! Period = 2*pi in radian w.r.t. the angle of the points 
  darg = pi2/np_interpol ! take this value to keep coherent with the original points 

  do i = 1, np_interpol, 1
    arg = (i-1) * darg
    pvin = ptnew_aux(i, : ) ! the sorted equispaced points
    
    ! subroutine gamm( pvin, n, ind_fun, h0, ind, p0, pv, isv2, cj2v) !--ckd
    call gamm(pvin, n, ind_fun, cj0, ind, p0, ptnew1, isv2, cj2v_plf)
    
    write(fcurve_interpol, *)  arg,  ptnew1 ! save all the four components? not necessary... 
    ptnew(i, :) = ptnew1
    
    !  ---  check gamm by the full state and energy    ! -- ckd, cj and gamm! ok!
!    call gr_cjplf(ptnew1, cj) ! to check  
  
!    print*, 'pv(ind_fun) = ', pvin  
!    print*, 'full state by gamm: ', ptnew1 ;   
!    print*, 'energy =', cj, 'H0 =', cj0;  print*; read*
    
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
!call system( 'gnuplot /home/yu/Dropbox/dipole/dat/curve_plf/curve_interpol.pl' )

! ---------------- Fourier analysis  ---------------------------
! Fourier coefficients + approximated curve, skip  the  component pos(ind) and vel(ind)

! the coefficients computed by general Fourier analysis:  gr_four.f90 
ffcs = 111;  fnfc = './dat/curve_plf/fcs.dat'
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
call init_curve(n, nitmax, tol, tol_err, opt) 

! -- save all the relevent parameter to a file 
open(222,  file = './dat/curve_plf/para_poinc.dat', status = 'replace') 
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
ffun = 222;  fnfun = './dat/curve_plf/fun.dat'
open(ffun,  file = fnfun) 
write(ffun, *) ' # Approximated Fourier representation, x(',ind_fun, ')' 

! check if \varphi(xi+rho) and  \varphi(xi) are the same curves in front of rho(radian)
ffun_rho = 25;  fnfun_rho  =  './dat/curve_plf/fun_rho.dat' 
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
 
  ! use gamma to get the full state   ! TODO: check
  call gamm(fun,     n, ind_fun, cj0, ind, p0, pv_fun, isv2, cj2v_plf)
  call gamm(fun_rho, n, ind_fun, cj0, ind, p0, pv_fun_rho, isv2, cj2v_plf)
  
!  call gr_cjplf(pv_fun, cj)
!  print*, 'fun ', fun  
!  print*, 'full state by gamm: ', pv_fun ;   
!  print*, 'energy =', cj, 'H0 =', cj0;  print*; read*

  ! check if pv_fun and pv_fun_rho are the same curve... also there is phase shift w.r.t. xi 
  ! --ckd!! the same curve, and the curve of y, y(xi+rho) in front of xi has a shift in xi 
  
  write(ffun, *)  arg,  pv_fun, pv_fun_rho 
 ! 
 ! TODO 
  ! check the error in the invariance equation of the first approximation of all the fourier coefficients
  ! do a poincare map to get P(\varphi(xi)) = ( \varphi(xi) ) 
  
!  subroutine poinc_tf_plf(pvin, tf, pvf)
!  call poinc_tf_plf(fun, tf, pvf )
!  tpc_fun(l) = tf 
  
!  ferr((l-1)*n+1 : l*n) = pvf(ind_fun)
  
!  write(ffun_rho, *) arg, pv_fun_rho,  pcf_rho(1:ndim), pvf

enddo   

 close(ffun)  
 
 
!  ---------- Refinement of the curve -------------
print*, 'before refine_curve_poinc, rho=', rho
read* 

!subroutine refine_curve_poinc(n, nf,  rho, c0, ck, sk, isref, gr_poinc)
call refine_curve_poinc(n, nf, rho, c0, ck, sk, isref, niter, PoincMap_plf)

! --- use this updated c0, ck, sk to compute refined curve ---- 
open(177, file = './dat/curve_plf/fcs_refn.dat', status='replace', access='append')
write(177, *) '# refined curve: c0 // ck // sk  format -- (10e24.14), rho(radian)=',rho
do i = 1, n, 1
  write(177, '(10e24.14)')  c0(i);     write(177,*) 
  write(177, '(10e24.14)')  ck(i, :);  write(177,*) 
  write(177, '(10e24.14)')  sk(i, :) 
  write(177,*);  write(177,*) ! two blank line2 to seperate components 
end do
 close(177)

! check the maximum norm of the refined Foureier coefficients 
 ckm = dlange('m', n, nf/2, ck(:, nf/2+1 : nf), n, 0.d0)
 skm = dlange('m', n, nf/2, sk(:, nf/2+1 : nf), n, 0.d0)
 
print*, 'maximum norm of CK and SK for the last', nf/2, 'Fourier harmonics'
print*,  ckm,  skm 
print*; read*
 
open(166, file = './dat/curve_plf/curve_refn.dat', access='append',status='replace')
write(166,*)  '# The refined curve:  xi     (x,y,vx,vy)    cj.  rho(radian) = ', rho 

darg = pi2/np_interpol 
! Evaluate the truncated Fourier series for \hat T(xi)
do l = 1,  np_interpol, 1  
  xi = (l-1) * darg
  call varphi( xi, n, nf, c0, ck, sk, pv_curv)
  pv_refn(l, :) = pv_curv
  
  call gamm(pv_curv, n, ind_fun, cj0, ind, p0, pv, isv2, cj2v_plf)
  
  ! compute the return time 
  call PoincMap_tf_plf(pv_curv, tpc, pvpc, 0, dp_dx, ispc)
  tpc_refn(i) =  tpc
  
  write(166, '(7e24.14)')  xi,  pv,  cj0, tpc 
    
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
 
! ------ compute the Fourier representation of \tau(xi)------------------
open(255, file = './dat/curve_plf/para_curve.dat', status='replace', access='append')
write(255, *) '# curve: H0, rho, T2, nf'

open(155, file = './dat/curve_plf/tpc_fun.dat', status='replace')
write(155,*) '# Return time: xi, tpc(from poinc), tpc_fun(Fourier), tau_fun(Fourier)'


call gr_four(tpc_refn, np_interpol, nf, csf, sif) 

! -- T2  is the constant term 
T2 = csf(0)

print*, 'We have rho =', rho, 'T2 = ', T2;  ! read*

! we need to save the energy, rho and T2 to the file... 
  
write(255, *)  cj0, rho, T2, nf
! -- compute the Fourier coefficients ctau, stau for the small time deviation 
!    from the Time-T map method  tau(xi)
  
call tau_fc(csf, sif, rho, nf, ctau, stau)
  
open(355, file = './dat/curve_plf/time_curve.dat', status='replace', access='append')
write(355, *) '# The deviated time-t2 curve from Poincare map curve: xi, tau, (x,y,vx, vy)'

open(455, file = './dat/curve_plf/tori_tau.dat', status='replace')
write(455,*) '# The orbits that connects the two curves)'

! Evaluate the truncated Fourier series for \hat T(xi)
do j = 1, np_interpol, 1  
  xi = (j-1) * darg 
  ! the fourier representation of the return time of the approximated curve 
  ! for 1d function, call the funtion four_seri is better    
  tpc = four_seri(xi, nf, csf,  sif)
  tau = four_seri(xi, nf, ctau, stau)
  write(155, *)  xi, tpc_refn(j), tpc, tau  
  
! for each point  integrate for time tau to obtain the time-T2 map, 
! check if it is a curve???  
!  subroutine plob(y0, t0, tf, ndim, nvar, tdir, ispl, ftag, deriv, gr_cj,  y) 
  pv_curv = pv_refn(j,:)
  call gamm(pv_curv, n, ind_fun, cj0, ind, p0, pv, isv2, cj2v_plf)
  call plob(pv, 0.d0, tau, 4, 4, 1, 1, 455, gr_plf, gr_cjplf,  pvf) 
  write(355, *) xi, tau, pvf   
enddo 
  

 close(255); close(155)
 

! ---- the final step, goes to the  torus!   
! 2017-03-16 10:17:47 --we skip this step, and put it into program tori_plf, just do separate routines for convenience.
 
!open(188, file = './dat/curve_plf/torus_refn2.dat',access='append', status='replace')
!! We have tried the approach of discretisized grid of theta1 and theta2
!! but it turns out a full expansion within [0,2pi] can not produce a complete torus 

!tf = 3.d2
!! Final approach: integrate from one point pv for a long enough time, say, 260 poincare maps
!call plob(pv, 0.d0, tf, 4, 4, 1, 1, 188, gr_plf, gr_cjplf,  pvf) 
!print*, 'pv, pvf'
!print*, pv; print*, pvf 
!print*, 'Finish the global torus!'


stop 

  
! TODO: fourier analysis for the return time \tilta T(xi)  
!  we have to record each tf, which are equispaced in rho(unit radian)
!  and store the coefficients as the ind column, just to save memory
call gr_four(tpc_fun, np_interpol, nf,  csf, sif) 
sifn(ind, :) = sif;  csfn(ind, :) = csf

open(155, file = './dat/curve_plf/tpc_fun.dat',access='append', status='replace')
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


 




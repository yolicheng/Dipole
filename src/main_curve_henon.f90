program curve_henon
! To modify for henon map from   curve_prtbp ! haven't started yet... 2017-03-15 18:41:02 

!  Study the planar Restricted Three Body Problem, and choose L4, which is of type 
!  center X center, so around L4 we have a plenty 2D tori. 
!    1. start from a point that is close to L4, and integrate the orbit 
!    2. compute the Poincare maps, npoinc = 10000, with Poincare section taken as x = p0 (or y= p0)
!    3. compute the rotation number of numerically computed invariant curve 
!    4. use the linear interpolation to obtain equally spaced points w.r.t. the angle between the relative position vector from L4 and the x-axis!
!    5. do a simple Fourier analysis to compute the Fourier coefficients as the initial guess to refine the curve 
!    6. Impose the invariance equations, and specify the energy level, to compute the invariant curve 
!    7. Globalize the curve to an invariant torus 


!  --- data saved in ./dat/curve_henon/  subfolder 
!  1. torus.dat         -- the original orbit of the 2D torus 
!  2. curve.dat         -- the invariant curve obtained by Poincare map 
!  3. ob_interpol. dat  -- linear interpolated points for Fourier analysis
!  4. fun.dat           -- approximated orbits by the Fourier representation 
!  
use dp_mod
use pi_mod
use emconst_mod 
use rtbpconst_mod
use sort_mod
use curve_poinc_mod

implicit none
!integer, parameter  ::  dp = kind(1.d0)
integer, parameter  ::  npoinc =  2000 , & !20 !1000  ! number of Poincare maps, to start test, 11 is enough?
                        ndim0   = 4, & ! for PRTBP
                        
                        npvar   = 20, & ! 20 for variational matrix, 4 for state vector 
                        np_interpol = 1000, &  ! in order to use fourier, we need odd number of points 
                        nf   =  64 ! 32  !8 
                        ! how many Fourier modes ???, Alex suggests of power 2, say, 32
                        ! Gerard said 5 is enough ....
                        ! use 2 to check... 
                         

! Local Variables
real(kind=dp) :: p0, tmax, xmax, pvi(ndim0), tf, pvf(ndim0),  cj0,   &  ! torus
                 xl4, yl4, & ! L4 equilibirum point
                 tpc(npoinc), xpc(npoinc, ndim0), cj, t, pv(ndim0), & !curve 
                 pt(npoinc, 2),  rho,  theta(npoinc), & ! rotnum
                 
                 ptold(npoinc, ndim0-2), ptnew_aux(np_interpol, ndim0-2),  & ! interpolation
                 ptnew(np_interpol, ndim0),  & ! interpolation 
                 
                 pv_in(np_interpol), fun(ndim0-2), pv_fun(ndim0), &  ! interpol + fourier
                 arg1, ptnew1(ndim0), pvin(ndim0-2),  darg, arg,  & ! gamm 
                 
                 pci(npvar),  pcf_rho(npvar), hminim, fun_rho(ndim0-2), & ! \varphi(xi+rho)
                 pv_fun_rho(ndim0),  tpc_fun(np_interpol), &              ! P(\varphi(xi))
                 c0(ndim0-2), ck(ndim0-2, nf), sk(ndim0-2, nf), &  ! coefs 
                 tol, tol_err, pv_curv(ndim0-2), xi  ! Netwon method  ! Netwon method
                 
! at this moment, we do not need vel(ind), it is just for test purpose.... 

! fourier coefs for the compoments skipping pos(ind) and vel(ind)
real(kind=dp), dimension(ndim0-2, 0:nf) ::  csfn, sifn 

real(kind=dp), dimension(0:nf)          ::  CSF, SIF   ! fourier coefs of 1 component
!                                           ctau, stau  ! fourier coefs for time deviation tau
                  
integer ::  i, k, l,  ncrs, ind_sorted(npoinc), &
            ftorus, ispl, fcurve, fcurve_interpol, ffcs, ffun, ffun_rho,  & 
            ind, dir, imax,  dircurve, & ! poinc 
            iscurve, isinterpol, isfourier, & !  control which subroutine to call
            ispc, & ! P(\varphi(xi))
            nitmax, isref, istmap
             
character(len=70) :: fntorus, fncurve, fncurve_interpol, fnfun_rho, fnfc, fnfun   

external :: poinc_prtbp  
external :: gr_cjprtbp, gr_prtbp, cj2v_prtbp, dvind_dx_prtbp  ! Jacobi constant + vector field 

real(kind=dp) ::  four_seri  !funcs

iscurve    = 0
isinterpol = 1
isfourier  = 1
istmap    = 0 ! do Fourier analysis for return time T
! ** NOTE ** we have to take a look at how the invariant curve moves, and set the value of dircurve
!            1: counterclockwise;  -1: clockwise
! we have -1 in this case... 
dircurve =  -1

! -- assign the value of mass ratio for EM system and put the values into the common module --- 
call init_emconst
print*, 'emrat = ', emrat 
call init_rtbpconst(emrat, ndim0)

print*, 'check the mass ration, mu =', mu, 'ndim=', ndim 
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
ftorus = 21;  fntorus    =  './dat/curve_henon/torus.dat'
fcurve   = 25;  fncurve  =  './dat/curve_henon/curve.dat'

if(iscurve == 1) then 
  open(fcurve  ,file = fncurve  , access ='append',status = 'replace' )
  write(fcurve  ,*)  ' # Invariant curve by Poincare map:  t      (x y vx vy)   cj  npc'
  
  ispl = 1;
  open(ftorus  ,file = fntorus, access ='append', status = 'replace' )
  write(ftorus ,*)  ' # Original 2D torus: t      (x y vx vy)   cj'
endif 


! -- compute the Poincare map   
! -- Initialization for the module  curve_poinc_mod 
nitmax  = 10 
tol     = 1.d-8
tol_err = 1.d-8
call init_curve_poinc(nf, ndim, nitmax, tol, tol_err) 

 
read*

! TODO: make a right call of poinc_n
if(iscurve == 1) then 
   call poinc_n( pvi,ndim, ndim,  1,  f npoinc, tf, pvf, ncrs, gr_prtbp, gr_cjprtbp)   
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
  
! -- compute the rotation number -- 
! be careful, which two columns to use to compute the rotation number...
! the indice are specified by init_curve_poinc for ind_fun in module curve_poinc_mod 
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
fcurve_interpol = 26;   fncurve_interpol  =  './dat/curve_henon/curve_interpol.dat'

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


! ---------------- Fourier analysis  ---------------------------
! Fourier coefficients + approximated curve, skip  the  component pos(ind) and vel(ind)

! the coefficients computed by general Fourier analysis:  gr_four.f90 
ffcs = 111;  fnfc = './dat/curve_henon/fcs.dat'
open(ffcs, file = fnfc) 
write(ffcs,*) '# For columns the components; for each line: Fourier Coef ck(i), sk(i), i=0,...,', nf 

sifn(:, ind) = 0.d0;   csfn(:, ind) = 0.d0 

do k = 1, n, 1 ! here we also do the Fourier analysis for vx, just for check purpose 
  
  pv_in =  ptnew(:, ind_fun(k) ) ! one component, the k-th column  
  call gr_foun( pv_in, np_interpol, nf, csf, sif)
  
  ! check gr_foun and gr_four(my routine) !--ckd, the same!
!  print*, 'csf and sif by gr_foun' !ckd
!  write(*,'(10f20.14)') csf;   write(*,'(10f20.14)') sif;   read* 
  
!  subroutine gr_four(f, n, m, csf, sif)
  call gr_four(pv_in, np_interpol, nf, csf, sif)
 
!  print*, 'csf and sif by my own routine gr_four'
!  write(*,'(10f20.14)') csf;   write(*,'(10f20.14)') sif;   read* 
!  
!  save the computed coefficients 
  write(ffcs,'(10f20.14)') csf ;  write(ffcs, *) !add a blank line  to seperate csf and sif
  write(ffcs,'(10f20.14)') sif 
  write(ffcs, *); write(ffcs, *)  !add two blank lines to seperate different components 
  read*   
  
  csfn(k, :) = csf
  sifn(k, :) = sif
end do
 close(ffcs)
! ------------ Approximated Fourier representation --------
! the approximated Fourier function, to check if they match the original data 
ffun = 222;  fnfun = './dat/curve_henon/fun.dat'
open(ffun,  file = fnfun) 
write(ffun, *) ' # Approximated Fourier representation, x(',ind_fun, ')' 

! can check if \varphi(xi+rho) and  \varphi(xi) are the same curves in front of rho(radian)
ffun_rho = 25;  fnfun_rho  =  './dat/curve_henon/fun_rho.dat' 
open(ffun_rho  ,file = fnfun_rho, access ='append',status = 'replace' )
write(ffun_rho  ,*)  '# rho  \varphi(xi+rho)  P( \varphi(xi) )'  
  
! for the call the refine_curve_poinc... 
  c0 = csfn(:, 0 ) ! the zero-th mode  
  ck = csfn(:, 1:nf )
  sk = sifn(:, 1:nf )
   
! Evaluate the truncated Fourier series with the above coefficients

! Period = 2*pi in radian w.r.t. the angle of the points 
!darg = pi2/np_interpol ! take this value to keep coherent with the original points 

do l = 1,  np_interpol, 1  
 
  arg = (l-1) * darg
    
  do k = 1, n ! only the components that we do Fourier analysis(y, vy)   
!   function four_seri( theta, m, cof, sif)  result(fun)
    csf = csfn(k, :); sif = sifn(k, :)
    fun(k)  = four_seri( arg,  nf, csf, sif)
    fun_rho(k) = four_seri(arg+rho, nf, csf, sif )
  end do
  
  ! compare four_seri and varphi(from module curve_poinc_mod) !--ckd ,they are the same 
!  call varphi( arg, c0, ck, sk, pv_curv)
!  print*, 'check four_seri and varphi,  only componet of index: ', ind_fun 
!  write(*, *)  arg, fun, pv_curv; read*
 
  ! use gamma to get the full state 
  call gamm(fun,     n, ind_fun, cj0, ind, p0, pv_fun,     cj2v_prtbp)
  call gamm(fun_rho, n, ind_fun, cj0, ind, p0, pv_fun_rho, cj2v_prtbp)
  
  write(ffun, *)  arg,  pv_fun, pv_fun_rho 
  
  ! ------   Poincare map and return time T, \tau ----- 
  ! do a poincare map to get P(\varphi(xi)) = \phi_t( \varphi(xi) ) 
  !subroutine poinc(sti, ndim, n, ind, p0, dir, imax, tmax, ispl, fob, tf,stf, hminim, ispc,deriv, gr_cj) 

 ! we need to compute the variational matrix, so put n=ndim*(ndim+1)
 ! initialize also the variational matrix
if(istmap == 1) then  
  pci = 0.d0
  pci(1:ndim) = pv_fun 
  if(npvar > ndim)   pci(ndim+1:npvar:ndim+1) = 1.d0 

  call poinc(pci, ndim, npvar, ind, p0,  1, 1, tmax, 0, 6, tf, pcf_rho, hminim, ispc, gr_prtbp,  gr_cjprtbp)
  tpc_fun(l) = tf 
  write(ffun_rho, *) arg, pv_fun_rho,  pcf_rho(1:ndim)
endif 

enddo   
  
 close(ffun)  
! stop 
 
 
!  ---------- Refinement of the curve -------------
print*, 'before refine_curve_poinc, rho=', rho
read* 

!subroutine refine_curve_poinc( h, rho, c0, ck, sk, isref, deriv, gr_cj, deriv_cj)
call refine_curve_poinc( cj0, rho, c0, ck, sk, isref, poinc_prtbp )

! --- use this updated c0, ck, sk to compute refined curve ---- 
open(177, file = './dat/curve_henon/fcs_refn.dat')
do i = 1, n, 1
  write(177, '(10e20.14)')  c0(i) 
  write(177, '(10e20.14)')  ck(i, :)
  write(177, '(10e20.14)')  sk(i, :) 
  write(177,*) ! a blank line to seperate the component 
end do
 close(177)
 
open(166, file = './dat/curve_henon/curve_refn.dat')

! Evaluate the truncated Fourier series for \hat T(xi)
do l = 1,  np_interpol, 1  
  xi = (l-1) * darg
  call varphi( xi, c0, ck, sk, pv_curv)
!  subroutine gamm( pvin, n, ind_fun, h0, ind, p0, pv, cj2v)
  call gamm(pv_curv, n, ind_fun, cj0, ind, p0, pv, cj2v_prtbp)
  write(166, '(10e20.14)')  xi,  pv
enddo 

stop 

  
! TODO: fourier analysis for the return time \tilta T(xi)  
!  we have to record each tf, which are equispaced in rho(unit radian)
!  and store the coefficients as the ind column, just to save memory
call gr_four(tpc_fun, np_interpol, nf,  csf, sif) 
sifn(ind, :) = sif;  csfn(ind, :) = csf

open(155, file = './dat/curve_henon/tpc_fun.dat')
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


 




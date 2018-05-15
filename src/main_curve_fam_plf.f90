program  curve_fam_plf 
! continuation along rho, fix the energy to obtain isoenergy family, which fits for a Poincare map similar to the famous phase portrait. 

!  TODO: main_curve_hfam_plf  does the continuation along energy, 

! 2017-03-16 23:38:56  fail to do along  rho  plus, try it again. 


! ** NOTE *
!  remember to copy the data files to corresponding subfolders 

! --- read the refined fourier coefficients from the file fcs_refn.dat 
!  Why don't we take the linear extrapolation for the prediction of new ck, sk 


use dp_mod
use pi_mod
use curve_mod  
use poinc_mod 
use plf_mod

implicit none

! ** NOTE **
! nf must keep coherent with the one used for curve_plf. 
! read from para_poinc.dat 
! nf = 512 for torus2 

integer, parameter ::  ndim0 = 4,       &    ! plf 
                       nf    = 512 ,    &    ! number of Fourier modes
                       np    = 2000,    &    ! number of points to test the new curve 
                       
                       ntori = 100           ! parameterization for tori  

real(kind=dp) ::  beta0 

integer  :: dircurve, i, j, k, niter, isv2                        
real(kind=dp) ::  pv_curv(2), pv(ndim0), pvf(ndim0), tf,   & ! poinc_tf_plf
                  c0(2), ck(2, nf), sk(2, nf),  pv_rho(2), &
                  c0pre(2), ckpre(2, nf), skpre(2, nf), rho_pre, &
                  c0_new(2), ck_new(2, nf), sk_new(2, nf),   & 
                  drho, drho_max, drho_min, nds  ! step control for continuaiton 
                  
 
integer    ::  isref, restart 

   
!  continuation 
real(kind=dp) ::  ckm, skm, darg,  tol, tol_err  
integer    ::  nitmax, opt
                              
!-- poinc_mod 
integer       :: dir0, imax0, ind0 
real(kind=dp) :: p00, cj0, xmax0, tmax0 

real(kind=dp) ::  rho, dxi, xi, rho_refn, pv_refn(np, 2)  ! for rho 

! for return time T and small time deviation \tau 
real(kind=dp) ::  tpc, pvpc(ndim0), tpc_refn(np), cj   
   
! for the torus  
real(kind=dp) ::  T2, theta1, theta2  
  
real(kind=dp), dimension(0:nf)  ::  csf,  sif, &    ! fourier coefs of 1 component
                                    ctau, stau    ! fourier coefs for time deviation tau
                  
real(kind=dp) ::  four_seri, dlange ! function type declaration
external ::   PoincMap_plf   


call init_plf
beta0 = 2.d0
call init_beta(beta0)



! ---- the final step, goes to the  torus!  
call system( 'gnuplot /home/yu/Dropbox/dipole/dat/curve_plf/curve_refn.pl' )
print*, 'Input the sense of the curve: 1: counterclockwise, -1: clockwise'; print* 
read*, dircurve

!dircurve = -1 ! TODO: keep coherent with the same curve 

! rough guess of poincare section, by observing directly at the plot of the tori 
ind0  = 2           ! x component 

! if ind = 2, we take y0 = 0, if ind = 1, we take x0 = 0.315 
if(ind0 == 1)  p00 = 0.315d0 
if(ind0 == 2)  p00 = 0.d0 

xmax0 =  2.d0 
tmax0 =  5.d1 ! a rough guess?? or not...

dir0  = 1           ! take the positive velocity 
imax0 = 1           ! consider every intersection, and choose by the direction of velocity

 cj0 = 4.3767487109222252d0  !--- the current studied one, h3+0.05
 
call init_poinc( ind0, p00, dir0, imax0, tmax0, xmax0,  ndim0, cj0) 


! Read a point from curve_refn, use as initial point and integrate, we will obtain a torus 
! The first line is comment, data structure: xi  (x,y,vx,vy)  cj
open(88, file = './dat/curve_plf/curve_refn.dat', status='old')
read(88,*);   read(88, *)  xi, pv, cj 

! Check the energy
call gr_cjplf(pv, cj)
print*, 'check the energy, cj == cj0?', cj, cj0;  print*; read*


! use the refined coefficients to compute a point with xi \in [0, pi2]
! -- this is read from the coefficients, but we only have x and vx, 
!    vy is computed as a function of (cj0, x, y=p0, vx)

! The first line is comment, data structure: c0 // ck // sk (10e24.14)

open(111, file = './dat/curve_plf/fcs_refn.dat', status='old')


read(111,*)
do i = 1, n, 1
  read(111, '(10e24.14)')  c0(i);    read(111, *)
  write(*,  '(10e24.14)')  c0(i) 
    
  read(111, '(10e24.14)')  ck(i, :); read(111, *)
  read(111, '(10e24.14)')  sk(i, :)
  read(111, *);   read(111, *)
  
!   -- check the read -- ! -- ckd
!  write(*, '(10e24.14)')  c0(i) 
!  write(*, '(10e24.14)')  ck(i, :)
!  write(*, '(10e24.14)')  sk(i, :) 
!  print*;  read*
end do

! check the rotation number in radian 
rho = 2.43609973300041d-2    ! torus1

!rho = 0.56911517174218806    ! torus2

! ---- compute the rotation number of the refined curve, check if it is 
!      the original fixed one...
! save the new evaluated curve, check if it overlap with curve_refn
open(12, file = './dat/curve_plf/curve_new.dat', status='replace')

dxi =  rho 
 
print*, 'c0 = ', c0 ; print*; read* 

do i = 1,  np , 1  
  xi = (i-1) * dxi
  call varphi( xi, n, nf, c0, ck, sk, pv_curv)
  write(12, *) xi,  pv_curv 
  pv_refn(i, :) = pv_curv
enddo 
 close(12)

! there is no point in computing rho using the approxiated Fourier functions, 
! since \xi is increasing monotonically. 
!call rotnum( pv_refn, np,  dircurve, rho_refn)
!rho_refn = rho_refn * pi2 

!print*, 'Original rho: ', rho 
!print*, 'rho of the refined curve: ', rho_refn 
!print*; read*


opt     = 3
nitmax  = 15 

! for continuation 1.d-5 is enough, if we want to look at one certain torus, do the refinement with higher precision again, say 1.d-10...
 
tol     = 1.d-8
tol_err = 1.d-8 

!tol     = 1.d-10
!tol_err = 1.d-10 
call init_curve(n, nf, nitmax, tol, tol_err, opt) 


!---- do the continuation w.r.t. rho
! we keep the energy, and change the value of rho

! --- use this updated c0, ck, sk to compute refined curve ---- 
open(166, file = './dat/curve_plf/curve_cont.dat', status='replace',access='append')
write(166,*)  '# The continuation of the curve family:  xi   x-y-vx-vy   cj'

open(177, file = './dat/curve_plf/fcs_cont.dat', status='replace',access='append')
write(177,*)  '# Coefs of the continued family :   (10e24.14)   c0,   ck,   sk'

open(188, file = './dat/curve_plf/fcs_rho_cont.dat', status='replace',access='append')
write(188, *)  '# For check:  rho,  c0,  ck(1, 1:10), ck(2, 1:10), sk(2, 1:10), sk(2, 1:10) '

! the original refined curve 
write(188, *)  rho,  c0,  ck(1, 1:10), ck(2, 1:10), sk(2, 1:10), sk(2, 1:10)
!write(188, *)  rho,  c0   
!write(188, '(10e24.14)')  ck(1, 1:10), ck(2, 1:10) 
!write(188, '(10e24.14)')  sk(2, 1:10), sk(2, 1:10); write(188,*)


! --------- Automatic step control --------------- 
! done -1.d-6, put into rho_famm/ subfolder 
drho  = 1.d-4 !
!drho = -1.d-6 ! -- have been computed already  
  
drho_max = 1.d-2 
drho_min = 1.d-8 


! the first curve 
restart  = 0
skpre = sk;  c0pre = c0;  ckpre = ck; rho_pre = rho

do  k = 1, 100, 1  
  
  nds = 1
  
  
  ! ------------- Predictor --------------- 
  if(restart == 1 .and. k>1 ) then 
  
    ! start with the previous one and use new dcj 
    sk = skpre;  c0 = c0pre;  ck = ckpre; rho = rho_pre 
  
  else 

    ! A linear extrapolation, introduce intermediate var c0_new, ck_new, sk_new 
    if(k > 1) then 
      c0_new =  c0 + (c0 - c0pre)  
      ck_new =  ck + (ck - ckpre)  
      sk_new =  sk + (sk - skpre) 
    endif 
  
    ! save the previous one  
    skpre = sk;  c0pre = c0;  ckpre = ck; rho_pre = rho  
  
    ! update c0, ck, sk as the guess for new curve 
    if(k > 1) then 
      c0 = c0_new; ck = ck_new; sk = sk_new 
    endif 
    
  endif 
    
  ! Update rho  
  rho  = rho + drho

  ! -- compute the linear prediction of CK SK for the new curve 
  ! put all the data in the same row for plot, other wise it's not easy to check the dependence on rho of ck,sk
  call refine_curve_poinc(n, nf, rho, c0, ck, sk, isref, niter, PoincMap_plf)
  
  ! 2017-03-17 10:05:14  -- if isref ==0, fail for refinement, get another chance with half the stepsize 
   
  ! -- modify stepsize drho based on iterations of Newton method 
  ! -- if we need to decrease drho, use the previous one as the new guess,
  !    and update to a smaller drho  
  if(niter > 5 .or. isref == 0)  then 
    drho = drho / 2 
    nds = .5d0
    restart = 1
    cycle 
  elseif(niter < 3)  then 
    nds = dmin1( dabs(2*drho), drho_max) / dabs(drho)
    drho = nds*drho 
    
  endif 
  
  if(dabs(drho) < drho_min) then 
    print*, 'Too small stepsize needed! Stop the continuation!'
    print*; read*
    exit 
  endif
  
  restart = 0 
  
  write(188, *) rho,  c0,  ck(1, 1:10), ck(2, 1:10), sk(2, 1:10), sk(2, 1:10)
!  write(188, '(10e24.14)') rho,  c0,  ck(1, 1:10), ck(2, 1:10), sk(2, 1:10), sk(2, 1:10)

  do i = 1, n, 1
    write(177, '(10e24.14)')  c0(i);     write(177,*) 
    write(177, '(10e24.14)')  ck(i, :);  write(177,*) 
    write(177, '(10e24.14)')  sk(i, :) 
    write(177,*); write(177,*) ! a blank line to seperate the component 
  end do

  ! check the maximum norm of the refined Foureier coefficients 
  ckm = dlange('m', n, nf/2, ck(:, nf/2+1 : nf), n, 0.d0)
  skm = dlange('m', n, nf/2, sk(:, nf/2+1 : nf), n, 0.d0)
 
  print*, 'maximum norm of CK and SK for the last', nf/2, 'Fourier harmonics'
  print*,  ckm,  skm 
  print*; read*
 
  darg = pi2/np  
  ! Evaluate the truncated Fourier series for \hat T(xi)
  do i = 1,  np, 1  
    xi = (i-1) * darg
    call varphi( xi, n, nf, c0, ck, sk, pv_curv)
    pv_refn(i, :) = pv_curv
  
    ! subroutine gamm( pvin, n, ind_fun, h0, ind, p0, pv, isv, cj2v)
    call gamm(pv_curv, n, ind_fun, cj0, ind, p0, pv, isv2, cj2v_plf)
    call gr_cjplf(pv, cj)
    write(166, '(6e24.14)')  xi,  pv, cj
  enddo 
 
  write(166, *);  write(166, *)
  write(177, *);  write(177, *)
 
  ! to check the result after each iteration
  close(166);   close(177);   close(188)
  open(166, file = './dat/curve_plf/curve_cont.dat', access='append', status='old')
  open(177, file = './dat/curve_plf/fcs_cont.dat', access='append', status='old')
  open(188, file = './dat/curve_plf/fcs_rho_cont.dat', access='append', status='old')

enddo  

 close(166);   close(177);   close(188)
 stop 


! ------ compute the Fourier representation of \tau(xi)------------------
! 
! Use the refined coefficients, Do a poincare map for the Fourier analysis for the return time 
! generate the points to do Fourier analysis for \tau(xi)

dxi  = pi2 / np   
do i = 1, np, 1
  xi = (i-1) * dxi 
  call varphi( xi, n, nf, c0, ck, sk, pv_curv)
  call varphi( xi+rho, n, nf, c0, ck, sk, pv_rho)
  
  ! subroutine poinc_tf_plf(pvin, tf, pvf)
  call poinc_tf_plf(pv_curv, tpc, pvpc)
  tpc_refn(i) =  tpc  
end do  

! Foureier analysis for the return time T

call gr_four(tpc_refn, np, nf, csf, sif) 

! -- T2  is the constant term 
T2 = csf(0)

print*, 'We have rho =', rho, 'T2 = ', T2;  read*

! we need to save the energy, rho and T2 to the file... 

print*, 'Only one step away from the torus! cheers!'; print*; read*
open(255, file = './dat/curve_plf/para_curve_poinc.dat', status='replace')
write(255, *) '# H0, rho, T2'
write(255, *) cj0, rho, T2
 close(255)

!TODO: 1- avoid the phase shift indetermination
!      2- continuaiton of the tori 
!      3- check why there is gap in the torus as a function of (theta1, theta2)  
!      4- use infinity norm to decide how many fourier modes to keep, 
!        -4.1 declare big enough array to store the CK, SK, and introduce a global parameter nf according to 
!             the precision we want to achieve, so 1.d-10 
!        -4.2 look at the error control strategy in Josep Maria's papar, Physic D

open(155, file = './dat/curve_plf/tpc_fun.dat', status='replace')
write(155,*) '# Return time: tpc, tpc_fun(Fourier)'

! Evaluate the truncated Fourier series for \hat T(xi)
do i = 1,  np, 1  
  xi = (i-1) * dxi
  
  ! the fourier representation of the return time of the approximated curve 
  ! for 1d function, call the funtion four_seri is better    
  tpc = four_seri( xi, nf, csf, sif)
   
  write(155, *)  xi, tpc_refn(i), tpc  
enddo 


! -- compute the Fourier coefficients ctau, stau for the small time deviation tau(xi)

!subroutine tau_fc( ct, st, rho,  nf, ctau, stau)
call tau_fc( csf, sif, rho, nf, ctau, stau)
  
  
! instead of the parameterization, just integrate one(any) point from the curve 
! for a long enough time 
!dxi =  pi2 / 1
!tf = 1.d5 

!open(33, file = './dat/curve_plf/torus_intg.dat')

!do i = 1, 1,  1

! ! evaluate a point on the curve 
! xi = (i-1)*dxi
! ! the point on the curve \varphi w.r.t. (theta1, theta2)
! call varphi( xi, n, nf, c0, ck, sk, pv_curv) 
! call gamm(pv_curv, n, ind_fun, h0, ind, p0, pv, isv2, cj2v_plf)
! 
! call gr_cjplf(pv, cj)
! print*, 'check the energy: ', cj; read*
! 
! call plob(pv, 0.d0, tf, ndim, 1, 33, 1, gr_plf, gr_cjplf, pvf)
! 
! call gr_cjplf(pvf, cj)
! print*, 'check the energy after integration: ', cj; read*
! 
!end do 


! - the parameterization of the torus, \psi
! \psi(theta1, theta2) = \Phi_tf  ( \varphi( theta1 - theta2/pi2*rho ) )
! where, tf = theta2/pi2*T2 + \tau( theta1 - theta2/pi2*rho )

! but in order to get a full curve, another way is 
! discretisize theta1 and xi in an uniformed grid, compute the corresponding theta2
! since a full period in xi for a full curve 

open(333, file = './dat/curve_plf/torus_refn_full.dat', status='replace')
open(444, file = './dat/curve_plf/curve_for_torus_refn_full.dat', status='replace')

dxi =  pi2 / ntori

do i = 1, ntori+1,  1

  xi = (i-1)*dxi
  
  ! the point on the curve \varphi w.r.t. (theta1, theta2)
  ! but \varphi essentially only depends on xi
  call varphi( xi, n, nf, c0, ck, sk, pv_curv) 
  call gamm(pv_curv, n, ind_fun, h0, ind, p0, pv, isv2, cj2v_plf)
    
  do j = 1, ntori + 1, 1
  
!    ! given theta1 and \vaphi,  compute theta 2 
!    ! xi = theta1 - theta2 / pi2 * rho

!    theta1 =  (j-1)*dxi 
!    theta2 = (theta1 - xi) * pi2 / rho 
!    
    theta2 =  (j-1)*dxi 
    theta1 =  xi - theta2 / pi2 * rho
    
    write(444, *) theta1, theta2, xi, pv
    
    ! integrate using the pv=\varphi(xi) for time tf
    tf = four_seri( xi, nf, ctau, stau) + theta2 / pi2 * T2

    call plob(pv, 0.d0, tf, 4, 1, 6, 0, gr_plf, gr_cjplf, pvf)
    
   ! check also the energy 
    call gr_cjplf(pvf, cj)
    write(333, *)  theta1, theta2, pvf,  cj 
    
  end do
  
  write(333, *); write(444,*)
end do 

stop

!open(111, file = './dat/curve_plf/torus_refn.dat')
!open(222, file = './dat/curve_plf/curve_for_torus_refn.dat')

!dxi =  pi2 / ntori

!do i = 1, ntori+1,  1

!  theta1 = (i-1)*dxi
!  
!  do j = 1, ntori+1, 1
!    theta2 = (j-1)*dxi 
!    
!    ! for each theta1 and theta 2, the new angle for \vaphi 
!    xi = theta1 - theta2 / pi2 * rho
!    
!    ! the point on the curve \varphi w.r.t. (theta1, theta2)
!    call varphi( xi, n, nf, c0, ck, sk, pv_curv) 
!    call gamm(pv_curv, n, ind_fun, h0, ind, p0, pv, isv, cj2v_plf)
!    write(222, *) theta1, theta2,  xi, pv
!    
!    ! integrate using the pv=\varphi(xi) for time tf
!    tf = four_seri( xi, nf, ctau, stau) + theta2 / pi2 * T2
!    
!!    subroutine plob(y0, t0, tf, n, tdir, ftag, ispl, deriv, gr_cj,  y) 
!!    print*, 'pv before plob:', pv; read*; 
!    call plob(pv, 0.d0, tf, ndim, 1, 6, 0, gr_plf, gr_cjplf, pvf)
!    
!   ! check also the energy 
!    call gr_cjplf(pvf, cj)
!    write(111, *) theta1, theta2, pvf, cj 
!    
!  end do
!  
!  write(111, *); write(222,*)
!end do 

end

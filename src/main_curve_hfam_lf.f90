program  curve_hfam_lf 
! continuation along energy, fix rho  for LF spatial problem

! ** NOTE *
!  remember to copy the data files to corresponding subfolders 
!  torus1/       torus2/


use dp_mod
use pi_mod
use curve_mod  
use poinc_mod 
use lf_mod

implicit none

! ** NOTE **
! nf must keep coherent with the one used for curve_lf. 
! read from para_poinc.dat 
! nf = 512 for torus2 

integer, parameter ::  ndim  = 6,   n0=4,  &    ! lf 
                       np    = 301,        &    ! number of points to test the new curve 
                       ntori = 100              ! parameterization for tori  

real(kind=dp) ::  beta0 

integer       ::  i, j, k, niter, isv2                         
real(kind=dp) ::  pv_curv(ndim-2), pv(ndim), pvf(ndim), tf  ! poinc_tf_lf

! Fourier coefficients 
real(kind=dp)              :: ctau0,  c0pre(n0), c0_new(n0)  
real(kind=dp), allocatable, dimension(:,:)  ::  ckpre, skpre, ck_new, sk_new  
real(kind=dp), allocatable, dimension(:)    ::  ctau, stau, ct, st
   
!  continuation 
integer       ::  isref, restart, nitmax, opt, nfmn, nf, nf_pre, nf_init
real(kind=dp) ::  cj_pre, dcj, dcj_max, dcj_min, nds,   &      ! step control for continuaiton 
                  tol, tol_err 

!-- poinc_mod 
integer       :: ind0, dir0, imax0,  ispc, ndim0  ! poinc 
real(kind=dp) :: p00, cj0, xmax0, tmax0  

real(kind=dp) ::  rho, dxi, xi !, pv_refn(np, ndim-2)  ! for rho 

! for return time T and small time deviation \tau 
real(kind=dp) ::  tpc, tau,  pvpc(ndim), tpc_refn(np), cj, xpc_all(np, ndim), hminim  
   
! for the torus  
real(kind=dp) ::  w1, w2, t2, theta1, theta2 
  
                  
real(kind=dp) ::  four_seri !, dlange ! function type declaration
external ::   PoincMap_lf   

call init_lf
beta0 = 2.d0
call init_beta(beta0)


open(88, file = './dat/curve_lf/para_poinc.dat', status='old')
read(88,*)
!  #  n, ind, p0, dir, cj0,  imax, tmax, ndim, nitmax, tol, tol_err
! take the same ind and p0 and nf as the previous refined curve. 

read(88, *)  n, ind0, p00, dir0, cj0, imax0, tmax0, ndim0, nitmax, tol, tol_err 
print*, 'Check poinc: n, ind0, p00, dir0, cj0, imax0, tmax0, ndim0, nitmax, tol, tol_err '
write(*, *)  n, ind0, p00, dir0, cj0, imax0, tmax0, ndim0, nitmax, tol, tol_err 
print*; read*

 cj = cj0 
call init_poinc( ind0, p00, dir0, imax0, tmax0, xmax0, ndim)
call init_poinc_h0(cj) 
 

print*, 'n, ind, ind_fun, cj0 from Poinc:', n, ind, ind_fun, cj0 
print*; read*


! ********* Initialize curve papameters ***************
opt = 3
call init_curve(ndim-2, nitmax, tol, tol_err, opt) 

open(112,  file = './dat/curve_lf/para_curve.dat', status = 'old') 
read(112,*) 

read(112, *)  nf_init, rho, cj0, w1, w2, t2 
print*, 'nf_init, rho, cj0 from para_curve.dat: ', nf_init, rho, cj0
print*;  read*

! use the refined coefficients to compute a point with xi \in [0, pi2]
! -- this is read from the coefficients, but we only have x and vx, 
!    vy is computed as a function of (cj0, x, y=p0, vx)

! The first line is comment, data structure: c0 // ck // sk (10e24.14)

call upd_four(0, nf_init)  
nf = nf_init
open(111, file = './dat/curve_lf/fcs_refn.dat', status='old')
read(111,*)
do i = 1, n, 1
  read(111, '(6e24.14)')  c0(i);    read(111, *)
  write(*,  '(10e24.14)')  c0(i) 
    
  read(111, '(10e24.14)')  ck(i, :); read(111, *)
  read(111, '(10e24.14)')  sk(i, :)
  read(111, *);  read(111, *)
end do

   print*, 'read  from the refined curve: n, nf,  size(c0, ck, sk):'
   print*, n, nf, size(c0), size(ck), size(sk); 
   print*;read*



! Update cj  
  cj  = cj0 
  dcj = -1.d-6
  cj  = cj + dcj
  call init_poinc_h0(cj)
  call upd_four(1, nf) 
  print*, 'cj0=', cj 
  call refine_curve_poinc(n, nf, rho, c0, ck, sk, isref, niter, PoincMap_lf, gr_cjlf, deriv_cjlf)
  
! before using the refined fourier coefs, update the value of nf if possible. -- DO not update!!
  call upd_four(1, nf) 

! *********** continuation for the family of tori *************************
! for continuation 1.d-5 is enough, if we want to look at one certain torus,
! do the refinement with higher precision again, say 1.d-10... 

!---**NOTE** Since we don't want to do the refinement agian, put a lower percision 
! 
!tol     = 1.d-8
!tol_err = 1.d-8 

!!tol     = 1.d-10
!!tol_err = 1.d-10 

!call init_curve(n, nitmax, tol, tol_err, opt) 


!---- do the continuation w.r.t. rho
! we keep the energy, and change the value of rho

! --- use this updated c0, ck, sk to compute refined curve ---- 
open(166, file = './dat/curve_lf/curve_cont.dat', status='replace',access='append')
write(166,*)  '# The continuation of the curve family:  xi  x-y-vx-vy   cj'

open(177, file = './dat/curve_lf/fcs_cont.dat', status='replace',access='append')
write(177,*)  '# Coefs of the continued family :   (10e24.14)   c0,   ck, sk'

open(255, file = './dat/curve_lf/para_curve_cont.dat', status='replace', access='append')
write(255, *) '# Continued family:  nf, rho, h0,  w1,  w2,  t2, c0, ck(:, 1), sk(:, 1) ...  ck(:, 4), sk(:, 4)'

open(155, file = './dat/curve_lf/tpc_fun_cont.dat', status='replace')
write(155,*) '# Return time: xi, tpc, tpc_fun(Fourier), tau_fun(Fourier)'

  
! --------- Automatic step control --------------- 
! done -1.d-6, put into rho_famm/ subfolder 
dcj     = -1.d-6 ! -1.d-6   
  
dcj_max = 1.d-3 
dcj_min = 1.d-8 


! the first curve, allocate ckpre, skpre
restart  = 0
nf = nf_init 

nf_pre = nf; cj_pre = cj; 
call alloc_arr_2d(n, nf_pre, ckpre); call alloc_arr_2d(n, nf_pre, skpre)  
skpre  = sk; c0pre = c0;  ckpre = ck; 

do  k = 1, ntori, 1  
  
  nds = 1
  
  ! ------------- Predictor --------------- 
  if(restart == 1 .and. k>1 ) then 
  
    ! start with the previous one and use new dcj 
    cj = cj_pre; nf = nf_pre
    call upd_four(0, nf_pre) ! declare empty array for c0, ck, sk 
    c0 = c0pre; ck = ckpre; sk = skpre   
 
  else 

    ! A linear extrapolation, introduce intermediate var c0_new, ck_new, sk_new
    ! ** NOTE ** Be careful with the dimension of each array  
    if(k > 1) then 
      call alloc_arr_2d(n, nf, ck_new); call alloc_arr_2d(n, nf, sk_new)
      nfmn = min0(nf_pre, nf)
      c0_new  =  2.d0*c0  - c0pre   
      ck_new(:,1:nfmn) =  2.d0*ck(:,1:nfmn)  - ckpre(:,1:nfmn)  
      sk_new(:,1:nfmn) =  2.d0*sk(:,1:nfmn)  - skpre(:,1:nfmn)
      
      if(nf_pre > nfmn) then     ! nf < nf_pre 
        ck_new(:, nfmn+1:nf_pre) = - ckpre(:, nfmn+1:nf_pre)
        sk_new(:, nfmn+1:nf_pre) = - skpre(:, nfmn+1:nf_pre)
      elseif(nf_pre < nfmn) then  ! nf > nf_pre
        ck_new(:, nfmn+1:nf_pre) =  2.d0*ck(:, nfmn+1:nf_pre)
        sk_new(:, nfmn+1:nf_pre) =  2.d0*sk(:, nfmn+1:nf_pre)
      endif 
    endif 
  
    ! save the previous one
    cj_pre = cj; nf_pre = nf
    call alloc_arr_2d(n, nf, ckpre); call alloc_arr_2d(n, nf, skpre)
    skpre = sk;  c0pre = c0;  ckpre = ck;   
  
    ! Update c0, ck, sk as the guess for new curve  for refinement 
    if(k > 1) then 
       call upd_four(0, nf) 
       c0 = c0_new; ck = ck_new; sk = sk_new 
    endif    
  endif 
    
  ! Update cj  
  cj  = cj + dcj
  call init_poinc_h0(cj)
  
  ! -- compute the linear prediction of CK SK for the new curve 
  ! put all the data in the same row for plot, other wise it's not easy to check the dependence on rho of ck,sk
  
    write(*, *) '# Initial guess for the new curve  c0, ck, sk : n, nf = ', n, nf, size(c0), size(ck), size(sk); 
    
    do i = 1, n, 1
      write(*, '(6e24.14)')   c0(i);       ! write(*,*) 
      write(*, '(10e24.14)')  ck(i, 1:10); ! write(*,*) 
      write(*, '(10e24.14)')  sk(i, 1:10) 
     write(*,*);   
    end do
    print*; read*
  
  print*, 'nf before init_nf: ', nf 
!  call init_nf(nf)  
  
  print*, 'read  the refined curve: n, nf,  size(c0, ck, sk):'
  print*, n, nf, size(c0), size(ck), size(sk); 
  print*;read*
   
  call refine_curve_poinc(n, nf, rho, c0, ck, sk, isref, niter, PoincMap_lf, gr_cjlf, deriv_cjlf)
  
  ! 2017-03-17 10:05:14  -- if isref ==0, fail for refinement, get another chance with half the stepsize 
   
  ! -- modify stepsize dcj based on iterations of Newton method 
  ! -- if we need to decrease dcj, use the previous one as the new guess,
  !    and update to a smaller dcj  
  if(niter > 5 .or. isref == 0)  then 
    dcj = dcj / 2 
    nds = .5d0
    restart = 1
    cycle 
  elseif(niter < 3)  then 
    nds = dmin1( dabs(2*dcj), dcj_max) / dabs(dcj)
    dcj = nds*dcj 
    
  endif 
  
  if(dabs(dcj) < dcj_min) then 
    print*, 'Too small stepsize needed! Stop the continuation!'
    print*; read*
    exit 
  endif
  
  restart = 0 
  
  ! *** suceed a new curve ***
  ! before using the refined fourier coefs, update the value of nf if possible.
  call upd_four(1, nf)  

  do i = 1, n, 1
    write(177, '(6e24.14)')   c0(i);     write(177,*) 
    write(177, '(10e24.14)')  ck(i, :);  write(177,*) 
    write(177, '(10e24.14)')  sk(i, :) 
    write(177,*); write(177,*) ! a blank line to seperate the component 
  end do

  dxi = pi2/(np-1)  
  
  ! Evaluate the truncated Fourier series for \hat T(xi)
  do i = 1,  np, 1  
    xi = (i-1) * dxi
    call varphi( xi, n, nf, c0, ck, sk, pv_curv)
    call gamm(pv_curv, n, ind_fun, cj, ind, p0, pv, isv2, cj2v_lf)
    
    call poinc(pv,  ndim, ndim, 1, tpc, pvpc, hminim, ispc, gr_lf, gr_cjlf) 
    tpc_refn(i)   = tpc
    xpc_all(i, :) = pv
    write(166, '(9e24.14)')  xi, pv, cj, tpc 
  enddo 
 
  write(166, *);  write(166, *)
  write(177, *);  write(177, *)
  
  close(166); close(177);   
  
  open(166, file = './dat/curve_lf/curve_cont.dat', access='append', status='old')
  open(177, file = './dat/curve_lf/fcs_cont.dat', access='append',   status='old')
   
  ! ------ compute the Fourier representation of \tau(xi)------------------
  call alloc_arr_1d(nf, ct); call alloc_arr_1d(nf, st)
  call gr_four(tpc_refn, np, nf, t2, ct, st) 


  ! --- save parameters for the family of curve.  
  write(255, '(I5, 9e24.14)', advance='no')  nf, rho, cj, pi2/t2, rho/t2, t2, c0 
  do i = 1, 4, 1
    write(255, '(8e24.14)', advance='no')  ck(:, i), sk(:, i)
  end do 
  write(255,*)
  
  ! ******** compute the Fourier coefficients ctau, stau for the small time deviation  ***********
  !    from the Time-T map method  tau(xi)
  ctau0 = 0.d0
  call alloc_arr_1d(nf, ctau); call alloc_arr_1d(nf, stau)
  call tau_fc( ct, st, rho, nf, ctau, stau)
  
  ! Evaluate the truncated Fourier series for \hat T(xi)
  do j = 1, np, 1  
    xi = (j-1) * dxi
  
    ! the fourier representation of the return time of the approximated curve 
    ! for 1d function, call the funtion four_seri is better    
    tpc = four_seri(xi, nf, t2,     ct,  st)
    tau = four_seri(xi, nf, ctau0,  ctau, stau)
    
    pv = xpc_all(j, :)
    call plob(pv, 0.d0, tau, ndim, ndim, 1, 0, 6, gr_lf, gr_cjlf, pvpc) 
    write(155, '(16E20.10)')  xi, tpc_refn(j), tpc, tau, pvpc, pv  
  enddo 
  
  write(155, *); write(155,*)
  close(255);   close(155) 
  
  
  open(155, file = './dat/curve_lf/tpc_fun_cont.dat', status='old', access='append')
  open(255, file = './dat/curve_lf/para_curve_cont.dat', status='old', access='append')
  
  print*, 'check para_curve_cont.dat'; print*; read*
enddo 


! check the time curve obtained from tau_fun

 
stop 






! 2017-03-18 09:44:08  finish here. 

print*, 'Only one step away from the torus! cheers!'; print*; read*
! TODO:  
!      2- compute the time deviation tau(xi), and integrate the Poincare map curve   
!         for tau(xi) to check if we get the Time-T map 
!      3- check why there is gap in the torus as a function of (theta1, theta2)  
!      4- use infinity norm to decide how many fourier modes to keep, 
!        -4.1 declare big enough array to store the CK, SK, and introduce a global parameter nf according to 
!             the precision we want to achieve, so 1.d-10 
!        -4.2 look at the error control strategy in Josep Maria's papar, Physic D

  
! instead of the parameterization, just integrate one(any) point from the curve 
! for a long enough time 
!dxi =  pi2 / 1
!tf = 1.d5 

!open(33, file = './dat/curve_lf/torus_intg.dat')

!do i = 1, 1,  1

! ! evaluate a point on the curve 
! xi = (i-1)*dxi
! ! the point on the curve \varphi w.r.t. (theta1, theta2)
! call varphi( xi, n, nf, c0, ck, sk, pv_curv) 
! call gamm(pv_curv, n, ind_fun, h0, ind, p0, pv, isv2, cj2v_lf)
! 
! call gr_cjlf(pv, cj)
! print*, 'check the energy: ', cj; read*
! 
! call plob(pv, 0.d0, tf, ndim, 1, 33, 1, gr_lf, gr_cjlf, pvf)
! 
! call gr_cjlf(pvf, cj)
! print*, 'check the energy after integration: ', cj; read*
! 
!end do 


! - the parameterization of the torus, \psi
! \psi(theta1, theta2) = \Phi_tf  ( \varphi( theta1 - theta2/pi2*rho ) )
! where, tf = theta2/pi2*T2 + \tau( theta1 - theta2/pi2*rho )

! but in order to get a full curve, another way is 
! discretisize theta1 and xi in an uniformed grid, compute the corresponding theta2
! since a full period in xi for a full curve 

open(333, file = './dat/curve_lf/torus_refn_full.dat', status='replace')
open(444, file = './dat/curve_lf/curve_for_torus_refn_full.dat', status='replace')

dxi =  pi2 / ntori

do i = 1, ntori+1,  1

  xi = (i-1)*dxi
  
  ! the point on the curve \varphi w.r.t. (theta1, theta2)
  ! but \varphi essentially only depends on xi
  call varphi( xi, n, nf, c0, ck, sk, pv_curv) 
  call gamm(pv_curv, n, ind_fun, h0, ind, p0, pv, isv2, cj2v_lf)
    
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

    call plob(pv, 0.d0, tf, 4, 1, 6, 0, gr_lf, gr_cjlf, pvf)
    
   ! check also the energy 
    call gr_cjlf(pvf, cj)
    write(333, *)  theta1, theta2, pvf,  cj 
    
  end do
  
  write(333, *); write(444,*)
end do 

stop

!open(111, file = './dat/curve_lf/torus_refn.dat')
!open(222, file = './dat/curve_lf/curve_for_torus_refn.dat')

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
!    call gamm(pv_curv, n, ind_fun, h0, ind, p0, pv, isv, cj2v_lf)
!    write(222, *) theta1, theta2,  xi, pv
!    
!    ! integrate using the pv=\varphi(xi) for time tf
!    tf = four_seri( xi, nf, ctau, stau) + theta2 / pi2 * T2
!    
!!    subroutine plob(y0, t0, tf, n, tdir, ftag, ispl, deriv, gr_cj,  y) 
!!    print*, 'pv before plob:', pv; read*; 
!    call plob(pv, 0.d0, tf, ndim, 1, 6, 0, gr_lf, gr_cjlf, pvf)
!    
!   ! check also the energy 
!    call gr_cjlf(pvf, cj)
!    write(111, *) theta1, theta2, pvf, cj 
!    
!  end do
!  
!  write(111, *); write(222,*)
!end do 

end

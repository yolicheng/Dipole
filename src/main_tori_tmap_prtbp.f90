! --- read the refined fourier coefficients from the file fcs_refn.dat 

program  tori_prtbp 

use dp_mod
use pi_mod
use emconst_mod
use rtbpconst_mod
use curve_poinc_mod  !, only : varphi
use poinc_mod 
implicit none

integer, parameter ::  ndim0 = 4,   &  ! prtbp 
!integer, parameter ::  n = 2, &
                       nf    = 256 , &          ! number of Fourier modes
                       np    = 2000, &          ! number of points to test the new curve 
                       ntori = 100   ! parameterization for tori  

integer :: dircurve , i, j                         
real(kind=dp) ::  pv_curv(2), pv(ndim0), pvf(ndim0), tf, & ! poinc_tf_prtbp
                  c0(2), ck(2, nf), sk(2, nf) 
 
! rotated coefficients 
real(kind=dp) ::  pv_curv2(2), ck2(2, nf), sk2(2, nf), tpc2, pvpc2(ndim0), pv_rho2(2), pv_rho(2)
                                    
!-- poinc_mod 
integer       :: dir0, imax0, ind0 
real(kind=dp) :: xl4, yl4, p00,  cj0, xmax, tmax0 

real(kind=dp) ::  rho, dxi, xi, rho_refn, pv_refn(np, 2)  ! for rho 

! for return time T and small time deviation \tau 
real(kind=dp) ::  tpc, pvpc(ndim0), tpc_refn(np), cj   
   
! for the torus  
real(kind=dp) ::  T2, theta1, theta2  
  
real(kind=dp), dimension(0:nf)  ::  csf,  sif, &    ! fourier coefs of 1 component
                                    ctau, stau    ! fourier coefs for time deviation tau
                  
                  
real(kind=dp) ::  four_seri ! function type declaration
external :: gr_prtbp, gr_cjprtbp, cj2v_prtbp  

! ---- the final step, goes to the  torus!  
!call system( 'gnuplot /home/yu/Dropbox/dipole/dat/curve_prtbp/curve_refn.pl' )
!print*, 'Input the sense of the curve: 1: counterclockwise, -1: clockwise'; print* 
!read*, dircurve

dircurve = -1 ! TODO: keep coherent with the same curve 

call init_emconst
call init_rtbpconst(emrat, ndim0)

! take a point that is close the L4 point,  with the big primary at (mu, 0, 0), 
! the triangular libration point L4 form a equilateral triangle with the two primaries. 
! so the location of L4 is ( -1/2+mu, sqrt(3)/2  ) and  L5 located at ( -1/2+mu, -sqrt(3)/2  )
xl4 = -.5d0 + mu 
yl4 = dsqrt(3.d0) / 2.d0

!outer bound for escape, nondimensinal unit 1 is big enough....
xmax =  2.d0 
tmax0 =  5.d1 ! a rough guess?? or not...

! Poincare section: x=-0.49 by observing the plot of the torus 
ind0  = 1           ! x component 
p00   = xl4   !yl4 !-0.49d0     ! rough guess, by observing directly at the plot of the tori  

dir0  = 1           ! take the positive velocity 
imax0 = 1           ! consider every intersection, and choose by the direction of velocity

!this is a test, to fix the energy level, and x0 = xl4,  and change the value of y, to have a look of a lot tori. 
 cj0 = 3.0000046704015602  
 
! 3.0000047683715843
 
 
call init_poinc( ind0, p00, dir0, imax0, tmax0, ndim0, cj0) 

!read a point from curve_refn, use as initial point and integrate, we will obtain a torus 
open(88, file = './dat/curve_prtbp/curve_refn.dat')
read(88, *)  xi, pv 

! check the energy
call gr_cjprtbp(pv, cj)
print*, 'check the energy, cj =', cj ; print*; read*


! use the refined coefficients to compute a point with xi \in [0, pi2]
! -- this is read from the coefficients, but we only have y and vy, we need to compute vx --- 
! n, ind_fun, cj0, ind, p0 -- all these parameters need to be initialized 

open(177, file = './dat/curve_prtbp/fcs_refn.dat')
do i = 1, n, 1
  read(177, '(10e24.14)')  c0(i);    read(177, *)
  read(177, '(10e24.14)')  ck(i, :); read(177, *)
  read(177, '(10e24.14)')  sk(i, :)
  read(177, *); read(177, *)
  
!   -- check the read -- ! -- ckd
  write(*, '(10e24.14)')  c0(i) 
  write(*, '(10e24.14)')  ck(i, :)
  write(*, '(10e24.14)')  sk(i, :) 
  print*;  read*
end do

! check the rotation number 
rho = 1.1140282949439548d0 ! the original one 

!Original rho:    1.1140282949439548     
!rho of the refined curve:    1.1140158607312580   
  
! ---- compute the rotation number of the refined curve, check if it is 
!      the original fixed one...
! save the new evaluated curve, check if it overlap with curve_refn
open(12, file = './dat/curve_prtbp/curve_new.dat', status='replace')

!dxi = rho / 100
dxi =  rho 
 
print*, 'c0 = ', c0 ; print*; read* 

do i = 1,  np , 1  
  xi = (i-1) * dxi
  call varphi( xi, n, nf, c0, ck, sk, pv_curv)
  write(12, *) xi,  pv_curv 
  pv_refn(i, :) = pv_curv
enddo 
 close(12)
 
call rotnum( pv_refn, np,  dircurve, rho_refn)
rho_refn = rho_refn * pi2 

print*, 'Original rho: ', rho 
print*, 'rho of the refined curve: ', rho_refn 
print*; read*
 

! ------ compute the Fourier representation of \tau(xi)------------------
! 
! Use the refined coefficients, Do a poincare map for the Fourier analysis for the return time 
! generate the points to do Fourier analysis for \tau(xi)

dxi  = pi2 / np   
do i = 1, np, 1
  xi = (i-1) * dxi 
  call varphi( xi, n, nf, c0, ck, sk, pv_curv)
  call varphi( xi+rho, n, nf, c0, ck, sk, pv_rho)
  
  ! subroutine poinc_tf_prtbp(pvin, tf, pvf)
  call poinc_tf_prtbp(pv_curv, tpc, pvpc)
  tpc_refn(i) =  tpc  
end do  

! Foureier analysis for the return time T

call gr_four(tpc_refn, np, nf, csf, sif) 

! -- T2  is the constant term 
T2 = csf(0)

print*, 'We have rho =', rho, 'T2 = ', T2; 

! we need to save the energy, rho and T2 to the file... 

print*, 'Only one step away from the torus! cheers!'; print*; read*
open(255, file = './dat/curve_prtbp/para_curve_poinc.dat', status='replace')
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

open(155, file = './dat/curve_prtbp/tpc_fun.dat', status='replace')
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

!open(33, file = './dat/curve_prtbp/torus_intg.dat')

!do i = 1, 1,  1

! ! evaluate a point on the curve 
! xi = (i-1)*dxi
! ! the point on the curve \varphi w.r.t. (theta1, theta2)
! call varphi( xi, n, nf, c0, ck, sk, pv_curv) 
! call gamm(pv_curv, n, ind_fun, h0, ind, p0, pv, cj2v_prtbp)
! 
! call gr_cjprtbp(pv, cj)
! print*, 'check the energy: ', cj; read*
! 
! call plob(pv, 0.d0, tf, ndim, 1, 33, 1, gr_prtbp, gr_cjprtbp, pvf)
! 
! call gr_cjprtbp(pvf, cj)
! print*, 'check the energy after integration: ', cj; read*
! 
!end do 


! - the parameterization of the torus, \psi
! \psi(theta1, theta2) = \Phi_tf  ( \varphi( theta1 - theta2/pi2*rho ) )
! where, tf = theta2/pi2*T2 + \tau( theta1 - theta2/pi2*rho )

! but in order to get a full curve, another way is 
! discretisize theta1 and xi in an uniformed grid, compute the corresponding theta2
! since a full period in xi for a full curve 

open(333, file = './dat/curve_prtbp/torus_refn_full.dat', status='replace')
open(444, file = './dat/curve_prtbp/curve_for_torus_refn_full.dat', status='replace')

dxi =  pi2 / ntori

do i = 1, ntori+1,  1

  xi = (i-1)*dxi
  
  ! the point on the curve \varphi w.r.t. (theta1, theta2)
  ! but \varphi essentially only depends on xi
  call varphi( xi, n, nf, c0, ck, sk, pv_curv) 
  call gamm(pv_curv, n, ind_fun, h0, ind, p0, pv, cj2v_prtbp)
    
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

    call plob(pv, 0.d0, tf, ndim, 1, 6, 0, gr_prtbp, gr_cjprtbp, pvf)
    
   ! check also the energy 
    call gr_cjprtbp(pvf, cj)
    write(333, *)  theta1, theta2, pvf,  cj 
    
  end do
  
  write(333, *); write(444,*)
end do 

stop

!open(111, file = './dat/curve_prtbp/torus_refn.dat')
!open(222, file = './dat/curve_prtbp/curve_for_torus_refn.dat')

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
!    call gamm(pv_curv, n, ind_fun, h0, ind, p0, pv, cj2v_prtbp)
!    write(222, *) theta1, theta2,  xi, pv
!    
!    ! integrate using the pv=\varphi(xi) for time tf
!    tf = four_seri( xi, nf, ctau, stau) + theta2 / pi2 * T2
!    
!!    subroutine plob(y0, t0, tf, n, tdir, ftag, ispl, deriv, gr_cj,  y) 
!!    print*, 'pv before plob:', pv; read*; 
!    call plob(pv, 0.d0, tf, ndim, 1, 6, 0, gr_prtbp, gr_cjprtbp, pvf)
!    
!   ! check also the energy 
!    call gr_cjprtbp(pvf, cj)
!    write(111, *) theta1, theta2, pvf, cj 
!    
!  end do
!  
!  write(111, *); write(222,*)
!end do 

end

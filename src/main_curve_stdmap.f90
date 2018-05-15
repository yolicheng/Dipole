program curve_stdmap

!  Study the planar Restricted Three Body Problem, and choose L4, which is of type 
!  center X center, so around L4 we have a plenty 2D tori. 
!    1. start from a point that is close to L4, and integrate the orbit 
!    2. compute the Poincare maps, nmap = 10000, with Poincare section taken as x = p0 (or y= p0)
!    3. compute the rotation number of numerically computed invariant curve 
!    4. use the linear interpolation to obtain equally spaced points w.r.t. the angle between the relative position vector from L4 and the x-axis!
!    5. do a simple Fourier analysis to compute the Fourier coefficients as the initial guess to refine the curve 
!    6. Impose the invariance equations, and specify the energy level, to compute the invariant curve 
!    7. Globalize the curve to an invariant torus 


!  --- data saved in ./dat/curve_stdmap/  subfolder 
!  1. torus.dat         -- the original orbit of the 2D torus 
!  2. curve.dat         -- the invariant curve obtained by Poincare map 
!  3. ob_interpol. dat  -- linear interpolated points for Fourier analysis
!  4. fun.dat           -- approximated orbits by the Fourier representation 
!  
use dp_mod
use pi_mod
use emconst_mod 
use sort_mod
use curve_poinc_mod

implicit none
!integer, parameter  ::  dp = kind(1.d0)
integer, parameter  ::  nmap =  2000 , & !20 !1000  ! number of Poincare maps, to start test, 11 is enough?
                        n  = 2, & ! dimension of standard map 
                        
                        np_interpol = 1000, &  ! in order to use fourier, we need odd number of points 
                        nf   =  2  !64 ! 32  !8 
                        ! how many Fourier modes ???, Alex suggests of power 2, say, 32
                        ! Gerard said 5 is enough ....
                        ! use 2 to check... 

! Local Variables
real(kind=dp) :: x0, p0, pti(n), ptf(n),  dpdx(n, n) , &  !curve by map 
                 ptall(nmap, n), rho,  theta(nmap),  & ! rotnum
                 ptold(nmap, n), ptnew_aux(np_interpol, n),  & ! interpolation
                 darg, arg, ptnew(np_interpol, n), pt1(n),   & ! interpolation 
                 pt_in(np_interpol), fun(n), fun2(n), &  ! Fourier + \varphi  
                 fun_rho(n),  fun_map(n), & ! \varphi(xi+rho)+ P(\varphi(xi))
                 
                 c0(n), ck(n, nf), sk(n, nf), &  ! coefs  
                 tol, tol_err, xi  ! Netwon method   
                 
! at this moment, we do not need vel(ind), it is just for test purpose.... 

! fourier coefs for the compoments skipping pos(ind) and vel(ind)
real(kind=dp), dimension(n, 0:nf)    ::  csfn, sifn 
real(kind=dp), dimension(0:nf)       ::  csf, sif    ! fourier coefs of 1 component
                  
integer ::  i, k, l, ind_sorted(nmap), ipt, &
            fcurve, fcurve_interpol, ffcs, ffun, ffun_rho,  dircurve, &  ! rotnum 
            iscurve, isinterpol, isfourier, & !  control which subroutine to call
            nitmax, isref, opt 
             
character(len=70) :: fncurve, fncurve_interpol, fnfun_rho, fnfc, fnfun   

real(kind=dp) ::  four_seri  !funcs

! which approach we want to take 
print*, 'Choose approach:   1: no constraint;  2: fix s1_1 = 0 by adding d s1_1 = 0;  '
read(*,*) opt 

iscurve    = 1
isinterpol = 1
isfourier  = 1

! (0, 0) is an elliptic point when k=-0.5
! take the initial point to be close to the Equilibrium, 
! a good start of the deviation  is of order 1.d-3
! NOTE: (0.5,0.5) is too much  
x0 = 0.d0 + 1.d-3
p0 = 0.d0 - 2.5d-3

! files to save  invariant curve + interpolated one + approximated one 
fcurve   = 25;  fncurve  =  './dat/curve_stdmap/curve.dat'

if(iscurve == 1) then 
  open(fcurve  ,file = fncurve  , access ='append',status = 'replace' )
  write(fcurve  ,*)  ' # Invariant curve by Standard map:  n      (x p)'
endif 

! ---------- Computation -----------------------------
! -- Initialization for the module  curve_poinc_mod 
nitmax  = 10 
tol     = 1.d-10
tol_err = 1.d-10
call init_curve_poinc(nf, n, nitmax, tol, tol_err, opt) 

pti = (/x0, p0/)
print*, 'before map, pti = ', pti
read*

if(iscurve == 1) then 
   do i = 1, nmap, 1
     write(fcurve, *) i,  pti 
!     write(fcurve, *) i, (/(dmod(pti(k), pi2), k=1,n)/)
     call std_map(pti, ptf, 0, dpdx)
     pti = ptf
   end do
   close(fcurve) ! alll the data is saved already inside poinc_n routine 
endif 

!stop ! -- checked the invariant curve 

! read to Poincare maps from the file  
open(fcurve, file = fncurve) 
read(fcurve, *) ! skip the first line 

! -- read the Poincare maps from the data file 
do i = 1, nmap, 1  
  read(fcurve, *) ipt, pt1 
  ptall(i, :) = pt1
enddo 

  ! ** NOTE ** we have to take a look at how the invariant curve moves, and set the value of dircurve
!            1: counterclockwise;  -1: clockwise
! we have -1 in this case...
! TODO, first plot the curve, and take a look at the direction.... 
! -- check ptnew by gnuplot --- discard at the monent 

!print*, 'sense of the curve: 1: counterclockwise, -1: clockwise'
!print*

!call system( 'gnuplot /home/yu/Dropbox/dipole/dat/curve_stdmap/curve.pl' )

!! dircurve =  -1
!print*; print*, 'sense of the curve: 1: counterclockwise, -1: clockwise'
!read*, dircurve 

dircurve = -1 ! by observation, clockwise


! -- compute the rotation number -- 
! be careful, which two columns to use to compute the rotation number...
! the indice are specified by init_curve_poinc for  in module curve_poinc_mod 
call rotnum( ptall, nmap, dircurve, rho )

print*, 'rho= (unit cycle)', rho  
rho = rho * pi2
print*, 'transfrom rho to radian: ', rho
read* 

! the new estimate angle of all the points, defined as theta(1) = 0, theta(2) = rho, ..., 
! given rho in radian
! theta = mod(rho*ind_sorted, 2pi) 
theta = (/ ( dmod(rho*i, pi2) , i = 0, nmap-1) /) ! implicit-do 
call sort_1d( theta, nmap, 1, ind_sorted)


! update all the points by an increasing order of theta 
ptall = ptall(ind_sorted, :) 

! -- interpolated curve (equally spaced in rho) ----
fcurve_interpol = 26;   fncurve_interpol  =  './dat/curve_stdmap/curve_interpol.dat'

! ----- linear interpolation to get equally spaced points in angle: rho  -----------
! rho in radian is the new estiamte of the curve

if( isinterpol == 1) then 
  open(fcurve_interpol,  file = fncurve_interpol)
  write(fcurve_interpol, *)  ' #  Interpolated curve: arg, (x, p)'
  
  ptold = ptall
  call  interpol( ptold, nmap, n, theta, np_interpol, ptnew_aux)
  
  ! Period = 2*pi in radian w.r.t. the angle of the points 
  darg = pi2 / np_interpol ! take this value to keep coherent with the original points 

  do i = 1, np_interpol, 1
    arg = (i-1) * darg
    pt1 = ptnew_aux(i, : ) ! the sorted equispaced points
    
    write(fcurve_interpol, *)  arg,  pt1  
    ptnew(i, :) = pt1
  end do
  
  close(fcurve_interpol) 
endif 

print*, 'Files names for data read and write'
print*, fncurve,  fncurve_interpol;   print*   !ck

! -- general Fourier analysis for nf Fourier modes, with np_interpol points  
! --  read from fcurve_interpol for the array ptnew --- 
if( isinterpol == 0) then 
  open(fcurve_interpol, file = fncurve_interpol) 
  read(fcurve_interpol,*) ! skip the comment line 
  do i = 1, np_interpol, 1
    read(fcurve_interpol,*) arg, pt1
    ptnew(i,:) = pt1
  end do
endif 

!stop !  -- checked the interpolation 

! -- check ptnew by gnuplot --- discard at the monent 
!call system( 'gnuplot /home/yu/Dropbox/dipole/dat/curve_stdmap/curve_interpol.pl' )

! ---------------- Fourier analysis  ---------------------------
! Fourier coefficients + approximated curve, skip  the  component pos(ind) and vel(ind)

! the coefficients computed by general Fourier analysis:  gr_four.f90 
ffcs = 111;  fnfc = './dat/curve_stdmap/fcs.dat'
open(ffcs, file = fnfc) 
write(ffcs,*) '# For columns the components; for each line: Fourier Coef ck(i), sk(i), i=0,...,', nf 


do k = 1, n, 1 ! here we also do the Fourier analysis for vx, just for check purpose 
  
  pt_in =  ptnew(:, k) ! one component, the k-th column  
  call gr_foun( pt_in, np_interpol, nf, csf, sif)
  
  ! check gr_foun and gr_four(my routine) !--ckd, the same!
!  print*, 'csf and sif by gr_foun' !ckd
!  write(*,'(10f20.14)') csf;   write(*,'(10f20.14)') sif;   read* 
  
!  subroutine gr_four(f, n, m, csf, sif)
  call gr_four(pt_in, np_interpol, nf, csf, sif)
 
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
ffun = 222;  fnfun = './dat/curve_stdmap/fun.dat'
open(ffun,  file = fnfun) 
write(ffun, *) ' # Approximated Fourier representation: rho    (x, p)' 

! can check if \varphi(xi+rho) and  \varphi(xi) are the same curves in front of rho(radian)
ffun_rho = 25;  fnfun_rho  =  './dat/curve_stdmap/fun_rho.dat' 
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
  call varphi(arg, c0, ck, sk,     fun)
  call varphi(arg+rho, c0, ck, sk, fun_rho )  
  
!  do k = 1, n ! only the components that we do Fourier analysis(y, vy)   
!!   function four_seri( theta, m, cof, sif)  result(fun)
!    csf = csfn(k, :) 
!    sif = sifn(k, :)
!    fun(k)     = four_seri(arg,  nf,  csf,  sif) ! for one component 
!    fun_rho(k) = four_seri(arg+rho, nf, csf, sif )
!  end do
  
!  ! check varphi ! --ckd, replace four_seri by varphi, which is more efficient. 
!  call varphi(arg, c0, ck, sk, fun2)
!  print*, 'xi  = ', arg
!  print*, 'check \varphi(xi) = ', fun2 
!  print*, 'check four_seri = ', fun 
!  print*; read* 
  
  ! compare four_seri and varphi(from module curve_poinc_mod) !--ckd ,they are the same 
  write(ffun, *)  arg,  fun, fun_rho 
  
  ! ------  Invariance equation  ----- 
  ! check the error in invariance equation  P(\varphi(xi)) = \phi_t( \varphi(xi) ) 
  call std_map(fun, fun_map,  0, dpdx)

  write(ffun_rho, *)  arg,  fun_rho, fun_map

enddo   
  
 close(ffun)  
! stop  !ckd
 
!  ---------- Refinement of the curve -------------
print*, 'before refine_curve_poinc, rho=', rho
read* 

!subroutine refine_curve_poinc( rho, c0, ck, sk, isref, gr_map)
call refine_curve_poinc(rho, c0, ck, sk, isref, std_map)

! --- use this updated c0, ck, sk to compute refined curve ---- 
open(177, file = './dat/curve_stdmap/fcs_refn.dat')
do i = 1, n, 1
  write(177, '(10e24.14)')  c0(i);      write(177, *)
  write(177, '(10e24.14)')  ck(i, :);   write(177, *)
  write(177, '(10e24.14)')  sk(i, :);   write(177, *)
  write(177,*) ! a blank line to seperate the component 
end do
 close(177)
 
open(166, file = './dat/curve_stdmap/curve_refn.dat')

! check the refined Fourier series  
do l = 1,  np_interpol, 1  
  xi = (l-1) * darg
  call varphi( xi, c0, ck, sk, pt1)
  write(166, '(10e24.14)') xi,  pt1
enddo 

stop 

  
end 


 




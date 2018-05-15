program tori_fam_plf
! the tori family from curve_cont.dat along the rotation number rho

! 2017-03-19 15:37:16  --- discard the refinement 

! the original integration works fine, index i in tori_cont.dat corresponds to index i+1 in curve_cont.dat 


! We have to refine very curve, to a smaller 
! read a point from file curve_cont.dat, and integrate for enough long time 

! to obtain a global torus 


use dp_mod
use pi_mod
use plf_mod 
use poinc_mod
use curve_mod 

implicit none

! Variables
integer, parameter    ::   nf = 512, np = 2000, ndim0 = 4 
integer ::  np_start, dircurve, i,  l, niter, isv2, isref

!-- curve 
integer       ::  nitmax, opt
real(kind=dp) :: x(6), pv(4), pvf(4), cj,  T2, & 
                 tf, beta0, pv_curv(2),        &
                 c0(2), ck(2, nf), sk(2, nf),  & !Foureier coef 
                 darg,  tol, tol_err   ! error control
              
                 
! -- Poinc 
integer       :: dir0, imax0, ind0 
real(kind=dp) :: p00, cj0, xmax0, tmax0 

real(kind=dp) ::  rho, xi, rho_refn, pv_refn(np, 2), pv4_refn(np, 4)  ! for rho 



character(len=70) ::  fncurve, fncurve_refn, fntori, fncurve_init, fncurve_para, fnfcs, fncj

external :: PoincMap_plf


call init_plf

beta0 = 2.d0
call init_beta(beta0)


! Pick one point among the np approximated points on the curve 
np_start = 166 

! -- torus1
tf = 50.d0  

!tf = 10.d2   ! torus2    the first 12 Poincare maps makes a closed loop

dircurve = -1 
 
 
ind0  = 2;  
if(ind0 == 1)  p00 = 0.315d0 
if(ind0 == 2)  p00 = 0.d0 

xmax0 =  2.d0; tmax0 =  5.d1  


dir0  = 1;  imax0 = 1           

 cj0 = 4.3767487109222252d0    !h3+0.05
 
call init_poinc( ind0, p00, dir0, imax0, tmax0, xmax0, ndim0, cj0) 

! check the rotation number in radian 
rho = 2.43609973300041d-2    ! torus1
!rho = 0.56911517174218806    ! torus2

opt     = 3
nitmax  = 15 

tol     = 1.d-10
tol_err = 1.d-10 

call init_curve(n, nf, nitmax, tol, tol_err, opt) 


! ---- the family of curve ------------ 
print*, '************ Original curve and torus ********************'

fncurve      = './dat/curve_plf/torus1/curve_cont.dat'
fncurve_refn = './dat/curve_plf/torus1/curve_cont_refn.dat'
fntori       = './dat/curve_plf/torus1/tori_cont.dat'
fncurve_init = './dat/curve_plf/torus1/curve_init.dat'
!fncurve_para = './dat/curve_plf/torus1/para_curve_cont.dat'
fncj  = './dat/curve_plf/torus1/cj_cont.dat'
fnfcs        = './dat/curve_plf/torus1/fcs_cont.dat'

!fncurve      = './dat/curve_plf/torus1/hfamp/curve_cont.dat'
!fnfcs        = './dat/curve_plf/torus1/hfamp/fcs_cont.dat'
!!fncurve_para  = './dat/curve_plf/torus1/hfamp/para_curve_cont.dat'
!fncj  = './dat/curve_plf/torus1/hfamp/cj_cont.dat'

!fntori       = './dat/curve_plf/torus1/hfamp/tori_cont2.dat'
!fncurve_init = './dat/curve_plf/torus1/hfamp/curve_init2.dat'
!fncurve_refn = './dat/curve_plf/torus1/hfamp/curve_cont_refn.dat'

!fncurve = './dat/curve_plf/torus1/rho_famm/curve_cont.dat'
!fntori = './dat/curve_plf/torus1/rho_famm/tori_cont.dat'
!fntori = './dat/curve_plf/torus1/rho_famm/curve_init.dat'

open(166, file = fncurve,       status='old') 
open(177, file = fnfcs,         status='old') 
!open(255, file = fncurve_para,  status='old') 
open(255, file = fncj,  status='old') 

open(200, file = fntori,  status='replace') 
open(300, file = fncurve_init,  status='replace') 
open(400, file = fncurve_refn,  status='replace') 
 
read(166, *); read(177,*); read(255,*) ! 1 comment line, 6 columns
write(200,*) '# Tori family.  t  (x,y,vx,vy)   cj';   close(200)
write(300,*) '# I.C. for curves:  (x,y,vx,vy)   cj';  close(300)
write(400,*) '# Refined curve family.  xi  (x,y,vx,vy)   cj';  close(400)

! do a loop for all the curves, read a point from curve_cont.dat ---
! together with fcs_cont.dat, and do the refinement again. 

 
do  
  
  ! the energy level, can get from cj_cont.dat 
!  read(255, *)  cj0, rho, T2   !nf0 
!  print*, 'cj0, rho, T2', cj0, rho, T2

  read(255, *)  cj0  !nf0 
  print*, 'cj0', cj0 
  
  
  call init_h0(cj0)
  
  ! read the Fourier coefficients, c0, ck, sk, and refine the curve   
  do i = 1, 2, 1
    read(177, '(10e24.14)')  c0(i);     read(177,*) 
    read(177, '(10e24.14)')  ck(i, :);  read(177,*) 
    read(177, '(10e24.14)')  sk(i, :) 
    read(177,*); read(177,*) ! a blank line to seperate the component 
  end do
  
  call refine_curve_poinc(2, nf, rho, c0, ck, sk, isref, niter, PoincMap_plf)
  
  ! The approximated curve 
  open(400,file=fncurve_refn, status='old')
  darg = pi2 / np 
  do l = 1, np, 1  
    xi = (l-1) * darg
    call varphi( xi, n, nf, c0, ck, sk, pv_curv)
    pv_refn(l, :) = pv_curv
    
    call gamm(pv_curv, n, ind_fun, cj0, ind, p0, pv, isv2, cj2v_plf)
    call gr_cjplf(pv, cj0)
    pv4_refn(l,:) = pv 
    write(400, '(10e24.14)')  xi,  pv,  cj0
    print*, xi, pv, cj0
!    print*; read*
  enddo 
  
  write(400, *);  write(400,*); close(400)
  
 
  print*,'Refine a curve! compare with the old one!'
  print*; read* 
  
  
  ! check the rotation number 
  call rotnum( pv_refn, np,  dircurve, rho_refn)
  rho_refn = rho_refn * pi2 

  print*, 'Original rho: ', rho 
  print*, 'rho of the refined curve: ', rho_refn 
  print*; read*

  do i = 1, np_start - 1   
    read(166, *, end = 99) 
  enddo    
   
  ! read for one point on one curve  
  read(166, *)  x 
  pv = pv4_refn(np_start, :)
  print*, 'Check pv: the continued curve, refined curve'
  print*, x(2:5); print*, pv;  print*; read* 
 
  call gr_cjplf(pv, cj)
  
  open(300,file=fncurve_init,status='old' )
  write(300, *)  rho_refn, pv, cj
  print*, 'pv0, cj0', rho_refn, pv, cj; 
  close(300); print*; read*
  
  
  ! Integrate from a point on the refined curve 
  open(200, file = fntori, status='old', access='append') 
  call plob(pv, 0.d0, tf, 4, 4, 1, 1, 200, gr_plf, gr_cjplf, pvf) 
  print*, 'pv, pvf'
  print*, pv; print*, pvf 
  print*, 'Finish one global torus!'; print*; read*
  close(200) 
  
  ! for the next curve 
  do i = 1, np+2    
    read(166, *, end = 99)  
  enddo   

enddo 

99 stop 
end program tori_fam_plf

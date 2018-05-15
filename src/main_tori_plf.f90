program tori_plf
! read a point from file curve.dat, and integrate for enough long time 
! to obtain a global torus 


use dp_mod
use pi_mod
use plf_mod 

implicit none

! Variables
integer :: i, np
real(kind=dp) :: x(5), pv(4), pvf(4), cj,    & 
                 x2(7), pv2(4), pvf2(4), cj2, &
                 tf, beta0 

call init_plf

beta0 = 2.d0
call init_beta(beta0)
np = 10 

!tf = 50.d0  ! torus1
tf = 10.d2   ! torus2    the first 12 Poincare maps makes a closed loop

! ---- the original curve and torus 
print*, '************ Original curve and torus ********************'
open(100, file = './dat/curve_plf/curve.dat', status='old')  ! no comment line, 7 columns
read(100, *)
do i = 1, np 
 read(100, *)  x2 
! print*, x2 
enddo 
pv2 = x2(2:5)
call gr_cjplf(pv2, cj2)
print*, 'pv0, cj0', pv2, cj2; print*; read*

open(200, file = './dat/curve_plf/torus_test.dat', status='replace') 
call plob(pv2, 0.d0, tf, 4, 4, 1, 1, 200, gr_plf, gr_cjplf, pvf2) 
print*, 'pv2, pvf2'
print*, pv2; print*, pvf2 
print*, 'Finish the global torus!'



! ---- the refined curve and torus
print*, '************ Refined curve and torus ********************' 
open(166, file = './dat/curve_plf/curve_refn.dat', status='old')  ! no comment line, 5 columns
read(166, *)
do i = 1, np*10
 read(166, *)  x 
! print*, x 
enddo 
pv = x(2:5)
call gr_cjplf(pv, cj)
print*, 'pv0, cj0', pv, cj; print*; read*


open(188, file = './dat/curve_plf/torus_refn.dat', status='replace') 
call plob(pv, 0.d0, tf, 4, 4, 1, 1, 188, gr_plf, gr_cjplf, pvf) 
print*, 'pv, pvf'
print*, pv; print*, pvf 
print*, 'Finish the global torus!'

stop 

end program tori_plf

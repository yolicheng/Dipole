program tori_fam_plf
! the tori family from curve_cont.dat along the rotation number rho

! We have to refine very curve, to a smaller 
! read a point from file curve_cont.dat, and integrate for enough long time 

! to obtain a global torus 


use dp_mod
use pi_mod
use plf_mod 

implicit none

! Variables
integer :: i, np,  np_start
real(kind=dp) :: x(6), pv(4), pvf(4), cj,     & 
                 tf, beta0 

character(len=70) ::  fncurve, fntori, fncurve_init

call init_plf

beta0 = 2.d0
call init_beta(beta0)


! Pick one point among the np approximated points on the curve 
np_start = 100 

! -- torus1
tf = 50.d0  
np = 2000 

!tf = 10.d2   ! torus2    the first 12 Poincare maps makes a closed loop

! ---- the family of curve ------------ 
print*, '************ Original curve and torus ********************'
 

fncurve      = './dat/curve_plf/torus1/curve_cont.dat'
fntori       = './dat/curve_plf/torus1/tori_cont.dat'
fncurve_init = './dat/curve_plf/torus1/curve_init.dat'


!fncurve      = './dat/curve_plf/torus1/hfamp/curve_cont.dat'
!fntori       = './dat/curve_plf/torus1/hfamp/tori_cont.dat'
!fncurve_init = './dat/curve_plf/torus1/hfamp/curve_init.dat'

!fncurve = './dat/curve_plf/torus1/rho_famm/curve_cont.dat'
!fntori = './dat/curve_plf/torus1/rho_famm/tori_cont.dat'
!fntori = './dat/curve_plf/torus1/rho_famm/curve_init.dat'

open(100, file = fncurve, status='old') 
open(200, file = fntori,  status='replace') 
open(300, file = fncurve_init,  status='replace') 
 
read(100, *)   ! 1 comment line, 6 columns
write(200,*) '# Tori family.  t  (x,y,vx,vy)   cj'; close(200)
write(300,*) '# I.C. for curves:  (x,y,vx,vy)   cj'


! do a loop for all the curves, read a point from curve_cont.dat ---
! together with fcs_cont.dat, and do the refinement again. 

 
do 
  do i = 1, np_start - 1   
    read(100, *, end = 99) 
  enddo    
   
  ! read for one point on one curve  
  read(100, *)  x  
  pv = x(2:5)
  call gr_cjplf(pv, cj)
  write(300, *)  pv, cj
  print*, 'pv0, cj0', pv, cj; print*; read*
  
  open(200, file = fntori, status='old', access='append') 
  call plob(pv, 0.d0, tf, 4, 4, 1, 1, 200, gr_plf, gr_cjplf, pvf) 
  print*, 'pv, pvf'
  print*, pv; print*, pvf 
  print*, 'Finish one global torus!'; print*; read*
  close(200) 
  
  ! for the next curve 
  do i = 1, np+2    
    read(100, *, end = 99)  
  enddo   

enddo 

99 stop 
end program tori_fam_plf

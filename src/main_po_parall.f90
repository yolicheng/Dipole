program main_po_parall 
! We want to get an orbit that is parallel to x-y plane, if it's circular, 
! it will be perfect for inspection in a hover orbit. 
! we deal with the one symmetric w.r.t. y-z plane... 
! this happens only in normal case

use dp_mod
use pi_mod
use lf_mod ! the system related parameters of lorentz force problem
use po_mod ! archived subroutines to compute families of symmetric p.o.
use poinc_mod

implicit none

integer, parameter ::  ndim = 6       ! the dimension of LF problem 
integer, parameter ::  neq  = 42,  &  ! compute also the variational matrix, 6+36=42
                       npo  = 1        

real(kind=dp) ::  beta0, p00, tmax0, xmax0  
 
integer      ::   ind0, dir0, imax0,  &  !poinc
                  dir_curve, idic, idfc, dsctr, anglectr, idcor, isarc, debug 

real(kind=dp) ::  tol, prsc       ! error control

! Local Variables
 
integer :: i,  ispl,   fpo, fpoinst,  fmmegp, fmmegv,  &
           tdir, ipo  

real(kind=dp) ::  poinst(6), ds,  ynew(npo,8), epsl_po, cj,  & !pofam
                  po0(6), tpo, pof(6), & ! plpo
                  dlf(6, 6), mmat(6, 6), wr_mm(6), wi_mm(6), vr_mm(6, 6), &  ! MM 
                  tp0, csangle_min, ds_max, ds_min, & 
                  po1(6), po2(6), po_new(6), tp1, tp2, tp_new, cj1, cj2, cj_new
                  
!                pom3(6), pom2(6), pom1(6), cj3(3), tp3(3), po_new(6) ! muller-discard 
                  
                  
                 
 character(len=70) :: fnpo, fnpoinst, fnmmegv, fnmmegp              

real(kind=dp), external :: dnrm2 ! from package... 

! debug or not 
debug = 0

! ---- for symmetric p.o. : 1-3 ------------ 
idcor = 1 ! for symmetric approach with poincare map  --- most used one...

!idcor = 2 ! symmetric + free time shooting 

! --- asymmetric p.o. 
! idcor = 4 ! Poincare map approach, section y(ind) = x0_fxd, 
! idcor = 5 ! additional constraint to fix one component, y(ind) = x0_fxd, the current used one
! idcor = 6 ! a 
!idcor = 7 ! improved version of idcor=5   

!! the minimum of the cosine of the angle between two consecutive vectors
! cos(0.1 rad) = 0.995;    cos(0.37 rad) = 0.95  
! csangle_min = 0.995d0  
    
call init_debug(debug, dsctr, anglectr, csangle_min,  isarc, ds_min, ds_max)

! Intialize the state vector of the equilibirum points 
call init_lf    

ind0 = 1; idic = ind0;   idfc = ind0


if(ind0 == 0) then 
print*, 'Input ind for Poincare section, X(ind) = p0, ind = ' 
read(*,*) ind0
endif

!print*, 'Please input beta (q/m):'
!read(*,*) beta0 
beta0 = 2.d0
call init_beta(beta0)   
    
print*, 'check, cs, ieq, beta', cs, ieq, beta 
read* 
    
! Jacobi matrix of the lorentz force with respective to the state 
call dflrtz(eq, dlf)
call gr_cjlf(eq, cj)
print*,'check energy!, cj, ieq,', cj, ieq, eq ; read*


! save the data of po and  its initial state 
fpo = 20;  fpoinst = 21;  fmmegv = 23;  fmmegp = 24 

! use new name, to avoid overwritten the good ones, and the type of approach selected is used as suffix
write(fnpo,    fmt='(a,i0,a,i0,a)')   './dat/po/case', cs, '/eq', ieq, '/po_parall.dat'  
write(fnpoinst,fmt='(a,i0,a,i0,a)')   './dat/po/case', cs, '/eq', ieq, '/poinst_parall.dat'
write(fnmmegv, fmt='(a,i0,a,i0,a)')   './dat/po/case', cs, '/eq', ieq, '/pommegv_parall.dat'
write(fnmmegp, fmt='(a,i0,a,i0,a)')   './dat/po/case', cs, '/eq', ieq, '/pommegp_parall.dat' 

print*, fnpo, fnpoinst, fnmmegv, fnmmegp

!open(fpo,file=fnpo, access ='append',status='replace')
!write(fpo,*)  '# t      (x y z vx vy vz)        cj'

!open(fpoinst,file=fnpoinst, access ='append',status='replace')
!write(fpoinst,*) '# TP  I.C.(x y z vx vy vz)    CJ    DCJ   ds'

!open(fmmegv, file=fnmmegv, access  ='append',status='replace') 
!open(fmmegp, file=fnmmegp, access ='append',status='replace')

open(fpo,file=fnpo, access ='append',status='old')
open(fpoinst,file=fnpoinst, access ='append',status='old')
open(fmmegv, file=fnmmegv, access  ='append',status='old') 
open(fmmegp, file=fnmmegp, access ='append',status='old')

! commomly used  varaibles
epsl_po = 1.d-3;   ds = 1.d-2;  tol  = 1.d-10 ;  prsc = 1.d-10

! initial guess by muller method ....  between pom3 and pom2, do the bisection ?  or continuation with very small step size? 

!  0.12930046758208E+01  0  -0.19474831094014E+00  0.10782828815608E+01 -0.10427410906443E+01  0   0  -0.23076621476532E+01  0.44408920985006E-15  0.04    1
!  0.12921416117285E+01  0  -0.23377288391382E+00  0.10702907983806E+01 -0.10392193141463E+01  0   0  -0.23086282504525E+01  0.22204460492503E-14  0.04   -1
!  0.12938136711696E+01  0  -0.27227309009366E+00  0.10611816279315E+01 -0.10333906221225E+01  0   0  -0.23067581394774E+01  0.26645352591004E-14  0.04   -1

!tp3  = (/0.12930046758208d1,   0.12921416117285d1, 0.12938136711696d1/)  ! by period
! cj3 = (/-0.23076621476532d1, -0.23086282504525d1, -0.23067581394774d1/) ! by energy 

tp1 =  0.12930046758208d1; tp2 = 0.12921416117285d1 
 cj1 =  -0.23076621476532d1; cj2 =  -0.23086282504525d1

! 1-the iterate: between pom3 and pom2, taken as po1, and po2 
po1 = (/0.d0, -0.19474831094014d0,  0.10782828815608d1,  -0.10427410906443d1,  0.d0,  0.d0 /)
po2 = (/0.d0, -0.23377288391382d0,  0.10702907983806d1,  -0.10392193141463d1,  0.d0,  0.d0 /)
!pom1 = (/0.d0, -0.27227309009366d0,  0.10611816279315d1,  -0.10333906221225d1,  0.d0,  0.d0 /)


! 2- iterate, between po2 and po_new 
! 0.12922539418621E+01     0.d0    -0.21432438620954E+00  0.10744286042532E+01 -0.10412717081158E+01     0.d0        0.d0    -0.23085024120330E+01 
!po_new  = (/0.d0, -0.21432438620954d0,  0.10744286042532d1,  -0.10412717081158d1, 0.d0, 0.d0/)
!tp_new  = 0.12922539418621d1; cj_new = -0.23085024120330d1
!po1 = po2; tp1 = tp2; cj1 = cj2; 
!po2 = po_new; tp2 = tp_new; cj2 = cj_new 


po1 = (/0.d0, -0.23377288391382d0,  0.10702907983806d1,  -0.10392193141463d1,  0.d0,  0.d0 /)
tp1 = 0.12921416117285d1 
 cj1 =  -0.23086282504525d1
! 

! 3-th iterate, between po1 and po_new, po_new is the new po2 
!  0.12921183253558E+01     0.d0    -0.22406571782006E+00  0.10723946893075E+01 -0.10403178977869E+01     0.d0        0.d0    -0.23086543410482E+01  0.88817841970013E-15   
tp2  = 0.12921183253558d1; cj2 = -0.23086543410482d1 
po2  = (/0.d0, -0.22406571782006d0,  0.10723946893075d1,  -0.10403178977869d1,  0.d0,  0.d0 /)

! 4-th iterate, between po2 and po_new, take po_new as the new po1 
!  0.12921101482912E+01     0.d0    -0.22892297846856E+00  0.10713515896572E+01 -0.10397867185515E+01     0.d0        0.d0    -0.23086635029787E+01  0.13322676295502E-14   
tp1 = 0.12921101482912d1 ; cj1  = -0.23086635029787d1 
po1 = (/0.d0, -0.22892297846856d0,  0.10713515896572d1,  -0.10397867185515d1,  0.d0,  0.d0 /)

! 5-th iterate, between po1 and po_new, po_new to be the new po2 
!  0.12921092899630E+01     0.d0    -0.22649534289459E+00  0.10718753423502E+01 -0.10400568382863E+01     0.d0        0.d0    -0.23086644986651E+01  0.35527136788005E-14  
tp2 = 0.12921092899630d1; cj2 = -0.23086644986651d1 
po2 = (/0.d0, -0.22649534289459d0, 0.10718753423502d1, -0.10400568382863d1,  0.d0,  0.d0 /)

! 6-th  iterate, between po2 and po_new, take po_new as the new po1
!  0.12921084707020E+01     0.d0    -0.22770939676731E+00  0.10716140140296E+01 -0.10399229058100E+01     0.d0        0.d0    -0.23086653725440E+01  0.22204460492503E-14 
tp1 =  0.12921084707020d1; cj1 = -0.23086653725440d1 
po1 = (/0.d0, -0.22770939676731d0, 0.10716140140296d1, -0.10399229058100d1,  0.d0,  0.d0 /)

! 7-th iterate, between po1 and po_new,  take po_new as the new po2 
!  0.12921085665096E+01     0.d0    -0.22710243036587E+00  0.10717448142419E+01 -0.10399901530170E+01     0.d0        0.d0    -0.23086652756509E+01 -0.13322676295502E-14  
tp2 =  0.12921085665096d1; cj2 = -0.23086652756509d1 
po2 = (/0.d0,  -0.22710243036587d0,  0.10717448142419d1, -0.10399901530170d1,  0.d0,  0.d0 /)

! 8-th iterate, between po1 and po_new, take po_new as the new po2 
!  0.12921084428626E+01     0.d0    -0.22740592910535E+00  0.10716794491583E+01 -0.10399566008910E+01     0.d0        0.d0    -0.23086654139165E+01  0.35527136788005E-14   
tp2 =  0.12921084428626d1; cj2 =  -0.23086654139165d1 
po2 = (/0.d0, -0.22740592910535d0,  0.10716794491583d1, -0.10399566008910d1,  0.d0,  0.d0 /)

! 9-th iterate, between po2 and po_new, take po_new as the new po1 
!  0.12921084391527E+01     0.d0    -0.22755766697553E+00  0.10716467408477E+01 -0.10399397718243E+01     0.d0        0.d0    -0.23086654180079E+01  0.13322676295502E-14   
tp1 =  0.12921084391527d1; cj1 =   -0.23086654180079d1 
po1 = (/0.d0,-0.22755766697553d0,  0.10716467408477d1, -0.10399397718243d1,  0.d0,  0.d0 /)



write(fpoinst, *) 
write(fpoinst, *) tp1, po1, cj1 
write(fpoinst, *) tp2, po2, cj2 

tdir  = 1; ispl = 1
call plob(po1, 0.d0,  tp1, 6, 6, tdir, ispl, fpo, gr_lf, gr_cjlf, pof) 
call plob(po2, 0.d0,  tp2, 6, 6, tdir, ispl, fpo, gr_lf, gr_cjlf, pof) 

! discard muller, which is to find the minimum... doesn't help a lot... 
!subroutine muller( n, x1, x2, x3, y(3),  x)
!call muller( 6, pom3, pom2, pom1, cj3,  po_new)


! Use bisection, and  do the approximation step by step ... 
po_new = (po1 + po2 ) / 2.d0 
tp_new = (tp1 + tp2) / 2.d0
print*, 'The new poinst = ', tp_new, po_new 
poinst = po_new 
print*; read*

!--------------------------------------------------------------------------------


! plot the initial guessed orbit  to see if it is a periodic orbit
open(333, file = './dat/po_guess.dat', access ='append',status='replace')
call plob(poinst, 0.d0,  tp_new,  6, 6, 1, 1, 333 , gr_lf, gr_cjlf, pof)  
print*, 'check the initial guessed orbit!'; close(333); read*

!  check which initializaiton routine to call  2017-02-19 16:03:36  
! checking idcor == 5... 
if (idcor  == 1  )  then !  for  idcor = 1, 2, 3, TODO 3
! ---- symmetric p.o. ----
  call init_symtpo(idic, idfc, idcor, ind0)    
  p0  = 0.d0;  imax0 = 1 
  xmax0 = 40.d0;  tmax0 = 100.d0 
endif 
 

if (idcor == 1 .or. idcor == 4 ) then  ! TODO 4 to check 
!  subroutine init_poinc( ind0, p00, dir0, imax0, tmax0, xmax0, ndim0, h00) 
  call init_poinc( ind0, p00, dir0, imax0, tmax0, xmax0, ndim, cj)
endif 

call init_errctr(tol, prsc)
call init_writedx 

! ds = 0.01d0
! subroutine pofam(yi, tp0, npo, dir,ds, fpoinst, ynew, ipo, deriv, gr_cj)
call pofam(poinst, tp0, npo, dir_curve,ds, fpoinst, ynew, ipo, gr_lf, gr_cjlf)
 
ipo = max0(ipo-1, 1) 

print*, 'PO finisned!'
print*, 'No of real p.o. computed = ', ipo 
read*

! --- plot the P.O ---
do i = 1, ipo
  
  po0 = ynew(i,2:7)
  tpo = ynew(i,1)

! ynew(i,:) = (/tp, yi, cj/)  from pofam
  print*, i, '-th P.O. TP: ', tpo
  
  write(*,'(8f12.8)') ynew(i,:) 
 
  ispl  = 1
  ! subroutine plob(y0, t0, tf, ndim, nvar, tdir, ispl, ftag, deriv, gr_cj,  y) 
  call plob(po0, 0.d0,  tpo, 6, 6, tdir, ispl, fpo, gr_lf, gr_cjlf, pof) 
  print*, 'Check the I.C. and F.C. of p.o.'
  print*, po0
  print*, pof  
  print*;  read*
  
  ! --- Monodramy matrix
  ! subroutine monomat(yi,tp, mmat, deriv, gr_cj)
  print*, 'refined initial state, tp, ynew'
  print*,  tpo, po0
  call monomat(po0, tpo, mmat, gr_lf, gr_cjlf)
  
! print mmat to file, mmat.dat  -- not necessary to save this 
!  do j = 1, n
!    write(*,'(6d20.10)') mmat(j,:)
!  enddo  
!  read*
   
!  write(fmmat, *)  ! add a blank line 
  
! analyze the stability of the monodramy matrix, a big step forward!
!subroutine eigrg(a,n,isv, wr,wi,vr)
  call eigrg(mmat, ndim, 1, wr_mm, wi_mm, vr_mm)
!  
!  print*, 'Eigenvalues and eigenvectors, mmat!!!'
!! it seems the real part of the eigenvectors of the monodramy matrix are nearly 1 
!! so the stability is not so straightforward
!! try the power method to see if we can get the dominant eigenvalue and eigenvector

!!  fmmegv = 6; fmmegp = 6 ! print to screen 
  print*, 'eigenvalues, real part'
  print*,  wr_mm 
  
  print*, 'eigenvalues, imaginary part'
  print*,  wi_mm 
  
  print*
  
  call prt_eigval( ndim, fmmegv, wr_mm, wi_mm )
  call prt_eigvec( ndim, fmmegp, wi_mm, vr_mm )
  
  write(fmmegp,*) ! add a blank line to seperate eigenvector matrix
enddo   

stop
end program main_po_parall



















  

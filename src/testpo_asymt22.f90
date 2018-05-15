program testpo_asymt ! for eq3 
! 20160406 
! test dfcr_asymt use the same data as eq3_bt10_ifam3, to test asymmetric approach, which do the refinement after 1 revolution w.r.t the initial condition

 
use lf_mod ! the system related parameters of lorentz force problem
use po_mod ! archived subroutines to compute families of symmetric p.o.

implicit none

integer, parameter ::  neq = 42, &  ! compute also the variational matrix, 6+36=42
		       npo = 50 ! 1 !35 ! NO of orbits in pofam 

real(kind=dp) :: pi = 4.d0*datan(1.d0)   
    
! the dimension of the problem is defined in lf_mod, n=6 

! lf_mod Module Variables
! Global :   beta, cs, sgn, eq

integer ::  cs, ieq, symt, asymt, ind, debug, dsctr
real(kind=dp) ::  beta, sec, hmax, hmin, e, tmax, tol, prsc ! error control
! Local Variables


integer :: i, j,  ivr,  dir, imax, fpo, fpoinst, ifam, fmmegp, fmmegv,  &
	   tdir, ipo

real(kind=dp) ::  dlf3(n,n),  & ! differential of vector field of lorentz force
                  wr(n),wi(n), vr(n,n), y0(neq), &   ! eigenspace of variational matrix  
                  poinst(6), ds, vrchs(6), ynew(npo,8), epsl_po, cj,  & !pofam
                  po0(6), tpo, pof(6), & ! plpo
                  dlf(n,n), mmat(n,n), wr_mm(n), wi_mm(n), vr_mm(n,n), &  ! MM 
                  yi(42), tf, yf(42), hminim, phi(6,6), yfck(42), phi1(6), & ! check phi 
                  tp0, csangle_min  

                  
integer ::  ispc, issymt, & !check phi
            ispo, niter 
            
                 
 character(len=70) :: fnpo, fnpoinst, fnmmegv, fnmmegp              

real(kind=dp), external :: dnrm2 ! from package... 


! debug or not 
debug =  1
dsctr =  0

!issymt = 1 ! for symmetric approach with poincare map 
!issymt = 0 ! for asymmetric approach with poincare map 
issymt = 2  ! the general approach 


! the minimum of the cosine of the angle between two consecutive vectors
! cos(0.1 rad) = 0.995;    cos(0.37 rad) = 0.95  
 csangle_min = 0.8d0  
 
call init_debug(debug, dsctr, csangle_min)

! Initialize the private variables eq1, sgn1 for  case 1
call init_lf_mod

! case 1: N=[1 0 0]; 2:  N=[0 1 0]; 3:  N=[0 0 1]; 
 cs = 1 ! use 2 to test the swap rule is coded correctly--ckd!

! the index of the equilibrium points
! 1:  q/m > 0,  x=0,y=0,z= \pm 1
! 2:  q/m < 0,  x = \pm (2 sqrt(9))^(1/3), y² = 2x²,z = 0
! 3:  q/m < 0   x = \pm (1/4/sqrt(2)^(1/3), y = 0,z² = x²  ! this the the case we study now 
ieq = 2 ! 3, x,0,z is the case that we study currently


! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy
! and  the index of column of vr to be used as the solution of the variational equatio, carefully checked by observing the eigenvectors

! for bt = 10, 3 families of p.o. (/1, 4, 6/) 
beta = 1.d0  

! --- 1st family  --- this one is hard to
!ivr = 1; ifam = 1; epsl_po = 1.d-4; ds = 1.d-3;  tol  = 1.d-10 ; prsc = 1.d-10


! --- 2nd family
! A little trick here, with greater values of error control tol = 1.d-9, we can continue further
!ivr = 4; ifam = 2; epsl_po = 1.d-6; ds = 1.d-5;  tol  = 1.d-9 ; prsc = 1.d-9 ! 1.d-10 we will only have 29 orbits, 1.d-9==>ipo=250 

! --- 3rd family
! this family is hard to continue, with a long period, tp = 44.5, ipo = 1
ivr = 6; ifam = 3; epsl_po = 1.d-6; ds = 1.d-3;  
!tol  = 1.d-9 ; prsc = 1.d-9 ! 1.d-9->ipo=163; 1.d-10->ipo=35

tol  = 1.d-10 ; prsc = 1.d-10

! for ifam = 1, ds = 1.d-3, npo = 40, espl_po = 1.d-5,  the continued family goes into a cylinder, 
!               also for ds = 5.d-4, npo = 80, the same 
!               test ds = 1.d-4 -- there isn't too much difference.

 

! Provided the case, beta and ieq,  beta, cs, sgn and eq are initialized by subroutine init_lf from module lf_mod
! subroutine init_lf(beta0, cs0, ieq)
call init_lf(beta, cs, ieq) 

print*, 'check the assignment with module' !--ckd
print*, 'beta, cs, sgn, eq', beta, cs, sgn, eq
read*

!  ***************** compute the 3 families of PO ***************************
! save the data of po and  its initial state 
fpo = 20;  fpoinst = 21;  fmmegv = 23; fmmegp = 24 

! remember to rename the data file for different families of po when beta = 10 
!write(fnpo,    fmt='(a,i0,a,i0,a)') './dat/eq3_bt', idint(beta), '_po', ifam, '.dat'  
!write(fnpoinst,fmt='(a,i0,a,i0,a)') './dat/eq3_bt', idint(beta), '_poinst', ifam, '.dat'
!write(fnmmegv, fmt='(a,i0,a,i0,a)') './dat/eq3_bt', idint(beta), '_mmegv', ifam, '.dat'
!write(fnmmegp, fmt='(a,i0,a,i0,a)') './dat/eq3_bt', idint(beta), '_mmegp', ifam, '.dat'

! use new name, to avoid overwritten the good ones
write(fnpo,    fmt='(a)') './dat/apo.dat'  
write(fnpoinst,fmt='(a)') './dat/apoinst.dat'
write(fnmmegv, fmt='(a)') './dat/apommegv.dat'
write(fnmmegp, fmt='(a)') './dat/apommegp.dat'

!fnpo = "'./dat/testpo.dat'"
!fnpoinst = "'.dat/testpoinst.dat'"
!fnmmegv = "'.dat/testmmegv.dat'"
!fnmmegp = "'.dat/testmmegp.dat'"

print*, fnpo, fnpoinst, fnmmegv, fnmmegp
read(*,*)

open(fpo,file=fnpo, access ='append',status='replace')
write(fpo,*)  '# t	(x y z vx vy vz)	cj'

open(fpoinst,file=fnpoinst, access ='append',status='replace')
write(fpoinst,*) '# TP	I.C.(x y z vx vy vz)	CJ	DCJ'

open(fmmegv, file=fnmmegv, access  ='append',status='replace') 
open(fmmegp, file=fnmmegp, access ='append',status='replace')
 
 
! Jacobi matrix of the lorentz force with respective to the state 
!subroutine dflrtz( x0, dlf)
call dflrtz(eq, dlf)

! check the energy 
!C1 = -3   
!C2 = 2.2894  
!C3 = 1.8899 

call gr_cjlf(eq, cj)
print*,'check energy!, cj, ieq,', cj, ieq, eq
!read(*,*)  

do i = 1, n
 write(*,'(6f8.4)') dlf(i,:) 
enddo

! compute the eigenvalues and eigenvectors  of dlf
!subroutine eigrg(a,n,isv, wr,wi,vr)
call eigrg(dlf,n,1, wr,wi, vr)
 

! Initialization of parameters to do numercial continuation
dir = 1
tdir = 1 ! integrate forward

! finally take 1.d-9 as the error control
!tol  = 1.d-10 ! err tolerance of target f to terminate Newton Method, from the Gerard's book, it should be 1.d-16?
!prsc = 1.d-10 ! For lf problem,  1.d-11 is too small for the second family, ok for the other 2 families

!tol  = 1.d-9    
!prsc = 1.d-9  

! Use the period as initial guess, TP = 2*pi/wi, wi is the pure imaginary eigenvalue
tp0 =  2*pi/dabs(wi(2*ifam)) ! the maximum time for the first return to poincare section

print*, 'tp0=', tp0, 'wi=', wi(2*ifam)
read*

print*, 'check', ivr,'-th column of vr to use', vr(:, ivr) 
vrchs = vr(:, ivr)
vrchs = vrchs/dnrm2(6,vrchs,1)
print*, dnrm2(6,vrchs, 1), vrchs

read* 
 
!poinst =  eq +  epsl_po * vr(:, ivr) ! the initial guess for po, move along the corresponding eigenvector for a small distance epsl_po 
poinst =  eq +  epsl_po * vrchs !  


if (issymt == 1) then 

!! --- Approach 1, by symmetric way 
  symt = 2 ! the 2nd symmetry: y=0 plane
  ind  = 2 ! y
  sec  = 0.d0 
  imax = 2
  call init_symtpo(symt)   ! forsymmetric p.o.

elseif(issymt == 0 ) then

! --- Approach 0, by general asymmetric way
  asymt = 1  ! still try y = 0 plane            !alternative: x=x0 plane
  ind   = 1  ! y
  sec   = poinst(ind) 
  imax  = 2
  
!  subroutine init_asympo(asymt ) ! for asymmetric p.o.
  call init_asymtpo(asymt)   ! for symmetric p.o.

elseif (issymt == 2) then 

! --- Approach 2, by general asymmetric way
  call init_grpo ! no input
endif 


! Step size and error control for rk78 integration 
hmin = 1.d-10
hmax = 1.d0
e    = 1.d-14! 1.d-13 is from Gerard's book 

call init_errctr(tol, prsc)

!  subroutine init_poinc(ind0, sec0, imax0,  hmin0, hmax0, e0)
call init_poinc(ind, sec, imax, hmin, hmax, e)
call init_writedx 


!----- test asymmetric!  Using the 4-th po, with the same initial guess 
!poinst = (/ 0.55183345638931047d0, 0.d0, 0.55128378693257030d0, 0.d0,-2.9512193334123396d-3, 0.d0/)  
 
!  I.C.           0.55989279109528511        0.0000000000000000       0.55981132730543903        0.0000000000000000       -4.3542408633655017E-004   0.0000000000000000     ! the second p.o. 
!poinst = (/ 0.55989279109528511d0, 0.d0, 0.55981132730543903d0, 0.d0, -4.3542408633655017d-4, 0.d0/)  ! 2nd p.o. 



!the good initial guess from po24, and by full asymmetry  
!poinst =  (/ 0.50436288926350648d0, 0.71327525584971330d0, -1.0452846326765423d-7,  &
!             1.000000000000000d-5, -1.000000000000000d-5, -2.1207904341154505d-6/)


!0.47738572264615D+02   
! 0.45329294567522D+00    0.00000000000000D+00    0.45112356361413D+00    0.00000000000000D+00   -0.22150193143966D-01    0.00000000000000D+00    0.19685896046402D+01    0.42188474935756D-14    0.10000000000000D-02

!poinst = (/ 0.45329294567522d0, 0.d0, 0.45112356361413d0, 0.d0, -0.22150193143966D-01, 0.d0/) 

!poinst = (/ 0.45499005098177D+00,    0.00000000000000D+00,    0.45215482243370D+00,    0.00000000000000D+00, & 
!  -0.21963320389781D-01,  0.00000000000000D+00 /)

!poinst = (/ 0.45329294567522d0, 0.d0, 0.45112356361413d0, 0.d0, -0.22150193143966D-01, 0.d0/) 

! test the continuation from the second families in the negative direction
!  0.4797289074490632D+02  0.4524274464352185D+00  0.0000000000000000D+00  0.4573275746506272D+00  0.0000000000000000D+00 -0.2226337907927061D-01  0.0000000000000000D+00  0.1975649365951931D+01  0.1998401444325282D-14  0.1000000000000000D-02
!poinst = (/0.4524274464352185d0, 0.d0,  0.4573275746506272d0,  0.d0, -0.2226337907927061d-1, 0.d0/)
!dir = -dir 
!ds = 2.d-3 ! try bigger step to see if is possible to obtain the second family of p.o.s with smaller size, without jumping to another family

!  subroutine pofam(yi, tp0, npo, dir,ds, fpo, ynew, ipo, deriv, gr_cj)

call pofam(poinst, tp0, npo, dir,ds, fpoinst, ynew, ipo, gr_lf, gr_cjlf)

!!subroutine plob(y0,t0,tf, n, tdir, ftag, deriv, gr_cj,  y) 
!po0 = poinst 
!print*, 'tp0, po0'
!print*, tp0, po0 
!read* 


!tdir = 1
!call plob(po0, 0.d0, 1.d0*tp0, 6, tdir, fpo, gr_lf, gr_cjlf, pof)
! close(fpo)
!read*
!stop 

!ds = 1.d-4

!! Initialize for mmat
!y0 = 0.d0
!y0(1:6) = poinst  
!y0(7:42:7) = 1.d0 !  the initial value of variational matrix if identity stored in y0(7:42)
! 



!!print*, 'before pofam'
!!  subroutine pofam(yi,npo,imax, dir,ds, fpo, ynew, i, deriv, gr_cj)
!call pofam(poinst, npo,imax, dir, ds, fpoinst, ynew, ipo, gr_lf, gr_cjlf)

ipo = max0(ipo-1, 1) 

print*, 'PO finisned!'
print*, 'No of real p.o. computed = ', ipo 
read*
! --- plot the P.O ---

do i = 1, ipo
  
  po0 = ynew(i,2:7)
  tpo = ynew(i,1)

!    ynew(i,:) = (/tp, yi, cj/) from pofam
  print*, i, '-th P.O. TP: ', tpo
  write(*,'(8f12.8)') ynew(i,:) 
  read(*,*) 

!subroutine plob(y0,t0,tf, n, tdir, ftag, ispl, deriv, gr_cj,  y) 
  call plob(po0, 0.d0, 1*tpo, 6, tdir, fpo, 1, gr_lf, gr_cjlf, pof)
    
  ! --- Monodramy matrix
!subroutine monomat(yi,tp, mmat, deriv, gr_cj)
  print*, 'refined initial state, tp, ynew', tpo, po0
  call monomat(po0, tpo, mmat, gr_lf, gr_cjlf)
  
! print mmat to file, mmat.dat  -- not necessary to save this 
!  do j = 1, n
!    write(*,'(6d20.10)') mmat(j,:)
!  enddo  
!  read*
   
!  write(fmmat, *)  ! add a blank line 
  
! analyze the stability of the monodramy matrix, a big step forward!
!subroutine eigrg(a,n,isv, wr,wi,vr)
  call eigrg(mmat,n,1, wr_mm, wi_mm, vr_mm)
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
  
  call prt_eigval( n, fmmegv, wr_mm, wi_mm )
  call prt_eigvec( n, fmmegp, wi_mm, vr_mm )
  
  write(fmmegp,*) ! add a blank line to seperate eigenvector matrix
enddo   


stop
end program testpo_asymt



















  

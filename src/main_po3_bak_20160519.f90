program po3 ! for eq3 
 
use lf_mod ! the system related parameters of lorentz force problem
use po_mod ! archived subroutines to compute families of symmetric p.o.

implicit none

integer, parameter ::  neq = 42, &  ! compute also the variational matrix, 6+36=42
		       npo = 60     !200 !35 ! NO of orbits in pofam 
		       
! the dimension of the problem is defined in lf_mod, n=6 

! lf_mod Module Variables
! Global :   beta, cs, sgn, eq

integer ::        cs, ieq, symt, asymt, ind, debug, dsctr
real(kind=dp) ::  beta, sec, hmax, hmin, e, tmax, tol, prsc ! error control
! Local Variables


integer :: i, j,  ivr,  dir, imax, fpo, fpoinst, ifam, fmmegp, fmmegv,  &
	   tdir, ipo

real(kind=dp) ::  dlf3(n,n),  & ! differential of vector field of lorentz force
                  wr(n),wi(n), vr(n,n), y0(neq), &   ! eigenspace of variational matrix  
                  poinst(6), ds, vrchs(6), ynew(npo,8), epsl_po, cj,  & !pofam
                  po0(6), tpo, pof(6), & ! plpo
                  dlf(n,n), mmat(n,n), wr_mm(n), wi_mm(n), vr_mm(n,n), &  ! MM 
                  yi(42), tf, yf(42), hminim, phi(6,6), yfck(42), phi1(6)  ! check phi 
                  
integer ::  ispc !check phi
                  
 character(len=70) :: fnpo, fnpoinst, fnmmegv, fnmmegp              

real(kind=dp), external :: dnrm2 ! from package... 


! debug or not 
debug = 0
dsctr = 0

call init_debug(debug, dsctr)

! Initialize the private variables eq1, sgn1 for  case 1
call init_lf_mod

! case 1: N=[1 0 0]; 2:  N=[0 1 0]; 3:  N=[0 0 1]; 
 cs = 1 ! use 2 to test the swap rule is coded correctly--ckd!

! the index of the equilibrium points
! 1:  q/m > 0,  x=0,y=0,z= \pm 1
! 2:  q/m < 0,  x = \pm (2 sqrt(9))^(1/3), y² = 2x²,z = 0
! 3:  q/m < 0   x = \pm (1/4/sqrt(2)^(1/3), y = 0,z² = x²  ! this the the case we study now 
ieq = 3 ! 3, x,0,z is the case that we study currently


! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy
! and  the index of column of vr to be used as the solution of the variational equatio, carefully checked by observing the eigenvectors

! for bt = 10, 3 families of p.o. (/1, 4, 6/) 
beta = 10.d0  

! --- 1st family  --- this one is hard to
!ivr = 1; ifam = 1; epsl_po = 1.d-4; ds = 1.d-3;  tol  = 1.d-10 ; prsc = 1.d-10


! --- 2nd family
! A little trick here, with greater values of error control tol = 1.d-9, we can continue further
!ivr = 4; ifam = 2; epsl_po = 1.d-6; ds = 1.d-5;  tol  = 1.d-9 ; prsc = 1.d-9 ! 1.d-10 we will only have 29 orbits, 1.d-9==>ipo=250 

! --- 3rd family
! this family is hard to continue, with a long period, tp = 44.5, ipo = 1
ivr = 6; ifam = 3; epsl_po = 1.d-6; ds = 1.d-3;  tol  = 1.d-9 ; prsc = 1.d-9 ! 1.d-9->ipo=163; 1.d-10->ipo=35

tol  = 1.d-11 ; prsc = 1.d-14

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
write(fnpo,    fmt='(a,i0,a,i0,a)') './dat/eq3_bt', idint(beta), '_po', ifam, '.dat'  
write(fnpoinst,fmt='(a,i0,a,i0,a)') './dat/eq3_bt', idint(beta), '_poinst', ifam, '.dat'
write(fnmmegv, fmt='(a,i0,a,i0,a)') './dat/eq3_bt', idint(beta), '_mmegv', ifam, '.dat'
write(fnmmegp, fmt='(a,i0,a,i0,a)') './dat/eq3_bt', idint(beta), '_mmegp', ifam, '.dat'


print*, fnpo, fnpoinst, fnmmegv, fnmmegp
read(*,*)

open(fpo,file=fnpo, access ='append',status='replace')
write(fpo,*)  't	(x y z vx vy vz)	cj'

open(fpoinst,file=fnpoinst, access ='append',status='replace')
write(fpoinst,*) 'TP	I.C.(x y z vx vy vz)	CJ	DCJ	Hminim'

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
read(*,*)  

do i = 1, n
  write(*,'(6f8.4)') dlf(i,:) 
enddo

! compute the eigenvalues and eigenvectors  of dlf
!subroutine eigrg(a,n,isv, wr,wi,vr)
call eigrg(dlf,n,1, wr,wi, vr)

! Initialization of parameters to do numercial continuation
dir = 1
 



! err tolerance of target f to terminate Newton Method, 1.d-16 from the Gerard's book????
! Question: with the precision in RK78 as 1.d-14, it is not possible to get smaller tol....
!tol  = 1.d-10 
!prsc = 1.d-10 ! For lf problem,  1.d-11 is too small for the second family, ok for the other 2 families

! finally take 1.d-9 as the error control 

! Use the period as initial guess, TP = 2*pi/wi, wi is the pure imaginary eigenvalue
tmax  = 2 * 2*3.15/dabs(wi(2*ifam)) ! the maximum time for the first return to poincare section

print*, 'tmax=', tmax, 'wi=', wi(2*ifam)
read*

! Step size and error control for rk78 integration 
hmin = 1.d-10
hmax = 1.d0
e    = 1.d-14! 1.d-13 is from Gerard's book 

symt = 2 ! the 2nd symmetry: y=0 plane
ind  = 2 ! y
sec  = 0.d0 
imax = 1
call init_sympo(symt)   ! forsymmetric p.o.
call init_errctr(tol, prsc)
call init_poinc(ind, sec, hmin, hmax, e, tmax)
call init_write 



print*, 'check', ivr,'-th column of vr to use', vr(:, ivr) 
vrchs = vr(:, ivr)
vrchs = vrchs/dnrm2(6,vrchs,1)
print*, dnrm2(6,vrchs, 1), vrchs
 
!poinst =  eq +  epsl_po * vr(:, ivr) ! the initial guess for po, move along the corresponding eigenvector for a small distance epsl_po 
poinst =  eq +  epsl_po * vrchs !  



! test the continuation from the second families in the negative direction
!  0.4797289074490632D+02  0.4524274464352185D+00  0.0000000000000000D+00  0.4573275746506272D+00  0.0000000000000000D+00 -0.2226337907927061D-01  0.0000000000000000D+00  0.1975649365951931D+01  0.1998401444325282D-14  0.1000000000000000D-02
!poinst = (/0.4524274464352185d0, 0.d0,  0.4573275746506272d0,  0.d0, -0.2226337907927061d-1, 0.d0/)
!dir = -dir 
ds = 2.d-3 ! try bigger step to see if is possible to obtain the second family of p.o.s with smaller size, without jumping to another family

! Initialize for mmat
y0 = 0.d0
y0(1:6) = poinst  
y0(7:42:7) = 1.d0 !  the initial value of variational matrix if identity stored in y0(7:42)
 
!print*, 'before pofam'
!  subroutine pofam(yi,npo,imax, dir,ds, fpo, ynew, i, deriv, gr_cj)
call pofam(poinst, npo,imax, dir, ds, fpoinst, ynew, ipo, gr_lf, gr_cjlf)

print*, 'PO finisned!'
print*, 'No of real p.o. computed = ', ipo-1
read*
! --- plot the P.O ---

do i = 1, ipo-1 
  
  po0 = ynew(i,2:7)
  tpo = ynew(i,1)

!    ynew(i,:) = (/tp, yi, cj/) from pofam
  print*, i, '-th P.O. TP: ', tpo
  write(*,'(8f12.8)') ynew(i,:) 
!  read(*,*) 

!  subroutine plob(y0,t0,tf,tdir, n, ftag, deriv, gr_cj,  y) 
  call plob(po0, 0.d0, 1*tpo, 6, tdir, fpo, gr_lf, gr_cjlf, pof)
    
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
end program po3



















  

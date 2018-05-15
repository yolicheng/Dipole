program po2 !-- the second equilibrium
!      Content 			e.g. of the name 	data structure 
!   1. eqmf...  		egmf3u_bt1  		 this is not done here, in eqmf subroutine...
!   2. poinst                   eq3_bt10_poinsti(i=1,2,3) 	  !write(fpo,'(10d24.16)') tp, yi, cj, dcj, hminim ! PoInSt.dat 
!   2. po 			eq3_bt10_poi(i=1,2,3)		  !write(ftag,'(7e20.10)') t, y, cj ! po.dat 
!   3. monodromy matrix  	eq3_bt10_mmegvi\egp(i=1,2,3) (the eigenvalues and respective eigenvectors)
  
! '(10d24.16)') tp, yi, cj, dcj, hminim  ---- PoInSt.dat 
! '(8e20.10)')  t, y, cj 		 ---- po.dat  for plob
!  real: real, complex: real+imaginary   ---- eq3_bt10_mmegvi.dat

! 20160308 -- finish po computaion for the 3rd equilibrium, except the mulitiple shooting method, which is not the focus currently.
!  For 2nd equilibrium, regardless of the value of beta, always have 4 pure imaginary +  2(real) eigenvalues

! 20160301
! add the monodromy matrix computaion, specify all the dimensionless unit 
! Question: for the computation, use the dimensionless unit, for the plot, use real data ! For EM system, we did the same
!           save all the unit in a file lf


! The important thing is not the continued family of p.o., learn the mulitiple shooting method later, but here we will not focus on this part

! 20160222 
! check the matrix norm associated to the vector norm

! to be modified later-- to debug pomod, something is definitely wrong!!!! cool! pomod seems fine! all I need to do is modify the other 
!  the tolerance and presion is to be careful assigned
! For Gerard's book, P96, the percision required for the modified Mewthon's method has been taken equal to 1.d-11 
!                         and the bound for local errors in RK78 routine was set to 1.d-13
!tol = 1.d-13; prsc =1.d-11 ! the suggested values
! because the error control for rk78 is 1.d-13, all the poincare map-related data has no smaller precision
!  so we cannot ask more precision in Newton method

! try the case with beta = 1, test if the center manifold can not very well be approximated by the linearized system

use lf_mod ! the system related parameters of lorentz force problem
use po_mod ! archived subroutines to compute families of symmetric p.o.

implicit none

integer, parameter ::  neq = 42, &  ! compute also the variational matrix, 6+36=42
		       npo = 400  ! NO of orbits in pofam 
		       
! the dimension of the problem is defined in lf_mod, n=6 

! lf_mod Module Variables
! Global :   beta, cs, sgn, eq

integer ::  cs0, ieq, symt
real(kind=dp) ::  beta0, sec0, hmax, hmin, e, tmax, tol, prsc ! error control


! Local Variables
integer :: i, j,  ivr,  dir, imax, fpo, fpoinst, ifam, fmmegp, fmmegv,  &
	   tdir, ipo

real(kind=dp) ::  dlf3(n,n),  & ! differential of vector field of lorentz force
                  wr(n),wi(n), vr(n,n), y0(neq), &   ! eigenspace of variational matrix  
                  poinst(6), ds, vrchs(6), ynew(npo,8), epsl_po, cj,  & !pofam
                  po0(6), tpo, pof(6), & ! plpo
                  dlf(n,n), mmat(n,n), wr_mm(n), wi_mm(n), vr_mm(n,n)  ! MM 
                 

 character(len=70) :: fnpo, fnpoinst, fnmmegv, fnmmegp              

!external :: gr_lf, gr_cjlf ! withou use statement, we don't need to declare external here
 

real(kind=dp), external :: dnrm2 ! from package NAG Fortran Library Routine
!double precision FUNCTION F06EJF (N, X, INCX)

! Initialize the private variables eq1, sgn1 for  case 1
call init_lf_mod

! case 1: N=[1 0 0]; 2:  N=[0 1 0]; 3:  N=[0 0 1]; 
 cs0 = 1 ! use 2 to test the swap rule is coded correctly--ckd!

! the index of the equilibrium points
! 1:  q/m > 0,  x=0,y=0,z= \pm 1
! 2:  q/m < 0,  x = \pm (2 sqrt(9))^(1/3), y² = 2x²,z = 0
! 3:  q/m < 0   x = \pm (1/4/sqrt(2)^(1/3), y = 0,z² = x²  ! this the the case we study now 
ieq = 2 ! 3, x,0,z is the case that we study currently


! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy
! and  the index of column of vr to be used as the solution of the variational equatio, carefully checked by observing the eigenvectors


!! for bt = 1, 1 family, choose the sixth column is the eigenvector. --- this one is hard to converge -done
beta0 = 1.d0; 

ifam = 1 ; ivr = 6; ! prsc=1.d-10 not working, 1.d-9 is ok! but the convergence is bad
epsl_po = 1.d-6; ds = 1.d-5!
   

!!!! for bt = 2, 1 family, choose the first column is the eigenvector  - done 
!beta0 = 2.d0; ivr = 2; ifam = 1! prsc=1.d-10, espl_po = 1.d-5, h0 = 1.d-3 
!epsl_po = 1.d-4; ds = 1.d-3 

! Provided the case, beta and ieq,  beta, cs, sgn and eq are initialized by subroutine init_lf from module lf_mod
! subroutine init_lf(beta0, cs0, ieq)
call init_lf(beta0, cs0, ieq) 

print*, 'check the assignment with module' !--ckd
print*, 'beta, cs, sgn, eq', beta, cs, sgn, eq


!  ***************** compute the 3 families of PO ***************************
! save the data of po and  its initial state 
fpo = 20;  fpoinst = 21;  fmmegv = 23; fmmegp = 24 

! the idea to write to mulitiple files, 
! ! build filename -- i.dat
!write(fn,fmt='(a,i0,a)') filenum, '.dat'

!! open it with a fixed unit number
!open(unit=outunit,file=fn, form='formatted')

! remember to rename the data file for different families of po when beta = 10, for ieq = 1, no need to put ifam suffix
write(fnpo,    fmt='(a,i0,a,i0,a)') './dat/eq',ieq,'_bt', idint(beta), '_po.dat'  
write(fnpoinst,fmt='(a,i0,a,i0,a)') './dat/eq',ieq,'_bt', idint(beta), '_poinst.dat'
write(fnmmegv, fmt='(a,i0,a,i0,a)') './dat/eq',ieq,'_bt', idint(beta), '_mmegv.dat'
write(fnmmegp, fmt='(a,i0,a,i0,a)') './dat/eq',ieq,'_bt', idint(beta), '_mmegp.dat'

!write(fnmmegp, fmt='(a,i0,a)') './dat/eq3_bt',idint(beta), '_mmegp.dat'!without specify the family

print*, fnpo, fnpoinst, fnmmegv, fnmmegp
read(*,*)

open(fpo,file=fnpo, access ='append',status='replace')
open(fpoinst,file=fnpoinst, access ='append',status='replace')
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
print*,'wi', wi
read*
! Pay attention to the distance move along the eigenvectors-- For EMRTBP, epsl_po= 1.d-3 is ok and ds =1.d-3 is also ok
!epsl_po = 1.d-6! the magnitude of the variation of the initial guess on the p.o. 

!epsl_po = 1.d-6! the magnitude of the variation of the initial guess on the p.o. 

! but for lf problem, if we take ds =1.d-3, the continued family goes to cylinder rather than in a plane 

! Initialization of parameters to do numercial continuation
dir = 1
!ds = 1.d-4 !5.d-4 !1.d-3 ! step size for the continuation
!ds = 1.d-4


!tol = 1.d-13; prsc =1.d-11 ! the suggested values

imax = 1      ! time of crossing through y=0 plane for the differential correction
tdir = 1 ! integrate forward

! finally take 1.d-9 as the error control
!tol  = 1.d-10 ! err tolerance of target f to terminate Newton Method, from the Gerard's book, it should be 1.d-16?
!prsc = 1.d-10 ! For lf problem,  1.d-11 is too small for the second family, ok for the other 2 families

tol  = 1.d-9 ! err tolerance of target f to terminate Newton Method, from the Gerard's book, it should be 1.d-16?
prsc = 1.d-9 ! For lf problem,  1.d-11 is too small for the second family, ok for the other 2 families


! Use the period as initial guess, TP = 2*pi/wi, wi is the pure imaginary eigenvalue
tmax  = 1.2 * 2*3.15/dabs(wi(2*ifam)) ! the maximum time for the first return to poincare section

print*, 'tmax=', tmax, 'wi=', wi(2*ifam)
read*

! Step size and error control for rk78 integration 
hmin = 1.d-10
hmax = 1.d0
e    = 1.d-14! 1.d-13 is from Gerard's book 

symt = 2
sec0 = 0.d0 

! subroutine init_po(symt, sec, hmin, hmax, e, tmax, tol, prsc) - for pomod
call  init_po(symt, sec0, hmin, hmax, e, tmax, tol, prsc)


! ifam is the family of periodic orbit to study, in total three families 
! instead of computing the 3 families in a loop,  try doing it seperately, such that we can forget about the failure of continuation due to the 
! the strong unstability of the p.o. with big eigenvalues


print*, 'check', ivr,'-th column of vr to use', vr(:, ivr) 
vrchs = vr(:, ivr)
vrchs = vrchs/dnrm2(6,vrchs,1)
print*, dnrm2(6,vrchs, 1), vrchs
read*
 
!poinst =  eq +  epsl_po * vr(:, ivr) ! the initial guess for po, move along the corresponding eigenvector for a small distance epsl_po 
poinst =  eq +  epsl_po * vrchs !  


! Initialize for mmat
y0 = 0.d0
y0(1:6) = poinst  
y0(7:42:7) = 1.d0 !  the initial value of variational matrix if identity stored in y0(7:42)
 
!print*, 'before pofam'
!  subroutine pofam(yi,npo,imax, dir,ds, fpo, ynew, i, deriv, gr_cj)
call pofam(poinst,npo,imax, dir, ds, fpoinst, ynew, ipo, gr_lf, gr_cjlf)

print*, 'PO finisned!'
print*, 'No of real p.o. computed = ', ipo-1
read*
! --- plot the P.O ---
!subroutine plob(y0,t0,tf,ftag, y) 

do i = 1, ipo-1 
  
  po0 = ynew(i,2:7)
  tpo = ynew(i,1)

!    ynew(i,:) = (/tp, yi, cj/) from pofam
  print*, i, '-th P.O. TP: ', tpo
  write(*,'(8f12.8)') ynew(i,:) 
!  read(*,*) 

!  subroutine plob(y0,t0,tf,tdir,ftag, deriv, gr_cj,  y) 
  call plob(po0, 0.d0, 1*tpo, tdir, fpo, gr_lf, gr_cjlf, pof)
    
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
end program po2



















  

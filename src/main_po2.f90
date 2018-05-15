program po2 !-- the second equilibrium
!      Content 			e.g. of the name 	data structure 
!   1. eqmf...  		egmf3u_bt1  		 this is not done here, in eqmf subroutine...
!   2. poinst                   eq3_bt10_poinsti(i=1,2,3) 	  !write(fpo,'(10d24.16)') tp, yi, cj, dcj, hminim ! PoInSt.dat 
!   2. po 			eq3_bt10_poi(i=1,2,3)		  !write(ftag,'(7e20.10)') t, y, cj ! po.dat 
!   3. monodromy matrix  	eq3_bt10_mmegvi\egp(i=1,2,3) (the eigenvalues and respective eigenvectors)
  
! '(10d24.16)') tp, yi, cj, dcj, hminim  ---- PoInSt.dat 
! '(8e20.10)')  t, y, cj 		 ---- po.dat  for plob
!  real: real, complex: real+imaginary   ---- eq3_bt10_mmegvi.dat

! 20160309 
! !  For 2nd equilibrium, 2 families of p.o. ., but only one form of eigenvalues for the equilibrium point 
!    consider   beta = 1 as example,  2 pure imaginary eigenvalues and 4 real ones
!  if the eigenvectors do not have the same kind of symmetry, we need to try the other method, 
!   option 1- look at the zero component. 

! 20160308 -- finish po computaion for the 3rd equilibrium, except the mulitiple shooting method, which is not the focus currently.
!  For 1st equilibrium,  only 1 family of p.o. . consider only beta = 1 (complex), and beta = 2(real)

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
                       npo = 40 ! 52... for x=z, ifam2  !   1 !35 ! NO of orbits in pofam 

real(kind=dp) :: pi = 4.d0*datan(1.d0)   
    
! the dimension of the problem is defined in lf_mod, n=6 

! lf_mod Module Variables
! Global :   beta, cs, sgn, eq

integer ::  cs, ieq, symt, asymt, ind, debug, dsctr, anglectr, issymt, isarc 
real(kind=dp) ::  beta, sec, hmax, hmin, e,  tol, prsc, tol_poinc, x0_fxd ! error control
! Local Variables


integer :: i, j,  ivr,  dir, imax, fpo, fpoinst, ifam, fmmegp, fmmegv,  &
           tdir, ipo

real(kind=dp) ::  dlf3(n,n),  & ! differential of vector field of lorentz force
                  wr(n),wi(n), vr(n,n),   &   ! eigenspace of variational matrix  
                  poinst(6), ds, vrchs(6), ynew(npo,8), epsl_po, cj,  & !pofam
                  po0(6), tpo, pof(6), & ! plpo
                  dlf(n,n), mmat(n,n), wr_mm(n), wi_mm(n), vr_mm(n,n), &  ! MM 
                  tp0, csangle_min, ds_max, ds_min 

                 
 character(len=70) :: fnpo, fnpoinst, fnmmegv, fnmmegp              

real(kind=dp), external :: dnrm2 ! from package... 

! debug or not 
debug =  0

! automatic step control depending on the number of iteration for Newton Method
dsctr =  1 

! contol the angle between the vector field of two latest consecutive points along the 
! characteristic curve
anglectr = 1

isarc = 1 ! along arc-length paramter
!isarc = 2 ! along period
!isarc = 3 ! along energy


!issymt  = 1 ! for symmetric approach with poincare map 
!issymt = 0 ! for asymmetric approach with poincare map 
issymt = 2  ! the general approach 


! the minimum of the cosine of the angle between two consecutive vectors
! cos(0.1 rad) = 0.995;    cos(0.37 rad) = 0.95  
 csangle_min = 0.995d0  

! bound vaules of stepsize ds
!  With the control of termination added when the arc step goes too smaller... 
!  better not to be too big, with 2.d-2, or 1.d-2, for now is good choice  

ds_max = 1.d-0
ds_min = 1.d-5 ! but this bound is not working for angle control....
    
call init_debug(debug, dsctr, anglectr, csangle_min,   isarc, ds_min, ds_max)

! Initialize the private variables eq1, sgn1 for  case 1
call init_lf_mod

! case 1: N=[1 0 0]; 2:  N=[0 1 0]; 3:  N=[0 0 1]; 
 cs = 1 ! use 2 to test the swap rule is coded correctly--ckd!

! the index of the equilibrium points
! 1:  q/m > 0,  x=0,y=0,z= \pm 1
! 2:  q/m < 0,  x = \pm (2 sqrt(9))^(1/3), y² = 2x²,z = 0
! 3:  q/m < 0   x = \pm (1/4/sqrt(2)^(1/3), y = 0,z² = x²  ! this the the case we study now 
ieq = 2 ! 


! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy
! and  the index of column of vr to be used as the solution of the variational equatio, carefully checked by observing the eigenvectors

!  for bt = 1, 1 family, choose the first column is the eigenvector  - done 
beta = 1.d0; ivr = 2;  ifam = 1 
epsl_po = 1.d-6; ds = 1.d-5 

tol  = 1.d-11! err tolerance of target f to terminate Newton Method, from the Gerard's book, it should be 1.d-16?
prsc = 1.d-14 ! For lf problem,  1.d-11 is too small for the second family, ok for the other 2 families


! finally take 1.d-9 as the error control
!tol  = 1.d-10 ! err tolerance of target f to terminate Newton Method, from the Gerard's book, it should be 1.d-16?
!prsc = 1.d-10 ! For lf problem,  1.d-11 is too small for the second family, ok for the other 2 families

! Provided the case, beta and ieq,  beta, cs, sgn and eq are initialized by subroutine init_lf from module lf_mod
! subroutine init_lf(beta0, cs0, ieq)
call init_lf(beta, cs, ieq) 

print*, 'check the assignment with module' !--ckd
print*, 'beta, cs, sgn, eq', beta, cs, sgn, eq
read*

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
read(*,*)  

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

!! Initialization of parameters to do numercial continuation
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


if (issymt == 1) then 

!! --- Approach 1, by symmetric way 
  symt = 2 ! the 2nd symmetry: y=0 plane
  ind  = 2 ! y
  sec  = 0.d0 
  imax = 1
  call init_symtpo(symt)   ! forsymmetric p.o.

elseif(issymt == 0 ) then

! --- Approach 0, by general asymmetric way
!  asymt = 1  ! still try y = 0 plane            !alternative: x=x0 plane
!  sec   = poinst(asymt) 
!  
  asymt = 2  ! still try y = 0 plane            !alternative: x=x0 plane
  sec = 0.d0 ! y = 0 still
  
  ind = asymt
  imax  = 2
  
!  subroutine init_asympo(asymt ) ! for asymmetric p.o.
  call init_asymtpo(asymt)   ! for symmetric p.o.

elseif (issymt == 2) then 

  call init_grpo 

!elseif (issymt == 3) then 

!! --- Approach 2, by general asymmetric way
!  ind = 1
!  x0_fxd =  eq(ind)
!  call init_grpo(ind, x0_fxd) ! no input
endif 


!  subroutine init_poinc(ind0, sec0, imax0,  hmin0, hmax0, e0)
if (issymt /= 2) then 
  ! Step size and error control for rk78 integration 
  hmin = 1.d-10
  hmax = 1.d0
  e    = 1.d-14! 1.d-13 is from Gerard's book 
  tol_poinc = 1.d-14 
  call init_poinc(ind, sec, tol_poinc, imax, hmin, hmax, e)
endif 


call init_errctr(tol, prsc)
call init_writedx 


!  subroutine pofam(yi, tp0, npo, dir,ds, fpo, ynew, ipo, deriv, gr_cj)
!ds = 0.01d0
call pofam(poinst, tp0, npo, dir,ds, fpoinst, ynew, ipo, gr_lf, gr_cjlf)


! Ask gerard about this, or check what happens for the first intersection with x-z plane 
! Symmetry :  x-z planar symmetry  +  y-axis symmetry

! Eigenvalues
!0.00000000  2.69451269      0.00000000  -2.69451269    1.41421356    -1.41421356     0.00000000  1.11337386      0.00000000  -1.11337386

! Eigenvectors:  real-imaginary  
!  0.1730 -0.0908    0.1730  0.0908    0.4082   0.2887  -0.0609  0.0000   -0.0609 -0.0000  
!  0.0000  0.2260    0.0000 -0.2260    0.0000   0.4082  -0.2264  0.1782   -0.2264 -0.1782  
!  0.0829  0.1579    0.0829 -0.1579   -0.4082   0.2887   0.0000 -0.5998    0.0000  0.5998  
!  0.2446  0.4661    0.2446 -0.4661    0.5774  -0.4082   0.0000 -0.0678    0.0000  0.0678  
! -0.6090  0.0000   -0.6090 -0.0000    0.0000  -0.5774  -0.1984 -0.2520   -0.1984  0.2520  
! -0.4255  0.2233   -0.4255 -0.2233   -0.5774  -0.4082   0.6679  0.0000    0.6679 -0.0000

! we have to explore how the orbits starting from nearby the equilibrium evolves, so take a neighborhood sphere of radius 1.d-7 around the equilibrium, take poincare map at time domain.

!After an approximate period, T = 2*pi/ w ( w is the pure imaginary eigenvalue)
 
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
  
  if( tpo < 1.d-3) return
  write(*,'(8f12.8)') ynew(i,:) 
!  read(*,*) 

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
end program po2









  

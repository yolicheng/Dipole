program po 

! this is using the initial states of halo and lyapunov orbits to test the refinement of po 

!      Content 			e.g. of the name 	data structure 
!   1. eqmf...  		egmf3u_bt1  		
!   2. poinst                   eq3_bt10_poinst 	  !write(fpo,'(10d24.16)') tp, yi, cj, dcj, hminim ! PoInSt.dat 
!   2. po 			eq3_bt10_po		  !write(ftag,'(7e20.10)') t, y, cj ! po.dat 
!   3. monodromy matrix  	eqmmat3_bt10 (the eigenvalues and respective eigenvectors)
  
! '(10d24.16)') tp, yi, cj, dcj, hminim  ---- PoInSt.dat 
! '(7e20.10)') t, y, cj 		 ---- po.dat 

! 20160303 
! check--!!!!  where is the problem, why after the call of gr_rk78, all the values in module pomod are lost and reassigned to strange values 

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
		       npo = 10, & ! NO of orbits in pofam 
		       nmf = 10     ! NO of orbits on the manifold

integer ::  cs0, ieq, vr_ind(3) , ivr, symt
real(kind=dp) ::  beta0, st0(n), x0, dlf(n,n)

! Local Variables
integer :: i, j, ifam, ffam, dir, imax, id, fpo, fpoinst, fmmat, fegp, fegv, tdir, ipo

real(kind=dp) ::  dlf3(n,n),  & ! differential of vector field of lorentz force
                  wr(n),wi(n), vr(n,n), y0(neq), &   ! eigenspace of variational matrix  
                  poinst(6), ds, vf(neq), ynew(npo,8), epsl_po, cj,  & !pofam
                  po0(6), tpo, pof(6), & ! plpo
                  mmat(n,n), wr_mm(n,n), wi_mm(n,n), vr_mm(n,n), vrchs(6), & ! MM 
                  ymfd(nmf,6),  epsl,  &! mfd 
                  sec, hmax, hmin, e, tmax, tol, prsc ! error control
                  
external ::   gr_rtbp, gr_cjrtbp
 
real(kind=dp), external :: dnrm2

! use this to initialize eq1, sgn1 
call init_lf_mod

! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy
! and  the index of column of vr to be used as the solution of the variational equatio, carefully checked by observing the eigenvectors

! for bt = 10, 3 families of p.o.
beta0 = 10.d0; vr_ind = (/1, 4, 6/); ffam = 1

!! for bt = 1, 1 family, ifam = 1, the sixth column is the eigenvector we choose, this one is hard to converge
!beta0 = 1.d0; vr_ind = (/6, 6, 6/); ffam = 1

!! for bt = 2, 1 family, ifam = 1, the first column is the eigenvector we choose
!beta0 = 2.d0; vr_ind = (/1, 1, 1/); ffam = 1

!  ***************** compute the 3 families of PO ***************************
! save the data of po and  its initial state 
fpo = 20;  fpoinst = 21 
if (beta0 == 10.d0) then 
  open(fpo,file='./dat/eq3_bt10_po.dat',access ='append',status='replace')  
  open(fpoinst,file='./dat/eq3_bt10_poinst.dat',access ='append',status='replace')  

else if (beta0 == 2.d0) then 
  open(fpo,file='./dat/eq3_bt2_po.dat',access ='append',status='replace')  
  open(fpoinst,file='./dat/eq3_bt2_poinst.dat',access ='append',status='replace')  

else if (beta0 == 1.d0) then 
  open(fpo,file='./dat/eq3_bt1_po.dat',access ='append',status='replace')  
  open(fpoinst,file='./dat/eq3_bt1_poinst.dat',access ='append',status='replace')  
endif 

! case 1: N=[1 0 0]; 2:  N=[0 1 0]; 3:  N=[0 0 1]; 
 cs0 = 1 ! use 2 to test the swap rule is coded correctly--ckd!

! the index of the equilibrium points
! 1:  q/m > 0,  x=0,y=0,z= \pm 1
! 2:  q/m < 0,  x = \pm (2 sqrt(9))^(1/3), y² = 2x²,z = 0
! 3:  q/m < 0   x = \pm (1/4/sqrt(2)^(1/3), y = 0,z² = x²  ! this the the case we study now 
  
ieq = 3 ! 3, x,0,z is the case that we study currently
 

! the state of the equilibrium point, positon+velocity
!eqst = (/1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0/)! test

! beta, cs, sgn, eq  are initialized by subroutine init_lf in  module lf_mod
! provided the case, beta, and ieq 

! subroutine init_lf(beta0, cs0, ieq)
call init_lf(beta0, cs0, ieq) 

print*, 'check the assignment with module' !--ckd
print*, 'beta, cs, sgn, eq', beta, cs, sgn, eq
read(*,*)

! the differential matrix of the lorentz force with respective to state 
!subroutine dflrtz( x0, beta, cs, sgn, dlf)
call  dflrtz( eq, dlf)
!call dflrtz( eq, beta, cs, sgn, dlf)

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

! compute the eigenvalues and eigenvectors 
!subroutine eigrg(a,n,isv, wr,wi,vr)
call eigrg(dlf,n,1, wr,wi,vr)

! check how wr is stored! -ckd, in column-wise order
! print*, 'vr', vr
! !read(*,*)  

epsl_po  = 1.d-6 ! the magnitude of the variation 

! initialization of parameters to do numercial continuation
dir = 1
ds =  1.d-3 !5.d-4 !1.d-3 ! step size for the continuation

!tol = 1.d-13; prsc =1.d-11 ! the suggested values

imax = 1      ! time of crossing through y=0 plane for the differential correction
tdir = 1 ! integrate forward

tol = 1.d-11 ! err tolerance of target f to terminate Newton Method, from the Gerard's book, it should be 1.d-16?
prsc = 1.d-11 ! For lf problem,  1.d-11 is too small for the second family, ok for the other 2 families

tmax = 30 
hmin = 1.d-10
hmax = 1.d0

e    = 1.d-15 ! 1.d-13 is from Gerard's book 
symt = 2
sec = 0.d0 

!  subroutine init_po(symt, sec, hmin, hmax, e, tmax, tol, prsc)
call  init_po(symt, sec, hmin, hmax, e, tmax, tol, prsc)



! ifam is the family of periodic orbit to study, in total three families 
! instead of computing the 3 families in a loop,  try doing it seperately, such that we can forget about the failure of continuation due to the 
! the strong unstability of the p.o. with big eigenvalues
do ifam = 1,  ffam!3 ! because the period is too big, and the accumulated numerical error is also big.

ivr = vr_ind(ifam) ! 

 y0 = 0.d0 
  y0(7:42:7) = 1.d0 !  the initial value of variational matrix if identity stored in y0(7:42)
 
!  print*, 'check', ivr,'-th column of vr to use', vr(:, ivr) 
!  vrchs = vr(:, ivr)
!  vrchs = vrchs/dnrm2(6,vrchs,1)
!  print*, dnrm2(6,vr(:, ivr),1), vrchs
 
  poinst =  eq +  epsl_po * vr(:, ivr) ! the ifam-th column of vr corresponds to the first eigenvalue

  print*, 'check', ifam,'-th poinst', poinst
  print*
  read(*,*)   

! check gr_lf
!  subroutine gr_lf(t,y,neq,f, beta,cs,sgn)
!print*, 'check variational matrix of lorentz force'


  do i = 1, 6, 2 ! Note: only works for kind 2 symmetry: 1 3 5, 
   y0(i) = poinst(i) ! the other three are zero!
  enddo

!print*, 'before gr_lf, check, y0,beta, cs,sgn',y0(1:6), beta,cs,sgn
 
!test halo, lyapunov and vertical lyapunov orbits --- ckd!
poinst = (/-0.8233804423d0, 0.d0, 0.0146213594d0, 0.d0, -0.1298107931d0, 0.d0/) ! halo 1
!poinst = (/-0.8233906879d0, 0.d0, 0.0016206515d0, 0.d0, -0.1263702002d0, 0.d0/) ! halo 2
!poinst = (/-0.8416266520d0, 0.d0, 0.d0, 0.d0, 0.0381494892d0, 0.d0/) ! lyap 1
!poinst = (/-0.8449604387d0, 0.d0, 0.d0, 0.d0, 0.0637254594d0, 0.d0/) ! lyap 1
!poinst = (/-0.8372892949d0, 0.d0, 0.0135816756d0,  0.d0, 0.0002773899d0, 0.d0/) ! vt_lyap 1
call pofam(poinst,npo,imax, dir, ds, fpoinst, ynew, ipo,  gr_rtbp, gr_cjrtbp)
!------------------------------------------------


!print*, 'before pofam'
!  subroutine pofam(yi,npo,imax, dir,ds, fpo, ynew, deriv, gr_cj)
!  call pofam(poinst,npo,imax, dir, ds, fpoinst, ynew, ipo, gr_lf, gr_cjlf)

  print*, 'PO finisned!'
 
! --- plot the P.O ---
!subroutine plob(y0,t0,tf,ftag, y) 
  do i = 1, ipo-1  
  !  write(fpoinst,'(7f12.8)') ynew(i,:) 
  
    po0 = ynew(i,2:7)
    tpo = ynew(i,1)
  
!    print*, i, '-th P.O. TP: ', tpo
!    write(*,'(8f12.8)') ynew(i,:) 
    !read(*,*) 

   call plob(po0, 0.d0, 2*tpo, tdir, fpo, gr_rtbp, gr_cjrtbp,  pof) ! test EM rtbp
!  subroutine plob(y0,t0,tf,tdir,ftag, deriv, gr_cj,  y) 
!    call plob(po0, 0.d0, 1*tpo, tdir, fpo, gr_lf, gr_cjlf,  pof)
  enddo   
  
  write(fpo,*) 
  write(fpoinst,*); write(fpoinst,*)
end do

  

!******************************************************************************
stop
end program po



















  

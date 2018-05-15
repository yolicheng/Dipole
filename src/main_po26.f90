program po26 !-- the second equilibrium -- test po , remember to change the name first to avoid overwritting
!!  suffix 6:  vx =  -1.d-5, vy =  1.d-5, vz = - ! -- done 

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
real(kind=dp), parameter :: pi = 4*datan(1.d0)

integer, parameter ::  neq = 42, &  ! compute also the variational matrix, 6+36=42
		       npo = 400  ! NO of orbits in pofam 
		       
! the dimension of the problem is defined in lf_mod, n=6 

! lf_mod Module Variables
! Global :   beta, cs, sgn, eq

integer ::  cs0, ieq, symt
real(kind=dp) ::  beta0, sec0, hmax, hmin, e, tmax, tol, prsc ! error control


! Local Variables
integer :: i, j,  ivr,  dir, imax, fpo, fpoinst, ifam, fmmegp, fmmegv,  &
	   tdir, iob, fpots, fpopcts, fpopotts, npt, iscj, ipo, k, isstep, ispo

real(kind=dp) ::  dlf3(n,n),  & ! differential of vector field of lorentz force
                  wr(n),wi(n), vr(n,n), y0(neq), &   ! eigenspace of variational matrix  
                  poinst(6), ds, vrchs(6), ynew(npo,8), epsl_po, cj,  & !pofam
                  po0(6), tpo, pof(6), & ! plpo
                  dlf(n,n), mmat(n,n), wr_mm(n), wi_mm(n), vr_mm(n,n), &  ! MM 
                  cj1, dt, epsl
                 

logical :: gnuplot_open = .false.

 character(len=70) :: fnpo, fnpoinst, fnmmegv, fnmmegp, fnpots, fnpopcts, fnpopotts            

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

!  for bt = 1, 1 family, choose the first column is the eigenvector  - done 
beta0 = 1.d0; ivr = 2; ifam = 1 
epsl_po = 1.d-6; ds = 1.d-5 
! finally take 1.d-9 as the error control
tol  = 1.d-10 ! err tolerance of target f to terminate Newton Method, from the Gerard's book, it should be 1.d-16?
prsc = 1.d-10 ! For lf problem,  1.d-11 is too small for the second family, ok for the other 2 families

! Provided the case, beta and ieq,  beta, cs, sgn and eq are initialized by subroutine init_lf from module lf_mod
! subroutine init_lf(beta0, cs0, ieq)
call init_lf(beta0, cs0, ieq) 

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

! Initialization of parameters to do numercial continuation
dir = 1
!ds = 1.d-4 !5.d-4 !1.d-3 ! step size for the continuation
!ds = 1.d-4


!tol = 1.d-13; prsc =1.d-11 ! the suggested values

imax = 1      ! time of crossing through y=0 plane for the differential correction
tdir = 1 ! integrate forward

! Use the period as initial guess, TP = 2*pi/wi, wi is the pure imaginary eigenvalue
!tmax  = 1.2 * 2*3.15/dabs(wi(2*ifam)) ! the maximum time for the first return to poincare section

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

epsl = 1.d-6
poinst = eq ! the equilibrium point 
npt =30
dt = 2*pi/npt
! take 30 points in a circle, this is for the position,  on the x-y plane
! for velocity, keep the vx and vy, modify vz to obtain the same energy ( the first one)

ifam = 1; 
fpots = 67; fpopcts = 68; fpopotts = 69
write(fnpots,  fmt='(a,i0,a,i0,a,i0,a)') './dat/eq',ieq,'_bt', idint(beta), '_pots',ifam,'.dat' 
write(fnpopcts,  fmt='(a,i0,a,i0,a,i0,a)') './dat/eq',ieq,'_bt', idint(beta), '_popcts',ifam,'.dat'   
write(fnpopotts,  fmt='(a,i0,a,i0,a,i0,a)') './dat/eq',ieq,'_bt', idint(beta), '_popotts',ifam,'.dat'   

isstep = 0
if( isstep == 0) then ! if we want to plot orbit one by one -- datafile with a 2 suffix
!  suffix 1:  vx =  1.d-5, vy = 1.d-5, vz = + ! to do 
!  suffix 2:  vx =  1.d-5, vy = 1.d-5, vz = - ! to do 

!  suffix 3:  vx =  1.d-5, vy = -1.d-5, vz = + !  
!  suffix 4:  vx =  1.d-5, vy = -1.d-5, vz = - ! -- done 

!  suffix 5:  vx = -1.d-5, vy = 1.d-5, vz = +
!  suffix 6:  vx = -1.d-5, vy = 1.d-5, vz = -

!  suffix 7:  vx = -1.d-5, vy = -1.d-5, vz = +
!  suffix 8:  vx = -1.d-5, vy = -1.d-5, vz = -


  write(fnpots,  fmt='(a,i0,a,i0,a,i0,a)') './dat/eq',ieq,'_bt', idint(beta), '_pots',ifam,'6.dat' 
  write(fnpopcts,  fmt='(a,i0,a,i0,a,i0,a)') './dat/eq',ieq,'_bt', idint(beta), '_popcts',ifam,'6.dat'   
! potential p.o.
  write(fnpopotts,  fmt='(a,i0,a,i0,a,i0,a)') './dat/eq',ieq,'_bt', idint(beta), '_popotts',ifam,'6.dat'   


  open(fpots,  file=fnpots, access ='append',status='replace')
  open(fpopcts,  file=fnpopcts, access ='append',status='replace')
  open(fpopotts,  file=fnpopotts, access ='append',status='replace')
  write(fpopotts, *) '# vx =  -1.d-5, 	vy =  1.d-5, 	vz = - ' ! suffix 1
  write(fpopotts, *) '# ds	ipo	i	j	 k	(x, y, z, vx, vy, vz)'
  
endif
 
tpo  = 1.0d0 * 2*pi/dabs(wi(2*ifam)) ! keep as 1T--a little bit greater than the T

!  subroutine plob(y0,t0,tf,tdir,ftag, deriv, gr_cj,  y) 
do i = 1, npt !  - pick the sphere -20160311
!do i = npt/4+1, npt/4+1  ! pick the x-y plane

do j = 1, npt
  ispo = 0
  poinst = eq ! the equilibrium point  

  poinst(1) = eq(1) + epsl  * dsin( (i-1)* dt ) * dcos( (j-1)* dt ) ! -x
  poinst(2) = eq(2) + epsl  * dsin( (i-1)* dt ) * dsin( (j-1)* dt ) ! -y 
  poinst(3) = eq(3) + epsl  * dcos( (i-1)* dt ) ! -z 
  
!  print*, dcos( (i-1)* dt )
!  read*

! randomly assign some small value for vx- vy, and the modify vz to obtain the fixed energy
  poinst(4) =  -1.d-5 ! vx
  poinst(5) =   1.d-5 ! vy
  
  if (i== 1 .and. j == 1) then  ! -- for the sphere
!  if (i==npt/4+1 .and. j == 1) then  ! -- for the x-y plane
    
    call gr_cjlf(poinst, cj1)
  else
    call lf2cj(poinst, 6, cj1, iscj) ! 6: vz to be modified
    if (iscj == 0) cycle
!    print*, 'new state', poinst
  endif 
  
  poinst(6) =  -poinst(6) ! vz take the negative value
  
  ipo = (i-1)*npt + j ! the index of the current orbit
!  print*, 'ipo, i, j, poinst'
!  print*,  ipo, i, j, poinst 

! intead of saving all the orbits, just save one, open and save inside one loop, this will make things easier 
  if( isstep == 1) then 
    open(fpots,  file=fnpots, access ='append',status='replace') 
    open(fpopcts,  file=fnpopcts, access ='append',status='replace')
  endif  
  
  do k = 1, 3
    call plob(poinst, (k-1)*tpo, k*tpo, 1, fpots, gr_lf, gr_cjlf, pof) !integrate forward
    write(fpopcts, *) k*tpo, pof(1:6) ! the poincare map at t=tpo in time domain
    write(fpots,*);
    
    ! the distance from the initial point 
!    ds = dsqrt( (pof(1)-0.504361916474301 )**2 + (pof(2) - 0.713275462622442)**2 + pof(3)**2 )
    
    ds = dsqrt( (pof(1)- poinst(1) )**2 + (pof(2) - poinst(2))**2 + ( pof(3)-poinst(3) )**2 )
    
    if( ds < 1.d-5) then 
      print*, 'ds, ipo, i, j, k, poinst'
      print*,  ds, ipo, i, j, k, poinst 
      ! save the initial state to file 
      write(fpopotts, *) ds, ipo, i, j, k, poinst 
      
      ispo = 1 ! keep this orbits
      if( isstep == 1)    read*
!      exit
    endif
    
    
!    write(*, *) k*tpo, pof(1:6) ! the poincare map at t=tpo in time domain
    poinst = pof
!    call system('gnuplot /home/yu/Desktop/dipole/dat/pots_step.pl')
!    read*
  enddo 
  
  write(fpots,*);  write(fpots,*);  write(fpopcts,*)
  
  if( isstep == 1) then 
    close(fpots)
    close(fpopcts)
    
    if(ispo == 1) then 
    call system('gnuplot /home/yu/Desktop/dipole/dat/pots_step.pl')
    read*
    endif
    
  endif
 
   
!  read*
enddo    
enddo 

stop     
 
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
call pofam(poinst,npo,imax, dir, ds, fpoinst, ynew, iob, gr_lf, gr_cjlf)

print*, 'PO finisned!'
print*, 'No of real p.o. computed = ', iob-1
read*
! --- plot the P.O ---
!subroutine plob(y0,t0,tf,ftag, y) 

do i = 1, iob-1 
  
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
end program po26









  

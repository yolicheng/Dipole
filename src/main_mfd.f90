program mfd


use lf_mod

implicit none

integer, parameter ::  neq = 42, & !
		       npo = 10, & ! NO of orbits in pofam
		       nmf = 10    ! NO of orbits on the manifold

integer ::  cs0, ieq, vr_ind(3), ind
real(kind=dp) ::  beta0, st0(n), x0, dlf(n,n)

! local variables
integer :: i, j, ifam, dir, imax, id, fpo, fpoinst, fmmat, fegp, fegv

real(kind=dp) ::  dlf3(n,n), dlfsswap(n),swap2(n), & ! differential of vector field of lorentz force
                  wr(n),wi(n), vr(n,n), y0(neq), &   ! stability of equilibrium points
                  poinst(6), ds,e2, vf(neq), ynew(npo,7), epsl_po, cj, aaa , & !pofam
                  po0(6), tpo, pof(6), & ! plpo
                  mmat(n,n), wr_mm(n), wi_mm(n), vr_mm(n,n), & ! MM 
                  ymfd(nmf,6), epsl ! mfd 
external :: gr_lf

! Initialize the private variables eq1, sgn1 
call init_lf_mod

! case 1: N=[1 0 0]; 2:  N=[0 1 0]; 3:  N=[0 0 1]; 
 cs0 = 1 ! use 2 to test the swap rule is coded correctly--ckd!

! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy
beta0 = 10.d0

! the index of the equilibrium points
! 1:  q/m > 0,  x=0,y=0,z= \pm 1
! 2:  q/m < 0,  x = \pm (2 sqrt(9))^(1/3), y² = 2x²,z = 0
! 3:  q/m < 0   x = \pm (1/4/sqrt(2)^(1/3), y = 0,z² = x²  ! this the the case we study now 
  
ieq = 3 ! 3, x,0,z is the case that we study currently

! Provided the case, beta and ieq,  beta, cs, sgn and eq are initialized by subroutine init_lf from module lf_mod
! subroutine init_lf(beta0, cs0, ieq)
call init_lf(beta0, cs0, ieq) 

print*, 'check the assignment with module' !--ckd
print*, 'beta, cs, sgn, eq', beta, cs, sgn, eq

! the differential matrix of the lorentz force with respective to state 
!subroutine dflrtz( x0, beta, cs, sgn, dlf)
call dflrtz( eq, beta, cs, sgn, dlf)

! check the energy 
!C1 = -3   
!C2 = 2.2894  
!C3 = 1.8899 

call gr_cjlf(eq, cj)
print*,'check energy!, cj, ieq,', cj, ieq, eq
read(*,*) aaa


! check the differential of the vector field
do i = 1, n
 write(*,'(6f8.4)') dlf(i,:) 
enddo

! The eigenvalues and eigenvectors of variatioal matrix at the equilibrium
!subroutine eigrg(a,n,isv, wr,wi,vr)
call eigrg(dlf,n,1, wr,wi,vr)
 
! the index of column of vr to be used as the solution of the variational equation
vr_ind = (/1, 4, 6/)

!  ***************** compute the 3 families of PO ***************************
! save the data of po and  its initial state 
!fpo = 20;  fpoinst = 21 
!open(fpo,file='./dat/po.dat',access ='append',status='replace')  
!open(fpoinst,file='./dat/poinst.dat',access ='append',status='replace')  

epsl_po  = 1.d-3 ! the magnitude of the variation 

! initialization of parameters to do numercial continuation
dir = 1
ds = 1.d-3 ! displacement of the continuation
id = 1   ! the type of symmetry,  1,2- xz plane, as in halo orbit, 
e2  = 1.d-10 !err tolerance for numerical continuation
imax = 1 ! time of crossing through y=0 plane for the differential correction

! save the eigenvalues and eigenvectors of monodromy
fmmat = 22; fegv = 23; fegp = 24 
open(fmmat,file='./dat/mmat.dat',access ='append',status='replace')  
open(fegv, file='./dat/mmegv.dat',access ='append',status='replace') 
open(fegp, file='./dat/mmegp.dat',access ='append',status='replace')
 
 
! ifam is the family of periodic orbit to study, in total three families 
do ifam = 1, 3

ind = vr_ind(ifam) ! 
y0 = 0.d0 
y0(7:42:7) = 1.d0 !  the initial value of variational matrix if identity stored in y0(7:42)


poinst =  eq +  epsl_po * vr(:, ind) ! the ifam-th column of vr corresponds to the first eigenvalue
print*, 'check', ind,'-th column of vr to use', vr(:, ind) 
print*
!read(*,*) aaa

print*, 'check', ifam,'-th poinst', poinst
print*
!read(*,*) aaa

! check gr_lf
!  subroutine gr_lf(t,y,neq,f, beta,cs,sgn)
!print*, 'check variational matrix of lorentz force'


do i = 1, 6, 2 ! Note: only works for kind 2 symmetry: 1 3 5, 
 y0(i) = poinst(i) ! the other three are zero!
enddo
 

!print*, 'before pofam'
call pofam(poinst, npo, imax, e2, dir, ds, id, 2, gr_lf, ynew)
print*, 'PO finisned!'



 
! --- Compute the Monodromy Matrix and Manifold------


do i = 1, npo !npo  ! -- test the first po 
  po0 = ynew(i,2:7)
  tpo = ynew(i,1)
  
! --- Monodramy matrix
!subroutine monomat(yi,tp, mmat, deriv)
  print*, 'refined initial state, tp, ynew', tpo, po0
  call monomat(po0, tpo, mmat, gr_lf)
  
! print mmat to file, mmat.dat  
  do j = 1, n
    write(fmmat,'(6d20.10)') mmat(j,:)
  enddo  
  
  write(fmmat, *)  ! add a blank line 
  
! analyze the stability of the monodramy matrix, a big step forward!
!subroutine eigrg(a,n,isv, wr,wi,vr)
  call eigrg(mmat,n,1, wr_mm, wi_mm, vr_mm)
  
  print*, 'Eigenvalues and eigenvectors, mmat!!!'
! it seems the real part of the eigenvectors of the monodramy matrix are nearly 1 
! so the stability is not so straightforward
! try the power method to see if we can get the dominant eigenvalue and eigenvector

!  fegv = 6; fegp = 6 ! print to screen 
  print*, 'eigenvalues, real part'
  print*,  wr_mm(1:6)
  
  print*, 'eigenvalues, imaginary part'
  print*,  wi_mm(1:6)
  
  print*
  
  call prt_eigval( n, fegv, wr_mm, wi_mm )
  call prt_eigvec( n, fegp, wi_mm, vr_mm )
  
  write(fegv,*)
  write(fegp,*)
  write(fegp,*)
!  
!  read(*,*) aaa
! 
!! mfd: unstable
!!subroutine mfdinst(ypo, tp, mmat, nmf, epsl, ymfd, deriv)
!  
!  epsl = 1.d-6
!  
!  call mfdinst(po0, tpo, mmat, nmf, epsl, ymfd, gr_lf) 
!       
!  do j = 1, nmf
!    write(*,'(6d20.10)') ymfd(j,:)
!  enddo 
!  
!  read(*,*) aaa
  
enddo

write(fegv,*)
write(fegp,*)

enddo


!******************************************************************************
stop
end program mfd 



















  

program vp1


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
                  mmat(n,n), wr_mm(n,n), wi_mm(n,n), vr_mm(n,n), & ! MM 
                  ymfd(nmf,6), epsl ! mfd 
external :: gr_lf


! use this to initialize eq1, sgn1 
call init_lf_mod


! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy
beta0 = 10.d0

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
!print*, 'beta, cs, sgn, eq', beta, cs, sgn, eq

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

do i = 1, n
 write(*,'(6f8.4)') dlf(i,:) 
enddo

! compute the eigenvalues and eigenvectors 
!subroutine eigrg(a,n,isv, wr,wi,vr)
call eigrg(dlf,n,1, wr,wi,vr)

! check how wr is stored! -ckd, in column-wise order
! print*, 'vr', vr
! read(*,*) aaa

! now we finish the eigenvalue and eigenvector part, goes to the computation of periodic orbit

!  Notice!!!!  - the general form of matrix differential equation
! x'(t)  = A x(t)
!  In the case where A has n distinct eigenvalues, this differential equation has the following general solution:
!  x(t) = c_1 * e^(lam1 * t) u_1 + c_2 * e^(lam2 * t) u_2 + ... 
! where lam1, lam2, ... are the eigenvalues of A; u_1, u_2, ... are the respective eigenvectors, c_1, c_2, ... are constants

! for the case complex eigenvalues, they always come in conjucated pair, as well as the eigenvectors.

!     we have the thoerem for the complex eigenvalues 
! Given a system x = Ax, where A is a real matrix. If x = x1 + i*x2 is 
! a complex solution, then its real and imaginary parts x1 , x2 are 
! also solutions to the system.
 
! 1 .the initial guess, linear approximation, which is the pertubation around the equilibrium point
! 2. refinement of the P.O. 
! 3. stability analysis of the P.O.
! 4. computation of invariant manifold (stable ones are to be interested)

! for the periodic orbits, the eigenvalues are pure imaginary

! ************** the first studied PO which is symmetric wrt yz(x=0)   plane 
! we are looking for the initial state of form (x,0,z,0,vy,0)

! The general form of the solution of the system equation is :

! X = c1 * e^(lam1 * i * t) * v1 + c2 * e^(-lam1 * i * t) * /bar v1 ! the first family 
!   + c3 * e^(lam2 * i * t) * v2 + c4 * e^(-lam2 * i * t) * /bar v2 ! the second family  
!   + c5 * e^(lam2 * i * t) * v3 + c6 * e^(-lam2 * i * t) * /bar v3 ! the third family 

! if the coefficients of the other two families are 0, then we obtain 1 family of periodic orbits

!  epsilon =   c1 * e^(lam1 * i * t) * v1 + c2 * e^(-lam1 * i * t) * /bar v1 

! Donate lam = lam_r + i * lam_i as a complex eigenvalue, the respective eigenvalue is v_r + i * v_i
!  consider the general form of the solution:
!      v * e^( (lam_r+ i * lam_i) * t ) 

! e^( i * lam_i * t) = cos( lam_i * t ) + i * sin( lam_i * t ) 

 
!  real part:      e^( lam_r * t) * ( v_r * cos( lam_i * t ) - v_i * sin( lam_i * t )  )
  
!  imaginary part: e^( lam_r * t) * ( v_r * sin( lam_i * t ) + v_i * cos( lam_i * t )  )

! so the the general form of the solution is : c1*Real + c2* Imag

! e^( lam_r * t) * (  c1 * ( v_r * cos( lam_i * t ) - v_i * sin( lam_i * t ) )  +
!                     c2 * ( v_r * sin( lam_i * t ) + v_i * cos( lam_i * t ) )    )

! where c1 and c2 are constants
! simplying 

! e^( lam_r * t) * ( v_r * ( c1 * cos( lam_i * t ) + c2 * sin( lam_i * t ) )  + 
!		     v_i * ( c2 * cos( lam_i * t ) - c1 * sin( lam_i * t ) )    )

! define:
!  sin( phi ) =  -c2 / sqrt( c1^2 + c2^2 )
!  cos( phi ) =   c1 / sqrt( c1^2 + c2^2 )

! further simplified:

! e^( lam_r * t) *  sqrt( c1^2 + c2^2 ) * (  v_r *   cos( lam_i * t  + phi )  
!		    			   - v_i *   sin( lam_i * t  + phi )  )

! (pay attention to the sign(minus -) before v_i)

! to obtain certain form of solution, we can take lam_i * t  + phi = 0 or pi/2 
!  to keep only the part of v_r or v_i



! Here, because the initial state should be of form  (x,0,z,0,vy,0)
! we need to look at the eigenvectors to decide which part to keep, the imaginary 
! or the real one 


! For the first family, pick v_r,  take lam_i * t  + phi = 0, t0 = -phi/lam_i
! we obtain

! epsilon_0 = C1 * v1^r, where C1 is any constant


! similarily, the initial state of the other 2 families are:
! take into account of the form of the eigenvectors, we take t0 = (pi/2-phi) / lam2 (or lam3)

! eta_0 = C2 * v2^i, where C2 is any constant
! xi_0  = C3 * v3^i, where C3 is any constant
 
! take C1(C2,C3) = 1.d-2, very small amplitude as the initial guess


! the solution is the small perturbation around the equilibrium point, add this 
! to the state of the equilibrium point, we obtain the initial guess for the P.O. 

! here is all the 3 pairs of eigenvectors of the chosen equilibrium point
! (  0.0224,  0.0000) (  0.0224, -0.0000) (  0.0000, -0.1926) (  0.0000,  0.1926) (  0.0000,  0.3575) (  0.0000, -0.3575)
! (  0.0000, -0.0239) (  0.0000,  0.0239) ( -0.0645,  0.0000) ( -0.0645, -0.0000) (  0.8424,  0.0000) (  0.8424, -0.0000)
! ( -0.0081,  0.0000) ( -0.0081, -0.0000) (  0.0000, -0.5255) (  0.0000,  0.5255) (  0.0000,  0.3793) (  0.0000, -0.3793)
! (  0.0000,  0.6634) (  0.0000, -0.6634) (  0.2825,  0.0000) (  0.2825, -0.0000) ( -0.0494,  0.0000) ( -0.0494, -0.0000)
! (  0.7075,  0.0000) (  0.7075, -0.0000) (  0.0000, -0.0945) (  0.0000,  0.0945) (  0.0000,  0.1163) (  0.0000, -0.1163)
! (  0.0000, -0.2412) (  0.0000,  0.2412) (  0.7707,  0.0000) (  0.7707, -0.0000) ( -0.0524,  0.0000) ( -0.0524, -0.0000)

! the index of column of vr to be used as the solution of the variational equation
vr_ind = (/1, 4, 6/)


!  ***************** compute the 3 families of PO ***************************
! save the data of po and  its initial state 
fpo = 20;  fpoinst = 21 
open(fpo,file='./dat/po.dat',access ='append',status='replace')  
open(fpoinst,file='./dat/poinst.dat',access ='append',status='replace')  

epsl_po  = 1.d-3 ! the magnitude of the variation 

! initialization of parameters to do numercial continuation
dir = 1
ds = 1.d-3 ! displacement of the continuation
id = 1   ! the type of symmetry,  1,2- xz plane, as in halo orbit, 
e2  = 1.d-10 !err tolerance for numerical continuation
imax = 1 ! time of crossing through y=0 plane for the differential correction

! ifam is the family of periodic orbit to study, in total three families 
do ifam = 1, 3
!ifam = 1 ! test the first one  -- 
!ifam = 2 ! test the second one 

ind = vr_ind(ifam) ! 
y0 = 0.d0 
y0(7:42:7) = 1.d0 !  the initial value of variational matrix if identity stored in y0(7:42)


poinst =  eq +  epsl_po * vr(:, ind) ! the ifam-th column of vr corresponds to the first eigenvalue
print*, 'check', ind,'-th column of vr to use', vr(:, ind) 
print*
!read(*,*) aaa

print*, 'check', ifam,'-th poinst', poinst
print*
read(*,*) aaa

! check gr_lf
!  subroutine gr_lf(t,y,neq,f, beta,cs,sgn)
!print*, 'check variational matrix of lorentz force'


do i = 1, 6, 2 ! Note: only works for kind 2 symmetry: 1 3 5, 
 y0(i) = poinst(i) ! the other three are zero!
enddo


!print*, 'before gr_lf, check, y0,beta, cs,sgn',y0(1:6), beta,cs,sgn


!subroutine gr_lf(t,y,neq,f)

!call gr_lf(0.d0, y0, 42, vf)

!print*, 'check gr_lf, vf', vf(1:6)

!print*, 'check the reshape of the variatioal matrix'
!mat =  reshape(fcopy(7:42), (/6,6/))
!do i = 1, 6
!  write(*,'(6f8.4)') mat(i,:) 
!enddo 



!print*, 'before pofam'
!subroutine pofam(yi,np,imax,e2,dir,ds,id,ntp, deriv, ynew)
call pofam(poinst, npo, imax, e2, dir, ds, id, 2, gr_lf, ynew)
print*, 'PO finisned!'
 
! --- plot the P.O ---
!subroutine plob(y0,t0,tf,ftag, y) 
do i = 1, npo  
  
  write(*,'(7f12.8)') ynew(i,:) 
  
  po0 = ynew(i,2:7)
  tpo = ynew(i,1)
  
  print*, 'TP: i', i, '-th P.O.', tpo
!  read(*, *) aaa
  call plob(po0, 0.d0, tpo, fpo, gr_lf, pof)

enddo   
  
  write(fpo,*) 
end do

! --------- --compute the monodromy matix of the period orbit-----------------

! save the eigenvalues and eigenvectors of monodromy
!fmmat = 22; fegv = 23; fegp = 24 
!open(fmmat,file='./dat/mmat.dat',access ='append',status='replace')  
!open(fegv, file='./dat/mmegv.dat',access ='append',status='replace') 
!open(fegp, file='./dat/mmegp.dat',access ='append',status='replace')
! 
! 
!! --- Compute the Monodromy Matrix and Manifold---
!do i = 1, npo  
!  po0 = ynew(i,2:7)
!  tpo = ynew(i,1)
!  
!!  print*, 'TP: i', i, '-th P.O.', tpo
!!  read(*, *) aaa
!!  call plob(po0, 0.d0, tpo, fpo, gr_lf, pof)
!!  
!! --- Monodramy matrix
!!subroutine monomat(yi,tp, mmat, deriv)
!  print*, 'refined initial state, tp, ynew', ynew
!  call monomat(po0, tpo, mmat, gr_lf)
!  
!! print mmat to file  

!  do j = 1, n
!    write(fmmat,'(6d20.10)') mmat(j,:)
!  enddo    
!  
!! analyze the stability of the monodramy matrix, a big step forward!
!!subroutine eigrg(a,n,isv, wr,wi,vr)
!  call eigrg(mmat,n,1, wr_mm, wi_mm, vr_mm)
!  
!  
!  print*, 'Eigenvalues and eigenvectors, mmat!!!'
!! it seems the real part of the eigenvectors of the monodramy matrix are nearly 1 
!! so the stability is not so straightforward
!! try the power method to see if we can get the dominant eigenvalue and eigenvector
!  fegv = 6; fegp = 6
!  call prt_eigval( n, fegv, wr_mm, wi_mm )
!  call prt_eigvec( n, fegp, wi_mm, vr_mm )
!  
!  write(fegv,*)
!  write(fegp,*)
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
!enddo




!******************************************************************************
stop
end program vp1



















  

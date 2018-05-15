program main_po
!2017-02-16 09:09:36    backup
!  --- and update the main_po.f90 to keep coherent with the updated subroutines.

!  2016-06-26 08:26:55  
! a lot of subroutines have been updated... we need do it also... 
! it is not so convenient to put everything in the main routine, better to seperate them to small but easily-called routines .

!2017-02-14 16:57:35 
! Try to update and refine the routines but give up, better to  finish the draft somehow with the old figures,  the replacement of finer figures  will be done later if we have enough time.

! Ask for the input to specify the case, equilibirum point and the family, in such a way that we can call the same main routine  

!  we now look at case 3, N = [0 0 1], which possesses the full symmetry? 
!  and is perfect for use with 4 symmetric p.o. around 4 symmetric equilibirum points 
!  we will look at the second kind of equilibirum point, the real part of while is also less than 1 in modulus for all beta, quite stable 
!  to start with, we take beta = 1 

!  the third one could also be good to use for formation flying, but the modulus of the real part is always great than 1.


!  take eq3, beta=10, ifam = 3 as example  --- to most studies case 

! 20160406  test asymmetric approach, which do the refinement after 1 revolution w.r.t the initial condition
! 2016018  add the anglectr to control the curvature of the characteristic curve, make sure it is a smooth curve... 
!          after the test, put the conclusion here, whether it is good to do the control on the curvature.
!       
 
use lf_mod ! the system related parameters of lorentz force problem
use po_mod ! archived subroutines to compute families of symmetric p.o.

implicit none

integer, parameter ::  neq = 42, &  ! compute also the variational matrix, 6+36=42
                       npo = 100 ! 60,  52... for x=z, ifam2  !   1 !35 ! NO of orbits in pofam 

real(kind=dp)      ::  pi = 4.d0*datan(1.d0)   
    
! the dimension of the problem is defined in lf_mod, n=6 

! lf_mod Module Variables
! Global :   beta, cs, sgn, eq

integer ::  cs, ind, idic, idfc, asymt,  debug, dsctr, anglectr, idcor, isarc 

! by input from module.... 
!integer :: ieq, 
!real(kind=dp) ::  beta

real(kind=dp) ::   sec, hmax, hmin, e,  tol, prsc, tol_poinc, x0_fxd ! error control
! Local Variables


integer :: i, j,  ivr,  dir, imax, fpo, fpoinst, ifam, fmmegp, fmmegv,  &
           tdir, ipo, isgn, idsgn

real(kind=dp) ::  dlf3(n,n),  & ! differential of vector field of lorentz force
                  wr(n),wi(n), vr(n,n),   &   ! eigenspace of variational matrix  
                  poinst(6), ds, vrchs(6), ynew(npo,8), epsl_po, cj,  & !pofam
                  po0(6), tpo, pof(6), & ! plpo
                  dlf(n,n), mmat(n,n), wr_mm(n), wi_mm(n), vr_mm(n,n), &  ! MM 
                  tp0, csangle_min, ds_max, ds_min 

                 
 character(len=70) :: fnpo, fnpoinst, fnmmegv, fnmmegp              

real(kind=dp), external :: dnrm2 ! from package... 

! TODO 
! ifam=3, stuck at 25th orbit, for the poincare map crossing, t=-2.d-4 fail to finish rk78, 

! ifam =1, stuck at 36th orbit,    poinc: t, y
! 1.8709195802188187E-004  -4.9322318364968389E-021   9.2589610217307418E-003  -1.7629589633198005E-004   1.0554145798777161       -2.8607775325232556E-003
! the orbit is too small.... stop the continuation.

! debug or not 
debug = 0

! automatic step control depending on the number of iteration for Newton Method
dsctr =  1 

! contol the angle between the vector field of two latest consecutive points along the 
! characteristic curve
anglectr = 1

! for the moment, we only use the first one 
isarc = 1 ! along arc-length paramter --- only this one is used so far 
!isarc = 2 ! along period
!isarc = 3 ! along energy

! Todo---- Check really carefully of the subroutine dsdecr.... this is where it is wrong? 
!  there is still something wrong here, in fact if anglectr = 0, it is supposed to obtain regular p.o.s with the possibility of jumping to another family ... 
! I'm not sure why it won't happen.....
! and I'm so regetful that I didn't save the version with the good result ... 

! --- the approaches for correction --- 
!  to determine the index of control and target variables

! ---- for symmetric p.o. : 1-3 ------------ 
idcor  = 1! for symmetric approach with poincare map  --- most used one...
!idcor  = 2 ! symmetric + free time shooting 

!idcor =  3 ! additional constraint to fix one component, y(ind) = x0_fxd, the current used one 

!idcor = 7 ! TODO  symmetric + fixed time shooting ! not-tested yet, could be use for one p.o. to obtain prescribed period 

!  ---- for symmetric p.o. : 4-6------------  
! idcor = 4 ! for asymmetric approach with poincare map--not tested...  
! idcor = 5 ! additional constraint to fix one component, y(ind) = x0_fxd, the current used one 

! idcor = 6 ! the general approach--- discard, all six state components as control variables  

 
! the minimum of the cosine of the angle between two consecutive vectors
! cos(0.1 rad) = 0.995;    cos(0.37 rad) = 0.95  
! csangle_min = 0.995d0  

 csangle_min = 0.95d0  

! bound vaules of stepsize ds
!  With the control of termination added when the arc step goes too smaller... 
!  better not to be too big, with 2.d-2, or 1.d-2, for now is good choice  

ds_max = 5.d-2
ds_min = 1.d-6 ! but this bound is not working for angle control....
    
call init_debug(debug, dsctr, anglectr, csangle_min,  isarc, ds_min, ds_max)

! Initialize the private variables eq1, sgn1 for case 1
call init_lf_mod

! case 1: N=[1 0 0];       2:  N=[0 1 0];    3:  N=[0 0 1]; 
! cs = 1 ! use 2 to test the swap rule is coded correctly--ckd!

print*,    'Please input the case: 1: N=[1 0 0];   2:  N=[0 1 0];  3:  N=[0 0 1]'
read(*,*)  cs 


! the index of the equilibrium points
! 1:  q/m > 0,  x=0,y=0,z= \pm 1
! 2:  q/m < 0,  x = \pm (2 sqrt(9))^(1/3), y² = 2x²,z = 0
! 3:  q/m < 0   x = \pm (1/4/sqrt(2)^(1/3), y = 0,z² = x²  ! this the the case we study now 

! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy
! and  the index of column of vr to be used as the solution of the variational equatio, carefully checked by observing the eigenvectors

! Intialize the state vector of the equilibirum points 
call init_lf(cs) 

!  Instead of pass the value of beta from routine... try read from the screen, which is more flexible    
if(cs == 1) then 
   ! if(ieq==1)  beta = 10.d0 ! TODO 
   if(ieq == 2)  beta = 1.d0
   if(ieq == 3)  beta = 4.d0 
  
 elseif(cs == 3) then 
   if(ieq == 2) beta = 1.d0
 endif  
 
call init_beta(beta)   
!    beta = beta0  ! still use this one, to save time...
    
!    print*, 'Please inpute beta:'
!    read(*,*) beta 
     
print*, 'check, cs, ieq,beta', cs, ieq, beta 
read* 
    
    
! Also compute the unit of distance, time and velocity, which are all function of beta
! but at this moment, we will not use the real quantity of the units. 
call lfunit(beta, runit, tunit, vunit)
    
print*, 'Do you want to change the sign of equilibirum point (Yes = 1) ?'
read(*,*)  isgn
    
! If the equilibirum points with different signs are not symmetric, the periodic orbits will also be different, we have to compute seperately 
     
if(isgn == 1) then 
  idsgn     = mod(ind,3)+1 
  eq(idsgn) = -eq(idsgn)
endif     



! try a different eq : we have four symmetric one, the first two where x=z are already done(x,z), (-x,-z)
! ---- for ieq = 3, --- the third one: x, -z
!eq(3) = -eq(3) 
!print*, 'New eq:', eq
!read* 
!read*
 
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
 

print*, 'Which column to use for initial guess of p.o.?'
read(*,*) ivr 

vrchs = vr(:, ivr)
print*, vrchs 


vrchs = vrchs/dnrm2(6,vrchs,1)
print*, dnrm2(6,vrchs, 1), vrchs

read* 


! commomly used  varaibles
epsl_po = 1.d-3;   ds = 1.d-2; tol  = 1.d-10 ;  prsc = 1.d-10
 
!tol  = 1.d-11;  prsc = 1.d-14

! the initial guess for po, move along the corresponding eigenvector for a small distance epsl_po 

poinst =  eq    +  epsl_po * vrchs !  
print*, 'poinst=', poinst 


if(cs == 1) then
  ! simply-symmtry: P2: x-z plane, ind = 2
  ! TODO : eq1: 0-0-z 
  if(ieq == 3 .or. ieq == 1 )  ind = 2 ! eq3: x-0-z, simple 
  
  if(ieq == 2 )  then 
     idcor = 5  ! eq2: x-y-0, general asymmetric approach and use 
  endif 
  
elseif(cs == 3) then 

  ind  = 1
  
endif 

! for planar symmetric approach, idic = idfc = ind 
if(idcor <= 3) then 
  idic = ind;   idfc = ind
endif 


!--------------------------------------------------------------------------------



!  ***************** compute the 3 families of PO ***************************
! save the data of po and  its initial state 
fpo = 20;  fpoinst = 21;  fmmegv = 23; fmmegp = 24 

! use new name, to avoid overwritten the good ones, and the type of approach selected is used as suffix
!write(fnpo,    fmt='(a,i0,a)') './dat/po',      idcor,'.dat'  
!write(fnpoinst,fmt='(a,i0,a)') './dat/poinst',  idcor,'.dat'
!write(fnmmegv, fmt='(a,i0,a)') './dat/pommegv', idcor,'.dat'
!write(fnmmegp, fmt='(a,i0,a)') './dat/pommegp', idcor,'.dat'

write(fnpo,    fmt='(a,i0,a,i0,a,i0,a)')   './dat/po/case', cs, '/eq', ieq, '/po', idcor,'.dat'  
write(fnpoinst,fmt='(a,i0,a,i0,a,i0,a)')   './dat/po/case', cs, '/eq', ieq, '/poinst',  idcor,'.dat'
write(fnmmegv, fmt='(a,i0,a,i0,a,i0,a)')   './dat/po/case', cs, '/eq', ieq, '/pommegv', idcor,'.dat'
write(fnmmegp, fmt='(a,i0,a,i0,a,i0,a)')   './dat/po/case', cs, '/eq', ieq, '/pommegp', idcor,'.dat' 

print*, fnpo, fnpoinst, fnmmegv, fnmmegp
!read(*,*)

open(fpo,file=fnpo, access ='append',status='replace')
write(fpo,*)  '# t	(x y z vx vy vz)	cj'

open(fpoinst,file=fnpoinst, access ='append',status='replace')
write(fpoinst,*) '# TP	I.C.(x y z vx vy vz)	CJ	DCJ	ds'

open(fmmegv, file=fnmmegv, access  ='append',status='replace') 
open(fmmegp, file=fnmmegp, access ='append',status='replace')



! Initialization of parameters to do numercial continuation
dir  = 1
tdir = 1 ! integrate forward

! finally take 1.d-9 as the error control
!tol  = 1.d-10 ! err tolerance of target f to terminate Newton Method, from the Gerard's book, it should be 1.d-16?
!prsc = 1.d-10 ! For lf problem,  1.d-11 is too small for the second family, ok for the other 2 families


if( mod(ivr,2) == 1) then 
  ifam = ivr/2+1
else 
  ifam = ivr/2
endif   


print*, 'check wi to determine ifam'
print*, 'wi=',  wi
print*, 'ifam=', ifam  
! Use the period as initial guess, TP = 2*pi/wi, wi is the pure imaginary eigenvalue
tp0 =  2*pi/dabs(wi(2*ifam)) ! the maximum time for the first return to poincare section
print*, 'tp0=', tp0, 'wi=', wi(2*ifam)
read*


! -------------------------------------------------------------------------------------------------
if (idcor <= 3 )  then !  for  idcor = 1, 2, 3, TODO 3
! ---- symmetric p.o. ----
  call init_symtpo2(idic, idfc, idcor, ind)    
  sec  = 0.d0;  imax = 1 
  
  
else! if ! idcor = 4, 5, 6. TODO 4 and 6, availale : 5
! ----- asymmetric p.o.

! --- Approach 0, by general asymmetric way -- discard at the moment
!  asymt = 1  ! still try y = 0 plane            !alternative: x=x0 plane
  
!  ind = 2 !? TODO this is the one to specify... 
!  sec = 0.d0 ! y = 0 still

   print*, 'General approach for p.o. computation! Input ind,  add constraint:  sec = poinst(ind)'
   read(*,*) ind 
   
   sec = eq(ind)
  
!  subroutine init_asymtpo2(idcor0, ind0, x0_fxd0)
  call init_asymtpo2(idcor, ind, sec)   

endif 


!  subroutine init_poinc(ind0, sec0, imax0,  hmin0, hmax0, e0)
if (idcor == 1 .or. idcor == 4 ) then 
  ! Step size and error control for rk78 integration 
  hmin = 1.d-10
  hmax = 5.d-1
  e    = 1.d-14  ! 1.d-13 is from Gerard's book 
  tol_poinc = 1.d-14 
  
  call init_poinc(ind, sec, tol_poinc, imax, hmin, hmax, e)
endif 


call init_errctr(tol, prsc)
call init_writedx 


!  subroutine pofam(yi, tp0, npo, dir,ds, fpo, ynew, ipo, deriv, gr_cj)
!ds = 0.01d0
call pofam(poinst, tp0, npo, dir,ds, fpoinst, ynew, ipo, gr_lf, gr_cjlf)
 
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
end program main_po



















  

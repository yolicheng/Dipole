program poformfly !  take eq3, beta=10, start with 2 p.o. of the same period
! look at the eigenvalues of the variational matrix of the equilibrium point, and select two with similar value 
! and do the differential correction to obtain p.o. with exactly the same period. check the evolution 

! Question: 
!  do we need to integrate the two orbits at the same time? no sure, but first, we need to compute it.

use lf_mod ! the system related parameters of lorentz force problem
use po_mod 

implicit none

! 	Module lf_mod based Parameter
!  dp, n, 		! the dimension of the problem is defined in lf_mod, n=6 
!  runit, vunit, tunit  ! unit of the fundamental quantities

integer, parameter  ::  neq = 6, &  ! we only need to focus on the position + velocity  
                        npo = 2! number of p.o.s to be integated simultaneously
  
real(kind=dp)       :: pi = 4.d0*datan(1.d0)   
   

integer ::        cs, debug 
real(kind=dp) ::  beta 

! Local Variables
real(kind=dp) ::  eqi(npo, 6), eqn(npo * 6), devn(npo*6), devi(npo,6), & ! the 4 symmetric equilibrium points 
                  tic, toc, devpos, devvel, mean, sd , & ! for normal distribution 
                  dlf(n,n), wr(n),wi(n), vr(n,n),  &    ! eigenspace of variational matrix 
                  cj, &  ! the energy 
                  tp0, poinstn(npo*n), poinsti(npo, n), poinst1(n), poinst2(n), poinst0(npo, n), & ! initial state 
                  sec, tol_poinc,  hmin, hmax, e, tf, yf(42), hminim  ! for poincare map to obtain the half evolution point
  
                  
integer :: i, j, tdir, isdev, iseed(4), idist,  fob, fdrv, fh, fegp, fegv, feq,  ind, imax, ispc, &
           id1, id2, symt  ! the index of the non-zero components in eq 
           
character(len=70) :: fnob, fndrv, fnh, fnegp, fnegv, fneq

!	subroutine from packages declared external here
real(kind=dp), external :: dnrm2 ! from BLAS-LEVEL1
external :: dlarnv ! from LAPACK to generate normal distribution random data 

! add tracking error or not  
isdev = 0 
idist = 3  ! normal [0, 1] 

! to obtain 2-sigma deviation, we set the standard deviation to be the maximum vaule/2
devpos =  1.d-1 / 2.d0 !  km     	100 m
devvel =  3.d-6 / 2.d0 !  km/s 		3 mm/s
  
  
! t=1.d4 ~= 800 days, for case with tracking error, the time is too long    
tf = 1.d3 ! integration time: adimensional

! debug or not 
debug =  0

! Initialize the private variables eq1, sgn1 for  case 1
call init_lf_mod

print*,    'Please input the case: 1: N=[1 0 0];   2:  N=[0 1 0];  3:  N=[0 0 1]'
read(*,*)  cs 
call init_lf(cs)
! ********************** case == 1 ******************* 
! case 1: N=[1 0 0]; 2:  N=[0 1 0]; 3:  N=[0 0 1]; 
! cs = 1 ! use 2 to test the swap rule is coded correctly--ckd!

! the index of the equilibrium points
! 1:  q/m > 0,  x=0,y=0,z= \pm 1
! 2:  q/m < 0,  x = \pm (2 sqrt(9))^(1/3), y² = 2x²,z = 0
! 3:  q/m < 0   x = \pm (1/4/sqrt(2)^(1/3), y = 0,z² = x²  ! this the the case we study now 

! ---------- ieq = 3 ------------ 2 p.o.s
!ieq = 3 ! 3, x,0,z is the case that we study currently
!id1 = 1;  id2 = 3; symt = 2
! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy
!beta = 10.d0   !  6 dimension center manifold 
! ------------------------------------------------

! ---------- ieq = 2 ------------ 4 p.o.s
!ieq = 2 ! 2, x,y,0, 0,0,0 , asymmetric case
!id1 = 1; id2 = 2 
!! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!!                     the angular velocity of the rotaion of the deputy
!beta = 1.d0   !  6 dimension center manifold 
! ------------------------------------------------ 


! ***************** case = 3 , N =[0 0 1] 
! -----  ieq = 2,   for +t direction, x-y plane symmetry, z-axis, origin 
! (x,y, z, vx,vy, vz ) ---> (x,y, -z, vx,vy, -vz ) 
!                           (-x,-y, z, -vx,-vy, vz) 
 cs = 3
 
! --- ieq =2 , beta = 1.d0 
!ieq = 2 ! 0, y,z, 0, 0,0 
id1 = 2; id2 = 3;    symt = 3 
beta  = 1.d0 

call init_beta(beta)
! Provided the case, beta and ieq,  beta, cs, sgn and eq are initialized by subroutine init_lf from module lf_mod

print*, 'check the assignment with module' !--ckd
print*, 'beta, cs, sgn, eq', beta, cs, sgn, eq
read*


! assignment of the equilibrium points around which we are going to compute the p.o.s

do i = 1, npo 
  eqi(i, : ) = eq
enddo 

! 1st eq: x = x_eq,  z = x 
! the second eq, x=-x_eq, z = x 
if(npo .ge. 2)  eqi(2, (/id1,id2/) )  = -eq( (/id1,id2/)  )


! comment these two, or check the eigenvalues of the Monodromy matrix to see if it is possible to get p.o. with the same period 

if (npo == 4 ) then 
! the third:  x = x_eq,  z = -x 
eqi(3, id2) = -eq(id2)

! the fourth eq, x= -x_eq, z = -x 
eqi(4, id1) = -eq(id1)
endif  

! ------------ to save the equilibirum points and eigenspace ----------------
write(fneq,  fmt='(a,i0,a)') './dat/po',npo, '_eq.dat'
write(fnegv, fmt='(a,i0,a)') './dat/po',npo, '_egv.dat'
write(fnegp, fmt='(a,i0,a)') './dat/po',npo, '_egp.dat'
 
print*,  fnegv, fnegp
read(*,*)

fegv = 30; fegp = 31; feq = 35
open(fegv, file=fnegv, access  ='append',status='replace') 
open(fegp, file=fnegp, access  ='append',status='replace')
open(feq,  file=fneq, access  ='append',status='replace')
 

print*, 'The points:'
do i = 1, npo 
  eqn((i-1)*6+1 : i*6) = eqi(i, :) 
  print*,  eqi(i, :)
  
  write(feq, '(6f20.14)') eqi(i,:) 
  read*
! check the eigenvaules and eigenvectors, possible similar values in eigenvalues that could lead to p.o. with similar period
   
! Jacobi matrix of the lorentz force with respective to the state 
  call dflrtz( eqi(i,:), dlf)
  call gr_cjlf(eqi(i,:), cj)
  print*,'check energy!, cj, ieq, eq: ', cj, ieq, eqi(i, :)
  read* 
  
  do j = 1, n
   write(*,'(6f8.4)') dlf(j,:) 
  enddo
  
  ! compute the eigenvalues and eigenvectors  of dlf
  !subroutine eigrg(a,n,isv, wr,wi,vr)
  call eigrg(dlf,n,1, wr,wi, vr)
  read* 
  
  ! print to files to save 
  call prt_eigval( n, fegv, wr, wi )
  call prt_eigvec( n, fegp, wi, vr )
  
  write(fegp,*) ! add a blank line to seperate eigenvector matrix
  
enddo 
print* ; ! read* 

!stop 

! the filename to save the useful data, the integration, relative position+velocity, angular momentum
fob = 20;  fdrv = 21;  fh = 22 

! use new name, to avoid overwritten the good ones, and the type of approach selected is used as suffix
write(fnob,   fmt='(a,i0,a,i0,a)') './dat/po',npo,'_ob', isdev,'.dat'  
write(fndrv,  fmt='(a,i0,a,i0,a)') './dat/po',npo,'_drv',isdev,'.dat'
write(fnh,    fmt='(a,i0,a,i0,a)') './dat/po',npo,'_h',  isdev,'.dat'


print*, fnob, fndrv, fnh; print*!ck
read(*,*)

open(fob,file=fnob, access ='append',status='replace')
write(fob,*)  '# Every 4 rows(p1,p2,p3,p4):  t	(x y z vx vy vz)	cj'

open(fdrv,file=fndrv, access ='append',status='replace')
write(fdrv,*) '# Every 6 rows(r12,r13,r14,r23,r24,r34): t	(dx dy dz dvx dvy dvz)'

open(fh, file=fnh, access  ='append',status='replace')
write(fh,*) '# t	angular momentum(hx hy hz)	< r14, h >' 

print*, 'tf = ', tf


! ************************** case 1 -------------------
! ieq = 3 
! Given the initial guess found by the continuation method, do we need to specify the period again by fixed-time shooting method?
!  My guess: I do not think there is any necessarity to recompute a p.o., use this one, and apply the symmetry to obtain the other one 
!            and integrate the 2 initial point simultaneously 

! ---------------  ifam = 1  --------------------


!---------------------------------------------------------------------




! ---------------  ifam = 2  --------------------
!  30-th p.o.  the size is relatively big  
!  0.42841573866374E+01  0.55905421061878E+00  0.00000000000000E+00  0.55534599921232E+00  0.00000000000000E+00 -0.10440740360204E-02  0.00000000000000E+00  0.18897977640739E+01  0.22204460492503E-15  0.00128000    1

! because the non-zero components of the initial state is (x,z,vy), so the one symmetric w.r.t x-z plane should be (-x, -z, -vy)
if( ieq  == 3 .and. cs == 1) then
   tp0 = 0.42841573866374d1
   poinst1 =  (/0.55905421061878d0,  0.d0,  0.55534599921232d0, 0.d0, -0.10440740360204d-2,  0.d0/)  
end if 
!-----------------------------------------------------------------------



! ----- ieq = 2 : the asymmetric case 

!0.56441453174155E+01  0.50436191647430E+00  0.70793627753308E+00  0.18379718925006E-01  0.19458540327179E-02  0.82046869444425E-02 -0.65200199829016E-03  0.22889666528315E+01  0.00000000000000E+00  0.00400000    1

if( ieq ==  2 .and. cs == 1 ) then
  tp0 = 0.56441453174155d1
  poinst1 =  (/ 0.50436191647430d0,  0.70793627753308d0,  0.18379718925006d-1,  & 
             0.19458540327179d-2,   0.82046869444425d-2, -0.65200199829016d-3 /)     
end if


! --------------------------------case 3 ----------------------------
if( ieq ==  2 .and. cs == 3 ) then
! 19-th no good....  0.17725193580209E+02  0.00000000000000E+00  0.99388826240297E+00  0.64342693948668E+00  0.25084528595441E-02  0.00000000000000E+00  0.00000000000000E+00 -0.16043428266082E+01  0.44408920985006E-15  0.01600000    1 
!  tp0 = 0.17725193580209d2
!  poinst1 =  (/ 0.d0,   0.99388826240297d0,  0.64342693948668d0, 0.25084528595441d-2,  & 
!               0.d0,   0.d0 /)   

! 4-th  0.17771483267230E+02  0.00000000000000E+00  0.10277155739168E+01  0.72458832059739E+00  0.18799848553220E-05  0.00000000000000E+00  0.00000000000000E+00 -0.15874190770946E+01  0.00000000000000E+00  0.00100000    1 

   tp0= 0.17771483267230d2 
   poinst1 = (/ 0.d0,  0.10277155739168d1,  0.72458832059739d0,  0.18799848553220d-5, 0.d0, 0.d0 /)
   
! 25-th -- not so good 
!  0.17430209063032E+02  0.00000000000000E+00  0.92659678196894E+00  0.52299380461568E+00  0.26349422738140E-01  0.00000000000000E+00  0.00000000000000E+00 -0.16997653751295E+01  0.15543122344752E-14  0.06400000    1   
!  tp0 = 0.17430209063032d2
!  poinst1 = (/ 0.d0,  0.92659678196894d0,  0.52299380461568d0,  0.26349422738140d-1 , 0.d0, 0.d0 /)

! 40-th  
!  0.11171321866461E+02  0.00000000000000E+00  0.93193400748506E+00  0.67855417943872E+00  0.51738711903930E+00  0.00000000000000E+00  0.00000000000000E+00 -0.18619444896883E+01 -0.22204460492503E-15  0.25600000   -1
  tp0 = 0.11171321866461d2 
  poinst1 = (/ 0.d0,  0.93193400748506d0, 0.67855417943872d0,  0.51738711903930d0 , 0.d0, 0.d0 /)
  
  tp0 = tp0 / 4.d0
  
  
  ! 40-th,  i=39,  -- the one that is parallel to x-y plane ! --- not good, in fact are two orbits.... almost overlap for 2 
!    0.18432164858035E+01  0.00000000000000E+00  0.36664748068655E+00  0.10117245290528E+01  0.10317592130296E+01  0.00000000000000E+00  0.00000000000000E+00 -0.23038657110273E+01  0.44408920985006E-15  0.25600000    1
  tp0 = 0.18432164858035d1
  poinst1 = (/ 0.d0,  0.36664748068655d0,  0.10117245290528d1,  0.10317592130296d1 , 0.d0,  0.d0 /)
  
  
  ! 46-th  i=45, the one before almost perpendicular to x-y plane-- turns out it is not almost perpendicular, but the size is small
!    0.35864311809759E+01  0.00000000000000E+00 -0.55680992271445E+00  0.11417028016985E+01  0.44424859433677E+00  0.00000000000000E+00  0.00000000000000E+00 -0.18033810944970E+01 -0.44408920985006E-15  0.25600000   -1

!  tp0 = 0.35864311809759d1
!  poinst1 = (/ 0.d0,  -0.55680992271445d0,  0.11417028016985d1,  0.44424859433677d0,  0.d0,  0.d0 /) 

! 47th, i = 46, almost perpendicular to x-y plane -- no big different with 46-th....
!  0.38029658821541E+01  0.00000000000000E+00 -0.70640514618501E+00  0.12077964137697E+01  0.26294530968888E+00  0.00000000000000E+00  0.00000000000000E+00 -0.18922397316259E+01  0.00000000000000E+00  0.25600000   -1
  tp0 = 0.38029658821541d1
  poinst1 = (/ 0.d0, -0.70640514618501d0,  0.12077964137697d1,  0.26294530968888d0, 0.d0, 0.d0 /)
  
! 32-th, i=31, a big one in size, just for testing   
!    0.34858334376319E+01  0.00000000000000E+00  0.95416926423625E+00  0.54917304800717E+00  0.34667273261016E+00  0.00000000000000E+00  0.00000000000000E+00 -0.17863900149851E+01  0.11102230246252E-14  0.12800000    1
  tp0 = 0.34858334376319d1 
  poinst1 = (/ 0.d0,  0.95416926423625d0,  0.54917304800717d0,  0.34667273261016d0, 0.d0, 0.d0 /)
end if

! the almost parallel one... 20160628 4iter-9th
! 0.18424900935116E+01  0.00000000000000E+00  0.34418186976432E+00  0.10183465512963E+01  0.10376341020454E+01  0.00000000000000E+00  0.00000000000000E+00 -0.23044605641102E+01 -0.13322676295502E-14  0.00000400    1
tp0 = 0.18424900935116d1
poinst1 = (/ 0.d0,  0.34418186976432d0,  0.10183465512963d1,  0.10376341020454d1, 0.d0, 0.d0 /)
  
! assign the first row 
poinsti(1, :) = poinst1  

!  --------------------- for case 1,N=[0,0,1] -------------------------------------------------------
! for the four symmetric orbits, they are the original soultion + three imapges for symmetric 4 point 
! --- p.o. 1 :  Original (x,   y,  z,  vx,  vy,  vz)  
! --- p.o. 2 :  Image 1  (-x, -y, -z, -vx, -vy, -vz)     -- symmetric w.r.t Origin      ----  '+' in time 

! --- p.o. 3 :  Image 2  ( x, -y,  z, -vx,  vy, -vz)     -- symmetric w.r.t  x-z plane  ----  '-' in time 
! --- p.o. 4 :  Image 3  (-x,  y, -z,  vx, -vy,  vz)     -- symmetric w.r.t  y-axis     ----  '-' in time 
 
! p.o. 3 and 4 are symmetric w.r.t the Origin in positive time direction 
 

! for ieq=3, the symmetric case, the p.o. is computed by half Origin + half image 2,  p.o.1 and p.o.3 is the same one 
! the problem is for 4 symmetric equilibrium point, with (x=z) the four solutions are actually 2 p.o.s around (x,0,z) and (-x,0,-z)

! While for (x=-z), the values of eigenvalues are different from the case where (x=z),  and the continuation fails because of some reason I donn't know...

! for ieq = 2, the asymmetric case,    
! given an p.o. ....  we only need to find the third image.... 


! add also the second group of p.o.s, first try to get p.o. with the same period. 

do i = 1, npo, 2 

! try the one with the same poinst1, but after exactly half revolution....
!subroutine plob(y0,t0,tf, n, tdir, ftag, ispl, deriv, gr_cj,  y), the point is, using plob to obtain the point at the half revolution is not so precise...
! poinc is better, if we ask the precision for the crossing to be 
  poinst1 = poinsti(i, :)
   
   ! --- for cs = 1, ieq=3---- half period difference in phase 
!  if( cs == 1 ) then 
    call plob(poinst1, 0.d0, tp0/2.d0, 6, 1,  6, 0, gr_lf,  gr_cjlf,  poinst2)! discard, not precise enough 
    print*, 'finish half period!', 'tf=', tp0/4.d0, 'yf=',  poinst2
    read* 

! try cs=3, ieq=2
!  elseif (cs == 3) then  
!    
!    poinst2 =  poinst1
!    
!  endif 
! use poincare map to make sure we find exactly the point  after half period...., but it only applies to symmetric p.o., 
! for unsymetric one, we have to use plob .... 

! ------------ poincare map for half period -----------------------   discard because it is inappropriate for asymmetric orbit
  
!  ind = 2 
!  sec = 0.d0
!  tol_poinc = 1.d-16
!  
!  imax = 1
!  hmin = 1.d-10
!  hmax = 1.d0
!  e    = 1.d-16
!  
!! -- initialize  for the second p.o.---------------  
!!  subroutine init_poinc(ind, sec, tol_poinc, imax,  hmin, hmax, e)
!  call  init_poinc(ind, sec, tol_poinc, imax, hmin, hmax, e)
!  
!!  subroutine poinc(yi,imax, tmax, neq, tf,yf, hminim, ispc, deriv, gr_cj) 
!  call poinc(-poinst1, 1, tp0, 42, tf,yf, hminim, ispc, gr_lf, gr_cjlf) 
!   
!  print*, 'finish poinc! The point after half revolution:'
!  print*, 'tf=', tf, 'yf=', yf(1:6)
!  read* 
!  
!! either with [-x, y, -z,  -vx, vy, -vz] or -poinst1, the step size for simultaneous integration is extremely small, < 1.d-6. so discard this approach.
!  poinst2 = yf(1:6)
! ------------ poincare map for half period -----------------------  

! the image 3 and poinst2 which is symmetric w.r.t the Origin  
  poinsti(i+1,:) = -poinst2 
  
  
  ! images 3 .... ! assign the third row ---- to keep in the same phase  
  if(i == 1 .and. npo == 4) then 
    if (cs == 1) then ! x-z plane, -t direction  -->(x, -y, z, -vx, vy, -vz)
      poinsti(3, (/1,3,5/) ) =  poinst1( (/1,3,5/) ) 
      poinsti(3, (/2,4,6/) ) = -poinst1( (/2,4,6/) ) 
      
    elseif (cs == 3) then ! x-y plane, t direction -->(x, y, -z, vx, vy, -vz)   
      poinsti(3, : ) =  poinst1
      poinsti(3, (/symt, symt+3/) ) =  -poinst1( (/symt, symt+3/) )
    endif
       
  endif 
 
  print*, 'I.C. for the second p.o.:'
  print*,  poinst2
  print*; read* 
enddo 


print* ,'check I.C. of', npo, 'p.o.s'


do i = 1, npo, 1
  poinstn((i-1)*6+1 : i*6) = poinsti(i, :) 
  print*,  poinsti(i, :)
enddo
print*; read* 


! ----- save the one without tracking error ------------ 
poinst0 = poinsti 


!! ********************* this is where we could add the noise   ********************
! we need to generate np * 6 random data 

if (isdev == 1) then 
! if we add the Gaussian random 2-sigma deviation as the tracking error on the initial point
! required magnitude of deviation in position: 100m, velocity: 3 mm/s
!  To obtain random data within 2-sigma deviation, we could use sigma = Max_dev / 2 

! we have to transform them to adimensional unit for computation 

! use Lapack Driver routine DLARNV to generate normal distribution random data
!  SUBROUTINE DLARNV( IDIST, ISEED, N, X )
! different seed generate different random data, but the same seed always produce the same one...

  iseed = (/29, 18, 3, 11/) ! use this one for poformfly
!  iseed = (/0, 0, 0, 1/) !  use this one for eq4
   
  call dlarnv( idist, iseed, npo*6,  devn)
  devi = reshape( devn, (/npo, 6/))

  print*, 'check the normal distribution deviation'
  do i = 1, npo 
    print*, devi(i, :)
  enddo 
  print* ; ! read*
  
  ! Check mean and standard deviation
  mean = SUM(devn)/ npo/ 6
  sd = DSQRT(SUM((devn - mean)**2)/ npo / 6)
  
  WRITE(*, "(A,F18.10)") "Mean = ", mean
  WRITE(*, "(A,F18.10)") "Standard Deviation = ", sd
  
  print*, 'Magnitude of deviation: devpos=', devpos, 'km,  devvel=', devvel, 'km/s'
  devpos = devpos / runit
  devvel = devvel / vunit
  
  print*, 'Adimensional Magnitude of deviation: devpos=', devpos, ',  devvel=', devvel
!  read* 

  ! add the deviation  -- for position 
  poinsti(:, 1:3)  = poinsti(:, 1:3)  + devpos * devi(:, 1:3) 
  ! -- for velocity
  poinsti(:, 4:6)  = poinsti(:, 4:6)  + devvel * devi(:, 4:6)  
  
endif


tf = 3.d0 * tp0 
!tf = 1.d3 
call cpu_time(tic)
!subroutine plob_eq_n(y, np, t0,tf, tdir,fob, fdrv, fh,   deriv_n, gr_cj) 
call plob_eq_n(poinstn, npo, 0.d0, tf, 1, fob, fdrv, fh, gr_lf_n, gr_cjlf) 

call cpu_time(toc)

print*, 'Elapsed time=', toc-tic


stop
end 





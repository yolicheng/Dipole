program curve_lf2
! deal with 6D phase space. check if the routine works.
! old version with an arbitraray torus is saved in main_curve_lf_pc.f90, and probably discard forever.

! GIVE UP THE POINCARE MAP APPROACH, AND TURN TO TIME_T2 MAP APPROACH.  2017-05-26 22:07:37 

! that takes the linear flow around an elliptic periodic orbit, more specific, centre X saddle (6D phase space) and centre (4D phase space), 
! as the initial guess of an invariant curve

! the eigenvalue (rho) of centre part (pure imaginary) serves as the rotation number (in radian), where T is the period.
! and the two frequencies are close to:  rho/T, 2pi/T. 
! TODO we have to check if we need to multiply the eigenvalue by 2pi as the rotation number 


! for Poincare map approach, we need one more step to integrate this curve for the next crossing with the Poincare section. 


! ** NOTE ** 
! DO not forget the initialization of the model parameters:   beta, sgn, ndim  

!  Study the planar norm case of lorentz force problem, take an elliptic periodic orbit, of type center X center
!     around it, we have plenty 2D tori, the one around x0 =0.35.  
 
!    1. take an unknown periodic orbits, compute the MM, and associated eigenvalues and eigenvectors
!    2. take the linear flow as the initial guess for invariant curve \phi: x0 + epsilon * (vr*cos(xi) - vi*sin(xi))
!       -- epsilon is the distance from the periodic orbit, vr+ivi is the associated eigenvector, 
!       -- x0 is an initial condition of the p.o., take the one on Poincare section just for convenience. 
!    3. Integrate this curve for the next crossing with the  Poincare section to obtain our initial curve \varphi 
!    5. do a simple Fourier analysis to compute the Fourier coefficients as the initial guess \bm X_0
!    6. Impose the invariance equations, and specify the energy level, to compute the invariant curve 
!    7. Globalize the curve to an invariant torus  -- actually, this is not well done, I just integrate for enough long time...


!  --- data saved in ./dat/curve_lf/  subfolder 
!  1. torus.dat            -- the original orbit of the 2D torus  ! not used at the moment 
!  2. curve_time.dat       -- the invariant curve obtained by Time map
!  3. curve_pc_org. dat    -- initial poincare curve obtained by integration from Time curve
!  4. curve_pc. dat        -- refined curve_approx
!  5. tori_pc.dat          -- the global tori by the Poincare curve 
!  4. fun.dat              -- approxiamted Fourier representation \varphi(xi+rho), the check \varphi(xi+rho), P( \varphi(xi) )
!  5. fcs.dat              -- fourier coefficients, For each coordinate, three parts: c0 // ck // sk

use dp_mod
use pi_mod
use sort_mod
use curve_mod
use poinc_mod 
use lf_mod 

implicit none

integer, parameter  ::  np    = 301    ! number of Poincare maps
 
integer, parameter  ::  ndim  = 6, &   ! dimension of phase sapce for spatial lf
                        n0    = ndim-2 ! for Poincare map space 
                        
integer, parameter  ::  nvar  = ndim*(ndim+1)                        

real(kind=dp) ::  beta0 

! --  ---- declare the size for Fourier coefficients of the original  Fourier analysis   ----
! how many Fourier modes ???, Alex suggests of power 2 
! Gerard said 60 something is enough ....

! use 2 to check... 
! TODO: apply Alex's strategy: claim a big enough array at first, then check the maximum norm 
!       of the last half, if it is greater than a given threshold, which should be a little bit smaller than 
!       the precision we ask for the Newton method,  we increase the size by 2

!       if it is smaller, then we decrease the size by half 
!       the arrays that store ck, sk are claimed dynamical arrays.                  

!       since they are just initial seed computed from a linear interpolated curve, not accurate enough as 
!       the estimate of the precision of the refined curve 
!       so we check again the refined Fourier coefficients by looking at the Maximum norm of the last half.

! TODO: JM's idea, paper  Physic D. P12. Section 2.2.6. Error estimate
!       increase the value of nf, until the error in Invariance equation is under a given tolerance. 1.d-10
!       where he uses multiple shooting with intermediate node M=50 


integer, parameter  ::  nf_init   =  64 ! 512 ! 128 !  !128 ! 64 !  !32 ! 32   !64 ! 32  !8 

! -- fourier coefs for Fourier analysis for the initial X_0  -- 

! We declare c00 allocatable since we cannot specify its dimension in the module as public array with n unknown.
real(kind=dp) :: csf0, ctau0 
real(kind=dp), allocatable,  dimension(:)   :: csf, sif, ctau, stau, ct, st! fourier coefs of 1 component
!real(kind=dp), allocatable                  :: ck0(:,:), sk0(:,:), c00(:) 


! -- fourier coefs for Refinement  -- 
! the number of Fourier modes we are going to use will be specified 
! after we check the maximum norm of ck and sk, all the following allocatable arrays have nf-columns  

integer           ::  nf, nf0,niter
                    
!  ------ Local Variables ----------- 
real(kind=dp) :: p00, tmax0, xmax0,    pvf(ndim),  cj0,     & ! torus
                 tp0, poinst(ndim),  mmat(ndim, ndim), wr(ndim), wi(ndim), vr(ndim, ndim), & ! MM 
                 v1(ndim), v2(ndim),  pv0(ndim),  xi, dxi, epsl,    & ! linear flow 
                 hminim, tpc, pvpc(ndim),  xpc_all(np,ndim),  & ! poinc 
                 rho,    pv_in(np),  & ! rotnum
                 tol, tol_err, pv_curv(n0),   &   ! Netwon method
                 pvf_rho(ndim),  pv_curv_rho(n0) 
                  
real(kind=dp) :: tpc_refn(np),  tau, T2                                  
                                                      
integer ::  i, j,   ivr, i_eig, ivr1, ivr2,  & 
            ind0, dir0, imax0, tdir, ispc, isv2,  & ! poinc 
            iscurve,  isfourier, & !  control which subroutine to call
            nitmax, isref, opt 
  
integer           ::  ftorus_time, fcurve_time, ftorus_pc, fcurve_pc_org, fcurve_pc,  & 
                      ftorus_conct, fcurve_pc_rho, fcurve_pc_fc           
character(len=70) ::  fntorus_time, fncurve_time, fntorus_pc,fncurve_pc_org, fncurve_pc, &
                      fntorus_conct, fncurve_pc_rho, fncurve_pc_fc   


!---------- interpolation --------------                   
real(kind=dp)     ::  pt(np, n0), rho2, theta(np), darg, arg, ptnew(np, ndim), &
                      ptnew1(ndim), pv(ndim), t, tf 
integer           ::  ind_sorted(np), dircurve, fcurve_interpol, ncrs
character(len=70) ::  fncurve_interpol


!---------- continuation --------------                   
! Fourier coefficients 
real(kind=dp)              ::  c0pre(n0), c0_new(n0)  
real(kind=dp), allocatable, dimension(:,:)  ::  ckpre, skpre, ck_new, sk_new  
   
!  continuation
integer, parameter :: ntori = 100 
integer       ::  k, restart,  nfmn,  nf_pre, itori 
real(kind=dp) ::  cj, cj_pre, dcj, dcj_max, dcj_min, nds    ! step control for continuaiton 
!---------- continuation --------------                   
                  
                  
                  
external :: PoincMap_lf ! general 2d map for lf  

real(kind=dp) :: four_seri,  datan_2pi  !funcs dlange, 
!real(kind=dp) :: MaxNorm   

call init_lf

!beta0 = 2.d0 

print*, 'Please input beta (q/m):'
read(*,*) beta0 
call init_beta(beta0) 

! which approach we want to take 
!!print*, 'Choose approach: (I prefer 3)'  
!!print*, '1: no constraint;  '
!!print*, '2: Avoid phase shift indetermination by adding one more equation'
!!print*, '3: Avoid phase shift indetermination by removing one cofficient C1_J = 0 or S1_J = 0'
!! 
!!read(*,*) opt 

opt = 3

! -- Initialization for the module  curve_mod, use poincare map approach, n=ndim-2
nitmax  = 20    
tol     = 1.d-10
tol_err = 1.d-10 

call init_curve(n0, nitmax, tol, tol_err, opt) 
  
  
iscurve    = 1   ! 1: Do poincare map of the torus to get the curve; other: curve already computed. 
isfourier  = 1   ! 1: Do fourier analysis


! files to save the original torus + invariant curve + interpolated one + approximated one 
ftorus_time    = 21;  fntorus_time    =  './dat/curve_lf/torus_time.dat'
fcurve_time    = 25;  fncurve_time    =  './dat/curve_lf/curve_time.dat'

ftorus_conct    = 22;  fntorus_conct    =  './dat/curve_lf/torus_conct.dat'
ftorus_pc      = 23;  fntorus_pc      =  './dat/curve_lf/torus_pc.dat'
fcurve_pc_org  = 26;  fncurve_pc_org  =  './dat/curve_lf/curve_pc_org.dat'
fcurve_pc      = 27;  fncurve_pc      =  './dat/curve_lf/curve_pc.dat'
fcurve_pc_rho  = 28;  fncurve_pc_rho  =  './dat/curve_lf/curve_pc_rho.dat'
fcurve_pc_fc   = 29;  fncurve_pc_fc   =  './dat/curve_lf/curve_pc_fc.dat'

open(ftorus_time,  file = fntorus_time, access ='append', status = 'replace')
write(ftorus_time,*)  ' # Torus from time curve:  xi   (x y z vx vy vz)   cj'
  
open(ftorus_pc,  file = fntorus_pc, access ='append', status = 'replace')
write(ftorus_pc,*)  ' # Torus from Poincare curve:  xi   (x y z vx vy vz)   cj  '

open(ftorus_conct,  file = fntorus_conct, access ='append', status = 'replace')
write(ftorus_conct,*)  ' # Torus from Poincare curve:  xi   (x y z vx vy vz)   cj  '
  
open(fcurve_pc,  file = fncurve_pc, access ='append', status = 'replace')
write(fcurve_pc,*)  ' # Refined Poinc curve:  xi   (x y z vx vy vz)   cj  tpc'

open(fcurve_pc_rho,  file = fncurve_pc_rho, access ='append', status = 'replace')
write(fcurve_pc_rho,*)  ' # The curve for (xi+rho):  xi   (x y z vx vy vz)_rho  (x y z vx vy vz)_PC   cj  '
  
    
! -- if we have computed the initial curve.  
if(iscurve == 1) then 
  open(fcurve_time, file = fncurve_time, access ='append', status = 'replace')
  write(fcurve_time,*)  ' # Initial Time + Poinc_org curve:  xi   (x y z vx vy vz)   cj'

  open(fcurve_pc_org, file = fncurve_pc_org  , access ='append', status = 'replace')
  write(fcurve_pc_org, *)  ' # Initial Poinc curve by NP crossing:  xi   (x y z vx vy vz)  cj  tpc'
  
  open(fcurve_pc_fc, file = fncurve_pc_fc  , access ='append', status = 'replace')
  write(fcurve_pc_fc, *)  ' # Fourier approximation of Poinc curve:  xi   (x y z vx vy vz)  cj'
endif 

! ---------- Computation -----------------------------
! Outer bound for escape, nondimensinal unit 1 is big enough....
xmax0 = 2.d0 
tmax0 = 5.d1        ! a rough guess, problem-dependent  
tdir  = 1           ! integrate forward.  
dir0  = 1           ! 1: positive velocity;   -1: v<0 
imax0 = 1           ! consider every intersection, and choose by the direction of velocity

! Poincare section:  PV(ind0) = p00, here we fix z=0 for the 
print*, 'Input the index of coordinate for Poincare section: X(ind)=0'
read*,  ind0;  print*

!ind0  = 2           ! 1: x,   2: y,  3: z

! rough guess of poincare section, by observing directly at the plot of the tori 
! if ind = 2, we take y0 = 0, if ind = 1, we take x0 = 0.315 
!if(ind0 == 1)  p00 = 0.315d0 
!if(ind0 == 2)  p00 = 0.d0 
!p00 = 0.d0 ! at the moment! 

! cj0 = 4.d0 
!! cj0 = 3.95d0
!call init_poinc_h0(cj0) 
     
!  ****   I.C. of torus 1 : x-y-z-vx-vy-vz-hn  normal case  *****

!    1. take an unknown periodic orbits, compute the MM, and associated eigenvalues and eigenvectors
! norm1, only one family, centre-saddle
!    0.219662920343E+01    0.000000000000E+00    0.932165693007E+00    0.701220365384E+00    0.460454350218E+00    0.000000000000E+00    0.000000000000E+00   -0.179869260122E+01    0.133226762955E-14  0.1600E-01   -1
!ind0 = 1
!tp0    = 0.219662920343d1
!poinst = (/0.d0, .932165693007d0,  0.701220365384d0, 0.460454350218d0, 0.d0, 0.d0/) ! centre X centre   
   
! saddle X centre h=-1.7    
!       0.264916780655E+01    0.000000000000E+00    0.995256291437E+00    0.650832514494E+00    0.322483730278E+00    0.000000000000E+00    0.000000000000E+00   -0.170564943658E+01    0.444089209850E-15  0.1600E-01   -1
!tp0    = 0.264916780655d1
!poinst = (/0.d0, .995256291437d0,  0.650832514494d0, 0.322483730278d0, 0.d0, 0.d0/) ! centre saddle to start with 
!   

! ---- norm3, ifam1, h=4
! take the energy close to 4.0, where we have 2 periodic orbits for the first planar  family, one centre X saddle, one centre X centre 
!    0.593491339810E+00     0.669782093436E+00    0.000000000000E+00    0.000000000000E+00    0.000000000000E+00   -0.576483748133E+00    0.000000000000E+00    0.399953643572E+01    0.444089209850E-15  0.8000E-02    1
!     0.203894618794E+00    0.449946679072E+00    0.000000000000E+00    0.000000000000E+00    0.000000000000E+00   -0.102850594859E+01    0.000000000000E+00    0.399450268905E+01    0.621724893790E-14  0.8000E-02    1

! --succeed with this one!
!ind0 = 2

!tp0    = 0.593491339810d0
!poinst = (/0.669782093436d0, 0.d0, 0.d0, 0.d0, -0.576483748133d0,  0.d0/) ! centre saddle to start with 

!tp0    = 0.203894618794d0
!poinst = (/0.449946679072d0, 0.d0, 0.d0, 0.d0, -0.102850594859d1,  0.d0/) ! centre X centre 

! -------norm3,  ifam2, h=4  2 orbits
! -- succeed with this one   - centre X saddle 
!    0.219773264803E+01    0.699131986266E+00    0.000000000000E+00   -0.277761694501E-12   -0.760963495201E-12   -0.158071934111E+00   -0.553048316996E+00    0.399619759445E+01    0.177635683940E-14  0.1600E-01   -1 # h=4

ind0 = 3    
tp0  = 0.219773264803d1  ! t2 from fft: 2.1977326477769235 
poinst = (/0.699131986266d0, 0.d0, 0.d0, 0.d0, -0.158071934111d0,  0.553048316996d0/) ! centre saddle to start with 


! fail this one no matter ind0=2 or 3  -- centre X centre
!    0.903410306172E+00    0.338496551294E+00    0.000000000000E+00    0.599636604742E-01   -0.596544282138E-01    0.135868534036E+01    0.350506167261E+00    0.400858698754E+01   -0.266453525910E-13  0.1600E-01   -1 # torus

!tp0    = 0.903410306172d0
!poinst = (/0.338496551294d0, 0.d0, .599636604742d-1, -0.596544282138d-1, 0.135868534036d1,0.350506167261d0/) ! centre saddle to start with 
!    

! *************** Initialize Poinc*************
p00 = poinst(ind0)
call init_poinc( ind0, p00, dir0, imax0, tmax0, xmax0, ndim) 

! Refine a periodic orbit to some specific value
call gr_cjlf(poinst, cj0)
call init_poinc_h0(cj0)

 
 ! -- save all the relevent parameter to a file 
open(222,  file = './dat/curve_lf/para_poinc.dat', status = 'replace') 
write(222, *) ' # n, ind, p0, dir, cj0,  imax, tmax, ndim, nitmax, tol, tol_err'
write(222, *)     n, ind, p0, dir, cj0,  imax, tmax, ndim, nitmax, tol, tol_err 
 close(222)
 
 
open(38, file = './dat/curve_lf/po_org1.dat', access ='append', status = 'replace' )
call plob(poinst, 0.d0, tp0, ndim, ndim, 1, 1, 38, gr_lf, gr_cjlf, pvf)

print*, 'check I.C. and F.C. of the periodic orbit:'
print*, poinst;  print*, pvf(1:ndim) 
print*; read* 
 
call monomat(poinst, ndim, tp0, mmat, gr_lf, gr_cjlf)


! Compute the eigenvalues and eigenvectors of dlf
call eigrg(mmat, ndim, 1, wr,wi, vr)
!wr0 = wr; wi0 = wi; vr0 = vr  
!call eig_sort( wr0, wi0, wr0, wi0, 6, 1, wr, wi, ind_st)
!vr = vr(ind_st, :)


print*,  'Which (pair) of eigvalue (central part) to use for initial guess of linear flow.?'
read(*,*)  ivr
i_eig = 2*ivr-1
print*,  wr(i_eig), '+i', wi(i_eig)
 
rho = datan_2pi(wr(i_eig), dabs(wi(i_eig)), 1) ! [0, 2pi]
print*, 'Initial guess of rho:', rho, dasin(dabs(wi(i_eig)))

print*, 'Input the columns of corresponding eigenvectors: '
read*, ivr1, ivr2
v1 = vr(:, ivr1); v2 = vr(:, ivr2) !

!v1 = vr(:, 2*ivr-1); v2 = vr(:, 2*ivr)
print*, v1;  print*, v2 ; 
print*; read*



!  2. take the linear flow as the initial guess for invariant curve \phi: x0 + epsilon * (vr*cos(xi) - vi*sin(xi))
!     -- epsilon is the distance from the periodic orbit, vr+ivi is the associated eigenvector, 
!     -- x0 is an initial condition of the p.o., take the one on Poincare section just for convenience.
 
! we can take: c0 = poinst, c1 = epsl * vr, s1 = - epsl * vi  
    
! 3. Integrate this curve (take NP points) for the next crossing with the  Poincare section to obtain our initial curve \varphi 

epsl = 1.d-5 
dcj  = 1.d-6 


dxi = pi2/(np-1)  ! take this value to keep coherent with the original points 

do i = 1,  np, 1  
 
  xi  = (i-1) * dxi
  
! x0 + epsilon * (vr*cos(xi) - vi*sin(xi))
  pv0 = poinst + epsl*( v1*dcos(xi) - v2*dsin(xi) ) 
  
  
  ! save the torus from the time curve
  if(i > 1) open(ftorus_time, file = fntorus_time, access='append', status='old')
  call plob(pv0, 0.d0, tp0, ndim, ndim, 1, 1, ftorus_time, gr_lf, gr_cjlf, pvf)
  
  
  ! compute the next crossing with poincare section
!  subroutine poinc(sti, ndim, nvar, tdir, tf, stf, hminim, ispc, gr_lf, gr_cjlf)
  call poinc(pv0, ndim, ndim, 1, tpc, pvpc, hminim, ispc, gr_lf, gr_cjlf) 
 
  ! the orbits that connects the poincare curve and time curve.
  if(i > 1) open(ftorus_conct, file=fntorus_conct, access='append', status='old')
  call plob(pv0, 0.d0, tpc, ndim, ndim, 1, 1, ftorus_conct, gr_lf, gr_cjlf, pvf)
  
  if(ispc == 0) then  ! just be sure we have the crossing 
    print*, 'Fail to reach the Poincare section!'
    print*; read*
    return
  
  else  
    xpc_all(i,:) = pvpc 
    
  endif  
  write(fcurve_time, *)  xi, pv0, pvpc, tpc ! save time_org together 
   
enddo    

 close(fcurve_time) 
 
print*, 'check time curve and initial poincare curve: they should be a curve!'
!print*, fncurve_time 
!print*; read*

! ********************* Interpolation **********************************
! ---- we have to parameterize the Poincare care, the angle parameter is not the same as the time curve------ 
!      keep that in mind!!!! 

pt = xpc_all(:,  ind_fun) 
do i = 1, np, 1
  if(ind == 2) then 
    theta = datan_2pi(pt(i, 1), pt(i, 2))
  else 
    theta = datan_2pi(pt(i, 2), pt(i, 1))
  endif 
end do

call sort_1d( theta, np, 1, ind_sorted)

! update the points by an increasing order of theta 
xpc_all = xpc_all(ind_sorted, :) 

! -- interpolated curve (equally spaced in rho) --
fcurve_interpol = 26;   fncurve_interpol  =  './dat/curve_lf/curve_interpol.dat'
open(fcurve_interpol,  file = fncurve_interpol, access='append', status='replace')
write(fcurve_interpol, *)  ' #  Interpolated curve: arg, PV,  the interpolated components (', ind_fun,'), '  
print*,'The index of components to do interpolation:', ind_fun
  
  pt  = xpc_all(:, ind_fun)
  call  interpol( pt, np, n, theta, np, ptnew)
  
  ! Period = 2*pi in radian w.r.t. the angle of the points 
  darg = pi2/(np-1) ! take this value to keep coherent with the original points 

  do i = 1, np, 1
    arg = (i-1) * darg
    pv  = ptnew(i, :) ! the sorted equispaced points
    
    call gamm(pv, n, ind_fun, cj0, ind, p0, ptnew1, isv2, cj2v_lf)
    
    write(fcurve_interpol, '(13e24.14)')  arg,  ptnew1, xpc_all(i,:)
    ptnew(i, :) = ptnew1
  end do
  
  close(fcurve_interpol) 

print*, '****** Finish interpolation ***********'; print*; read*

! -- For curve initialization, c0, ck, sk are global variable from curve_mod module 
nf = nf_init 
call upd_four(0, nf) 
call alloc_arr_1d(nf, csf);  call alloc_arr_1d(nf, sif)
! ------------- Do Fourier Analysis using the interpolated points --------
do k = 1, n, 1 
  pv_in =  ptnew(:, ind_fun(k) ) ! one component, the k-th column  
  call gr_four(pv_in, np, nf, csf0, csf, sif)
  c0(k)    = csf0
  ck(k, :) = csf
  sk(k, :) = sif
end do


open(14, file='ck_fcs.dat', access='append', status='replace')
write(14,*) '# Check Fourie Coefficients: c0, ck, sk.  format -- (10e24.14)'
write(14, *) '# Initialization: n, nf = ', n, nf 
do i = 1, n, 1
  write(14, '(6e24.14)')   c0(i);     write(14,*) 
  write(14, '(10e24.14)')  ck(i, :);  write(14,*) 
  write(14, '(10e24.14)')  sk(i, :) 
  write(14,*);  write(14,*) ! two blank line2 to seperate components 
end do
  close(14); open(14, file='ck_fcs.dat', access='append', status='old')
  
  
! check if the initial Fourier coefficients are good enough? 
do i = 1,  np, 1  
  xi = (i-1) * dxi
  call varphi( xi,     n, nf, c0, ck, sk, pv_curv)
  call gamm(pv_curv,   n, ind_fun, cj0, ind, p0, pvf, isv2, cj2v_lf)
  write(fcurve_pc_fc,  '(8e24.14)')   xi, pvf, cj0 
enddo  
 close(fcurve_pc_fc)
 
 print*, '****** Finish Fourier Approximation ***********'
 
! TODO: Debug!!! The results seem wrong after we call upd_four to do the 
!       assignment and update of the Fourier Coefficients 
 
 print*, 'From FFT: n, nf=', n, nf 
 
 ! update the number of harmonics to use 
call upd_four(1, nf) !, c0, ck, sk)
!print*, 'afte upd_four, nf =', nf; ! print*; read* 

print*, 'before the refinement! n, nf:', n, nf;! print*; read*  

open(177, file = './dat/curve_lf/fcs_refn.dat', status='replace', access='append')
write(177, *) '# refined curve: c0 // ck // sk  format -- (10e24.14), rho(radian)=',rho

open(133,  file = './dat/curve_lf/para_curve.dat', status = 'replace', access='append') 
write(133, *) ' # nf, rho, h0,  w1,  w2,  t2, c0, ck(:, 1), sk(:, 1) ...  ck(:, 4), sk(:, 4)'



! ***************** Continutation ****************************
do itori = 1, 5, 1 
   
   if(itori > 1) cj0 =  cj0 + dcj 
print*, 'Refine a new curve!, cj0 = ', cj0
call init_poinc_h0(cj0)

!subroutine refine_curve_poinc(n, nf, rho, c0, ck, sk, isref, iter, gr_poinc)
call refine_curve_poinc(n, nf, rho, c0, ck, sk, isref, niter, PoincMap_lf, gr_cjlf, deriv_cjlf)

if (isref == 0) then 
  print*, 'Fail for the refinement of the curve! stop!'
  stop 
  
else   
  !  update the number of harmonics
  nf0 = nf
  call upd_four(1, nf)!, c0, ck, sk)
!  print*, 'nf0, nf: ', nf0, nf; print*; read*
endif   

! --- use this updated c0, ck, sk to compute refined curve ---- 
do i = 1, n, 1
  write(177, '(6e24.14)')   c0(i);     write(177,*) 
  write(177, '(10e24.14)')  ck(i, :);  write(177,*) 
  write(177, '(10e24.14)')  sk(i, :) 
  write(177,*);  write(177,*) ! two blank line2 to seperate components 
end do
! close(177)

 
! Evaluate the truncated Fourier series for \hat T(xi)
do i = 1,  np, 1  
  xi = (i-1) * dxi
     
  call varphi( xi,     n, nf, c0, ck, sk, pv_curv)
  call varphi( xi+rho, n, nf, c0, ck, sk, pv_curv_rho)
  
!  print*, 'n, ind, ind_fun', n, ind, ind_fun; print*; read*
  call gamm(pv_curv,     n, ind_fun, cj0, ind, p0, pvf, isv2, cj2v_lf)
  call gamm(pv_curv_rho, n, ind_fun, cj0, ind, p0, pvf_rho, isv2, cj2v_lf)

  ! compute the return time 
!  call PoincMap_tf_lf(pv_curv, tpc, pvpc, 0, dp_dx, ispc) ! discard this one 
  call poinc(pvf, ndim, ndim, 1, tpc, pvpc, hminim, ispc, gr_lf, gr_cjlf) 
  tpc_refn(i)   = tpc
  xpc_all(i, :) = pvf

  write(fcurve_pc,     '(9e24.14)')   xi, pvf, cj0, tpc 
  write(fcurve_pc_rho, '(14e24.14)')  xi, pvf_rho, pvpc, cj0 
    
  if(i > 1) open(ftorus_pc, file=fntorus_pc, access='append', status='old')
  call plob(pvf, 0.d0, tpc, ndim, ndim, 1, 1, ftorus_pc, gr_lf, gr_cjlf, pvpc)
enddo 

write(fcurve_pc,*);  write(fcurve_pc,*);  
write(fcurve_pc_rho,*);  write(fcurve_pc_rho,*);

! close(fcurve_pc); close(fcurve_pc_rho)

! --- for the computation of the two basic frequencies 
! Do a FFT to obtain X_0,  one component - by - one component
call alloc_arr_1d(nf, ct); call alloc_arr_1d(nf, st)
call gr_four(tpc_refn, np, nf, t2, ct, st) 
 
print*, 'Compare T2 and the Tp:', t2, tp0; print*; read* 
 
! -- save all the relevent parameter of the curve to a file 
write(133, '(I5, 9e24.14)', advance='no')  nf, rho, cj0, pi2/t2, rho/t2, t2, c0 
do i = 1, 4, 1
  write(133, '(8e24.14)', advance='no')  ck(:, i), sk(:, i)
end do 
write(133,*); ! close(133)

print*, 'Finish para_curve!'; !print*; read*
 
enddo 
 close(133); close(fcurve_pc); close(fcurve_pc_rho); close(177)

stop 



! ******** compute the Fourier coefficients ctau, stau for the small time deviation  ***********
!    from the Time-T map method  tau(xi)

!open(112, file = './dat/curve_lf/tpc_fun.dat', status='replace', access='append')
!write(112,*) ' # xi, tpc, tpc_fun(Fourier), tau_fun(Fourier), curve_time_fun'
! ctau0 = 0.d0
!call alloc_arr_1d(nf, ctau); call alloc_arr_1d(nf, stau)
!call tau_fc( ct, st, rho, nf, ctau, stau)

!dxi = pi2/(np-1)  
!! Evaluate the truncated Fourier series for \hat T(xi)
!do j = 1, np, 1  
!  xi = (j-1) * dxi
!  
!  ! the fourier representation of the return time of the approximated curve 
!  ! for 1d function, call the funtion four_seri is better    
!  tpc = four_seri(xi, nf, t2,     ct,  st)
!  tau = four_seri(xi, nf, ctau0,  ctau, stau)
!  
!  pvf = xpc_all(j, :)
!  call plob(pvf, 0.d0, tau, ndim, ndim, 1, 0, 6, gr_lf, gr_cjlf, pvpc) 
!  write(112, '(10E20.10)')  xi, tpc_refn(j), tpc, tau, pvpc   
!enddo 
!  
! close(112) 
!print*, 'Finish tau and Tpc_fun'; print*; read*


 cj0 =  cj0 + dcj 
print*, 'Refine a new curve!, cj0 = ', cj0
call init_poinc_h0(cj0)
call refine_curve_poinc(n, nf, rho, c0, ck, sk, isref, niter, PoincMap_lf, gr_cjlf, deriv_cjlf)
print*, 'isref:', isref 
print*, 'Refine the new curve!'; print*; read*

dxi = pi2/(np-1)
! Evaluate the truncated Fourier series for \hat T(xi)
do i = 1,  np, 1  
  xi = (i-1) * dxi
  call varphi( xi,     n, nf, c0, ck, sk, pv_curv)
  call varphi( xi+rho, n, nf, c0, ck, sk, pv_curv_rho)
  
  call gamm(pv_curv,     n, ind_fun, cj0, ind, p0, pvf, isv2, cj2v_lf)
  call gamm(pv_curv_rho, n, ind_fun, cj0, ind, p0, pvf_rho, isv2, cj2v_lf)

  write(fcurve_pc,     '(9e24.14)')   xi, pvf, cj0, tpc 
  write(fcurve_pc_rho, '(14e24.14)')  xi, pvf_rho, pvpc, cj0 
  
  open(ftorus_pc, file=fntorus_pc, access='append', status='old')
  call plob(pvf, 0.d0, tpc, ndim, ndim, 1, 1, ftorus_pc, gr_lf, gr_cjlf, pvpc)
enddo 
write(fcurve_pc,*);  write(fcurve_pc,*);  
write(fcurve_pc_rho,*);  write(fcurve_pc_rho,*);



 cj0 =  cj0 + dcj 
print*, 'Refine a new curve!, cj0 = ', cj0
call init_poinc_h0(cj0)
call refine_curve_poinc(n, nf, rho, c0, ck, sk, isref, niter, PoincMap_lf, gr_cjlf, deriv_cjlf)
print*, 'isref:', isref 
print*, 'Refine the new curve!'; print*; read*

dxi = pi2/(np-1)
! Evaluate the truncated Fourier series for \hat T(xi)
do i = 1,  np, 1  
  xi = (i-1) * dxi
  call varphi( xi,     n, nf, c0, ck, sk, pv_curv)
  call varphi( xi+rho, n, nf, c0, ck, sk, pv_curv_rho)
  
  call gamm(pv_curv,     n, ind_fun, cj0, ind, p0, pvf, isv2, cj2v_lf)
  call gamm(pv_curv_rho, n, ind_fun, cj0, ind, p0, pvf_rho, isv2, cj2v_lf)

  write(fcurve_pc,     '(9e24.14)')   xi, pvf, cj0, tpc 
  write(fcurve_pc_rho, '(14e24.14)')  xi, pvf_rho, pvpc, cj0 
  
  open(ftorus_pc, file=fntorus_pc, access='append', status='old')
  call plob(pvf, 0.d0, tpc, ndim, ndim, 1, 1, ftorus_pc, gr_lf, gr_cjlf, pvpc)
enddo 


 close(fcurve_pc);  close(fcurve_pc_rho); !close(ftorus_pc)

stop 


end 


 




program curve_time_lf 
!2017-05-26 22:09:27  TURN TO THIS APPROACH, GIVE UP THE POINCARE MAP APPROACH 
 
! check only the first one, and then check the continuation.... please 

! deal with 6D phase space. check if the routine works.
! Instead of Poincare map method, we took the time T_2 method 

! the eigenvalue (rho) of centre part (pure imaginary) serves as the rotation number (in radian),
! the two frequencies are close to:  rho/T, 2pi/T. 

!  --- for the transversal case: 
! the new T and rho are: δ = 2pi / (±ν +2pi j) T,  rho = δ * 2pi / T =  (2pi)^2 / (±ν +2pi j) * T 
! and the two basic frequencies are close to:  rho/δ,  2pi/δ. 

! TODO we have to check if we need to multiply the eigenvalue by 2pi as the rotation number 

! ADD the multiple shooting, since the saddle part makes it hard to converge if its modulus is much greater than 1.


! ** NOTE ** 
! DO not forget the initialization of the model parameters:   beta, sgn, ndim  

!  1. take an unknown periodic orbits, compute the MM, and associated eigenvalues and eigenvectors
!  2. take the linear flow as the initial guess for invariant curve \phi: x0 + epsilon * (vr*cos(xi) - vi*sin(xi))
!     -- epsilon is the distance from the periodic orbit, vr+ivi is the associated eigenvector, 
!  3. do a simple Fourier analysis to compute the Fourier coefficients as the initial guess \bm X_0
!     it is only necessary when istrans=1
 
!  4. Impose the invariance equations to compute the invariant curve
!     the Jacobi matrix and system equations are to be specified according to opt 

!  5. Globalize the curve to an invariant torus, by integrating from one period
!

!  --- data saved in ./dat/curve_lf/time/  subfolder 
!  1. torus.dat            -- the original orbit of the 2D torus  ! not used at the moment 
!  2. curve_time.dat       -- the invariant curve obtained by Time map
!  4. fun.dat              -- approxiamted Fourier representation \varphi(xi+rho), the check \varphi(xi+rho), P( \varphi(xi) )
!  5. fcs.dat              -- fourier coefficients, For each coordinate, three parts: c0 // ck // sk


use dp_mod
use pi_mod
use sort_mod
use lf_mod 
use curve_time_mod


implicit none

! perhaps we need to do multiple shooting for big TP, for example, to start transverally to the planar p.o. norm3, fam1 

integer, parameter  ::  np    = 101 ! 65    !  number of points on the curve 
 
integer, parameter  ::  ndim  = 6,  n=ndim  ! dimension of phase sapce for spatial lf
                        
integer, parameter  ::  nvar  = ndim*(ndim+1)                        

real(kind=dp) ::  beta0 

! --  ---- declare the size for Fourier coefficients of the original  Fourier analysis   ----
! how many Fourier modes ???, Alex suggests of power 2 
! Gerard said 60 something is enough ....

! nf: apply Alex's strategy: claim a big enough array at first, then check the maximum norm 
!     of the last half, if it is greater than a given threshold, which should be a little bit smaller than 
!     the precision we ask for the Newton method,  we increase the size by 2
!     if it is smaller, then we decrease the size by half 

!     the arrays c0,  ck, sk are claimed allocatable arrays.                  

!     since they are just initial seed computed from a linear interpolated curve, not accurate enough as 
!     the estimate of the precision of the refined curve 
!     so we check again the refined Fourier coefficients by looking at the Maximum norm of the last half.

! --- this is harder to implement 
! TODO: JM's idea, paper  Physic D. P12. Section 2.2.6. Error estimate
!       increase the value of nf, until the error in Invariance equation is under a given tolerance. 1.d-10
!       where he uses multiple shooting with intermediate node M=50 


!integer, parameter  ::  nf_init   =  2 ! to debug, check every entry of Jacobi Matrix  

integer, parameter  :: nf_init   = 32 ! 16  ! 32 ! 32   ! 16 !8! 16 !64 !8 !32  !64 ! 512 ! 128 !   !64 ! 32  !8 
integer, parameter  :: m0        = 1   ! test with 2 

! -- fourier coefs for Fourier analysis for the initial X_0  -- 
real(kind=dp) :: csf0 !,   dc0(n) !, rho_old, tp_old  
real(kind=dp), allocatable,  dimension(:)   :: csf, sif  ! fourier coefs of 1 component
!real(kind=dp), allocatable,  dimension(:,:) ::  sk_copy,   sk_old 

! -- fourier  coefs for Refinement  -- 
! the number of Fourier modes we are going to use will be specified 
! after we check the maximum norm of ck and sk, all the following allocatable arrays have nf-columns  

integer       :: nf ! , nf0,niter !, isrestart, nf_copy, nf_old, ncol
                    
!  ------ Local Variables ----------- 
real(kind=dp) :: pv0(ndim), pv1(ndim), pvf(ndim), cj0,     & ! torus
                 tp0, poinst(ndim),  mmat(ndim, ndim), wr(ndim), wi(ndim), vr(ndim, ndim), & ! MM 
                 v1(ndim), v2(ndim), v10(ndim), v20(ndim),  xi, dxi, epsl,    & ! linear flow 
                 xpc_all(np,ndim),    & ! poinc 
                 pv_in(np), tol, tol_err   ! Netwon method
!                  pv_curv(ndim),  pv_curv_rho(ndim) ,  pvpc(ndim), 
               
                  
                                                      
integer ::  i, ivr, i_eig, ivr1, ivr2,  ind0,    iscurve,  nitmax, opt, isincr

  
integer           ::  ftorus_time, fcurve_time, fcurve_time_refn, ftori_time_refn         
character(len=70) ::  fntorus_time, fncurve_time, fncurve_time_refn, fntori_time_refn   


!---------- Transversal curve --------------           
integer :: istrans        
real(kind=dp)     ::  pv0_new(ndim), dt, tp1, rho0,  eta,  tdir, & 
                      pvf2(ndim),  xpc_all_trans(np, ndim), phi(ndim, ndim), pvil(nvar), pvfl(nvar)


! -- Fourier coefficients 
!real(kind=dp)              ::     step
integer       :: nf_min, nf_max, ntori 
real(kind=dp) :: ds, ds_min, ds_max, dir 
   
!-------------------------- continuation -------------------------------
integer       ::  k,  tori, itori, im, nrow, row_st !   ind_c00,
real(kind=dp) ::  cj  !,   dcj         ! step control for continuaiton 
!-------------------------- continuation -------------------------------
                  
                  
external :: TimeMap_lf ! general 2d map for lf  

real(kind=dp) ::  datan_2pi !, dnrm2  ! external functions 



call init_lf


tori = 0 

print*, 'Please input beta (q/m):'
!read(*,*) beta0 

beta0 = 2.d0 
call init_beta(beta0) 

! Time-T map method 
opt     = 10
call init_opt(opt)
 
!opt = 100  ! fail 


istrans = 0 

! -- Initialization for the module  curve_mod, use poincare map approach, n=ndim-2
nitmax  = 20    
tol     = 1.d-10
tol_err = 1.d-10 

call init_curve(ndim, m0, nitmax, tol, tol_err) 

  
iscurve    = 1   ! 1: Do poincare map of the torus to get the curve; other: curve already computed. 

! files to save the original torus + invariant curve + interpolated one + approximated one 
ftorus_time      = 21;  fntorus_time      =  './dat/curve_lf/time/torus_time.dat'
fcurve_time      = 25;  fncurve_time      =  './dat/curve_lf/time/curve_time.dat'
fcurve_time_refn = 28;  fncurve_time_refn =  './dat/curve_lf/time/curve_time_refn.dat'
ftori_time_refn  = 29;  fntori_time_refn  =  './dat/curve_lf/time/tori_time_refn.dat'

open(ftorus_time,  file = fntorus_time, access ='append', status = 'replace')
write(ftorus_time,*)  ' # Initial torus:  xi   (x y z vx vy vz)   cj'

open(fcurve_time, file = fncurve_time, access ='append', status = 'replace')
write(fcurve_time,*)  ' # Initial Time + rho  curve:  xi   (x y z vx vy vz)   cj'
    
open(fcurve_time_refn,  file = fncurve_time_refn, access ='append', status = 'replace')
write(fcurve_time_refn,*)  ' # Refined Time + rho  curve:  xi   (x y z vx vy vz)   cj'
 
open(ftori_time_refn,  file = fntori_time_refn, access ='append', status = 'replace')
write(ftori_time_refn,*)  ' # Refined Time + rho  tori:  xi   (x y z vx vy vz)   cj'
       
    
! ---------- Computation -----------------------------
! Outer bound for escape, nondimensinal unit 1 is big enough....
!  ****   I.C. of torus 1 : x-y-z-vx-vy-vz-hn  normal case  *****

!    1. take an unknown periodic orbits, compute the MM, and associated eigenvalues and eigenvectors --- not converge
! norm1, only one family, centre-saddle
!    0.219662920343E+01    0.000000000000E+00    0.932165693007E+00    0.701220365384E+00    0.460454350218E+00    0.000000000000E+00    0.000000000000E+00   -0.179869260122E+01    0.133226762955E-14  0.1600E-01   -1
tori    = 10
ind0    = 1
tp0     = 0.219662920343d1
poinst  = (/0.d0, .932165693007d0,  0.701220365384d0, 0.460454350218d0, 0.d0, 0.d0/) ! centre X centre   
!istrans = 1   

!ind0    = 2
   
! saddle X centre h=-1.7    
!       0.264916780655E+01    0.000000000000E+00    0.995256291437E+00    0.650832514494E+00    0.322483730278E+00    0.000000000000E+00    0.000000000000E+00   -0.170564943658E+01    0.444089209850E-15  0.1600E-01   -1
!tp0    = 0.264916780655d1
!poinst = (/0.d0, .995256291437d0,  0.650832514494d0, 0.322483730278d0, 0.d0, 0.d0/) ! centre saddle to start with 
!   

! ---- norm3, ifam1, h=4
! take the energy close to 4.0, where we have 2 periodic orbits for the first planar  family, one centre X saddle, one centre X centre 
!    0.593491339810E+00     0.669782093436E+00    0.000000000000E+00    0.000000000000E+00    0.000000000000E+00   -0.576483748133E+00    0.000000000000E+00    0.399953643572E+01    0.444089209850E-15  0.8000E-02    1
!     0.203894618794E+00    0.449946679072E+00    0.000000000000E+00    0.000000000000E+00    0.000000000000E+00   -0.102850594859E+01    0.000000000000E+00    0.399450268905E+01    0.621724893790E-14  0.8000E-02    1

! -- planar one is difficult, we succeed with istrans = 1 ! 
!tori    = 11 
!ind0    = 3
!tp0     = 0.593491339810d0
!poinst  = (/0.669782093436d0, 0.d0, 0.d0, 0.d0, -0.576483748133d0,  0.d0/) ! centre saddle to start with 
!istrans = 1 

!tp0    = 0.203894618794d0
!poinst = (/0.449946679072d0, 0.d0, 0.d0, 0.d0, -0.102850594859d1,  0.d0/) ! centre X centre 

! -------norm3,  ifam2, h=4  2 orbits
! -- succeed with this one   - centre X saddle 
!    0.219773264803E+01    0.699131986266E+00    0.000000000000E+00   -0.277761694501E-12   -0.760963495201E-12   -0.158071934111E+00   -0.553048316996E+00    0.399619759445E+01    0.177635683940E-14  0.1600E-01   -1 # h=4

! up to w1/w2 = 8.99  # N_3 fam2  tori1 
!tori = 21     
!ind0 = 3    
!tp0  = 0.219773264803d1  ! t2 from fft: 2.1977326477769235 
!poinst = (/0.699131986266d0, 0.d0, 0.d0, 0.d0, -0.158071934111d0,  0.553048316996d0/) ! centre saddle to start with 

! 0.69913198626600E+00    0.21869538702075E-18   -0.26463317253161E-20   -0.33891657647829E-20   -0.15807193411100E+00    0.55304831699600E+00

! good this one no matter ind0=2 or 3  -- centre X centre 
!    0.903410306172E+00    0.338496551294E+00    0.000000000000E+00    0.599636604742E-01   -0.596544282138E-01    0.135868534036E+01    0.350506167261E+00    0.400858698754E+01   -0.266453525910E-13  0.1600E-01   -1 # torus

! up to w1/w2 = 7.99 
!tori   = 22
!ind0   = 2
!tp0    = 0.903410306172d0
!poinst = (/0.338496551294d0, 0.d0, .599636604742d-1, -0.596544282138d-1, 0.135868534036d1,0.350506167261d0/) 


! *************** Initialize Poinc*************
open(38, file = './dat/curve_lf/time/po_org1.dat', access ='append', status = 'replace' )
call plob(poinst, 0.d0, tp0, ndim, ndim, 1.d0, 1, 38, gr_lf, gr_cjlf, pvf)
 close(38)
 
call gr_cjlf(poinst, cj0)
call init_time_h0(cj0)

print*, 'h0', cj0; print*; read*

print*, 'check I.C. and F.C. of the periodic orbit:'
print*, poinst;  print*, pvf(1:ndim) 
print*; read* 
 
 
call monomat(poinst, ndim, tp0, mmat, gr_lf, gr_cjlf)

! Compute the eigenvalues and eigenvectors of dlf
call eigrg(mmat, ndim, 1, wr,wi, vr)

print*,  'Which (pair) of eigvalue (central part) to use for initial guess of linear flow.?'
if(tori == 0)  read(*,*)  ivr

if(tori == 10)   ivr = 2
if(tori == 11)   ivr = 3
if(tori == 21)   ivr = 2 
if(tori == 22)   ivr = 1 

i_eig = 2*ivr-1
print*,  'Take ', i_eig, '-th eigenvalue!'
print*,  wr(i_eig), '+i', wi(i_eig)
 
rho0 = datan_2pi(wr(i_eig), dabs(wi(i_eig)), 1) ! [0, 2pi]
print*, 'Initial guess of rho:', rho0,  dasin(dabs(wi(i_eig)))

print*, 'Input the columns of corresponding eigenvectors: '
if(tori == 21) then 
   ivr1 = 3; ivr2 = 4 
elseif(tori == 11) then 
   ivr1 = 5; ivr2 = 6 
elseif(tori == 10) then 
   ivr1 = 2; ivr2 = 3 
elseif(tori == 22) then 
   ivr1 = 1; ivr2 = 2
endif    

if(tori == 0)  read*, ivr1, ivr2


v1  = vr(:, ivr1); v2 = vr(:, ivr2) !
print*, v1;  print*, v2 ; 

! normalize the two eigenvector 
v10 = v1; v20 = v2 
!v10 = v1 / dnrm2(ndim, v1, 1);  v20 = v2 / dnrm2(ndim, v2, 1)
print*; read*


!print*, 'Input the index of component to set to 0  to avoid curve indetermination(1-3): '
!read*, ind_c00 

call init_ind_c0(ind0)


!  2. take the linear flow as the initial guess for invariant curve \phi: x0 + epsilon * (vr*cos(xi) - vi*sin(xi))
!     -- epsilon is the distance from the periodic orbit, vr+ivi is the associated eigenvector, 
!     -- x0 is an initial condition of the p.o., take the one on Poincare section just for convenience.
 
! we can take: c0 = poinst, c1 = epsl * vr, s1 = - epsl * vi  
    
! 3. Integrate this curve (take NP points) for the next crossing with the  Poincare section to obtain our initial curve \varphi 
 

epsl = 1.d-6

if(istrans == 1) then 
  open(411, file='./dat/curve_lf/time/torus_trans.dat', status='replace', access='append')
  write(411,*) '# The transverse torus! t, x,y,z,vx,vy,vz'

  open(412, file='./dat/curve_lf/time/curve_trans.dat', status='replace', access='append')
  write(412,*) '# The transverse curve. curve, curve_rho,  \varphi(eta) ! eta, x,y,z,vx,vy,vz'
endif 


dxi = pi2/np  ! take this value to keep coherent with the original points 

do i = 1,  np, 1  
 
  xi  = (i-1) * dxi
  
  pv0 = poinst + epsl*( v10*dcos(xi) - v20*dsin(xi) ) 
  
  if(i == 1) then 
    rho = rho0
    call gr_cjlf(pv0, cj0)
    call init_time_h0(cj0)
    call init_time_tp(tp0)
    print*, 'tp0, tp',  tp0, tp ; print*; read* 
  endif 
  
  ! save the torus from the time curve
  if(i > 1)  open(ftorus_time, file = fntorus_time, access='append', status='old')
  call plob(pv0, 0.d0, tp0, ndim, ndim, 1.d0, 1, ftorus_time, gr_lf, gr_cjlf, pvf)
  call gr_cjlf(pv0, cj)
  
  xpc_all(i,:) = pv0 
    
  write(fcurve_time, *)  xi, pv0, pvf, cj   
   
   
  if (istrans == 0)  cycle 
   
   
 ! ****** transversal curve ****************   TO CHECK
!  dt = \eta/pi2 *T
  eta = xi/pi2
  dt = eta * tp0;  eta = - eta * rho0 
 
  pv0_new = poinst + epsl*( v10*dcos(eta) - v20*dsin(eta) ) 
  
  ! Initialize the energy 
  if(i == 1) then 
    tp1  = pi2/rho0*tp0  ! for the transverse period ; the new rho is rho_trans = tp1*pi2/tp0 
    call init_time_tp(tp1)
  
    rho = tp1*pi2/tp0  ! - 2.d0*pi2 
    print*, 'original rho, tp, transveral rho, tp:', rho0, tp0,  rho, tp1
    print*; read*
     
    call gr_cjlf(pv0_new, cj0)
    call init_time_h0(cj0)
    print*, 'tp1, tp',  tp1, tp ; print*; read* 
  endif 
  
  tdir = dsign(1.d0, dt)
  
  if(i == 1) then 
    pvf = pv0_new
  else   
    call plob(pv0_new, 0.d0, dabs(dt), ndim, ndim, tdir, 0, 6, gr_lf, gr_cjlf, pvf)
  endif 
  
  xpc_all_trans(i, :) = pvf 
    
  ! save the torus from the transveral curve_approx
  if(i > 1) open(411, file = './dat/curve_lf/time/torus_trans.dat', access='append', status='old')
  call plob(pvf, 0.d0, tp1, ndim, ndim, 1.d0, 1, 411, gr_lf, gr_cjlf, pvf2)
  
  write(412, *)  xi, pvf,  pvf2 , pv0_new !  
enddo  
write(412, *);  write(412, *) ! start a new block 
   

! -- Initialization of c0, ck, sk, taking into consideration of multiple shooting ---- 
!  Starting “longitudinally” to the periodic orbit: 
!    
!  Starting “transversally” to the periodic orbit: 
!    we take the Fourier coefficients of ϕj from a Fourier expansion of φj(δ/m)(Lψ(0,η)).

nf = nf_init 
call upd_four(0, nf)   ! c0, ck, sk be initialized as zeros 

if(istrans == 0) then 
  ! im = 0, A^0 
  c0(1:n) = poinst; ck(1:n,1) = epsl*v10; sk(1:n,1) = -epsl*v20
  
  ! for the rest im = 1, m-1, the Fourier coefficents can be approximated by the evolution of X_0 and eigenvector 
  ! \phi(X_0) and Phi_(TP/m) * v1, v2 , A^0, A^1, ..., A^m-1 
  dt   = tp / m 
  tdir = dsign(1.d0, dt)
  
  pvil = 0.d0 
  pvil(1:ndim) = poinst 
  pvil(ndim+1:nvar:ndim+1) = 1.d0 
  
  
  do i = 1, m-1
    row_st = i*n 
    
    call plob(pvil, 0.d0, dt, ndim, nvar, tdir,  0, 0, gr_lf, gr_cjlf, pvfl)
!    write(fcurve_time, * ); write(fcurve_time, *)

    pvil = pvfl
    
    c0(row_st+1 : row_st+n) = pvfl(1:ndim) 
    phi = reshape(pvfl(ndim+1 : nvar), (/ndim, ndim/) ) ! the STM 
    
    print*, 'STM for ', i, '-th segment'; print*
    do k = 1, ndim 
      print*, phi(k,:)
    enddo
    print*; read*
    
    v1 = matmul(phi, v10);         v2 = matmul(phi, v20)
!    v1 = v1 / dnrm2(ndim, v1, 1);  v2 = v2 / dnrm2(ndim, v2, 1);
    
    ck(row_st+1: row_st+n, 1) =  epsl*v1  ! c1^im
    sk(row_st+1: row_st+n, 1) = -epsl*v2  ! s1^im
     
  enddo 
  
  close(fcurve_time); close(ftorus_time)
  
elseif( istrans == 1) then  
  
  ! do  fourier analysis on the transversal curve to get the initial guess of C0^0, CK^, SK^0
 
  ! for the rest A^im, im = 1, m-1, we integrate for TP/m  the intermediate curve, and do FFT on these curves 

  xpc_all = xpc_all_trans ! for im = 0 
  
  dt   = tp / m 
  tdir = dsign(1.d0, dt)
  
  call alloc_arr_1d(nf, csf);  call alloc_arr_1d(nf, sif)
  
     
  do im = 0, m-1   
   
    row_st = im * ndim
  
    ! ------------- Do Fourier Analysis for  C0^im, CK^im, SK^im --------
    ! for im=0, the first segment 
    do k = 1, ndim, 1
    
      nrow = row_st + k 
       
      pv_in =  xpc_all(:, k) ! one component, the k-th column  
      call gr_four(pv_in, np, nf, csf0, csf, sif)
      c0(nrow)    = csf0
      ck(nrow, :) = csf
      sk(nrow, :) = sif
    end do
    
    
   
    if(im == m-1)  exit 
    
    ! Integrate TP/m  to obtain the intermediate curve, saved in xpc_all, and then do FFT
    do i = 1, np 
      pv0 = xpc_all(i, :)
      call plob(pv0, 0.d0, dt, ndim, ndim, tdir, 0, 0, gr_lf, gr_cjlf, pvf)
      xpc_all(i, :) = pvf 
    enddo 
    
  enddo 
  
  
  ! ------------- approximated fourier representation ------------
  do im = 1, m
    do i = 1,  np, 1  
      xi = (i-1) * dxi
      call varphi( xi, n, im-1, nf, c0, ck, sk, pv0)
     
      if(im == m)  then  
         ! the last segment 
        call varphi( xi+rho, n, 0,  nf, c0, ck, sk, pv1)
      else
        call varphi( xi, n, im, nf, c0, ck, sk, pv1)
      endif    
    
      call gr_cjlf(pv0, cj) !;  call gr_cjlf(pv1, cj1)
    
      ! check the initial guess of the transversal curve, the error in the invariance equation 
      call plob(pv0, 0.d0, dt, ndim, ndim, tdir, 0, 0, gr_lf, gr_cjlf, pvf2) 
      write(412,  '(21e24.14)')   xi, pv0, pv1, pvf2, maxval(dabs(pvf2-pv1)),  cj     
    enddo  
    
    print*,' check the Foureier approximation of transvesal curve!'; print*; read*
    write(412,*); write(412, *)
  enddo 
     
  ! update the number of harmonics to use   
  call check_tail(isincr);   call upd_four(isincr, nf)
   
  print*, '****** Finish Fourier Approximation, check rho and tp ***********'; print*; ! read*; read*; read*
  
  close(412);  close(411) 
endif
 
 
! --- for the computation of the two basic frequencies 
open(133,  file = './dat/curve_lf/time/para_curve_time.dat', status = 'replace', access='append') 
write(133, *) ' # nf,  rho, tp, h0,  w2,  w1,  c0, ck(:, 1), sk(:, 1) ...  ck(:, 4), sk(:, 4)'

open(177, file = './dat/curve_lf/time/fcs_refn.dat', status='replace', access='append')
write(177, *) '# refined curve: c0 // ck // sk  format -- (10e24.14), rho(radian)=',rho


!subroutine tori_cont_time( ds, ntori, opt0, opt1, fcs_refn,  fcurve_time_refn,  & 
!       fpara_curve, ftori_time_refn, gr_map, deriv, gr_cj, deriv_cj)
ds     = 1.d-3
ds_max = 1.d-1
ds_min = 1.d-6

ds     = 1.d-4
ds_max = 2.d-2

!ds_max = 1.d0

!if(tori == 22) then 
!  ds     = 1.d-4
!  ds_min = 1.d-6
!  ds_max = 1.d-3
!endif 

nf_max = 128
nf_min = 16
call init_cont(nf_min, nf_max, ds_min, ds_max)

ntori = 100  
dir   = -1.d0 

call tori_cont_time( ds, dir, ntori, 10, 11, itori, 177,  fcurve_time_refn,  & 
       133, ftori_time_refn, TimeMap_lf, gr_lf, gr_cjlf, deriv_cjlf)


 close(177); close(133); close(fcurve_time_refn); close(ftori_time_refn)
 stop 





end 


 




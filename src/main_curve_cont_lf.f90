program curve_time_cont_lf
! 2017-05-26 22:09:27  try to do the continuation further using the last curve of the obtained results 


! the eigenvalue (rho) of centre part (pure imaginary) serves as the rotation number (in radian),
! the two frequencies are close to:  rho/T, 2pi/T. 

!  --- for the transversal case: 
! the new T and rho are: δ = 2pi / (±ν +2pi j) T,  rho = δ * 2pi / T =  (2pi)^2 / (±ν +2pi j) * T 
! and the two basic frequencies are close to:  rho/δ,  2pi/δ. 

! TODO we have to check if we need to multiply the eigenvalue by 2pi as the rotation number 

! ADD the multiple shooting, since the saddle part makes it hard to converge if its modulus is much greater than 1.


!  1. read from para_curve_time to get nf need for each curve, better use the last one, or we can take what ever we want 
!  2. do the continuation from the last obtained curve 
!  3. Globalize the curve to an invariant torus, by integrating NP = 2nf +1 points for one period

!  --- data saved in ./dat/curve_lf/time/  subfolder 
!  1. tori_cont2.dat            -- the original orbit of the 2D torus  ! not used at the moment 
!  2. curve_cont2.dat       -- the invariant curve obtained by Time map
!  3. fcs_cont2.dat              -- fourier coefficients, For each coordinate, three parts: c0 // ck // sk


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


integer, parameter  :: nf_init   = 32 ! 16  ! 32 ! 32   ! 16 !8! 16 !64 !8 !32  !64 ! 512 ! 128 !   !64 ! 32  !8 
integer, parameter  :: m0        = 1   ! test with 2 

! -- fourier coefs for Fourier analysis for the initial X_0  -- 
real(kind=dp) :: csf0  
real(kind=dp), allocatable,  dimension(:)   :: csf, sif  ! fourier coefs of 1 component

! -- fourier  coefs for Refinement  -- 
! the number of Fourier modes we are going to use will be specified 
! after we check the maximum norm of ck and sk, all the following allocatable arrays have nf-columns  

integer       :: nf  
                    
!  ------ Local Variables ----------- 
real(kind=dp) :: pv0(ndim), pv1(ndim), pvf(ndim), cj0,     & ! torus
                 tp0, poinst(ndim),  mmat(ndim, ndim), wr(ndim), wi(ndim), vr(ndim, ndim), & ! MM 
                 v1(ndim), v2(ndim), v10(ndim), v20(ndim),  xi, dxi, epsl,    & ! linear flow 
                 xpc_all(np,ndim),    & ! poinc 
                 pv_in(np), tol, tol_err   ! Netwon method
               
                  
                                                      
integer ::  i, ivr, i_eig, ivr1, ivr2,  ind0,   iscurve,  nitmax, opt, isincr

  
integer           ::  ftori_time, fcurve_time, fcurve_time_refn, ftori_time_refn, fpara, fcs        
character(len=70) ::  fntori_time, fncurve_time, fncurve_time_refn, fntori_time_refn, fnpara, fncs   


!---------- Transversal curve --------------           
integer :: istrans        
real(kind=dp)     ::  pv0_new(ndim), dt, tp1, rho0,  eta,  tdir, & 
                      pvf2(ndim),  xpc_all_trans(np, ndim), phi(ndim, ndim), pvil(nvar), pvfl(nvar)


! -- Fourier coefficients 
!real(kind=dp)              ::     step
integer       :: nf_min, nf_max, ntori 
real(kind=dp) :: ds, ds_min, ds_max, dir 
   
!-------------------------- continuation -------------------------------
integer       ::  j, k,  tori, itori, im, nrow, row_st !   ind_c00,
real(kind=dp) ::  cj, w1, w2           
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
 
istrans = 0 

! -- Initialization for the module  curve_mod, use poincare map approach, n=ndim-2
nitmax  = 20    
tol     = 1.d-10
tol_err = 1.d-10 

call init_curve(ndim, m0, nitmax, tol, tol_err) 

! open the old file, pay attention to the corresponding folder 
fcurve_time_refn = 28;   fpara = 30; fcs = 31

!fncurve_time_refn =  './dat/curve_lf/time/curve_time_refn.dat'
!fnpara            =  './dat/curve_lf/time/para_curve_time.dat'
!fncs              =  './dat/curve_lf/time/fcs_refn.dat'

!fncurve_time_refn =  './dat/curve_lf/time/torus22/curve_time_refn.dat'
!fnpara            =  './dat/curve_lf/time/torus22/para_curve_time.dat'
!fncs              =  './dat/curve_lf/time/torus22/fcs_refn.dat'

fncurve_time_refn =  './dat/curve_lf/time/torus11/curve_time_refn.dat'
fnpara            =  './dat/curve_lf/time/torus11/para_curve_time.dat'
fncs              =  './dat/curve_lf/time/torus11/fcs_refn.dat'

!fncurve_time_refn =  './dat/curve_lf/time/torus11_trans/curve_time_refn.dat'
!fnpara            =  './dat/curve_lf/time/torus11_trans/para_curve_time.dat'
!fncs              =  './dat/curve_lf/time/torus11_trans/fcs_refn.dat'


!fncurve_time_refn =  './dat/curve_lf/time/torus22/cont1/curve_cont2.dat'
!fnpara            =  './dat/curve_lf/time/torus22/cont1/para_cont2.dat'
!fncs              =  './dat/curve_lf/time/torus22/cont1/fcs_cont2.dat'



open(fcurve_time_refn, file = fncurve_time_refn, status = 'old')
read(fcurve_time_refn, *);

open(fpara,  file = fnpara,   status = 'old')
read(fpara, *)

open(fcs,   file = fncs,   status = 'old')
read(fcs,*) 
 
    
    
! files to save the original torus + invariant curve + interpolated one + approximated one 
ftori_time       = 21;  fntori_time       =  './dat/curve_lf/time/tori_cont2.dat'
fcurve_time      = 25;  fncurve_time      =  './dat/curve_lf/time/curve_cont2.dat'
open(ftori_time,  file = fntori_time, access ='append', status = 'replace')
write(ftori_time,*)  ' # Continued tori:  xi   (x y z vx vy vz)   cj'

open(fcurve_time, file = fncurve_time, access ='append', status = 'replace')
write(fcurve_time,*)  ' # Continued curve: Time + rho  curve:  xi   (x y z vx vy vz)   cj'
    
 
! ---------- Computation -----------------------------
! Outer bound for escape, nondimensinal unit 1 is big enough....
!  ****   I.C. of torus 1 : x-y-z-vx-vy-vz-hn  normal case  *****

!    1. take an unknown periodic orbits, compute the MM, and associated eigenvalues and eigenvectors --- not converge
! norm1, only one family, centre-saddle
!    0.219662920343E+01    0.000000000000E+00    0.932165693007E+00    0.701220365384E+00    0.460454350218E+00    0.000000000000E+00    0.000000000000E+00   -0.179869260122E+01    0.133226762955E-14  0.1600E-01   -1
!tori    = 10
ind0    = 1
tp0     = 0.219662920343d1
poinst  = (/0.d0, .932165693007d0,  0.701220365384d0, 0.460454350218d0, 0.d0, 0.d0/) ! centre X centre  
!istrans = 1   
   
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

! up to w1/w2 = 8.99
!tori = 21    
!ind0 = 3    
!tp0  = 0.219773264803d1  ! t2 from fft: 2.1977326477769235 
!poinst = (/0.699131986266d0, 0.d0, 0.d0, 0.d0, -0.158071934111d0,  0.553048316996d0/) ! centre saddle to start with 

! 0.69913198626600E+00    0.21869538702075E-18   -0.26463317253161E-20   -0.33891657647829E-20   -0.15807193411100E+00    0.55304831699600E+00

! good this one no matter ind0=2 or 3  -- centre X centre 
!    0.903410306172E+00    0.338496551294E+00    0.000000000000E+00    0.599636604742E-01   -0.596544282138E-01    0.135868534036E+01    0.350506167261E+00    0.400858698754E+01   -0.266453525910E-13  0.1600E-01   -1 # torus

! up to w1/w2 = 7.99 --- norm, eq3,  1,3
!tori   = 22 
!ind0   = 2
!tp0    = 0.903410306172d0
!poinst = (/0.338496551294d0, 0.d0, .599636604742d-1, -0.596544282138d-1, 0.135868534036d1,0.350506167261d0/) 

if(istrans == 1) then 
  open(411, file='./dat/curve_lf/time/torus_trans2.dat', status='replace', access='append')
  write(411,*) '# The transverse torus! t, x,y,z,vx,vy,vz'

  open(412, file='./dat/curve_lf/time/curve_trans2.dat', status='replace', access='append')
  write(412,*) '# The transverse curve. curve, curve_rho,  \varphi(eta) ! eta, x,y,z,vx,vy,vz'
endif 


! *************** Initialize  C0,CK,SK *************
call init_ind_c0(ind0)

dxi = pi2/np  ! take this value to keep coherent with the original points 

! read from the old para_curve_cont.dat file for nf, and curve_time_refn.dat for c0,ck,sk 

! write(fpara_curve, '(I5, 11e24.14, 1f5.0)', advance='no')  nf, rho, tp, h0, pi2/tp, rho/tp,  c0(1:n), dir
!  do i = 1, n, 1
!    write(fpara_curve, '(12e24.14)', advance='no')  ck(:, i),  sk(:, i)
!  end do 
 
!  do i = 1, n, 1
!    write(fcs, '(6e24.14)')   c0(i);     write(fcs,*) 
!    write(fcs, '(10e24.14)')  ck(i, :);  write(fcs,*) 
!    write(fcs, '(10e24.14)')  sk(i, :) 
!    write(fcs,*);    write(fcs,*) ! two blank lines to seperate the components in phase vector
!  end do

! first do awk to get para.dat for only 12 parameters
call alloc_arr_1d(n, c0)

itori = 1 
do 
!  if(itori > 62) goto 99 
!  if(itori > 30) goto 99  ! tori =  11
  read(fpara,*, end=99)  nf, rho, tp, h0, w1, w2,  c0(1:n), dir
  print*, 'itori, nf, rho, tp, h0', itori, nf, rho, tp, h0
  call upd_four(0, nf)
  
  do i = 1, n, 1
    read(fcs, '(6e24.14)')   c0(i);     read(fcs,*) 
    read(fcs, '(10e24.14)')  ck(i, :);  read(fcs,*) 
    read(fcs, '(10e24.14)')  sk(i, :) 
    read(fcs,*);     read(fcs,*) ! two blank lines to seperate the components in phase vector
  end do
  
  itori = itori + 1
enddo

99  itori = itori - 1 

print*, 'nf, itori =', nf, itori 
print*,'check c0: ', c0 ; read*
 
 
! --- for the computation of the two basic frequencies 
open(133,  file = './dat/curve_lf/time/para_cont2.dat', status = 'replace', access='append') 
write(133, *) ' #Continued family:  nf,  rho, tp, h0,  w1,  w2,  c0, ck(:, 1), sk(:, 1) ...  ck(:, 4), sk(:, 4)'

open(177, file = './dat/curve_lf/time/fcs_cont2.dat', status='replace', access='append')
write(177, *) '# Continued family: c0 // ck // sk  format -- (10e24.14), rho(radian)=',rho


ds     = 1.d-3
ds_max = 1.d-1
ds_min = 1.d-5

!ds     = 1.d-3
!ds_max = 1.d0
ds_min = 1.d-6
nf_max = 256


!ds     = 1.d-2
!!ds_max = 1.d0
!if(tori == 22) then 
!  ds     = 1.d-4
!  ds_min = 1.d-6
!  ds_max = 1.d-3
!endif 

!nf_max = 128
nf_min = 16
call init_cont(nf_min, nf_max, ds_min, ds_max)

ntori = 100  
dir   = -1.d0 

call tori_cont_time( ds, dir, ntori, 10, 11, itori, 177,  fcurve_time,  & 
       133, ftori_time, TimeMap_lf, gr_lf, gr_cjlf, deriv_cjlf)


 close(177); close(133); close(fcurve_time); close(ftori_time)
 stop 




end 


 




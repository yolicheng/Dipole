program main_po_ck
! ** NOTE ** 
! For each family, the parameter for error control and continuation steps should be different 

! 2017-05-16 18:31:30  
! Compute the bifurcated family from norm3_ifam2   first branch. 
!    0.216263759969E+01    0.701006714321E+00    0.000000000000E+00    0.000000000000E+00    0.000000000000E+00   -0.148214724511E+00    0.481970030333E+00    0.407300824595E+01    0.000000000000E+00  0.6400E-01   -1
 

! 2017-05-07 00:41:25 
! check if the det(M) = 1 for Monodromy matrix of the periodic orbit. 


! 2017-03-29 09:26:50 
! Add the trace computation and sort part, or go to main_trace.f90 to do another step 

! 2017-02-16   10:36:20 
! better to do small routines, that this main program only deal with the symmetric p.o.s
! compute all the families, and their ending, do also the characteristic curve. 

!
! 2017-02-14 16:57:35 
! Try to update and refine the routines but give up, better to  finish the draft somehow with the old figures,  
! the replacement of finer figures will be done later if we have enough time.
! do small routines, here, we only deal with symmetric p.o.s

! Ask for the input to specify the case, equilibirum point and the family, in such a way that we can call the same main routine  
 
! The normal case with N = [0 0 1] possesses the full symmetry, 
! and is perfect for use with 4 symmetric p.o. around 4 symmetric equilibirum points 
! we will look at the second kind of equilibirum point, the real part of while is also less than 1 in modulus for all beta, quite stable 
! start with, we take beta = 1 

!  the third one could also be good to use for formation flying, but the modulus of the real part is always great than 1.


! 20160406  test asymmetric approach, which do the refinement after 1 revolution w.r.t the initial condition

! 2016018  add the anglectr to control the curvature of the characteristic curve, make sure it is a smooth curve... 
!          after the test, put the conclusion here, whether it is good to do the control on the curvature.
!       
use dp_mod
use lf_mod ! the system related parameters of lorentz force problem
use po_mod ! archived subroutines to compute families of symmetric p.o.
use pi_mod
use poinc_mod

implicit none

integer, parameter ::  ndim = 6       ! the dimension of LF problem 
integer, parameter ::  neq  = 42,  &  ! compute also the variational matrix, 6+36=42
                       npo  = 100 ! 200 !600 !1000 !1000  !600    ! 60,  52... for x=z, ifam2  !   1 !35 ! NO of orbits in pofam 

! Global  Variables, the assignment in the main routine can use an appendix 0 to avoid confliction. 
!lf_mod:      beta, cs, sgn, eq    
!poinc_mod:   p0, tmax, xmax, h0 
!             n,  ind_vel, ind, dir, imax  
   

real(kind=dp) ::  beta0, p00, tmax0, xmax0  
 
integer      ::   ind0, dir0, imax0,  &  !poinc
                  dir_curve, idic, idfc, dsctr, anglectr, idcor, isarc, debug 

real(kind=dp) ::  tol, prsc       ! error control

! Local Variables
 
integer :: i, j, ispl, ivr, fpo, fpoinst, ifam, fmmegv, &
           tdir, ipo, isgn, idsgn, & 
           ftrace, ftrace_st, find_st, ind_st(6) !sort MM_eig

real(kind=dp) ::  wr(6), wi(6), vr(6, 6),   &   ! eigenspace of variational matrix  
                  poinst(6), ds, vrchs(6), ynew(npo,8), epsl_po, cj,  & !pofam
                  po0(6), tpo, pof(6), & ! plpo
                  dlf(6, 6), mmat(6, 6), wr_mm(6), wi_mm(6), vr_mm(6, 6),  det, &  ! MM 
                  tr_re(6), tr_im(6),  wr0(6), wi0(6), tr_str(6), tr_sti(6),  & ! sort trace
                  tp0, csangle_min, ds_max, ds_min, poinst1(6), poinst2(6) 

                 
 character(len=70) :: fnpo, fnpoinst, fnmmegv,  fntrace, fntrace_st, fnind_st             

real(kind=dp), external :: dnrm2 ! from package... 

! TODO 
! ifam = 3, stuck at 25th orbit, for the poincare map crossing, t=-2.d-4 fail to finish rk78, 

! ifam = 1, stuck at 36th orbit,    poinc: t, y
! 1.8709195802188187E-004  -4.9322318364968389E-021   9.2589610217307418E-003  -1.7629589633198005E-004   1.0554145798777161       -2.8607775325232556E-003
! the orbit is too small.... stop the continuation.

! debug or not 
debug = 0

dir_curve = -1  ! plus family 

!test dir_curve =1 for the orbits after the second bifurcation for norm2-ifam2, bif3 
!dir_curve = 1  

! automatic step control depending on the number of iteration for Newton Method
dsctr = 1!  1

! contol the angle between the vector field of two latest consecutive points along the 
! characteristic curve
anglectr = 1! for n3, fam2, no control is better  anglectr = 0, goes to the first bifurcated family. 

! rad3, fam1, smooth, anglectr = 1, 1.d-11, dir_curve = -1
! rad3, fam2, dir_curve  = 1, terminate to a termination ifam2_end/ ??? or continue further to strange orbits ifam2/ 

! rad3, fam3, dir_curve = -1
!epsl_po = 1.d-4;      ds = 1.d-4; ds_max = 1.d-3 

! rad2, fam1, dir_curve = 1

! for the moment, we only use the first one 
isarc = 1 ! along arc-length paramter --- only this one is used so far 
!isarc = 2 ! along period
!isarc = 3 ! along energy

! Todo---- Check really carefully of the subroutine dsdecr.... this is where it is wrong? 
! there is still something wrong here, in fact if anglectr = 0, it is supposed to obtain regular p.o.s with the possibility of jumping to another family ... 
! I'm not sure why it won't happen.....
! and I'm so regetful that I didn't save the version with the good result ... 


! ---- for symmetric p.o. : 1-3 ------------ 
idcor = 1 ! for symmetric approach with poincare map  --- most used one...


!idcor = 2 ! symmetric + free time shooting 

! rad2, x-z, y=0, only one symmetry wrt y=0plane can be used, 
! --- asymmetric p.o. 
! idcor = 4 ! Poincare map approach, section y(ind) = x0_fxd, 
! idcor = 5 ! additional constraint to fix one component, y(ind) = x0_fxd, the current used one
 
! idcor = 6 !   
!idcor = 7 ! improved version of idcor=5   


!! the minimum of the cosine of the angle between two consecutive vectors
! cos(0.1 rad) = 0.995;    cos(0.37 rad) = 0.95  
 csangle_min = 0.995d0  

! csangle_min = 0.95d0  

! bound vaules of stepsize ds
! With the control of termination added when the arc step goes too smaller... 
! better not to be too big, with 2.d-2, or 1.d-2, for now is good choice  

! this boundary matters a lot  --- non-symmetric approach
! for case1, eq3, ifam2 we do ds_max = 1.d-2; ds_min = 1.d-5; anglectr=0, isarc=1, idcor=5, ind =2, po=0 
! ds_max   = 1.d-2;  ds_min = 1.d-5;  
! anglectr = 0;      isarc = 1;  idcor = 5;  ind0 = 2;  p00 = 0.d0;  ! dir_curve = -1; ! npo = 600
!! **NOTE** only the non-symmetric approach with big ds_max as 1,d-1 or 2.d-2 could jump the bifurcation 

! eq3, norm, ifam2, symt, indc=indf=4, ind=3, only the first family, fail to arrive to the second family.
!ds_max = 1.d-2
!ds_min = 1.d-6 ! but this bound is not working for angle control....
!   

!   ..............ifam1, ds_max = 2.d-2; ds_min = 2.d-5        ! npo = 500 , dir_curve = -1  
! rad3, ifam1, x=-z, very far... 
!ds_max = 4.d-2
!ds_min = 2.d-5 ! but this bound is not working for angle control....

! rad3,  ifam2, x=-z, dir_curve = 1 
!ds_max = 1.d-3
!ds_min = 2.d-5 ! but this bound is not working for angle control....
!!    

! rad3, dir_curve = -1, ifam1, ds_max = 2.d-2 

! rad2, 2 family1, one smooth,connections,  npo = 500, enough, ds_max2.d-1, dir =1  
! rad2, fam2, end nowhere.  npo = 1500, enough, ds_max 8.d-1 , dir = -1

! tan 1 
ds_max = 1.d-1
ds_min = 1.d-5 ! but this bound is not working for angle control....
 
! tan_2 -ifam1  1.d-2 
ds_max = 4.d-2
!ds_min = 1.d-5 ! but this bound is not working for angle control....

!norm2, ifam2, detect the termination, we have to use very small steps 
ds_min = 1.d-5 

ds_max = 1.d-2

!ds_max = 4.d-3  

!ds_max = 1.d-3  
ds_max = 2.d-3 
 
!ds_max = 4.d-2 ! for postbif3, but actually stop before bif3

! if we take ds_max = 1.d-3, we go to bifurcation at ipo=55, the first crossing with Re(tr)=2 
! if we tkae ds_max = 4.d-3, we jump the first bifurcation and goes to the next one ...
! now, let's try bigger values to see if we can pass the second bifurcation too 

! ds_max = 2.d-2, we manage to continue close to the second bifurcation. with bifurcation... (from the very beginning guess)
! ds_max = 1.d-2, goes to the second bifurcated family (from the very beginning guess)


! case1, eq1, the same as ifam1, dir_curve = 1, ds_max = 2.d-1; ds_min = 2.d-5

! case1, eq2  ! ivr = 1, ivr = 5
!if(cs == 1 .and. ieq == 2) then 
!  dir_curve = -1;  anglectr = 1;         ! ivr = 1, ivr = 5
!  ds_max = 1.d-3
!  ds_min = 1.d-6  
!endif  


if(cs == 2 ) then 
  if(ieq == 3 .or. ieq == 1) then
     idcor = 1; anglectr = 1
  else
     idcor = 5; anglectr = 0
  endif  
endif 

print*, 'idcor,anglectr', idcor,anglectr
  
call init_debug(debug, dsctr, anglectr, csangle_min,  isarc, ds_min, ds_max)

! beta, the ratio between the angular velocity of mean motion of the chief around the earth and 
!       the angular velocity of the rotaion of the deputy

! Intialize the state vector of the equilibirum points 
call init_lf  

ind0 = 0; ivr = 0

!if(cs == 1) then   ! normal case... 

!  if(ieq == 2)  then 
!    ind0  = 1;  ivr = 6 
!  endif    
!  
!else if(cs == 2) then  ! radial case 
!! TODO 2017-02-21 11:18:27 
!  ! simply-symmtry: P2: x-z plane, ind = 2, y = 0
!  ! TODO : eq1: 0-0-z 
!  if(ieq == 3 .or. ieq == 1 )  ind0 = 2 ! eq3: x-0-z, simple 
!  
!  if(ieq == 2 )  then  ! no symmetry in the eigenvectors
!     idcor = 5  ! eq2: x-y-0, general asymmetric approach and use 
!  endif 
!  
!endif

if(ind0 == 0) then 
  print*, 'Input ind for Poincare section, X(ind) = p0, ind = ' 
  read(*,*) ind0
endif


! Instead of pass the value of beta from routine... 
! try read from the screen, which is more flexible, and we will not mix different cases.   

print*, 'Please input beta (q/m):'
read(*,*) beta0 
call init_beta(beta0)   
    
print*, 'check, cs, ieq, beta', cs, ieq, beta 
read* 
    
! At the moment --- skip this real unit --- TODO  
! Also compute the unit of distance, time and velocity, which are all function of beta
! but at this moment, we will not use the real quantity of the units. 
! call lfunit(beta, runit, tunit, vunit)
    
print*, 'Do you want to change the sign of equilibirum point (Yes = 1) ?'
read(*,*)  isgn
    
! If the equilibirum points with different signs are not symmetric,  
! the periodic orbits will also be different, we have to compute seperately 
if(isgn == 1) then 
  idsgn     = mod(ind0, 3) + 1 
  eq(idsgn) = -eq(idsgn)
endif     
print*, 'New eq:', eq
read* 

! try a different eq : we have four symmetric one, the first two where x=z are already done(x,z), (-x,-z)
! ---- for ieq = 3, --- the third one: x, -z
!eq(3) = -eq(3) 

 
! Jacobi matrix of the lorentz force with respective to the state 
call dflrtz(eq, dlf)

! check the energy 
call gr_cjlf(eq, cj)
print*,'check energy!, cj, ieq,', cj, ieq, eq 
read*

!print*, 'The Jacobi Matrix'
!do i = 1, ndim
! write(*,'(6f8.4)') dlf(i,:) 
!enddo
!print*; read*

! Compute the eigenvalues and eigenvectors of dlf
call eigrg(dlf, ndim, 1, wr,wi, vr)

!if(ivr == 0) then 
  print*, 'Which column to use for initial guess of p.o.?'
  read(*,*) ivr 
!endif 

vrchs = vr(:, ivr)
print*, vrchs 

vrchs = vrchs/dnrm2(6,vrchs,1)
print*, dnrm2(6,vrchs, 1), vrchs

read* 


! commomly used  varaibles-- general case
!epsl_po  =  1.d-3;      ds = 1.d-3;
tol     = 1.d-10 ;  prsc = 1.d-10

epsl_po = 1.d-3;      ds = 1.d-3;
!tol  = 1.d-11 ;  prsc = 1.d-11  
!tol  = 1.d-12 ;  prsc = 1.d-12 

! for  ! case1, eq2, ivr = 5, that is ifam2, we have to take higher precision to detect the termination of the family.
if(cs == 1 .and. ieq == 2 .and. ivr == 5  ) then 
!  epsl_po = 1.d-4;   ds = 1.d-4;
  epsl_po = 1.d-3;   ds = 1.d-3;     ! test for postbif2, norm2 ifam2
!  epsl_po = 1.d-2;   ds = 1.d-2;    ! test for postbif3, norm2 ifam2    
!  tol  = 1.d-11 ;  prsc = 1.d-11
endif 


! the initial guess for po, move along the corresponding eigenvector for a small distance epsl_po 
poinst =  eq   +  epsl_po * vrchs !  
print*, 'poinst=', poinst 


! For the simply-symmetric periodic orbits, the  index of zero-component 
! in initial and final configurations are the same.  idic = idfc = ind 
if(idcor <= 3) then 
  idic = ind0;   idfc = ind0
endif 

!if(cs == 1 .and. ieq == 3 .and. ivr <3  ) then 
!  idic = 4; idfc = 4  !  
!endif 

!idic = 4; idfc = 4 ! normal, eq3   ! ---discard.... 
!ind0 = 2
!--------------------------------------------------------------------------------

!  ***************** compute the 3 families of PO ***************************
! Initialization of parameters to do numercial continuation
dir0  = 0  ! the first crossing regardless of the direction of the velocity 
tdir  = 1 ! integrate forward

! finally take 1.d-9 as the error control
! tol  = 1.d-10 ! err tolerance of target f to terminate Newton Method, from the Gerard's book, it should be 1.d-16?
! prsc = 1.d-10 ! For lf problem,  1.d-11 is too small for the second family, ok for the other 2 families


if( mod(ivr,2) == 1) then 
  ifam = ivr/2+1
else 
  ifam = ivr/2
endif   


print*, 'check wi to determine which family to study'
print*, 'wi=',  wi
print*, 'ifam=', ifam  
! Use the period as initial guess, TP = 2*pi/wi, wi is the pure imaginary eigenvalue
tp0 =  2*pi/dabs(wi(2*ifam)) ! the maximum time for the first return to poincare section
print*, 'tp0=', tp0, 'wi=', wi(2*ifam)
read*
 
! do a trick for norm3, fam2, complicated family, with a bifurcation 
! use a refined p.o.  
if(cs == 1 .and. ieq == 3 .and. ivr >=5 .and. idcor == 1) then  ! symmetric approach

!    0.198691879704E+01    0.693361395984E+00    0.000000000000E+00    0.790548520772E-09    0.641695286943E-12   -0.128583970319E-05   -0.999992904571E-03    0.432674771093E+01    0.888178419700E-15  0.1000E-02   -1
  tp0 =   0.198691879704d1
  poinst = (/0.693361395984d0, 0.d0, 0.d0, 0.d0, -0.128583970319d-5,  -0.999992904571d-3/)

!    0.198724757056E+01    0.693396317399E+00    0.000000000000E+00   -0.345649574372E-13   -0.668451317761E-14   -0.370741431022E-03   -0.169902896835E-01    0.432645991458E+01    0.000000000000E+00  0.4000E-02   -1
!tp0 =    0.198724757056d1
!poinst = (/ 0.693396317399d0, 0.d0, 0.d0, 0.d0, -0.370741431022d-3,  -0.169902896835d-1/)

  idic  = 4;  idfc = 4 ! normal, eq3, fam2
  imax0 = 1
endif  
 

!if(cs == 2 .and. ieq == 2 .and. ivr > 3 )then 
!if(cs == 2 .and. ieq == 2   )then 
!!  0.56641973494575E+01  0.50779432783733E+00  0.71327546262244E+00  0.78057249709902E-01  0.95872798304379E-02  0.54062183698304E-01 -0.59123787741211E-01  0.22769301518562E+01  0.44408920985006E-15  0.00500000    1 

!tp0   = 5.6641973494575d0
!poinst = (/0.50779432783733d0, 0.71327546262244d0, 0.78057249709902d-1, & 
!           0.95872798304379d-2, 0.54062183698304d-1,  -0.59123787741211d-1/)

! ---check the termination of the last orbit...
!    0.494069644968E+01    0.429219426309E+00    0.713275462622E+00    0.422292367972E-01   -0.215974185230E+00   -0.674220592223E+00    0.379991680857E-03    0.181287323131E+01   -0.222044604925E-14  0.1250E-03    1

! tried, cannot continue further... do not know why...
!tp0 = 0.494069644968d1
!poinst = (/0.429219426309d0, 0.713275462622d0, 0.422292367972d-1,  & 
!          -0.21597418523d0,-0.674220592223d0, 0.379991680857d-3/)
!  
!          
!call gr_cjlf(poinst, cj)
!print*, 'Check the energy: ', cj, 'cj0=0.22769301518562E+01'
!print*; read*
!endif 
 
! plot the initial guessed orbit  to see if it is a periodic orbit
open(333, file = './dat/po_guess.dat', access ='append',status='replace')
call plob(poinst, 0.d0,  tp0,  6, 6, 1, 1, 333 , gr_lf, gr_cjlf, pof)  
print*, 'check the initial guessed orbit!'; close(333); read*


! We have to add the general approach here, since not all the periodic orbits can be computed by applying the symmetry.

!  check which initializaiton routine to call  2017-02-19 16:03:36  
! checking idcor == 5... 
if (idcor  == 1  )  then !  for  idcor = 1, 2, 3, TODO 3
! ---- symmetric p.o. ----
  call init_symtpo(idic, idfc, idcor, ind0)    
  p0  = 0.d0;     imax0 = 1 
  xmax0 = 40.d0;  tmax0 = 100.d0 
  
elseif(idcor == 5 .or. idcor == 7 .or. idcor == 6)  then 

! ----- asymmetric p.o. ---
   print*, 'General approach for p.o. computation! Input ind,  add constraint:  p0 = poinst(ind)'
   read(*,*) ind0 
   
   p0 = eq(ind0)

   call init_asymtpo2(idcor, ind0, p0)   

endif 

if (idcor == 1 .or. idcor == 4 ) then  ! TODO 4 to check 
  call init_poinc( ind0, p00, dir0, imax0, tmax0, xmax0, ndim, cj)
endif 

call init_errctr(tol, prsc)
call init_writedx 



! save the data of po and  its initial state 
fpo = 20;  fpoinst = 21;  fmmegv = 33; 

! use new name, to avoid overwritten the good ones, and the type of approach selected is used as suffix
write(fnpo,    fmt='(a,i0,a,i0,a)')   './dat/po/case', cs, '/eq', ieq, '/po.dat'  
write(fnpoinst,fmt='(a,i0,a,i0,a)')   './dat/po/case', cs, '/eq', ieq, '/poinst.dat'
write(fnmmegv, fmt='(a,i0,a,i0,a)')      './dat/po/case', cs, '/eq', ieq, '/pommegv.dat'
print*, fnpo, fnpoinst, fnmmegv  

open(fpo,file = fnpo, access = 'append', status='replace')
write(fpo,*)  '# t      (x y z vx vy vz)        cj'

open(fpoinst,file=fnpoinst, access ='append',status='replace')
write(fpoinst,*) '# TP  I.C.(x y z vx vy vz)    CJ    DCJ   ds'
write(fpoinst,*) '# beta, eq: ', beta,  eq 
write(fpoinst,*) '# ivr,  wi: ', ivr, wi ; write(fpoinst,*)

open(fmmegv,    file = fnmmegv,    access  = 'append', status='replace') 
write(fmmegv, *) '# Eig of MM:   real/imag'


!----- compute and sort trace  ---------------
ftrace = 34; ftrace_st = 35; find_st = 36
write(fntrace,    fmt='(a,i0,a,i0,a)')   './dat/po/case', cs, '/eq', ieq, '/trace.dat'  
write(fntrace_st,fmt='(a,i0,a,i0,a)')   './dat/po/case', cs, '/eq', ieq, '/trace_st.dat'
write(fnind_st, fmt='(a,i0,a,i0,a)')      './dat/po/case', cs, '/eq', ieq, '/trace_ind_st.dat'
print*, fntrace, fntrace_st, fnind_st
open(ftrace, file = fntrace, access  ='append', status='replace') 
write(ftrace, *) '# Trace:',  n, '(real / imag)'  

open(ftrace_st, file = fntrace_st, access  ='append', status='replace') 
write(ftrace_st, *) '# Sorted ', n, 'trace (real / imag)'  

open(find_st, file = fnind_st, access  ='append', status='replace') 
write(find_st, *) '# Sorted index of trace, original: [1, 2, 3, 4, 5, 6]'

!test for norm2/ifam2/ doesn't make any sense'

! continue after the second bifurcation.  
! -- the last one,   
!    0.715798941242E+01    0.242931723939E+00    0.000000000000E+00    0.521686327546E+00    0.000000000000E+00    0.694907245389E-01    0.000000000000E+00   -0.719261054605E+00   -0.177635683940E-14  0.1280E-01    1

!poinst = (/0.242931723939d0, 0.d0,  0.521686327546d0, 0.d0,  0.694907245389d-1, 0.d0/) 

! -- the second last one.. 
!    0.712873929449E+01    0.247725648607E+00    0.000000000000E+00    0.533042967665E+00    0.000000000000E+00    0.662892435444E-01    0.000000000000E+00   -0.708783505886E+00   -0.133226762955E-14  0.1280E-01   -1

!poinst = (/0.247725648607d0, 0.d0,  0.533042967665d0, 0.d0,  0.662892435444d-1, 0.d0/) 
 ! start from this point using ds = 1.d-2, and ds_max = 1.d-1, finally, we obtain families after the third bifurcation... but we have to use dir_curve =1 
 ! and first go back to the equilibirum, then back again. 
 
! --- finally the one returned to the equilibrium from the above p.o., is given in below  -- fail to start from here.  
!     0.722415492654E+01    0.232967499759E+00    0.000000000000E+00    0.498526460005E+00    0.000000000000E+00    0.777286380215E-01    0.000000000000E+00   -0.743199445954E+00   -0.333066907388E-14  0.8000E-03    1
! poinst = (/0.232967499759d0, 0.d0,  0.498526460005d0, 0.d0,  0.777286380215d-1, 0.d0/) 
 

! --- after the second bifurcation  ----  start from the last three lines at bt2_postbif2/poinst.dat 
!      0.719771459954E+01    0.236963614117E+00    0.000000000000E+00    0.507390540845E+00    0.000000000000E+00    0.710498705680E-01    0.000000000000E+00   -0.733530771854E+00    0.244249065418E-14  0.3200E-02    1
!    0.720621073538E+01    0.235733623012E+00    0.000000000000E+00    0.504492455903E+00    0.000000000000E+00    0.716260470212E-01    0.000000000000E+00   -0.736591379419E+00    0.888178419700E-15  0.3200E-02    1
!    0.721464556884E+01    0.234518515026E+00    0.000000000000E+00    0.501668476865E+00    0.000000000000E+00    0.725032805718E-01    0.000000000000E+00   -0.739635683681E+00   -0.888178419700E-15  0.3200E-02    1
! 
!poinst = (/0.236963614117d0, 0.d0,  0.507390540845d0, 0.d0,  0.710498705680d-1, 0.d0/) 
!poinst = (/0.234518515026d0, 0.d0,  0.501668476865d0, 0.d0,  0.725032805718d-1, 0.d0/)      

! -- before the third bifurcation... 
!   0.722303914699E+01    0.233222505860E+00    0.000000000000E+00    0.498911053569E+00    0.000000000000E+00    0.757143810237E-01    0.000000000000E+00   -0.742722497891E+00    0.310862446895E-14  0.5000E-03    1
!    0.722340492488E+01    0.233149244021E+00    0.000000000000E+00    0.498788840371E+00    0.000000000000E+00    0.761937322401E-01    0.000000000000E+00   -0.742869565016E+00   -0.288657986403E-14  0.5000E-03    1
!poinst = (/0.233149244021d0, 0.d0,  0.498788840371d0, 0.d0,  0.761937322401d-1, 0.d0/)      

! after the third bifurcation, before the fourth one. the last two 
!   0.699124351779E+01    0.373598691035E+00    0.000000000000E+00    0.851356709806E+00    0.000000000000E+00   -0.804651960928E-01    0.000000000000E+00   -0.659915452844E+00    0.133226762955E-14  0.4000E-01   -1
!    0.703539716892E+01    0.385656998514E+00    0.000000000000E+00    0.881443883720E+00    0.000000000000E+00   -0.103902177494E+00    0.000000000000E+00   -0.675543800951E+00    0.666133814775E-15  0.4000E-01   -1
!    0.708680156807E+01    0.397597278605E+00    0.000000000000E+00    0.910681968980E+00    0.000000000000E+00   -0.128443246901E+00    0.000000000000E+00   -0.693814441251E+00    0.144328993201E-14  0.8000E-01    1

!poinst = (/0.373598691035d0, 0.d0,  0.851356709806d0, 0.d0,  -0.804651960928d-1, 0.d0/)   ! -1
!poinst = (/0.385656998514d0, 0.d0,   0.881443883720d0, 0.d0,  -0.103902177494d0, 0.d0/)    ! -1
!poinst = (/0.397597278605d0, 0.d0,  0.910681968980d0, 0.d0,  -0.128443246901d0, 0.d0/)     !  1 -fail to pass the fourth bifurcation...


!    0.720874979171E+01    0.421248626454E+00    0.000000000000E+00    0.966908481309E+00    0.000000000000E+00   -0.180104513560E+00    0.000000000000E+00   -0.737506895488E+00    0.222044604925E-15  0.8000E-01   -1
!    0.735806311348E+01    0.444997975629E+00    0.000000000000E+00    0.102439984829E+01    0.000000000000E+00   -0.231559700090E+00    0.000000000000E+00   -0.793213318327E+00   -0.222044604925E-15  0.8000E-01    1

!poinst = (/0.421248626454d0, 0.d0,  0.966908481309d0, 0.d0,  -0.180104513560d0, 0.d0/)      
!poinst = (/0.444997975629d0, 0.d0,  0.102439984829d1, 0.d0,  -0.231559700090d0, 0.d0/)   ! too close to a bifurcation..  1 gives bif4, -1 gives bif5  
!  
   
! use the last two to interpolate for a middle point that might stay between the fourth-fifth bifurcations 
!poinst2 = (/0.421248626454d0, 0.d0,  0.966908481309d0, 0.d0,  -0.180104513560d0, 0.d0/)      
!poinst1 = (/0.444997975629d0, 0.d0,  0.102439984829d1, 0.d0,  -0.231559700090d0, 0.d0/)   ! too close to a bifurcation..    
!poinst =  (poinst1 + poinst2) / 2.d0     
      
! for the last point, we are able to continue for orbits in between bif4-bif5 
!    0.723618704537E+01    0.426239887951E+00    0.000000000000E+00    0.977778298056E+00    0.000000000000E+00   -0.192327479701E+00    0.000000000000E+00   -0.747417603158E+00    0.355271367880E-14  0.2000E-02   -1
!    0.723520735293E+01    0.426091822930E+00    0.000000000000E+00    0.977383076171E+00    0.000000000000E+00   -0.192059464133E+00    0.000000000000E+00   -0.747067436394E+00   -0.277555756156E-14  0.5000E-03   -1
!    0.723419491292E+01    0.425944486357E+00    0.000000000000E+00    0.976971985081E+00    0.000000000000E+00   -0.191816073433E+00    0.000000000000E+00   -0.746707362805E+00    0.000000000000E+00  0.5000E-03   -1
!    0.723313124502E+01    0.425799188847E+00    0.000000000000E+00    0.976536507251E+00    0.000000000000E+00   -0.191615497830E+00    0.000000000000E+00   -0.746332408535E+00   -0.122124532709E-14  0.5000E-03   -1

!poinst = (/0.426239887951d0, 0.d0,  0.977778298056d0, 0.d0,  -0.192327479701d0, 0.d0/)   !  use 1 to continue forward (reverse the original sense)
!poinst = (/0.425799188847d0, 0.d0,  0.976536507251d0, 0.d0,  -0.191615497830d0, 0.d0/)   !  use 1 to continue forward (reverse the original sense)
 
 
! for bif1 

!0.697516735976E+01    0.281164903912E+00    0.000000000000E+00    0.615262887524E+00    0.000000000000E+00    0.457060469438E-01    0.000000000000E+00   -0.654240221855E+00   -0.155431223448E-14  0.2500E-04   -1 ! still in old family 
!poinst = (/0.281164903912d0, 0.d0,  0.615262887524d0, 0.d0,  0.457060469438d-1, 0.d0/)   !  for bif1

!    0.697522798385E+01    0.281211912379E+00    0.000000000000E+00    0.615144950499E+00    0.000000000000E+00    0.444644376869E-01    0.000000000000E+00   -0.654273473166E+00    0.444089209850E-15  0.2000E-03   -1
!poinst = (/0.281211912379d0, 0.d0,  0.615144950499d0, 0.d0,  0.444644376869d-1, 0.d0/)   !  for bif1

!    0.695320070531E+01    0.282208853974E+00    0.000000000000E+00    0.611126576864E+00    0.000000000000E+00    0.122990532732E-05    0.000000000000E+00   -0.656771947253E+00    0.199840144433E-14  0.8000E-03   -1
! poinst = (/0.282208853974d0, 0.d0,  0.611126576864d0, 0.d0,  0.122990532732d-5, 0.d0/)   !  for bif1 --- fail to go further. 
 

! bif4-bif5, take inverse sense, dir_curve = -1 , 40-th orbit in /bt2_bif4/  and save data in bt2_po4_5, 
!    0.730806594360E+01    0.438270057500E+00    0.000000000000E+00    0.100458898457E+01    0.000000000000E+00   -0.220844063541E+00    0.000000000000E+00   -0.773503746863E+00    0.111022302463E-15  0.1250E-03    1

poinst = (/0.438270057500d0, 0.d0,  0.100458898457d1, 0.d0,  -0.220844063541d0, 0.d0/)   ! take dir_curve = -1

! bif3_4, 1-th in  ./bt2_postbif3_4/ save data in bt2_po3_4
!   0.723313124508E+01    0.425799188847E+00    0.000000000000E+00    0.976536507251E+00    0.000000000000E+00   -0.191615497830E+00    0.000000000000E+00   -0.746332408535E+00    0.266453525910E-14  0.1000E-02    1

poinst = (/0.425799188847d0, 0.d0,  0.976536507251d0, 0.d0,  -0.191615497830d0, 0.d0/)   ! take dir_curve = 1

! bif2_3, 1-th in  ./bt2_bif3/ save data in bt2_po2_3
!    0.715798941257E+01    0.242931723939E+00    0.000000000000E+00    0.521686327546E+00    0.000000000000E+00    0.694907245389E-01    0.000000000000E+00   -0.719261054605E+00   -0.299760216649E-14  0.1000E-02    1
poinst = (/0.242931723939d0, 0.d0,  0.521686327546d0, 0.d0,  0.694907245389d-1, 0.d0/)   ! take dir_curve = 1  

! bif1_2, last 4-th in  ./bt2_postbif1/ save data in bt2_po1_2
!    0.707507201735E+01    0.257373849133E+00    0.000000000000E+00    0.556435795031E+00    0.000000000000E+00    0.617238160434E-01    0.000000000000E+00   -0.689638145002E+00    0.777156117238E-15  0.1280E-01   -1  
! poinst = (/0.257373849133d0, 0.d0,  0.556435795031d0, 0.d0,  0.617238160434d-1, 0.d0/)   ! take dir_curve = 1  
 
! bif0_1, last 9-th in  ./bt2_postbif1/ save data in bt2_po0_1
!    0.697686637673E+01    0.280623767667E+00    0.000000000000E+00    0.613979314944E+00    0.000000000000E+00    0.465095616952E-01    0.000000000000E+00   -0.654839110962E+00    0.133226762955E-14  0.1280E-01    1   
!  poinst = (/0.280623767667d0, 0.d0,  0.613979314944d0, 0.d0,  0.465095616952d-1, 0.d0/)   ! take dir_curve = 1     
 
 
 
! check norm3, ifam2, the termination point A is a bifurcation point? and the following orbits are another family? 205-th
! I say this because there is a sharpe turn at A  
!    0.109100753893E+01    0.411995608951E+00    0.000000000000E+00    0.773465353168E-11    0.294653114168E-12    0.122286670473E+01   -0.153899844781E-01    0.386800198613E+01    0.488498130835E-14  0.1600E-01    1
tp0 = 0.109100753893d1
poinst = (/0.411995608951d0, 0.d0,  0.773465353168d-11,0.294653114168d-12,  0.122286670473d1, -0.153899844781d-1/)   ! take dir_curve = 1     

!# the first  one after the termination A, check G if small also?
!    0.109128971924E+01    0.412031631511E+00    0.000000000000E+00    0.149233990567E-02   -0.373697625411E-04    0.122273840739E+01    0.165740842135E-01    0.386784480642E+01    0.124344978758E-13  0.1600E-01   -1

tp0 = 0.109128971924d1
poinst = (/0.412031631511d0, 0.d0,  0.149233990567d-2, -0.373697625411d-4,  0.122273840739d1, 0.165740842135d-1/)      



! norm3_ifam2, bifurcated family 
!    0.216263759969E+01    0.701006714321E+00    0.000000000000E+00    0.000000000000E+00    0.000000000000E+00   -0.148214724511E+00    0.481970030333E+00    0.407300824595E+01    0.000000000000E+00  0.6400E-01   -1

poinst = (/0.701006714321d0, 0.d0,  0.d0, 0.d0, -0.148214724511d0, 0.481970030333d0/)         

! rad3_dif , ifam3 
!    0.391584934196E+02    0.417872128549E+00    0.000000000000E+00   -0.423850393832E+00    0.000000000000E+00    0.394703447367E-01    0.000000000000E+00    0.204660848512E+01    0.710542735760E-14  0.1600E-02    1
! #  next branch same family     
!0.404254693762E+02    0.404669042618E+00    0.000000000000E+00   -0.411601900332E+00    0.000000000000E+00    0.409369930174E-01    0.000000000000E+00    0.208207294830E+01    0.488498130835E-14  0.3200E-01    1
 
poinst = (/0.404669042618d0, 0.d0,  -0.411601900332d0, 0.d0, 0.409369930174d-1, 0.d0/)  
       
 
open(333, file = './dat/po_guess.dat', access ='append',status='replace')
call plob(poinst, 0.d0,  tp0,  6, 6, 1, 1, 333 , gr_lf, gr_cjlf, pof)  
print*, 'check the initial guessed orbit!'; close(333); read*

          
! ----------- do the refinement and continuation ---------------
! ds = 0.01d0
! subroutine pofam(yi, tp0, npo, dir,ds, fpoinst, ynew, ipo, deriv, gr_cj)
call pofam(poinst, tp0, npo, dir_curve, ds, fpoinst, ynew, ipo, gr_lf, gr_cjlf)
 
ipo = max0(ipo-1, 1) 

print*, 'PO finisned!'
print*, 'No of real p.o. computed = ', ipo 
read*

open(133, file = './dat/detm.dat', access='append',status='replace')


! --- plot the P.O ---
do i = 1, ipo
  
  po0 = ynew(i,2:7)
  tpo = ynew(i,1)

! ynew(i,:) = (/tp, yi, cj/)  from pofam
  print*; print*, '************', i, '-th P.O. TP: ', tpo
  
  write(*,'(8f12.8)') ynew(i,:) 
 
  ispl  = 1
  ! subroutine plob(y0, t0, tf, ndim, nvar, tdir, ispl, ftag, deriv, gr_cj,  y) 
  call plob(po0, 0.d0,  tpo, 6, 6, tdir, ispl, fpo, gr_lf, gr_cjlf, pof) 
  print*, 'Check the I.C. and F.C. of p.o.'
  print*, po0
  print*, pof  
  print*;  !read*
  
  ! --- Monodramy matrix
  ! subroutine monomat(yi,tp, mmat, deriv, gr_cj)
  print*, 'refined initial state, tp, ynew'
  print*,  tpo, po0
  call monomat(po0, tpo, mmat, gr_lf, gr_cjlf)
  
!  subroutine detmat( a, n, det)
  call detmat( mmat, 6, det)
  write(133, *) det 
  
  
! analyze the stability of the monodramy matrix, only eigenvalues are needed 
  call eigrg(mmat, ndim, 0, wr_mm, wi_mm, vr_mm)
  
  ! save the pommegv.dat 
  do j = 1, ndim, 1
    write(fmmegv,    '(2e20.10)', advance='no')  wr_mm(j),  wi_mm(j)
  end do
  write(fmmegv,*)
  
  !-- compute the trace and sort them   
  call trace( wr_mm, wi_mm, 6, tr_re, tr_im)
  
  if(i == 1) then 
    wr0 = tr_re;    wi0 = tr_im
    call mmeig_sort( wr0, wi0,  tr_re, tr_im, 6, 1,  tr_str, tr_sti, ind_st)
  else 
    call mmeig_sort( wr0, wi0, tr_re, tr_im, 6, 0,  tr_str, tr_sti, ind_st) 
  endif   
  wr0 = tr_str; wi0 = tr_sti 
   
  ! save the ordered index and trace 
  write(find_st, *)   ind_st
  
  do j = 1, 6, 1
    write(ftrace,    '(2f20.10)', advance='no')  tr_re(j),  tr_im(j)
    write(ftrace_st, '(2f20.10)', advance='no')  tr_str(j),  tr_sti(j)
  end do
  write(ftrace, *);  write(ftrace_st, *)  
  
enddo   

 close(fmmegv); close(ftrace);  close(ftrace_st)
 
 
 
stop
end program main_po_ck



















  

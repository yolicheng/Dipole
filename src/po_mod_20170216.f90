module po_mod 

! 2016-08-07 17:17:28 
!   Add the computation of the differential of the Poincare map 

! 20160411 
!  Combine the two approaches: Poincare map + General into this routine 

!    we have 2 possible approaches to deal with the asymmetric p.o. 
!    1. fix one component of the state vector, use this as the poincare section, and call poinc_pomod 
!    2. do the refinement w.r.t X0(6)+TP, altogether 7 unknowns, and try to minimize the difference between the initial condition and the final
!       condtion, call gr_pomod 

  
! For the continuation: Adams predictor + champ(vector field update) 

! user-defined subroutine, make sure to be compiled before the callee module 
! 1.  detmat   	the determinant of a square matrix 



! 20160407
! For the asymmetric approach, please note that because we are refine w.r.t the initial condition
!   So the Jacobi Matrix w.r.t the correction on the initial state is actually Phi - I  

! 20160317
!   Remember! You shouldn't (generally) try to implement a very efficient matrix multiplication yourself and instead consider using a library such as e.g. atlas, which can optimize for cache size, or a vendor library such as MKL or ACML.
!  So, if possible, try to use the intrinsic function or use a library rather than do it by yourself

! 20160316
!  Modified to general case for both symmetric and asymmetric p.o., initialize by init_asympo and init_asympo respectively.
! ---todo: the continue subroutines: adams and   champ

! 20160315
! this is for symmetric p.o. computation, in this case, we have 3 control variables and 2 target varaibles. 
! but for asymmetric orbit,  it is no longer available....

! 20160306 -- No more computer games! modify the module to work perfectly, and pass the name of  deriv and gr_cj as parameters 
! - For the external subroutines the module calls, they should be complied before the module, or the module will be linked to empty external ones, without complain  and error.

! 20160304
! there is some problems here, tha call of gr_rk78 makes all the values from module losing their values... 
!  the save statement doesn't help
!  it seems I have to include all the routines inside this module, it turns out that if there 

! Archive all the subroutines to compute periodic orbit in this module, save the time to copy the related subroutines for every problem
! use f90 free format, instead of fixed f77 format
! 1- Advantage:  don't need to pass some common parameters through all the subroutines, declared as public to share among all the subroutines

! 2- Idea:  use poinc to compute the intersecion with the section y(ind)=sec, use modified Newton method to refine p.o.
! defaultly we deal with 2 equations(target varaibles) with three unknowns(control variables)

! 3- Error and termination control 
!    Newton method is so efficient that, after 5 iterations, if it is still not convergent, we need to check what happens
!    The tolerance(error in target variables) and presion(of the control variables) is to be careful assigned 
!    For Gerard's book, P96, the percision required for the modified Mewthon's method has been taken equal to 1.d-11 
!                         and the bound for local errors in RK78 routine was set to 1.d-13

!    Error in fm = dsqrt( f(1)*f(1) + f(2)*f(2)) ! Gerard uses tol = 1.d-16, quite confusing, should be greater than 1.d-13  

!    tol = 1.d-13; prsc =1.d-11 ! the suggested values
!    but for lf problem, we take  tol = 1.d-11; prsc =1.d-11  

!    because the error control for rk78 is 1.d-13, all the poincare map-related data has no smaller precision
!    so we cannot ask more precision in Newton method

! 4- Comments: do not use the obsolete syntax goto in f90 

! 5- Fortran has 15 significant digits for  double real numbers 

! 6- Symmetry - initial condition - final condtion
!    For the discussion of how to utilize the symmetry to compute the p.o., please refer to Rusell, global search for p.o. P6-symmetries
!    For x-z planar symmetry,       given initial condition(x0,0,z0,0 0,dy0,0) --> T/2 --> (x0,0,z0, 0,-dy0,0)
!    For x-axis symmetry(vt_lyap),  given initial condition(x0,0,z0,0 0,dy0,0) --> T/4 --> (x0,0,0,  0, dyf,dzy) 

! 7- Termination control for the poincare map 
!    treat as lose of convergence if it spends too long for the next interation with the poincare section

! 8- -- discard this item!! Check the infinity norm of matrix phi and g - from Lapac!!!!! use LLS solver instead
!     DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
!     phinormi = dlange('I', 6, 6, phi, 6, phinm)
!     subroutine dgecon(norm,n,a,lda,anorm,rcond,work,iwork,info) -- work(4*n), iwork(n)
!     call dgecon('I',6, phi, 6, phinormi,rcond, work,iwork,info)

!  	Contains Subroutines  		  		function used:
!  1. init_po	initialize public varaibles	   	none
!  2. pofam  	main subroutine				dfcr, champ, adams
!         general routine to integrate the orbit in time interval [t0, t0+tf] 

!  3. dfcr   						poinc, deriv, gr_cj , fctn, deltx
!  4. fctn  						poinc, deriv, gr_cj 
!  5. poinc 						sect, gr_rk78, deriv, gr_cj 
!  6. champ, adams, deltx_lns, sect, deriv,gr_cj 	 	none ---- basic subroutines  

! Note: gr_rk78 is outside this module, that is a general integration routine

!   **********  Subroutine  Declaration ************ 

! ----------------  Module Variable Initialization-------------------
!  subroutine init_poinc(ind0, sec0, hmin0, hmax0, e0, tmax0)
!  subroutine init_errctr(tol0, prsc0)
!  subroutine init_asympo(asymt ) ! for symmetric p.o.
!  subroutine init_sympo(symt )   ! for asymmetric p.o.

 
!  subroutine pofam(yi,npo,imax, dir,ds, fpo, ynew, i, deriv, gr_cj)
!        write(fpoinst,'(10d24.16)') tp, yi, cj, dcj ! PoInSt.dat 
!        ynew(i,:) = (/tp, yi, cj/)
       
!  subroutine dfcr(y0,imax, tf,yf,g, cj,dcj, ispo, deriv, gr_cj)
!  subroutine fctn(x,init,imax, f,g,y,vf,tf, ispc, deriv, gr_cj)
!  subroutine poinc(yi,imax, tf,yf, hminim, ispc, deriv, gr_cj) 
!  subroutine sect(y, neq, g, dg)  
 

!   *****************   Public  variables   ***************************** 
!    Assigned by subroutine init_lf before call any function from this module tine fctn
!  sec, ind	 	the value to specify the Poincare section y(ind) = sec 
!  ntp 			the real period is TP = ntp * tp(to go to the first crossing)
!  tari, (i=1,2)	target varaibles
!  ctri, (i=1,2,3) 	control varaibles
! ************************************************************************************************

implicit none

save
integer, parameter, private :: dp = kind(1.d0) !to avoid compliance with dp defined in other modules
!  ( this module is supposed to be  self-contained, and for general use, define dp for private use)


! ******************** Declare Public *****************
! public by default, the most commonly used varaibles which are also accessible by any subroutine that uses this module
real(kind=dp), parameter ::  mu = 0.12150586d-1   !mu = 0.012150585609624  ! this is just to check halo orbit in other polar orbit exploration

integer, private :: debug, dsctr, anglectr,  isarc, issymt   ! do the assignment from main routine, to aviod the compilation of link before every excution
 
integer  :: ntp,  &  ! ntp is the be multiplied to the tf by poinc subroutine to obtain the full period 
            ntar, nctr  ! number of target variables(equations) and control variables(free variables)

! declare  allocatable arrays to store the index of target and control variables  -- to check
integer, allocatable :: tar(:), ctr(:)
            

! Declare external for the routines or functions to be called the routines from library lapack  

external  deriv,  gr_cj ! i don't know if this helps

! -------------- declaration of subroutines from lapack or nag( external package)
!external dgecon ! subroutine-- not inside this module
real(kind=dp), external :: dlange ! compute the norm of array,  dlange is a function, declare the type
real(kind=dp), external, private :: dnrm2 ! compute the norm of array, from package NAG Fortran Library Routine

! the files to check the iteration, the correction
integer, private ::  fdx 
    
!  private - only shared by the subroutines, cannot be accessable from outside the module 
!  can avoid passing the shared parameters level-by-level
integer, private       ::  ind, imax, idcor,  tfree
real(kind=dp), private ::  csangle_min,  sec, hmin, hmax, e, & ! Bound values for rk78, for the integration of orbit
                     tmax, tol, prsc, & ! tolerance control for the termination of Poincare map and Newton method
                     ds_max, ds_min, &
                     x0_fxd, tol_poinc ! cor_gr 

contains  
  include 'cor_poinc.f90' ! Correction by poincare map approach, not reliale for asymmetric case
  include 'cor_gr.f90'    ! correction by general approach 
  
  include 'curv_angle.f90' ! compute the anlge between 3 consecutive points
  
  ! ***********  prediction *****************
  include 'pre_arclen.f90' 
!  include 'pre_prd.f90'   !TODO: un-debugged
!  include 'pre_cj.f90'  
  
!******************** Module-related Parameters ********************************
  subroutine init_poinc(ind0, sec0, tol_poinc0, imax0,  hmin0, hmax0, e0)
! Assignment of shared variables for poincare map-related subroutine  + rk78 (for integration)
! imax is explicitly assigned here, which is also implicitly assigned in init_asymtpo and init_symtpo
! so if we want to modify a little bit how the routine works, remember to call this one after these two subroutines init_asymtpo and init_symtpo

  integer, intent(in) :: ind0, imax0
  real(kind=dp), intent(in) :: sec0, hmax0, hmin0, e0, tol_poinc0  

! speicially for only call of poinc, ind may be reassigned later by init_asym or init_symt subroutine
  ind  = ind0 
  sec  = sec0
  imax = imax0

  hmin = hmin0
  hmax = hmax0
  e    = e0
  tol_poinc  = tol_poinc0
  
  if (debug == 1) then          
    print*, 'check init_poinc: ind, sec, imax, hmin, hmax, e'
    print*, ind, sec, imax,  hmin, hmax, e
    read(*,*)
  endif
   
!  tmax = tmax0 
  
  end subroutine init_poinc
  
! ******************************************************************************  
   ! contol the working mode: To debug, or to execute
! ******************************************************************************  
  subroutine init_debug( debug0, dsctr0, anglectr0, csangle_min0,  isarc0, ds_min0, ds_max0 )
  integer, intent(in) :: debug0, dsctr0, isarc0, anglectr0 
  real(kind=dp), intent(in) :: csangle_min0, ds_min0, ds_max0
  
  debug = debug0
   
 !  utomatic control of the stepsize for continuation, with dsctr = 0, we obtain 2 families of p.o.s, which i showed to Gerard.
  dsctr =  dsctr0

! control the angle between 2 consecutive vectors.
  anglectr = anglectr0

  csangle_min = csangle_min0 !csangle_min
  ds_min = ds_min0
  ds_max = ds_max0
  isarc = isarc0
  
  end subroutine init_debug
  
  
! ******************************************************************************  
  ! the error control for termination of Newton Method
! ******************************************************************************  
  subroutine init_errctr(tol0, prsc0 )

  real(kind=dp), intent(in) :: tol0, prsc0
  
  tol = tol0 
  prsc = prsc0
    
  end subroutine init_errctr
  

! ****************************************************************************** 
! Idea: Fix the value a component of the initial state as the poincare section x = x0, 
!       and generally an orbit will pass any orbit on it twice, so we take the second intersection with the x = x0 
!       with the velocity in x component to be in the same direction with the initial one. And then use differential corrector method
!       to make sure the final point overlaps with the initial one . 

! A potential bug: what if the initial point with x=x0 is exactly the point that only has one intersection with the x=x0 plane.... 
!                  It rarely happens, but there is possibility.

! Solution: Take the most secure approach, use gr_po instead....

!  The energy h0 and x0 is fixed, look for the first return to x=x0 plane, we only need to ask four final component in state to be as the 
!  initial one, the fifth will be satisfied simultaneously.  And we have 5 initial components to modify(except x0)

! constraint funcion:   ds^2 = ( (yf-y0)^2 + (zf-y0)^2 + (vxf-vx0)^2 + (vyf-vy0)^2 ) to be minimum

! asym   		contol  		target 
! 1   	 (y0,z0, vx0,vy0,vz0) - (2,3, 4,5,6)   (yf, zf, vxf, vyf) 

!  To ensure that enough memory is available to allocate space for your array, make use of the STAT option of the ALLOCATE command:

!   ALLOCATE ( A(N), STAT = AllocateStatus)
!   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!   Here, AllocateStatus is an integer variable. AllocateStatus takes the value 0 if allocation is successful or some other machine dependent value of there is insufficient memory.

!   An array can be released from memory by using the DEALLOCATE command
!   DEALLOCATE (A, STAT = DeAllocateStatus)
! ******************************************************************************  

  subroutine init_asymtpo(asymt)
  integer, intent(in) :: asymt  ! index of symmetry, used to assign the public values accordingly
  integer :: allocatestatus
 
! flag of symmetric p.o. 
  issymt = 0
  ntp = 1    
  imax = 2
  
! number of control and target variables
  ntar = 4 
  nctr = 5 
  
! fix the value of x(ind) =sec , and try to modify the other component to make the first return to x=x0 plane goes to the same point
  ind = asymt   
  
! allocate memory for ctr and tar   
  allocate ( ctr(nctr), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"
    
  allocate ( tar(ntar), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"

 
  if ( asymt == 1 ) then 
    tar = (/2,3,4,5/)   
    ctr = (/2,3,4,5,6/) 
  
  elseif (asymt == 2) then 
    tar = (/1,3,4,5/)   
    ctr = (/1,3,4,5,6/) 
    
  elseif (asymt == 3) then  
    tar = (/1,3,4,5/)   
    ctr = (/1,3,4,5,6/) 
       
  endif 
 
  if (debug == 1) then          
    print*, 'check init_asympo: asymt, ind,  control and target variables'
    print*, asymt, ',', ind, ',', ctr, ',', tar 
    read(*,*)
  endif
  
  end  subroutine init_asymtpo 


! ****************************************************************************** 
! For asymmetric p.o.:  Initialize the index of the control variables and target variables as in the state vector
!   The Poincare section is specified as:  y(ind) = sec
!   Note: to use these symmetries, the invariant transformation in time should be t-->-t
! question: symt = ind? NO! for vertical lyapunov orbit, we take y=0 as Poincare section, but use symmetry 2, so cannot join as one parameter

! It is not ok to just pass the parameter as input, without assignment by h=h0
! ****************************************************************************** 
  subroutine init_symtpo(symt )
  
  integer, intent(in) :: symt  ! index of symmetry, used to assign the public values accordingly

  integer ::  allocatestatus
  
! flag of symmetric p.o. 
  issymt = 1
  imax   = 1 
   
! --- allocate memory for ctr and tar 
  nctr = 3
  ntar = 2  
  
  allocate ( tar(ntar), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"

  allocate ( ctr(nctr), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"
 
 
! ***************************  planar symmetry  *******************************
! 1: x=0(y-z) plane, 2: y=0(x-z) plane, 3: z=0(x-y) plane
 
! ----- Take symt = 2 for example   
!  initial state:  	(x,0,z, 0,vy,0)
!  stop condition: 	Poincare map to get the first crossing with y=0 plane
!  control variables:	x, z, vy  	-- 1,3 (1,2,3 except ind ); 3+symt
!  target variables: 	vx = 0, vz = 0 	-- 3  +   1,3(1,2,3 except symt )

! tha State Transit Matrix:   6-6 
! if symt == 2, ctr1 = 1, ctr2 = 3, ctr3=5;  tar1 = 4, tar2 = 6
! g (the variational matrix) : 
! | phi(tar1, ctr1), phi(tar1, ctr2), phi(tar1, ctr3) |
! | phi(tar2, ctr1), phi(tar2, ctr2), phi(tar2, ctr3) |
! --------------------------------------------------------

! symmetry:   Initial condition	        index of control var   	Past time      target var
!    1:	    	x=0, yz-plane,	   	 2,3,3+1 (y, z, vx)   	  T/2	 3+2, 3+3 (vy, vz) (yz plane)     
!    2: 	y=0, xz-plane,   	 1,3,3+2 (x, z, vy)   	  T/2	 3+1, 3+3 (vx, vz) (xz plane)     
!    3: 	z=0, xy-plane,   	 1,2,3+3 (x, y, vz)   	  T/2	 3+1, 3+2 (vx, vy) (xy plane)  

! ************************** axis symmetry ************************** 
!    4:x-axis	y=0, xz-plane  		 1,3,3+2  (x, z, vy)  	  T/4    3, 4 (z, vx) (perpendicular to x-axis) -- vt_lyap case 
! to be added the other two

! in fact, now we only need to deal with the case 2,  test vt_lyap later----to be modified
! ! symmetry:  index of control var   of target var
! 5: y=0, xz-plane,  1,3,3+2 (x, z, vy)     3+1, 3+3 (vx, vz)      
! 6: z=0, xy-plane,  1,2,3+3 (x, y, vz)     3+1, 3+2 (vx, vy)   
   
  if (symt == 1) then  ! 
    ind = 1;   ntp = 2 ! the Poincare sect y(ind) = sec 
    ctr = (/2, 3, 4/)  ! control variables
    tar = (/5, 6/) ! target variables
      
  elseif (symt == 2)  then 
    ind = 2;   ntp = 2 
      
    ctr  = (/1, 3, 5/)  
    tar  = (/4, 6/)  
      
  elseif (symt == 3)  then 
    ind = 3;   ntp = 2 
      
    ctr = (/1, 2, 6/)  
    tar = (/4, 5/) 

! axis symmetry      
  elseif (symt == 4)  then  ! x-axis + x-z planar symmetry
    ind =  2;  ntp = 4  
    ctr  = (/1, 3, 5/)   
    tar  = (/3, 4/)  
      
! 20160318 todo -- the other cases are to be added, not needed at this moment         
  endif
  
  if (debug == 1) then          
    print*, 'check init_po: symt, ind, y0, control and target var'
    print*, symt, ',', ind,sec, ',', ctr, ',', tar 
    read(*,*)
  endif   
  
  end subroutine init_symtpo


!! ****************************************************************************** 
!! The initializaiton for general approach to refine p.o.: 
!!  Initialize the index of the control variables and target variables as in the state vector
!! ****************************************************************************** 
!  subroutine init_grpo

!  integer ::  allocatestatus
!  
!  issymt = 2
!  ntp = 1
!! --- allocate memory for ctr and tar 
!  nctr = 7
!  ntar = 6 
!  
!  allocate ( tar(ntar), stat = allocatestatus)
!  if (allocatestatus /= 0) stop "*** not enough memory ***"

!  allocate ( ctr(nctr-1), stat = allocatestatus)
!  if (allocatestatus /= 0) stop "*** not enough memory ***"
!  
!  tar = (/1,2,3,4,5,6/) 
!  
!  ! the control variables in fact have 7 components, X0(6) + TP , here we only need the index of the ones that are the component of the state  
!  ctr = (/1,2,3,4,5,6/) 
!  
!  end subroutine init_grpo 
!  
  
! ****************************************************************************** 
!  This  subroutine is to open the file to save the intermediate result, including
!  the files to check the correction,  matrix norm  for every iteration of Newton method 
!  Turns out the matrix is in good condition, because gnewnorm is small, that means the error will lead to a small correction(multiplied by the matrix),
!  but pay attention to the numerical error, which is 1.d-13 previously, but I ask the tol<1.d-14, stupid!!!!
 
! ****************************************************************************** 
  subroutine init_writedx  

!! Bound values for RK78, shared by all subroutines that calls rk78 (poinc, fctn)
  if(debug == 1) then 
    print*, 'hmin,hmax, e, tmax,tol,prsc', hmin,hmax, e,tmax,tol,prsc ;  read(*,*)
  endif 
  
  fdx = 32 
  
  if (issymt == 0) then 
    open(fdx,file='./dat/adx.dat',access ='append',status='replace')  
  else if (issymt == 1) then 
    open(fdx,file='./dat/dx.dat',access ='append',status='replace') 
  elseif (issymt == 2) then 
    open(fdx,file='./dat/grdx.dat',access ='append',status='replace') 
  endif 
  
       
  write(fdx,*) '# dxm	 errfm	 dx	index of correction(', ctr, ')'
    
  end subroutine init_writedx
  
  
!******************************* Subroutines ***********************************
  
  
!************************  Numerical Continuation ****************************** 
! Main subroutine to compute np new p.o.s starting from one initial guess yi(6), 
! and do the numerical continuation along the arc-length parameter.

!  To make sure the initial state is periodic, refine it without check for the first one.
!  the initial condition for all p.o.s are stored in file fpo

! 	Input Varialbles
!  yi(6)	initial condition (ie, x,y=0,z, vx=0, vy, vz=0)
!  tp0		the initial guess for tpo( use as tmax for poincare map approach, or just tp0 for general case)
!  npo		number of new p.o.  if np=1, just do differential correction, without continuation

!  dir		dirction along the vector field, +1:increase,-1:decrease
!  ds		step size for the continuation (adaptive control), lower and upper bound [ds_min, ds_max] 
!  fpo		file tag for the initial states of the p.o.s

! 	Output Variables
!  ynew(npo,8)	the initial condition for new p.o.s, npo*8 ( TP, XPO, CJ )
!    		save the initial p.o in the first row
!  i 		the real number of available p.o.s is (i-1)

!  	Private Module-base Varaibles
! tol		tolerance for vx,vz to be zero (ie,1.d-13 or 1.d-16 in gerard's routine)
! prsc 		presion to terminate the iteration(ie, 1.d-13)
! nctr, ntar, tar, ctr, ntp   
  

!  Finally revised : 20160411 

! todo: 1. chaml is only used for check,  delete this one later after a thorough debug of the routine 
!       2. 

! Routine used: dfcr, fctn, champ, adams
! --------------------------------------------------------------------------
 
  subroutine pofam(yi, tp0, npo, dir,ds, fpoinst, ynew, ipo, deriv, gr_cj)
  
  implicit none
  integer, parameter :: n = 42  !  both the state and the variational equations are required, so n=42


  integer, intent(in)      :: npo,  fpoinst 
  integer, intent(out)     :: ipo
  integer, intent(inout)   :: dir  ! may be reversed
  
  real(kind=dp), intent(in)     :: yi(6)  ! initial guess is updated
  real(kind=dp), intent(inout)  :: ds, tp0  ! adjust based on the number of iteration in Newton Method
  real(kind=dp), intent(out)    :: ynew(npo,8) 


  external :: deriv, gr_cj  ! vector field  & jacobi constant
  
! Local Variables  -- be careful with the dimension
  integer :: ispo, niter, ncol, nrow, & 
             nds, nds_new, incrds, con_stop, dir0  ! the counter of the available vectors with the previous step size  
             
  real(kind=dp) ::  y0(6), yf(n), g(ntar, nctr), yctr(nctr), &
                cj, cj2, dcj, &
                curv(5, nctr), curv3(3, nctr), cham(4, nctr), cham5(5,nctr), &
                curv_prd(5, nctr+1), curv_cj(5, nctr+1), csangle, vf(nctr), vf6(6) 
! f(ntar),
 
 !  adaptive control of the stepsize for continuation, with dsctr = 0, we obtain 2 families of p.o.s, which i showed to Gerard.
  if(imax == 2  )  ntp = 1 ! is it really the case? ?? for the doubly-symmetric case ???

  y0   = yi   ! copy the initial state, and keep yi unchanged
  ynew = 0.d0 ! assign to 0 for the initial state of new p.o.s 
  
  print*, '************************ start pofam*******************'; print*

  if (debug == 1) then 
    call gr_cj(yi, cj) !ck
    print*, 'cj, yi', cj, yi ; !read(*,*) !ck
  endif
  

  ipo = 1 
  nds = 0 ! the available vectors with the current stepsize in champ, use this we can generate cham for nds <= 4
  nds_new = 0 
  
  csangle   = 0.d0 ! the angle between two consecutive vectors 
  
!  ------------------------- numercial continuation --------------------------- 
  do 

  ! reinitialize the flags for control of stepsize to be zeros
    incrds = 0   ! according to the iterations
    con_stop = 0 ! too small stepsize required
    
  ! if we want to change the value of the iteration counter ipo, instead of the automatic way do i = 1,n 
  ! we use if syntax to detect the termination of the do loop
    if ( ipo > npo ) exit 
    
 
 ! y0 is only effectively updated if Newton Method succeeds(ispo = 1) 
 ! if it fails, call dsdecr to update a new yctr for the next p.o., and if the step size required is too small, terminate the continuation 
    
 ! --- write the correction onto screen and into file dx.dat to check-------------- 
    print*  
    print* ,ipo,'-th p.o.', '    dir=', dir, '        ds=', ds, '    nds_new=', nds_new, '       <v1,v2>=', csangle;  print*
    
    write(fdx,*); write(fdx,*) '****************************************'
    write(fdx,*)   ipo,'-th p.o.', '       dir=', dir, '        ds=', ds, '    nds_new=', nds_new, '       <v1,v2>=',csangle
    
    
    ! initial condition before the refinement
    write(fdx, *)  'I.C.         ', y0
    write(*, *)    'I.C.         ', y0  ! print to screen   
 ! --------------- ------------------------------- --------------- ---------------   
 
    
!---- refine the p.o. using differential correction -------------- 
! --  with initial and updated refined condition as y0, together with period tp0 -----------

!   subroutine dfcr(y0, tp0, yf, g, cj,dcj, ispo, niter, deriv, gr_cj) 
    call dfcr(y0, tp0, yf, g, cj,dcj, ispo, niter, deriv, gr_cj)  
    

    if( ispo == 1) then 
! check the angle between the last 3 consecutive points, if not in the admissable bound, decrease the arc step without update curv
      print*; print*, ' ************** Refined initial state: tp, y0 ******************'
      print*, tp0, y0; print*
      
      ! -------------------- check curve ------------------
      ! if nds > 2, at least 2 available rows in cham and curv , and 3 computed p.o.s
      ! compute the angle between the 3 consecutive points, and decide here if we need to decrease ds or not
     
      ncol = min0(nctr,6)  ! available rows in curv 
      nrow = min0(nds,5)   ! available column in curv
      
      if (nds .ge. 2 ) then  
        curv3(1:2, :) = curv( nrow-1 : nrow,  : )
        curv3(3, 1:ncol) =  y0(ctr)
        if (issymt == 2) curv(3, 7) = tp0
       
        call curv_angle(curv3, nctr, csangle)
        
        ! need to decrease the arc step to satisfy the angle constraint  
        if ( dabs(csangle) < csangle_min .and. anglectr == 1)  then 
          if ( dsctr == 0 ) then
            print*, '|cos<v1,v2>| =', csangle, '<', csangle_min,'Terminate the continuation!'
            read*; return 
          else ! dsctr == 1  
            incrds = -1  
          endif   
        endif 
          
      endif   
      !--------------------------------------------------------
      
      ! -------------- update curv ---------------------       
      if (incrds /= -1) then  ! update curv only if we do not need to decrease ds

        if(nds .le. 4) then 
          curv(nds+1, 1:ncol )  = y0(ctr) 
          if (issymt == 2) curv(nds+1,7) = tp0
    
        else 
          curv(1:4,:) = curv(2:5,:)
          
          curv(5, 1:ncol ) = y0(ctr)
          if (issymt == 2) curv(5, 7) = tp0
        endif 
        
        
        if (isarc == 2) then  ! along the period
          curv_prd( min0(nds+1, 5), nctr+1) = tp0 

        elseif(isarc == 3) then ! along the energy...
          curv_cj( min0(nds+1, 5), nctr+1) = cj
        endif 
        
         
        !---------- evaluate if we need to double the arcstep ---------- 
        if (niter < 3 .and. dsctr == 1 .and. nds_new .ge. 4 .and. 2.d0*ds .le. ds_max )  incrds = 1 

        if(debug == 1) then 
          print*, 'ck! incrds = 1 if and only if nds_new >=4', nds_new, incrds !ck
          read*  
        endif 
          
      endif  
      ! ----------------------------------------------------     
    
    else 
    
! ************************   ispo = 0  ****************************************
!  the only case that we still have a chance to start again with a smaller stepsize with ispo = 0 is ds>ds_min 
!  because of the round error, ds is not equal to ds_min exactly, using the system precision as the criteria of equivalence.
        
      if (dsctr == 1 .and. dabs(ds - ds_min) .ge. 1.d-14 ) then 
        incrds = -1 
      else 
        print*, 'Terminate the continuation process!' 
        read*;   return 
        
      endif  
      
    endif ! ispo    
   

! there are only two possible situations once we arrive here :
!  --1.  ispo = 1 
!        if niter < 3,  ds = ds * 2.   (condition: ispo == 1 )
  
!  --2.  ispo = 0 .and. dsctr = 1 .and. niter > 5 .and. ds < ds_min, decrease the step size and start again the correction later ....
!        if niter > 5,  ds = ds / 2; 

! bound interval [ds_min, ds_max],  for ipo > 1
! for each value of ds, do at least 4 orbits, in order to provide cham(1:4, ctr) with the same step size for adams
! after each change of value of ds, we do more 4 orbits to avoid the previous points to be too far from the present one. 

    if(debug == 1)  then  !ck
      print*, 'ds, niter, incrds, nds_new', ds, niter, incrds, nds_new; read*
    endif 


    if ( npo > 1 ) then   
 ! ------------------------------  update yctr  --------------------------------   
      if ( incrds /= -1 ) then 
!  use the initial state of current good p.o. as the guess for next p.o., and then add a small variation 
!  when no decrease is needed, yctr is assigned as the current y0 and tp0(if necessary)
!  with incrds = -1, yctr is updated by routine dsdecr 
    
        yctr( 1 : ncol ) = y0(ctr)  
        if (issymt == 2) yctr(7) = tp0
      endif 

!----------------  use arc-length parameter  -----------------------------------
      if ( isarc == 1 .or. ipo < 3) then 
        dir0 = dir 
!  subroutine pre_arclen(ipo, ds, nds, dir, cham, cham5, curv, yctr, c, incrds,  con_stop) 
        call  pre_arclen(ipo, ds, nds, dir, g, cham, cham5, curv, yctr, csangle, incrds, con_stop)  
        if ( isarc /= 1 ) dir = dir0 
      else
        call deriv(0.d0,yf(1:6),6,vf6) ! compute the vector field
        print*, 'yf:',  yf(1:6)
        print*, 'vf6:', vf6
        vf = vf6(ctr)
        
        print*, 'vf=', vf
        read*

! TODO: at this moment, haven't make the idea clear...        
!        if(isarc == 2) then  ! along the period   --TODO 
!          curv_prd(: , 1:nctr) = curv
!          call pre_prd( ipo, curv_prd, vf, g,  nds, ds, dir, incrds, csangle, yctr, con_stop)
!        
!        elseif(isarc == 3) then! along the energy -- TODO
!          curv_cj(: , 1:nctr) = curv
!          call pre_prd( ipo, curv_cj, vf, g,  nds, ds, dir, incrds, csangle, yctr, con_stop)
!        endif 
        
      endif 
      
     ! ---  stop condtion check ----------
      if ( con_stop == 1)  return 
    endif    


! update the counter for each vaule of arc step ds (reset after each change)   
    if(incrds == 0 ) then 
      nds_new = nds_new + 1
    else 
      nds_new = 1 ! reset to 1
    endif   


! ------------ save the refined 'good' initial state of ipo-th p.o. --------------------  

    if( ispo == 1 .and. incrds /= -1 ) then  
      
      tp0 = ntp * tp0  ! the full period
      ynew(ipo,:) = (/tp0, y0, cj/)
    
  ! write to file labeled by fpo passed from the main routine, PoInSt.dat  
      write(fpoinst,'(9e22.14, 1f12.8, i5)')  tp0, y0, cj, dcj, ds, dir  ! fpoinst - PoInSt.dat 
      
 ! ----------- check the energy ---------------------------   
      if(debug == 1)  then 
        call gr_cj(y0, cj2)
      
! check to see if there is any point to keep cj, dcj as the output! ckd--cj = cj2 anyway
        print*, 'check energy-- cj, dcj, cj2:', cj, dcj, cj2  
        write(*,'(9e22.14, 1f12.8, i5)')  tp0, y0, cj, dcj, ds, dir  !  PoInSt.dat   !  !ck print to screen
        read(*,*)
      endif 
!--------------------------------------------
  
      if(npo .eq. 1)  return  !  refinement for only one p.o.
   
    endif   
 
    
!  we succeed in orbaining an initial guess for the new p.o, update the value of counter and  state + period(if necessary)     
    ipo = ipo + 1 
    y0(ctr) =  yctr(1: ncol)
    if (issymt == 2)  tp0 = yctr(7)
    
    if (debug == 1) then   
      print*, 'Numerical continuation, new state!' !ck
      print*,  y0(1:6);  print*;  read* 
    endif 
  enddo


  close(fpoinst) ! close poinst.dat

  return
  end subroutine pofam


!**************************   Differential Correction  ************************* 
!  Using modified Newton method to refine the p.o,  provided with the initial conditon.

!  We have nctr unkowns and ntar equations(nctr>ntar, inderterminant equations),  using LLS(Linear least square problem) solve 
!  and take least-square solution.

!  The point is get the derivative of the target variables with respect to the control variables


! The cases when this refinement fails to reach convergence
!   1. too many iterations, niter>5.  based on the quadratic convergence of Newton method, 5 iteration is enough!
!   2. For the poincare map approach, too long time to arrive to the next intersection acrossing poincare section,  ispc = 0


! for Newton method, alway check two things:
!   1. the error in target varaibles    --- stop condition, satisfied means finish the correction
!   2. the modulus of the correction    --- if too small, no point to continue ...

!  For the general approach, instead of Poincare map approach, we use the whole variational matrix to do the refinement w.r.t X0 + TP
!  so we 6 equations and 7 unkonwns(period included) 

! the 7 unkowns are (xi, yi, zi, vxi, vyi, vzi, T) 
! and the 6 constraint equations are satisfied to return to the same initial point 
!  (xfi-xi=0, yf-yi=0, zf-zi=0, vxf-vxi=0, vyf-vyi=0, vzf-vzi=0)
 
 
!******************************************************************************* 
! 	Input 
!  y0 	initial state
! tp0	inital guess of the period(depend on which approach to apply)

! 	Output
!  y0(tp0)  	the updated state and period
!  tp0		elapsed time to reach to the next mirror configuration,  which is 1 period for asymmetric p.o. 
!  yf		final state at the epoch of the mirror configuration
!  g 		Jacobi Matrix of the F(the periodicity equation ), will be used for the continuation
!  cj		Jacobi constant
! dcj		difference of the Jacobi constant between the initial and final point along the p.o.
! ispo		flag to indicate if the correction succeeds(1) or fails(0)
! niter		number of iteration for the Newton Method, used as criteria to adjust the stepsize for continuation

! deriv		the subroutine to compute vector field of the equation of motion
! gr_cj		the subroutine to compute the Jacobi Constant(or energy)
! 
  subroutine dfcr(y0, tp0, yf, g, cj,dcj, ispo, niter, deriv, gr_cj)
  
  implicit none  
  integer, parameter :: dp = kind(1.d0), n = 42 ! both the state and the variational equations are desired, so n=42
  
  
  integer, intent(out) :: ispo, niter
  real(kind=dp), intent(out) :: yf(n), g(ntar, nctr), cj, dcj 
  real(kind=dp), intent(inout)  :: y0(6), tp0
  
  external ::  deriv, gr_cj
  
  
! Local Varaibles
  integer        :: i, debug, ispc
  real(kind=dp)  :: f(ntar), pf(n), dx(nctr), dxm, fm, cj2 
                    ! yit(6) !  yi(n), y0copy(6) ! , dx2(nctr)  
  logical :: ok
  
  ! Saved the intermediate orbits. if file poit.dat is opened, close it, use inquire to obtain the status, 
  inquire( unit = 12, opened = ok )
  if(ok)  close(12)
  open(12, file='./dat/poit.dat',access ='append',status='replace')
 
  if (debug == 1)   print*, 'start dfcr, y0!', y0   

  niter = -1 
  ispo  = 1 ! defaultly set Newton Method as successful
  dx = 0.d0;   dxm = 0.d0
  
! initial condition before the refinement
  if (issymt == 2) then 
    write(fdx, *)  'I.C.         ', y0, tp0
    write(*, *)    'I.C.         ', y0, tp0  ! print to screen   
    
  else
    tmax = 1.5d0 * tp0 ! the maximum possible time for poincare approach to reach to the imax-th poincare section
  endif   
  
 ! --------------- --------------- --------------- --------------- ---------------   
    
! iteration to do the refinement, iter is the counter
  do  
    
    niter = niter+1 ! counter starts from 0
    
    if(debug == 1) then 
      print*, niter, '-th iteration'   
      read(*,*)
    endif
 
      
! Check the error in target function before the correction
! if niter > 6, too many iterations, bad initital guess, treat as a failure of the refinement, no matter ispo = 1 or not

! From Carles' class: For Newton method, which is quadratically convergent, if the iteration is great than 5 or 6, then there is something wrong, provided that the initial guess is a good one, which is  close to the solution.

    if (niter > 5) then 
      print*, 'Newton method iterates > 6, dxm, fm',  dxm,  fm
      ispo = 0
      read*
      return
    endif   
 
    call gr_cj(y0, cj)

    if (issymt == 2) then 
! General approach: Integrate to the next mirror configuration to compute G and F

!   subroutine gr_fctn(y0, tp0, f, g, yf, deriv, gr_cj )
      call gr_fctn(y0,  tp0,  f, g, yf, deriv, gr_cj)


    else  !! Poincare map approach   

!  subroutine fctn(x,init, f,g, tf,y,vf, ispc, deriv, gr_cj)
      call fctn(y0, 0, tmax,  f, g, tp0, yf, pf, ispc, deriv, gr_cj)  
      
    endif 
    
    close(12)  ! close the file poit.dat to be able to plot for check
    
!   check the variation in energy for the poincare map here, instead of inside the fctn subroutine    
    call gr_cj(yf(1:6), cj2)
    dcj = cj2 - cj  

    if(debug == 1) then
      print*, 'dcj, cj0, cjf', dcj, cj, cj2   !ckd ! very small difference
      print* 
      print*, 'Finish fctn, g='
      do i = 1,  ntar
        print*,  g(i,:)
      enddo  
      
      print*;       read(*,*)
    endif 

    fm = dnrm2(ntar, f, 1)
    
! save the correction in file dx.dat--- need to figure out what is the best option..., for all orbit I think...
! this fm is corresponding to the dx and dxm from the last iteration, start from the 1st iteration
   
    write(fdx,*) 'F.C.         ', yf(1:6)!finial condition of the poincare map point
    ! difference between the final and the initial condition, check only after a full revolution
    if(ntp==1) write(fdx,*) 'F.C.-I.C.    ', yf(1:6)-y0(1:6)
    
    !print to screen
    write(*,*) 'F.C.         ', yf(1:6)
    if(ntp==1)  write(*,*) 'F.C.-I.C.    ', yf(1:6)-y0(1:6) 
    
    if(debug == 1) read*
   
 ! precision to terminate the correction, 1.d-16 from Gerard, which is of non-sense, because it beyonds the precision of double real
    if ( fm .le. tol)  then 
      print*, 'fm < tol', fm,  tol
      return   ! try 1.d-11  >should be no less than the numerical precision we are using: e=1.d-13 (error in rk78)
    endif 
    
! solve the LLS problem to obtain the correction     
!   subroutine deltx_lns(f, g,  dx)

    call deltx_lns(f, g,  dx)
    
    
    if (debug == 1) then 
      print*, 'By Lapack LLS solver: dx'
      print*, dx 
      print*;    read*
    endif
    
    ! the  precision is the modulus of the correction, or dxm = dnrm2(nctr, dx, 1)
    dxm = dnrm2(nctr, dx, 1)
 
    ! if the correction is below the precision, regard as success of the refinement
    ! if it is too large, treat as a failure...
    
    if ( dxm .lt. prsc )  then 
      print*, '|dx|=', dxm, '<', prsc, 'stop iteration!'
      if(debug == 1)   read(*,*) ! in debug mode, pause to check the result 
      return
    elseif (dxm .gt. 2.d0) then 
      print*, 'The correction needed is too big! dxm = ', dxm
      ispo = 0
      return
    endif   

   ! update the initial state, if the precision is not satisfied 
    if(issymt == 2 ) then 
   ! general case: 7 control varaibles
      y0(ctr) = y0(ctr) + dx(1:nctr-1)
      tp0 = tp0 + dx(nctr)
      
    else 
      y0(ctr) = y0(ctr) + dx 
    endif 
    
    
    write(fdx,*) ! add a blank line to seperate iterations
    write(fdx,*) niter+1, '-th iteration'
    write(fdx,*) '|DX0|,|f|,DX0', dxm, fm, dx 
    
    if(issymt == 2 ) then 
      write(fdx,*) 'U.I.C.       ', y0, tp0 ! the updated initial conditon
    else 
      write(fdx,*) 'U.I.C.       ', y0 ! the updated initial conditon   
    endif 
    
    ! print to screen 
    write(*,*)
    write(*,*) niter+1, '-th iteration'
    write(*, *)  '|DX0|,|f|,DX0', dxm, fm, dx
    if(issymt == 2 ) then 
      write(fdx,*) 'U.I.C.       ', y0, tp0 ! the updated initial conditon
    else 
      write(fdx,*) 'U.I.C.       ', y0 ! the updated initial conditon   
    endif 
    
    if(debug == 1)  read*
    
  enddo 
   
  return
  end subroutine dfcr
  
!******************************* deltx_lns  ************************************
! Use Lapack driver dgelsy to solve the least norm problem  Ax = b,
! which is better than computing the inverse of the matrix followed by the matrix multiplication 
!   here it is G * X = -F   (A=G : 6 X 7, B = -F)

! Finally we choose the solver dgelsd, which is more reliable, meanwhile the time consumed is more or less the same

!- general solution for the Modified Newton method 
!  We have X(k) = X(k-1) + deltX(k-1)
!  where we require :   G ( deltX(k-1) ) = - F( X(k-1) )

!  Minimize the norm
!     (deltX)T * Q * deltX
!      Q is a given diagonal positive definite weight matrix,  here,just set it to identity

!  So, the solution of this problem of minima is given by
!     (deltXk-1)= -Q-1 (G)T * (G Q-1 (G)T)-1 * F( X(k-1) ) 
!     (deltXk-1)= -(G)T* (G*(G)T)-1 * F( X(k-1) ) 

! We have compared dgelsy and dgelsd for sloving this problem, the diffrenence in results is of order 1.d-15. 
! Here we choose dgelsd, which is faster according to the Lapack Official Manuel, while dgelsd is more efficient in the rank deficient LSP

! TODO: we have already had a general routine deltx, consider using that one... 

!	Input 
!  ntar		number of rows in G 
!  nctr		number of columns in G
!  g		Jacobi Matrix of F w.r.t (X0 + T) 
!  f		Target function
!   
!	Output 
!  dx 		the least-norm correction by sloving G * DX = -F 

! module po_mod based parameter 
! integer:: ntar, nctr 

!  Finaly revised by Yu 20160411
!*****************************************************************************

  subroutine deltx_lns(f, g,  dx)
  
  implicit none 
  integer, parameter :: dp = kind(1.d0)
  
  real(kind=dp), intent(in)  :: g(ntar, nctr), f(ntar)
  real(kind=dp), intent(out) :: dx(nctr)

    
! Local Varaibles
! for linear equation solver
  integer, parameter :: lwork = 1600 !  3*nctr + 64*(ntar+1), use 1600 to allocate enough memory for work 
  integer :: jpvt(nctr), rank, info, nrhs 
  
  real(kind=dp)  ::  a(ntar, nctr), b(nctr), rcond, work(lwork) !dgelsy
  
  !dgelsd
  integer :: iwork(1600) 
  real(kind=dp) :: s(ntar) !, a2(ntar, nctr), b2(nctr)
  
  external  ::  dgelsd  ! lapack lls solver
  
  nrhs = 1
 
 ! copy g and f, because after the call of dgelsy, a and b will be modified
  a = g 
  b(1:ntar) =  -f
 
 ! check this least-norm solver with the symmetric approach, 2 equations and 3 unknowns
  
! Initialize JPVT to be zero so that all columns are free
  jpvt = 0
 
! Choose RCOND to reflect the relative accuracy of the input data, for double precision, 
! use a small one, a good option is one smaller than the cpu precision. 
  RCOND = 1.d-20
 
  CALL DGELSD(ntar, nctr,nrhs, a, ntar, b, nctr, S, RCOND, RANK, WORK, LWORK, IWORK, INFO)
  
  
! Print solution to screen to check 
  if( info == 0 ) then 
 
!    WRITE (*,*) 'DGELSD-- Least norm solution' ! this one is more aliable, so better to keep this one maybe 
!    WRITE (*, '(1X,7e24.14)')  (B(I),I=1, nctr)

    ! assign the solution to dx  
    dx = b
    
  else 
    WRITE (*,*) 'SSL Solver Failed to converge!'
    read*
    return
  endif
  
  end subroutine deltx_lns
  
   
end module po_mod

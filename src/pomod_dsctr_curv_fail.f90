module pomod 
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
!  6. champ, adams, deltx, sect, deriv,gr_cj 	 	none ---- basic subroutines  

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

integer, private :: debug, dsctr, issymt   ! do the assignment from main routine, to aviod the compilation of link before every excution
 
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
integer, private       ::  ind, imax
real(kind=dp), private ::  sec, hmin, hmax, e, & ! Bound values for rk78, for the integration of orbit
			   tmax, tol, prsc ! tolerance control for the termination of Poincare map and Newton method


contains 

  include 'poinc_pomod.f90'
  include 'gr_pomod.f90'
  
!******************** Module-related Parameters ********************************
  subroutine init_poinc(ind0, sec0, imax0,  hmin0, hmax0, e0)
! Assignment of shared variables for poincare map-related subroutine  + rk78 (for integration)
! imax is explicitly assigned here, which is also implicitly assigned in init_asymtpo and init_symtpo
! so if we want to modify a little bit how the routine works, remember to call this one after these two subroutines init_asymtpo and init_symtpo

  integer, intent(in) :: ind0, imax0
  real(kind=dp), intent(in) :: sec0, hmax0, hmin0, e0  

! speicially for only call of poinc, ind may be reassigned later by init_asym or init_symt subroutine
  ind  = ind0 
  sec  = sec0
  imax = imax0

  hmin = hmin0
  hmax = hmax0
  e    = e0
  
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
  subroutine init_debug( debug0, dsctr0 )
  integer, intent(in) :: debug0, dsctr0
  
  debug = debug0
   
 !  adaptive control of the stepsize for continuation, with dsctr = 0, we obtain 2 families of p.o.s, which i showed to Gerard.
  dsctr =  dsctr0
  
  end subroutine init_debug
  
  
! ******************************************************************************  
  ! the error control for termination of Newton Method
! ******************************************************************************  
  subroutine init_errctr(tol0, prsc0)

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


! ****************************************************************************** 
! The initializaiton for general approach to refine p.o.: 
!  Initialize the index of the control variables and target variables as in the state vector
! ****************************************************************************** 
  subroutine init_grpo

  integer ::  allocatestatus
  
  issymt = 2
  ntp = 1
! --- allocate memory for ctr and tar 
  nctr = 7
  ntar = 6 
  
  allocate ( tar(ntar), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"

  allocate ( ctr(nctr-1), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"
  
  tar = (/1,2,3,4,5,6/) 
  
  ! the control variables in fact have 7 components, X0(6) + TP , here we only need the index of the ones that are the component of the state  
  ctr = (/1,2,3,4,5,6/) 
  
  end subroutine init_grpo 
  
  
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
    open(fdx,file='./dat/dx.dat',access ='append',status='replace')  
  else if (issymt == 1) then 
    open(fdx,file='./dat/adx.dat',access ='append',status='replace') 
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
!   yi(6)	initial condition (ie, x,y=0,z, vx=0, vy, vz=0)
!   npo		number of new p.o.  if np=1, just do differential correction, without continuation

!   dir		dirction along the vector field, +1:increase,-1:decrease
!   ds		step size for the continuation (adaptive control), lower and upper bound [ds_min, ds_max] 
!   fpo		file tag for the initial states of the p.o.s

! 	Output Variables
!   ynew(npo,8)	the initial condition for new p.o.s, npo*8 ( TP, XPO, CJ )
!     		save the initial p.o in the first row
!   i 		the real number of available p.o.s is (i-1)

!   	Private Module-base Varaibles
!  tol		tolerance for vx,vz to be zero (ie,1.d-13 or 1.d-16 in gerard's routine)
!  prsc 	presion to terminate the iteration(ie, 1.d-13)
!  nctr, ntar, tar, ctr, ntp   
  
!  Routine used: dfcr, fctn, champ, adams
!  Finally revised : 20160411 
! --------------------------------------------------------------------------
 
  subroutine pofam(yi, tp0, npo, dir,ds, fpo, ynew, ipo, deriv, gr_cj)
  
  implicit none
  integer, parameter :: n = 42  !  both the state and the variational equations are required, so n=42


  integer, intent(in) 	   :: npo,  fpo 
  integer, intent(out) 	   :: ipo
  integer, intent(inout)   :: dir  ! may be reversed
  
  real(kind=dp), intent(in)     :: yi(6)  ! initial guess is updated
  real(kind=dp), intent(inout)  :: ds, tp0  ! adjust based on the number of iteration in Newton Method
  real(kind=dp), intent(out)    :: ynew(npo,8) 


  external :: deriv, gr_cj  ! vector field  & jacobi constant
  
! Local Variables  -- be careful with the dimension
  integer :: fout, ispo, i, j, niter, & 
             nds, incrds  ! the counter of the available vectors with the previous step size  
             
  real(kind=dp) ::  y0(6),  y0_org(6), yf(n), f(ntar), g(ntar, nctr), cham(4, nctr), tp, yctr(nctr), &
                    ds_max, ds_min, & ! bound for stepsize of the continuation
   		    cj, cj2, dcj, &
   		    chaml(npo, nctr), curv(2, nctr), cham5(npo,nctr) ! too keep all the vectro field along the characteristic curve
 
 
 !  adaptive control of the stepsize for continuation, with dsctr = 0, we obtain 2 families of p.o.s, which i showed to Gerard.
  if(imax == 2  )  ntp = 1 ! is it really the case? ?? for the doubly-symmetric case ???

  y0   = yi  ! copy the initial state, and keep yi unchanged
  ynew = 0.d0 ! assign to 0 for the initial state of new p.o.s 
  
  ! decide if ds should be increased by 2 (incrds = 1), decreased by half(incrds=-1), or keep the current value
  incrds = 0  
  
  ! bound vaules of stepsize ds
  ds_max = 1.d-1
  ds_min = 5.d-4
    
  if (debug == 1) then 
    print*, 'hmin,hmax,e, tol, prsc ', hmin,hmax,e, tol, prsc 
    read(*,*)
  endif
  
  if (debug == 1) then 
    print*, 'start pofam!' ; !read(*,*)!ck
    call gr_cj(yi, cj) !ck
    print*, 'cj, yi', cj, yi ; !read(*,*) !ck
  endif


  ipo = 1 
  nds = 0 ! the available vectors with the current stepsize in champ 
  
  do 
  ! if we want to change the value of the iteration counter ipo, instead of the automatic way do i = 1,n 
  ! we use if syntax to detect the termination of the do loop
    
    if ( ipo > npo ) exit 
 
 ! back up the initial guess as yi , which will only be update if Newton Method succeeds(ispo =1) 
 ! if it fails, we may start again with ynew(ipo-2) and ds/2,  only for the continuation 
    
 ! --- write the correction onto screen and into file dx.dat to check 
    print*  
    print* , ipo, '-th p.o.','***********', '  dir=', dir, 'ds=', ds ;  print*
    
    write(fdx,*); write(fdx,*) '****************************************'
    write(fdx,*)   ipo,'-th p.o.', '       dir=', dir, '        ds=', ds, '    nds=', nds
    
    ! initial condition before the refinement
    write(fdx, *)  'I.C.         ', y0
    write(*, *)    'I.C.         ', y0  ! print to screen   
 ! --------------- --------------- --------------- --------------- ---------------   
 
    
!---- refine the p.o. with initial and updated refined condition as y0
!   subroutine dfcr(y0, tp0, yf, g, cj,dcj, ispo, niter, deriv, gr_cj)
    call dfcr(y0, tp0, yf, g, cj,dcj, ispo, niter, deriv, gr_cj)  

 
    ! Add the control of stepsize ds here, bound [ds_min, ds_max] for ipo > 1
    ! for each value of ds, do at least 4 orbits, in order to provide cham(4, ctr) for adams
  
    ! if niter > 6,  ds = ds / 2;   
    ! if niter < 3,  ds = ds * 2.
    
    if ( ipo > 1 .and. dsctr == 1)  then 
      
      if(niter .gt. 6) then 
      
        incrds = -1 ! we need to decrease the ds 
        
        ! check if it is possible that ispo = 1 here? if not, delete this part... 
        if(ispo == 1 ) then 
          print*, 'niter>6, but Newton successes!',  'niter=', niter, 'ispo=', ispo
          read*
        endif   
        
       
        if( ispo == 0 ) then 
        ! Deal with the case, where ds=ds_min, but still fails to reach convergent after 6 iterates of Newton Method
        ! regard this case as failure of the continuation and stop
          if (ds .eq. ds_min) then  
            print*, 'P.O. not obtained with the smallest stepsize! Terminate the continuation!'
            read*
            return
         
          else 
          
  !  ----- decrease the stepsize by half, this part has not been checked
  ! since all the computations are not very accurate, in this case we do the linear interpolation to add an additional 'middle' point, 
  ! but keep the number of available p.o.s, without doing the refinement, treat the new 'middile' point as a good one. 
            if (ipo == 2) then 
            ! if the initial guess for ds is too big for the second p.o.. It rarely happens, but we deal with it here just in case
              nds = 1 
              y0 = ynew(ipo-1, 2:7)
              cham(1,:) = chaml(ipo-1, :)
              curv(1,:) = y0(ctr) 
              ds = dmax1(ds_min,  ds / 2)
              
            elseif (ipo > 2) then  
            ! we have at least 2 'good' p.o.s
            !  if ds/2 > ds_min, we can add a middle point bewteen the last 2 good curv,  
            !  else, we only have 1 available good point 
              
              if( ds/2 >= ds_min) then 
                cham(1,:) =  cham(nds-1,:)
                cham(2,:) = ( cham(nds,:) + cham(nds-1,:) ) / 2.d0
                cham(3,:) =  cham(nds,:)
                
                curv(1,:) = ( curv(1,:) + curv(2,:) ) / 2.d0
                nds = 3
                
              else 
                ds = ds_min 
                cham(1,:) = cham(nds,:)
                curv(1,:) = curv(2,: )
                curv(2,:) = 0.d0
                
                nds = 1
              endif
              
            endif 
            ! fail with cham(np,:), try with the new guess of cham, ds and curv 
            ipo = ipo -1
            yctr = curv(min0(nds,2), :) 
            
            call adams(nds, yctr, ds, cham, curv) 
            y0( ctr ) = yctr
            
            ! check y0 
            print*, 'Start again with new y0 :', y0 
            read*
            cycle
          endif
          
        endif   
          
      elseif (niter < 3) then 
       !--------  increase the step size, without exceeding ds_max
        if( nds .ge. 4 ) then 
        
           incrds = 1
 
        endif
             
      endif 
       
    endif 

 ! if there is no control of the step size, we need to detect if the dfcr succeeds or not 
    if (dsctr == 0) then 
       if(ispo == 0) return 
    endif 
 
    print*, 'ds, niter', ds, niter 
!    read*
    
    ! update to save the last 2 points on the characteristic curve
    if(ipo .le. 2) then 
      curv(ipo,:) = y0(ctr)
    else 
      curv(1,:) = curv(2,:)
      curv(2,:) = y0(ctr)
    endif 
    
    
    tp = ntp * tp0  ! the full period
    ynew(ipo,:) = (/tp, y0, cj/)
    
  ! write to file labeled by fpo passed from the main routine, PoInSt.dat  
    write(fpo,'(10e22.14)')  tp, y0, cj, dcj  ! fout - PoInSt.dat 
 
 ! ----------- check the energy ---------------------------   
    if(debug == 1)  then 
      print*, 'finish dfcr! g='
      do i = 1, ntar 
        print*, g(i,:)
      enddo  
      
      print*
      
      print*, 'ds, niter', ds, niter 
      read*
      call gr_cj(y0, cj2)
      
! check to see if there is any point to keep cj, dcj as the output! ckd--cj = cj2 anyway
      print*, 'cj, dcj:', cj, dcj   
      write(*,'(10e22.14)')  tp, y0, cj, dcj  !  !ck print to screen
      read(*,*)
    endif 
!--------------------------------------------
  
    if(npo .eq. 1)  return  !  refinement for only  one p.o.
   
    
! Numerical continuation along arc-length paramter
! Using Adams predictor to get the new initial condition

! be careful if we change the stepsize, the cham here, should have vectors of the same step  

    nds = nds + 1
    
    call champ(nds, g, dir, cham)     ! to modify to general form 
    
    ! For the current 'good' p.o, vector field  is kept as the last available(nds-th) line in cham 
    chaml(ipo,:) = cham(min0(nds, 4), : )
    cham5(nds,:) = cham(min0(nds, 4), : )
    
    if (incrds == 1) then 
       ds = dmin1(ds_max, ds*2) 
     ! we have the previous 5 points to provide us 2 points with stepsize as 2*ds 
       cham(1:3,:) = cham5(nds-4: nds: 2, : ) 
       nds = 3
       cham5 = cham(1:nds,:) !update cham5 to be the one always have the same step size
       incrds = 0
    endif 
    
    ! check the stepsize, available vectors and so on 
!    write(fdx, *) 'nds=', nds 
    write(*, *) 'nds=', nds, 'ipo=',ipo, 'ds=',  ds 
    
    if (debug == 1 ) then
      print*, 'cham'
      do i = 1, 4
        print*, cham(i,:)
      enddo 
      print*
    
      print*, 'chaml'
      do i = 1, ipo
        print*, chaml(i,:)
      enddo 
      print*;     read*
    endif 
    
    yctr = y0(ctr)  
    
    if(debug == 1)  then 
      print*, 'Before adams, check!'
      print*, 'dir, nds, ds, ipo', dir, nds, ds, ipo 
      
      print*, 'yctr',  yctr 
      print*, 'y0',    y0 
      print*, 'cham', cham(nds,:); print* 
      
      print*, 'curv'  
      do i = 1, 2 
        print*, curv(i,:)
      enddo  
       
      read(*,*)   
    endif 

!  subroutine adams(np,x, hh, cam, curv)
    call adams(nds, yctr, ds, cham, curv) 
    
    y0( ctr ) = yctr  

    if (debug == 1) then   
      print*, 'Numerical continuation, new state!' !ck
      print*, y0(1:6)
      print*; read* 
    endif 
    
! we succeed in orbaining an initial guess for the new p.o, update the value of counter      
    ipo = ipo + 1 

  enddo

  close(fout) ! close poin.dat

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
  integer        ::    i, j, debug, ispc
  real(kind=dp)  :: f(ntar), phi(6,6), yi(n), pf(n), vf(6), dx(nctr), dxm, fm, cj2, &
                    yit(6), y0copy(6), dx2(nctr)  
  
  logical :: ok  ! status of the file poit.dat, opened or not, which saves the orbit of all the iterations                  
 
  if (debug == 1)   print*, 'start dfcr, y0!', y0   

  niter = -1 
  ispo  = 1 ! defaultly set Newton Method as successful
  dx = 0.d0;   dxm = 0.d0
  
! Saved the intermediate orbits. if file poit.dat is opened, close it, use inquire to obtain the status, 
  inquire( unit = 12, opened = ok )
  if(ok)  close(12)
  open(12,file='./dat/poit.dat',access ='append',status='replace')
 
 
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
 
      
!for Newton method, which is quadratically convergent. if the iteration is great than 5 or 6, then there is something wrong, provided that the initial guess is a good one, which is  close to the solution.

    if (niter > 6) then 
      print*, 'Newton method iterates > 6, dxm, fm',  dxm,  fm
      if(fm .gt. tol) ispo = 0 ! tol for fm, and prsc for dxm
      read*
      return
    endif   
 
    call gr_cj(y0, cj)

    if (issymt == 2) then 
! General approach: Integrate to the next mirror configuration to compute G and F

!   subroutine gr_fctn(y0, tp0, f, g, yf, deriv, gr_cj )
      call gr_fctn(y0,  tp0,  f, g, yf, deriv, gr_cj)
    
    else 
    
! Poincare map approach
!  subroutine fctn(x,init, f,g, tf,y,vf, ispc, deriv, gr_cj)
      call fctn(y0, 0,  f, g, tp0, yf, pf, ispc, deriv, gr_cj)  

      ! avoid too small p.o.    
      if( tp0  < 1.d-3)  then 
        print*, 'Too small period. tp =', ntp * tp0
        return 
      endif
      
    endif 

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
    read*
   
 ! precision to terminate the correction, 1.d-16 from Gerard, which is of non-sense, because it beyonds the precision of double real
    if ( fm .le. tol)  then 
      print*, 'fm < tol', fm,  tol
      return   ! try 1.d-11  >should be no less than the numerical precision we are using: e=1.d-13 (error in rk78)
    endif 
    
! solve the LLS problem to obtain the correction     
!   subroutine deltx_lns(f, g,  dx)
    call deltx_lns(f, g,  dx)
    print*, 'By Lapack LLS solver: dx'
    print*, dx 
    print*;    read*
    
    
    ! the  precision is the modulus of the correction, or dxm = dnrm2(nctr, dx, 1)
    dxm = dnrm2(nctr, dx, 1)
 
    ! if the correction is below the precision, regard as success of the refinement
    if ( dxm .lt. prsc )  then 
      print*, '|dx|=', dxm, '<', prsc, 'stop iteration!'
      if(debug == 1)   read(*,*) ! in debug mode, pause to check the result 
      return
    endif 
    
   ! update the initial state, if the precision is not satisfied 
    
    
    if(issymt == 2 ) then 
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
  
!***********************************************************************
! Adams-bashforth predictor method, Linear multistep method is to computationally 
! solve the ordinary differential equations

! 20160329  Multistep methods refer to several previous points and derivative values. 
!           Adams method uses fixed step size for all the previous points.

!    So be careful with the input cam 

!          Input parameters:
!  np         	number of computed periodic orbits
!  x(*)       	initial conditions of the last p.o.
!  nctr 	number of components in X
!  hh         	integration step along arc length of solution locus
!  cam(*,*)   	the vector field of the latest four points 
!  curv(2,*)	the latest 2 points along the characteristic curve

!        Output parameters:
!  X(*)       new aproximated initial conditions for p.o.

!  Revised by Yu, 20160411
! --------------------------------------------------------------------------
  subroutine adams(np,x, hh, cam, curv)
  implicit none 

! Input and Output declaration 
  integer, intent(inout) :: np 
!  real(kind=dp), intent(in) ::  cam(4, nctr)
  real(kind=dp), intent(inout) :: hh, x(nctr), curv(2, nctr), cam(4, nctr)
  
! Local Variables 
  integer :: i, incrds
  real(kind=dp) :: v1(nctr), v2(nctr), v12
  
  incrds = 0 
  
    ! -------- check before the detection ------
    print*, 'Check  ds, np, before adams'
    print*,   hh, np  
    print*, 'cam, curv'
    do i = 1, 4
      print*, cam(i,:) 
    enddo  
    
    print*
    do i = 1, 2
      print*,  curv(i,:)
    enddo 
    print*
     

! Adams 4-steps perdictor method for each component of X  
  do i = 1, nctr
   
    if (np .ge. 4) then 
      x(i) = x(i)+ hh*( 55.d0*cam(4,i) - 59.d0*cam(3,i) + 37.d0*cam(2,i) - 9.d0*cam(1,i) ) / 24.d0

    elseif(np == 1) then 
      x(i) = x(i) + hh*cam(1,i)
    
    elseif(np == 2) then 
      x(i) = x(i) + .5d0*hh*( 3.d0*cam(2,i) - cam(1,i) )
    
    elseif(np == 3) then 
      x(i) = x(i) + hh*( 23.d0*cam(3,i) - 16.d0*cam(2,i) + 5.d0*cam(1,i) ) / 12.d0
      
    endif 
  enddo
  
  ! This proves not a good option to adjust the stepsize...... if < cham(np,:), cham(np-1) > is approximately equal to 0, it is not possbile to obtain a good guess given the last two points curv.... 
  ! so instead, we should do the adjustment on cham, instead of on curv 
  
  ! detect if the guess for the third points, based on the vector field of the previous np curv, that would make the angle between 
  ! the lastest 3 points be less than 0.1 rad, i.e., cos( v2, v1) > 0.995
  if (np .ge. 2 .and. dsctr == 1) then 
    v1 = curv(2,:) - curv(1,:)
    
    ! do interpolation to find a good intermediate point that makes smooth continuation
    do 
    
      v2 = x - curv(2,:)
      v12 = dot_product(v1, v2)/ dnrm2(nctr, v1, 1) / dnrm2(nctr, v2, 1)
    
      if (v12 .ge. 0.995d0)  exit
      
      incrds = -1
      x = ( x + curv(2,:) ) / 2  ! linear interpolation 
      
      hh =  dnrm2( nctr, x - curv(2,:), 1)
      print*, 'v1, v2, x'
      print*, v1
      print*, v2
      print*, x 
      read* 
      
      print*, 'Interpolation! ds, v12, X3'
      print*, hh, v12, x 
      read*
    enddo  
    
    ! update the new stepsize ds and the new available curv and the corresponding cham(only 1 point in this case)
    if (incrds == -1) then
      cam(1,:) = cam(np,:)
      curv(1,:) = curv(2,:)
      curv(2,:) = 0.d0
      hh =  dnrm2( nctr, x - curv(2,:), 1)
      np = 1  
    endif 
    
    
   ! -------- check v12 here, to see if we really get a good guess for x ------
    print*, 'Check v12,  ds, np'
    print*, v12, hh, np  
    print*, 'cam, curv'
    do i = 1, 4
      print*, cam(i,:)
    enddo  
    
    print*
    do i = 1, 2
      print*, curv(i,:)
    enddo 
    !--------------------------------------------------------------------------
    
  endif 
 
  return
  end subroutine adams
 
  
!*******************************************************************************
  subroutine champ(np,g, dir, cham)
!  Modified for the asymmetric p.o. case, compute the vector field along the characteristic curve by sloving the Kernel(DG)
!  The kenel V is a combination of Vi, which is the determinant of submatrix DG discarding i-th column
  
!     COMPUTATION OF THE VECTOR FIELD GIVING THE CHARACTERISTIC CURVE
!     OF THE FAMILY

!	  INPUT PARAMETERS:
!     NP	     	NUMBER OF THE LAST COMPUTED P.O.
!     G(*,*)         	JACOBIAN MATRIX OF F(*) (SEE SUBROUTINE TRACA)
!     ntar, nctr	The dimension of matrix G( ntar X nctr)
!     dir	     	SENSE ON THE CHARACTERISTIC CURVE
!	   
! 	  OUTPUT PARAMETERS:
!     cham(*,*)  	VECTOR FIELD ON THE CHARACTERISTIC CURVE AT THE LAST 4
!		     	POINT, ONLY THE LAST ROW OF cham(*,*) IS COMPUTED
!		     	IN THIS SUBROUTINE 
! 	Module-based Parameters
! 	 ntar, nctr 
! 
! --------------------------------------------------------------------------
 
  implicit none 
 
  integer, intent(in) :: np
  integer, intent(inout) ::  dir
  real(kind=dp), intent(in) :: g(ntar, nctr) 
  real(kind=dp), intent(out) :: cham(4, nctr)  
  
! Local Variables 
  integer :: i, j, nnp
  real(kind=dp) :: a(nctr), sm,  c, &
  		   gi(ntar, ntar), det 
 
  ! --- this maybe the key part to be modified
  
  ! the vector field of the characteristic curve. which is its integral curve. 
  ! the characteristic curve is obtained by integrating the vector field a using Adams method
  
  ! the following is the computation of the vector field a, the components of a are the derivative of the control variable w.r.t the arc-length paramter s.
  
! d ctr(i) / d s, where s = sum( s_i**2) is the arc length
! the vector field along the characteristic curve, which is also Ker(G), such that G * v = 0 
! 1. v = | G1 G2 ... Gn |  
!    where Gi (i=1,..,n) is the determinant of the submatrix of G discarding i-th column
! 	 vi = (-1)^(i+1) * det(Gi) 

! 2. Then normalize v to obtain a unit vector field
! Check 20160407  error is 1.d-15 compared with the direct element-by-element computation, discard the latter approach    
 
! compute the component of kernel v  
!   subroutine detmat( a, n, det)
  
  ! check the original g 
!  print*, 'Original G:'
!  do i = 1, ntar
!    print*, g(i,:)
!  enddo 
  
  ! The vector field for the current point, of which the Jacobi Matrix is G
  do i = 1, nctr
   
   ! construct Gi, which is the submatrix of G discarding i-th column
    if ( i== 1) then 
      gi = g(:, 2 : nctr) 
       
    elseif (i < nctr) then  ! discard i-th column
      gi(:, 1 : i-1)    = g(:, 1   : i-1) 
      gi(:, i : nctr-1) = g(:, i+1 : nctr) 
     
    else 
      gi(:, 1 : i-1)    = g(:, 1   : i-1)  
        
    endif 
    
    !  compute the determinant of gi  
    call detmat(gi, ntar, det)
    
    a(i) = det * (-1)**(i+1)
  
  enddo 
  
  ! normalize
  a = a / dnrm2(nctr, a, 1) 
 
 
 ! update the matrix cham, which contain the vector filed of the last four curv
  if(np .gt. 4)  then 
    cham(1:3, :) = cham(2:4, :)  
    cham(4, :)   = a 
    nnp = 4
    
  else 
    cham(np, :) = a  
    if (np .eq. 1)   return
    nnp = np
    
  endif 
 
! detect whether the new point obtained is a bifurcating point. ---- Modify this part 
! If the vector product of the new vector field and the previous one is < 0, change the sign of the vector field
  c = dot_product( cham(nnp,:), cham(nnp-1, :) )
  
!  c = 0.d0
!  do  i = 1,3
!    c = c + cham(nnp,i)*cham(nnp-1,i)
!  enddo

  if( c .gt. 0.d0 )  return 
  
!  this part for the reverse of the sense, detect bifurcation? -- need to be understood later
!  comment this part... to see what happens, the comparision of the result shows what???? 
!  Add the conclusion here! 

  dir = -dir
  cham(nnp, :) = - cham(nnp, :)
  write(*,*)  'Reversal in the Sense along the Characteristic Curve!'
  print*, '<v1, v2>', c 
  print*, cham(nnp-1, :)
  print*, cham(nnp, :)
  read*
  
  return
  end subroutine champ
  
  

!******************************* deltx_lns  ************************************
! Use Lapack driver dgelsy to solve the least norm problem  Ax = b,
! which is better than computing the inverse of the matrix followed by the matrix multiplication 
!   here it is G * X = -F   (A=G : 6 X 7, B = -F)

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

!	Input 
!  ntar		number of rows in G 
!  nctr		number of columns in G
!  g		Jacobi Matrix of F w.r.t (X0 + T) 
!  f		Target function
!   
!	Output 
!  dx 		the least-norm correction by sloving G * DX = -F 

! module pomod based parameter 
! integer:: ntar, nctr 

!  Finaly revised by Yu 20160411
!*****************************************************************************

  subroutine deltx_lns(f, g,  dx)
  
  implicit none 
  integer, parameter :: dp = kind(1.d0)
  
  real(kind=dp), intent(in) :: g(ntar, nctr), f(ntar)
  real(kind=dp), intent(out) :: dx(nctr)

    
! Local Varaibles
! for linear equation solver
  integer, parameter :: lwork = 1600 !  3*nctr + 64*(ntar+1), use 1600 to allocate enough memory for work 
  integer :: jpvt(nctr), rank, info, i, nrhs 
  
  real(kind=dp)  ::  a(ntar, nctr), b(nctr), rcond, work(lwork) !dgelsy
  
  
  !dgelsd
  integer :: iwork(1600) 
  real(kind=dp) :: s(ntar), a2(ntar, nctr), b2(nctr)
  
  external  ::  dgelsy, dgelsd  ! lapack lls solver
  
  ! test the CPU time elapsed for dgelsd and deglsy
  real(kind=dp) :: ti, tf, dt1, dt2
  
  
  nrhs = 1
 
 ! copy g and f, because after the call of dgelsy, a and b will be modified
  a = g 
  b(1:ntar) =  -f
 
 ! check this least-norm solver with the symmetric approach, 2 equations and 3 unknowns
  
! Initialize JPVT to be zero so that all columns are free
  jpvt = 0
 
! Choose RCOND to reflect the relative accuracy of the input data, for double precision, 
! use a small one... 
  RCOND = 1.d-20
  
  call cpu_time(ti)   
! Solve the least squares problem min( norm2(b - Ax) ) for the x of minimum norm.
!  SUBROUTINE DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,  WORK, LWORK, INFO )  
  CALL DGELSY(ntar, nctr, nrhs, A, ntar, B, nctr, JPVT,RCOND, RANK, WORK, LWORK, INFO)
   
  call cpu_time(tf)  
  dt1 = tf - ti
  
  a2 = g 
  b2(1:ntar) =  -f
  CALL DGELSD(ntar, nctr,nrhs, A2, ntar, B2, nctr, S, RCOND, RANK, WORK, LWORK, IWORK, INFO)
  
  
  call cpu_time(ti)
  dt2 = ti - tf 
  
  print*, 'Elapsed time for dgelsy and  dgelsd:', dt1, dt2 
  read*
  
  
! Print solution
  if( info == 0 ) then 
    WRITE (*,*) 'DGELSY-- Least norm solution'
    WRITE (*, '(1X,7e24.14)')  (B(I),I=1, nctr)
  
    WRITE (*,*) 'DGELSD-- Least norm solution' ! this one is more aliable, so better to keep this one maybe 
    WRITE (*, '(1X,7e24.14)')  (B2(I),I=1, nctr)
  
  
    ! assign the solution to dx  
    dx = b
    read*
    
! Print the effective rank of A
!    WRITE (*,*)
!    WRITE (*,*) 'Tolerance used to estimate the rank of A'
!    
!    WRITE (*,'(3X,1P,E11.2)')  RCOND
!    
!    WRITE (*,*) 'Estimated rank of A'
!    WRITE (*,'(1X,I6)') RANK

  else 
    WRITE (*,*) 'Failed to converge!'
    read*
    return
  endif
  
  end subroutine deltx_lns
  
   
end module pomod

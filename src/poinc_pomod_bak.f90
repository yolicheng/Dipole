module symt_pomod 

! 20160411 
!  Final update:  This module is to refine p.o. with poincare map approach, for either symmetric one or asymmetric one. 
!                 The general approach in which the correction is all the whole state vector + the peroid is done is gr_pomod

!  Note: 
!    For the asymmetric approach, note that we are refining w.r.t the initial condition
!    So the Jacobi Matrix w.r.t the variation on the initial state is actually Phi - I  


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

! 8- Check the infinity norm of matrix phi and g - from Lapack
! DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
!    phinormi = dlange('I', 6, 6, phi, 6, phinm)
!   subroutine dgecon(norm,n,a,lda,anorm,rcond,work,iwork,info) -- work(4*n), iwork(n)
!    call dgecon('I',6, phi, 6, phinormi,rcond, work,iwork,info)

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
!        write(fpoinst,'(10d24.16)') tp, yi, cj, dcj, hminim ! PoInSt.dat 
!        ynew(i,:) = (/tp, yi, cj/)
       
!  subroutine dfcr(y0,imax, tf,yf,g, cj,dcj,hminim, ispo, deriv, gr_cj)
!  subroutine fctn(x,init,imax, f,g,y,vf,tf, hminim, ispc, deriv, gr_cj)
!  subroutine poinc(yi,imax, tf,yf, hminim, ispc, deriv, gr_cj) 
!  subroutine sect(y, neq, g, dg)  
 

!   *****************   Public  variables   ***************************** 
!    Assigned by subroutine init_lf before call any function from this module 
!  sec, ind	 	the value to specify the Poincare section y(ind) = sec 
!  ntp 			the real period is TP = ntp * tp(to go to the first crossing)
!  tari, (i=1,2)	target varaibles
!  ctri, (i=1,2,3) 	control varaibles
! ************************************************************************************************

implicit none

save
integer, parameter, private :: dp = kind(1.d0) !to avoid compliance with dp defined in other modules
!( this pofam is supposed to be  self-contained, and for general use,

! ******************** Declare Public *****************
! public by default, the most commonly used varaibles which are also accessible by any subroutine that uses this module
real(kind=dp), parameter ::  mu = 0.12150586d-1   !mu = 0.012150585609624  ! this is just to check halo orbit in other polar orbit exploration

integer, private :: debug, dsctr, issymt   ! do the assignment from main routine, to aviod the compilation of link before every excution
 
integer  :: ntp,  &  ! ntp is the be multiplied to the tf by poinc subroutine to obtain the full period 
            ntar, nctr  ! number of target variables(equations) and control variables(free variables)


! declare  allocatable arrays to store the index of target and control variables  -- to check
integer, allocatable :: tar(:), ctr(:)
            

! Declare external for the routines or functions to be called the routines from library lapack  
!external dgecon ! subroutine-- not inside this module
 external  deriv,  gr_cj ! i don't know if this helps

! -------------- declaration of subroutines from lapack or nag( external package)
real(kind=dp), external :: dlange ! compute the norm of array,  dlange is a function, declare the type
real(kind=dp), external, private :: dnrm2 ! compute the norm of array, from package NAG Fortran Library Routine


 ! the files to check the iteration, matrix norm --discard, mostly caused by numerical error
integer, private :: fphi, fgnew, fdx 
    
!  private - only shared by the subroutines, cannot be accessable from outside the module 
!  can avoid passing the shared parameters level-by-level
integer, private       ::   ind
real(kind=dp), private ::  sec, hmin, hmax, e, & ! Bound values for rk78, for the integration of orbit
			   tmax, tol, prsc ! tolerance control for the termination of Poincare map and Newton method

contains 

!******************** Module-related Parameters ********************************
  subroutine init_poinc(ind0, sec0,  hmin0, hmax0, e0, tmax0)
  ! Assignment of shared variables for the poincare map  + rk78 (for integration)
  integer, intent(in) :: ind0 
  real(kind=dp), intent(in) :: sec0, hmax0, hmin0, e0, tmax0 

! speicially for only call of poinc, ind may be reassigned later by init_asympo or init_sympo subroutine
  ind = ind0 
  sec = sec0

  hmin = hmin0
  hmax = hmax0
  e = e0
  tmax = tmax0 
  
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
! For asymmetric p.o.:  Initialize the index of the control variables and target variables  in the state vector
! asym   		contol  		target 
! The energy h0 and x0 is fixed, look for the first return to x=x0 plane, we only need to ask four final component in state to be as the initial one, the fifth will be satisfied simultaneously.  And we have 5 initial components to modify(except x0)
! constraint funcion:   ds^2 = ( (yf-y0)^2 + (zf-y0)^2 + (vxf-vx0)^2 + (vyf-vy0)^2 ) to be minimum
! 1   	 (y0,z0, vx0,vy0,vz0) - (2,3, 4,5,6)   (yf, zf, vxf, vyf) 

!  To ensure that enough memory is available to allocate space for your array, make use of the STAT option of the ALLOCATE command:

!   ALLOCATE ( A(N), STAT = AllocateStatus)
!   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!Here, AllocateStatus is an integer variable. AllocateStatus takes the value 0 if allocation is successful or some other machine dependent value of there is insufficient memory.

!An array can be released from memory by using the DEALLOCATE command
!   DEALLOCATE (A, STAT = DeAllocateStatus)
! ******************************************************************************    
  subroutine init_asympo(asymt)

  integer, intent(in) :: asymt  ! index of symmetry, used to assign the public values accordingly
  
  integer :: allocatestatus
 
! not symmetric, so we need to return to the initial point...
!  this means, we need to take the second intersection with the 
! poincare section x = x0, with the velocity in x component to be in the same direction with the initial one?

! flag of symmetric p.o. 
  issymt = 0
  
! need to check .....
  ntp = 1 !  
   
  if ( asymt == 1 ) then 
  ! fix the value of x, and try to modify the other component to make the first return to x=x0 plane goes to the same point
    ind = 1 
    
  ! number of control and target variables
    ntar = 4 
    nctr = 5 
  
  ! allocate memory for ctr and tar   
    allocate ( ctr(nctr), stat = allocatestatus)
    if (allocatestatus /= 0) stop "*** not enough memory ***"
    
    allocate ( tar(ntar), stat = allocatestatus)
    if (allocatestatus /= 0) stop "*** not enough memory ***"

    tar = (/2,3,4,5/)   
    ctr = (/2,3,4,5,6/) 
  
  elseif (asymt == 2) then 
    ! fix the value of y, and try to modify the other component to make the first return to x=x0 plane goes to the same point
    ind = 2
    
  ! number of control and target variables
    ntar = 4 
    nctr = 5 
  
  ! allocate memory for ctr and tar   
    allocate ( ctr(nctr), stat = allocatestatus)
    if (allocatestatus /= 0) stop "*** not enough memory ***"
    
    allocate ( tar(ntar), stat = allocatestatus)
    if (allocatestatus /= 0) stop "*** not enough memory ***"

    tar = (/1,3,4,5/)   
    ctr = (/1,3,4,5,6/) 
     
  endif 
 
  if (debug == 1) then          
    print*, 'check init_asympo: asymt, ind, y0, control and target var'
    print*, asymt, ',', ind,sec, ',', ctr, ',', tar 
    read(*,*)
  endif
  
  end  subroutine init_asympo 


! ****************************************************************************** 
! For asymmetric p.o.:  Initialize the index of the control variables and target variables as in the state vector
!   The Poincare section is specified as:  y(ind) = sec
!   Note: to use these symmetries, the invariant transformation in time should be t-->-t
! question: symt = ind? NO! for vertical lyapunov orbit, we take y=0 as Poincare section, but use symmetry 2, so cannot join as one parameter

! It is not ok to just pass the parameter as input, without assignment by h=h0
! ****************************************************************************** 
  subroutine init_sympo(symt )
  
  integer, intent(in) :: symt  ! index of symmetry, used to assign the public values accordingly

  integer ::  allocatestatus
  
  ! flag of symmetric p.o. 
  issymt = 1
  
  
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
  
  end subroutine init_sympo
  
! ****************************************************************************** 
!  This  subroutine is to open the file to save the intermediate result, including
!  the files to check the correction,  matrix norm  for every iteration of Newton method 
!  Turns out the matrix is in good condition, because gnewnorm is small, that means the error will lead to a small correction(multiplied by the matrix),
!  but pay attention to the numerical error, which is 1.d-13 previously, but I ask the tol<1.d-14, stupid!!!!
 
! ****************************************************************************** 
  subroutine init_write  

!! Bound values for RK78, shared by all subroutines that calls rk78 (poinc, fctn)
  if(debug == 1) then 
    print*, 'hmin,hmax, e, tmax,tol,prsc', hmin,hmax, e,tmax,tol,prsc ;  read(*,*)
  endif 
  
  fphi = 30;  fgnew = 31;  
  open(fphi,file='./dat/phinorm.dat',access ='append',status='replace')
  open(fgnew,file='./dat/gnewnorm.dat',access ='append',status='replace')   

  write(fphi,*) '# Phi :    norm_Inf      norm_row (6)' 
  write(fgnew,*) '# g^T(g*g^T)^(-1) :     norm_inf      norm_row(', nctr, ')' 
  
  fdx = 32 
  open(fdx,file='./dat/dx.dat',access ='append',status='replace')  
  write(fdx,*) '# dxm	 errfm	 dx	index of correction(', ctr, ')'
    
  end subroutine init_write
  
  

!**************************   Differential Correction  ************************* 
!  Using modified Newton method to refine the p.o, provided with the initial conditon.
!  We have nctr unkowns and ntar equations(nctr>ntar, inderterminant equations), we take least-square solution.
!  The point is get the derivative of the target variables with respect to the control variables

! The cases when this refinement fails to reach convergence
!   1. too many iterations, niter>6.  based on the quadratic convergence of Newton method, 6 iteration is enough!
!   2. too long time to arrive to the next intersection acrossing poincare section, ispc = 0


! for Newton method, alway check two things:
!   1. the error in target varaibles    --- stop condition, satisfied means finish the correction
!   2. the modulus of the correction    --- if too small, no point to continue ...

!  	Input Variables:
!  y0(6)	initial condition(x,y=0,z, vx=0, vy, vz=0)
!  imax 	the times of the crossing through the section(y=0)

!  	Output Variables:
!  y0, yf	the refined initial and final conditon
!  tf		half period or one period, depends on the number of intsec-imax
!  g(ntar, nctr)	JACOBIAN MATRIX OF F(*), ntar by  nctr,  for symmetric p.o. :  3 variables and 2 equations
!  cj, dcj	the initial and variation of the jacobi constant
!  hminim	minimum step used in the integration of the p.o.

!   	Private Module-base Varaibles
!  tol		tolerance for vx,vz to be zero (ie,1.d-13 or 1.d-16 in gerard's routine)
!  prsc 	presion to terminate the iteration(ie, 1.d-13)
  
!     the equations: f1= vx = 0, f2 = vz = 0 symmetric wrt xz plane : halo and planar lyapunov orbit, 1st crossing in T/2
!                    f1= z = 0,  f2 = vx = 0 symmetric wrt x-aixs : vertical lyapunov orbit, first crossing in T/4
!  x-axis: (x,y,z,t) --> (x, -y,-z, -t)
  

!  function used: poinc, deriv, gr_cj, fctn, deltx
! --------------------------------------------------------------------------
  subroutine dfcr(y0,imax, tf,yf,g, cj,dcj,hminim, ispo, niter, deriv, gr_cj)
   
  implicit none  
  integer, parameter :: n = 42 ! both the state and the variational equations are desired, so n=42
  
  integer, intent(in)  :: imax ! number of crossings of the Poincare section
  integer, intent(out) :: ispo, niter
  real(kind=dp), intent(out) :: tf, yf(n), g(ntar, nctr), cj, dcj, hminim 
  real(kind=dp), intent(inout)  :: y0(6) 
  external :: deriv, gr_cj
  
! Local Varaibles
  integer ::   ispc, i, j
  real(kind=dp)  :: f(ntar), phi(6,6), yi(n),pf(n), dx(nctr), dxm, fm, cj2, &
                    yit(6), y0copy(6), dx2(nctr)
                    
  logical :: ok  ! status of the file poit.dat, opened or not, which saves the orbit of all the iterations                  
  
  if (debug == 1)   print*, 'start dfcr, y0!', y0   

  niter = -1 
  ispo = 1 ! defaultly set Newton Method as successful
  dx = 0.d0; dxm = 0.d0
  
! Saved the intermediate orbits. if file poit.dat is opened, close it, use inquire to obtain the status, 
  inquire( unit = 12, opened = ok )
  if(ok)  close(12)
  open(12,file='./dat/poit.dat',access ='append',status='replace')
  
  
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
    
!  subroutine fctn(x,init,imax, f,g,y,vf,tf, hminim, ispc, deriv, gr_cj)
    call fctn(y0, 0, imax, f, g, yf, pf, tf, hminim, ispc, deriv, gr_cj)  
    
    if( tf  < 1.d-3)  then ! avoid too small p.o. 
      print*, 'Too small period. tp =', ntp * tf
      return 
    endif
     
    
    if(debug == 1) print*, 'After fctn: tf, yf', tf, yf(1:6)
    
! ---- plot the intermediate orbit to check -------------   
! possible debug here, if tf < 1.d-3, which is the stepsize of the first step of integration
! it is not possible to start here....

!    subroutine plob(y0,t0,tf, n, tdir, ftag, deriv, gr_cj,  y)    
!    call plob(y0,0.d0, tf, 6, 1, 12, deriv, gr_cj,  yit)   ! there is some problem here?
    
    if ( ispc == 0 ) then
      ispo = 0 
      close(12)
      return 
    endif 
! ------------------------------------------
   
!   check the variation in energy for the poincare map here, instead of inside the fctn subroutine    
    call gr_cj(yf(1:6), cj2)
    dcj = cj2 - cj  

    if(debug == 1) then
      print*, 'ck cj! dcj, cj0, cjf', dcj, cj, cj2   !ckd ! very small difference
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
   
 ! precision to terminate the correction, 1.d-16 from Gerard, which is of non-sense, because it beyonds the precision of double real
    if ( fm .le. tol)  then 
      print*, 'fm < tol', fm,  tol
      return   ! try 1.d-11  >should be no less than the numerical precision we are using: e=1.d-13 (error in rk78)
    endif 
    
    ! compare the two approach of solving dx
    ! compute the correction, using the least-norm solution 
    call deltx_lns(f, g, ntar, nctr, dx)
    print*, 'lns: dx', dx
    read*
    
    call deltx(f,g, dx2)
     
    print*, 'inv: dx', dx2 
    read*
       
    ! the  precision is the modulus of the correction, or dxm = dnrm2(nctr, dx, 1)
    dxm = dnrm2(nctr, dx, 1)
 
 
    ! if the correction is below the precision, regard as success of the refinement
    if ( dxm .lt. prsc )  then 
      print*, '|dx|<', prsc, 'stop iteration!'
      if(debug == 1)   read(*,*) ! in debug mode, pause to check the result 
!      return
    endif 
    
   ! update the initial state, if the precision is not satisfied 
    y0(ctr) = y0(ctr) + dx
    
    write(fdx,*) ! add a blank line to seperate iterations
    write(fdx,*) niter+1, '-th iteration'
    write(fdx,*) '|DX0|,|f|,DX0', dxm, fm, dx 
    write(fdx,*) 'U.I.C.       ', y0 ! the updated 
    
    ! print to screen 
    write(*,*)
    write(*,*) niter+1, '-th iteration'
    write(*, *)  '|DX0|,|f|,DX0', dxm, fm, dx
    write(*,*)   'U.I.C.       ', y0 
    if(debug == 1)  read*
    
  enddo 
   
  return
  end subroutine dfcr


! ***************************** fctn *****************************************
!  by Gerard, without any modification, Only available for planar lyapunov and halo orbit 
!   modified by Yu for general purpose, 20160218
  
!     AUXILIARY SUBROUTINE FOR THE COMPUTATION OF THE FUNCTION F(*) 
!	 AND ITS JACOBIAN MATRIX AFTER A HALF OR ONE REVOLUTION

!          INPUT PARAMETERS:
!     X(*)     	 	INITIAL POINT (x,y=0,z,xdot=0,ydot,zdot=0)ON y=0 WITH xdot=0,zdot=0
!     y(*)     	 	INITIAL POINT IF IND=1
!     init      	Flag FOR DETERMINATION OF INITIAL POINT, 1: y<-vf(initialized by vf) 
!	 	 	0: initialize the first 6 components to be the state x, 7-42: identity matrix for phi
!     imax       	number of crossing through the Poincare section 
!     deriv/gr_cj 	 The subroutine to compute vector field: gr_rtbp or gr_lf

! 		OUTPUT PARAMETERS:
!     F(*)       The target variables(f1=0, f2=0 )
!     G(*,*)     JACOBIAN MATRIX OF F(*) (w.r.t the control variables)
!     y(*)       FINAL POINT UNDER THE Poincare MAP
!     vf(*)       VECTOR FIELD AT Z(*), also used to pass the initial value of y if (init=1)
!     tf        time to go from the initial point to z(*)

!! discard!   cj, dcj     the initial and variation of the jacobi constant
!!     hminim     minimum step used in the integration of the p.o.

! 		Global Parameters  from module pofam
!   ind, tar, ctr, ntar, nctr 

! period:  symmetric wrt xz plane, 1st crossing at T/2
!          symmetric wrt x-aixs : vertical lyapunov orbit, first crossing in T/4

!  subroutine used: poinc, deriv
!********************************************************************************
  subroutine fctn(x,init, f,g, tf,y,vf, hminim, ispc, deriv, gr_cj)
  
  implicit none  
  integer, intent(in)  :: init 
  integer, intent(out)  ::  ispc
  real(kind=dp), intent(in)  :: x(6)
  real(kind=dp), intent(out) :: f(ntar), g(ntar, nctr), vf(42), tf, hminim  
  real(kind=dp), intent(inout)  :: y(42)
  
! local varaibles
  integer        ::  i, j 
  real(kind=dp)  ::  phi(6,6), yi(42), f1, f2, & 
  		     phinm(6), phinormi, gnm(ntar), gnormi, cj, &
  		     vftar(ntar,1)  
    
  external deriv, gr_cj ! the vector field :  derive(t,x,n,f)
 
  
  if(debug == 1)  then 
    print*, 'start fctn, ind=', ind, 'x0=', x
    read(*,*)
  endif      
   
  if (init .eq. 0) then  ! not initialized
  ! y0 and variational matrix need to be initialized here 
    yi = 0.d0
    yi(1:6)= x  
    yi(7:42:7) = 1.d0
  else 
    yi = y  
  endif  

!  subroutine poinc(yi,imax, neq, tf,yf, hminim, ispc, deriv, gr_cj) 
  call poinc(yi,imax, 42, tf,y,  hminim, ispc, deriv, gr_cj) ! imax determines the number of crossing we are considering as
  
  if (ispc == 0) return !fail to reach to Poincare section, stop  the correction 
  
! the target variables : is the difference between the imax-return to poincare section and the initial point
  f = y( tar ) - x( tar )

! check if the assignment by indice is ok  ! ckd  
  if (debug == 1) then  
    print*, 'tar, f(', ntar,')', tar, f  
    print*, 'y(tar) =', y(tar) 
    read* 
    print*
  endif
  
  call deriv(0.D0,y, 42, vf)  ! compute the vector field
 
! G is the Jacobi matrix of f( ntar target variables) with respect to nctr control variables
! and Phi is the variatioal matrix computed as y(7-42)
  phi = reshape(y(7:42), (/6,6/))   ! reshape in column wise 
 
 ! For the asymmetric case, G = Phi - I  - vf( tar(i) ) / vf(ind) * Phi
 ! remember not to forget to abstract the identity matrix 
  if (issymt == 0 ) then 
    print*, 'Asymmteric case! Phi= Phi-I'
    read*
    do i = 1, 6
      phi(i,i) = phi(i,i) - 1.d0
    enddo
  endif 
  
 
! G = Phi - vf( tar(i) ) / vf(ind) * Phi, the coefficient matrix of target error to obtain correction
  do i = 1, ntar 
    do j = 1, nctr 
      g(i,j) = phi( tar(i), ctr(j) ) - phi(ind, ctr(j) ) * vf( tar(i) ) / vf(ind) 
    enddo 
  enddo  

 !------------------------------------------------------
 ! for debug, to check the key variables
  if (debug == 1) then 
    print*, 'phi'
    do i = 1, 6
      write(*,'(6f20.14)') phi(i,:)
    enddo 
    
    print*, 'ind, tar, ctr'
    print*, ind, ',',  tar, ',', ctr ; print*
    
!    print*, 'beta4=v(4)/v(2), beta6=v(6)/v(2)'
!    print*, vf(4)/vf(2), vf(6)/vf(2), vf(4), vf(6), vf(2)
! 
    print*, 'Jacobi Matrix of target variable w.r.t control variables' 
    do i = 1, ntar 
      write(*,*)  phi(tar(i), ctr)
    enddo 
    read*
    
!    print*, 'Jacobi Matrix for the correction'
!    do i = 1, ntar
!      print*, g(i,:)
!    enddo 

  ! --- check the infinity norm of matrix phi and g
! DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
    phinormi = dlange('I', 6, 6, phi, 6, phinm)
    print*
    print*, 'Phi: norm_inf, norm_row' 
    print*, phinormi, phinm
    write(fphi,*) phinormi, phinm
    read*
  endif 
  
  return
  end subroutine fctn
  

!*******************************************************************************    
  subroutine poinc(yi,imax, neq, tf,yf, hminim, ispc, deriv, gr_cj) 
! Determination of the imax-th passage of the orbit through the Poincare section defined by subroutine sect(here,y=0) 
!  Note:  only one intersection is record here!
! to make this subroutine more general, use deriv as an argument as  the function to compute the vector field
!   of the common form:  deriv(T,X,N,F)
     
! 	Input Variables
!  yi(*)	initial point 
!  imax		the number of time for intersecting the section
!  neq 		dimension of the state: 6--position+velocity  42: ...+ variational matrix

! 	Private Module-based Varaible
!  tmax, hmax, hmin, e, dp

!  	Output Variables
!  tf  	 	time spent by the orbit to go from yi to yf
!  yf(*)  	first cut of the orbit which passes by yi with the surface of section
!  hminim	minimum step used in the integration of the p.o.
!  ispc 	flag of the success to return to the poincare section

! function used: gr_rk78, sect  
! --------------------------------------------------------------------------

  implicit none  
!  integer, parameter  ::  neq = 42 ! dimension of the state: 6-state vector, 42-also the variatioal matrix
  
  integer, intent(in)  ::  imax, neq
  integer, intent(out)  ::  ispc
  real(kind=dp), intent(in)  :: yi(neq) 
  real(kind=dp), intent(out) :: tf, yf(neq), hminim
  external deriv, gr_cj
  
! Local varaibles
  integer :: i, iter
  real(kind=dp)  ::  g, gi,  dg(neq), t, dh, dy, &
  	  	     y(neq), r(13,neq),b(neq),f(neq), h  ! gr_rk78 
  
  ispc = 1 ! default value is 1, if fails, set to 0
  call sect(yi, neq, g, dg)
  
!! A little trick to avoid the mis-detection of sign change at the start point, ie, g0 = -1.d-9, g1 = 1.d-9, that is not what we want  
!  if(dabs(g) .lt. 1.d-9)  g = 0.d0   ! cancel numberical error, if the compoent of the state is less than 1.e-9, treat as 0 
! 
! initial state
  h = 1.d-3
  t = 0.d0 ! initial time
  hminim = h
  y = yi 
  
  if (debug == 1) then 
    print*, 'Poinc- hmin,hmax,e, h, t, tmax, imax', hmin,hmax,e, h, t, tmax, imax
!    print*, 'y0', yi(1:6)
    print*
    read*
  endif
  
      
  do i = 1, imax
    
    if(dabs(g) .lt. 1.d-9)  g = 0.d0   ! cancel numberical error, if the compoent of the state is less than 1.e-9, treat as 0 
     
    do  !look for the next intersecion across the Poincare section 
      gi = g ! the previous value of g
    
      call gr_rk78(t,y,neq, h,hmin,hmax,e, r,b,f, deriv )
      
      if (debug == 1)   write(*,'(8f18.14)') t, y(1:6), h ! print the correction process on the screen
      
      call sect(y, neq, g, dg)
      
!if it spends too much time to go to the next intersecion across y(ind)=sec plane, seen as failed 
      if (dabs(t) > dabs(tmax) ) then 
        ispc = 0
        print*, 'Maximum time exceeded!, t>tmax', t, tmax
        return
      endif
      
      hminim = dmin1(hminim, h) 
      if (gi*g < 0.d0 ) exit ! terminate of different sign
      
    enddo
  enddo  
 
  if(debug == 1) then 
    print*, 'ck crossing the Poincare section', gi, g;     
    read(*,*) !ck
  endif  

!   REFINEMENT OF THE INTERSECTION POINT YF(*) USING THE NEWTON'S METHOD
!   TO GET A ZERO OF THE FUNCTION G (SEE SUBROUTINE SECCIO), precision is 1.d-13 ?
! pay attention here, there is possibility that it falls into endless iteration here
  iter = 0
  if (debug == 1) print*, 'Iteration for poincare map! t - g - dh'
  do  
    iter = iter + 1
!    if (debug == 1) print*, iter,'-th iteration for intersection with poincare section'
    if (iter > 10 ) then 
     print*, 'Iteration exceeds 10!'; read*
     ispc = 0
     return
    endif 
    
    if(dabs(g) .le. 1.d-13)  exit  ! detect the y(ind)-sec is within the tolerance, values 1.d-13 from Gerard's routine
    
    call deriv(t, y, neq, f)

    dy = 0
    do i = 1, neq
      dy = dy + f(i)*dg(i)
    enddo
    dh = - g/dy
    
! print*, 'before gr_rk78, dh', dh
! to make sure this step works, we need to set hmin to the required stepsize 
! Need to debug here, for possible stuck here....

    call gr_rk78(t,y, neq, dh, dabs(dh),hmax,e, r,b,f, deriv) ! remember that gr_rk78 is only 1 step
    
!    call gr_rk78(t,y, neq, dh, dabs(dh), dabs(dh),e, r,b,f, deriv) ! remember that gr_rk78 is only 1 step
    
!    print*, 'after gr_rk78'
    call sect(y, neq, g, dg)
    
    if(debug == 1) then 
      print*, iter, ':', t, g, dh    
      read*   !ck
    endif
      
  enddo 
 
! we get y=0, t is tp, yi is the intersection	
  yf = y
  tf = t
  
  if(debug == 1) then 
    write(*,*)'Poinc finished, g, tf, yf:', g, tf
    print*, yf(1:6); print* 
  endif 

  return
  end subroutine poinc 
 
!*****************************************************************************
  subroutine sect(y, neq, g, dg)
  ! the surface of section, defined by g =  y(ind) - y0
! 	Input parameters:
! y(*)      	the state vector 
! neq 		the dimension of the state y, 6: position+velocity; 42: + variational equation

!  	Output parameters:
! g 		funciton that equated to 0 gives the surface of section
! dg(*) 	gradient of function g


!  Global Variables from module  
! ind  		the index of the component of y to be used as Poincare section  
! sec 		the Poincare section to be specified as y(ind) - sec

! 20160218 -by Yu
! --------------------------------------------------------------------------

  implicit none 
  
  integer, intent(in) :: neq
  real(kind=dp), intent(in) :: y(neq) 
  real(kind=dp), intent(out) :: g, dg(neq)

  integer :: i

  g = y(ind) - sec
   
  do  i = 1, neq
    dg(i) = 0
  enddo 

  dg(ind) = 1.d0

  return       
  end subroutine sect


! ****************************  deltx  *****************************************
!  20160407  
!  discard this approach by computing the inverse of the matrix followed by matrix multiplication
!  replace with the least-norm solution by Lapack  Linear least square (LLS) solve deglsy  
!  The latter one can deal with ill-conditinal matrix 

!  Nov.26,2014 -- general solution for the Modified Newton method with 3 free variables(control), 2 target equations
!  We have Xk = X(k-1) + deltX(k-1)
!  where we require :   G ( deltX(k-1) ) = - F( X(k-1) )

!  It is better to just do the computation here, and check if we reach to the tolerance in the upper subroutine

!  Minimize the norm
!     (deltX)T * Q * deltX
!      Q is a given diagonal positive definite weight matrix,  here,just set it to identity

!  So, the solution of this problem of minima is given by
!     (deltXk-1)= -Q-1 (G)T * (G Q-1 (G)T)-1 * F( X(k-1) ) 
!     (deltXk-1)= -(G)T* (G*(G)T)-1 * F( X(k-1) ) 

!	Input Variables
!  f(*)		Target Variables to be set 0, (here vx,vz )
!  g(*,*) 	Jacobian Matrix of f(*) 

!	Output Variables
!  dx   	the difference for the initial state to meet the final target

!	Module-based Varaible
!  ntar, nctr 

!  subroutine used: none
! ------------------------------------------------------------------------------
 
  subroutine deltx(f,g, dx)
  implicit none 
  
  real(kind=dp), intent(in)  :: f(ntar), g(ntar, nctr) 
  real(kind=dp), intent(out) :: dx(nctr)

! Local Variables  
  integer :: i, j, k, &
  	     n, info ! for the norm computed by lapack
  real(kind=dp) :: g2(ntar, ntar), g2inv(ntar,ntar), gnew(nctr, ntar),  & 
                   gnm(nctr),  gnormi, g2m, g2inv2(ntar, ntar), gnew2(nctr, ntar), dx2(nctr)! compute the matrix norm to see if the equation is in well condition

!  checked minv, the inv using lapack is not accurate enough, use direct computation subroutine minv instead

!  g2 = g*g^T   
  g2 = matmul( g, transpose(g) )
 
  call minv(g2, g2inv, ntar)
!  read*
  
  gnew  = matmul( transpose(g), g2inv )

! Dx = - Gnew * F :   nctr-by-ntar  X  ntar-by-1 ===> nctr-by-1
  dx = - matmul(gnew, f)
  
  if(debug == 1) then ! check the matrix norm of g^T(g*g^T)^(-1)
!   DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )

    gnormi = dlange('I', nctr,  ntar, gnew,  nctr, gnm)
    print*, 'g^T(g*g^T)^(-1), : norm_inf,  norm_row(3)'
    print*,  gnormi, gnm
    write(fgnew,*) gnormi, gnm
    print*
    
!    print*, 'g^T(g*g^T)^(-1)'
!    do i = 1, nctr
!      print*, gnew(i,:)
!    enddo

    read*
    print*, 'dx', dx  ! ckd, exactly the same
  endif 
     
  return      
  end subroutine deltx  
end module symt_pomod

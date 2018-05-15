module asym_pomod 
! 20160315
! this is for asymmetric p.o. computation, in this case,  we have 4 equations and 5 unkowns
! unlike symmetric case, where we have 2 equations and 3 unkowns: 2 target varaibles and 3 control variables, but for asymmetric orbit,  it is no longer available....


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
!  subroutine init_po(sym, sec0, hmin, hmax, e, tmax, tol, prsc)
!  subroutine pofam(yi,npo,imax, dir,ds, fpo, ynew, i, deriv, gr_cj)
!    write(fpoinst,'(10d24.16)') tp, yi, cj, dcj, hminim ! PoInSt.dat 
!    ynew(i,:) = (/tp, yi, cj/)
       
!  subroutine dfcr(y0,imax, tf,yf,g, cj,dcj,hminim, ispo, deriv, gr_cj)
!  subroutine fctn(x,init,imax, f,g,y,vf,tf, hminim, ispc, deriv, gr_cj)
!  subroutine poinc(yi,imax, tf,yf, hminim, ispc, deriv, gr_cj) 
!  subroutine sect(y, neq, g, dg)   
!  subroutine deriv(t,x,n,f)
!  subroutine gr_cj(x, cj)

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
 
integer  :: ind, ntp, &  ! 
            tar1, tar2, ctr1, ctr2, ctr3  ! index of the control and target variables

! Declare external for the routines or functions to be called the routines from library lapack  
!external dgecon ! subroutine-- not inside this module
 external  deriv,  gr_cj ! i don't know if this helps

real(kind=dp), external :: dlange ! dlange is a function, declare the type

 ! the files to check the iteration, matrix norm --discard, mostly caused by numerical error
integer, private :: fphi, fgnew, fdx 
  
!  private - only shared by the subroutines, cannot be accessable from outside the module 
!  can avoid passing the shared parameters level-by-level
integer, private       ::  debug
real(kind=dp), private ::  sec, hmin, hmax, e, & ! Bound values for rk78, for the integration of orbit
			   tmax, tol, prsc ! tolerance control for the termination of Poincare map and Newton method

contains 

!  include 'gr_rk78.f90' ! should be a stand alone subroutine
      
!******************** Module-related Parameters ********************************
  subroutine init_poinc(ind0, sec0, hmin0, hmax0, e0, tmax0)
  ! for the poincare map 
  real(kind=dp), intent(in) :: ind0, sec0, hmax0, hmin0, e0, tmax0 
  
  ind = ind0
  sec = sec0
  
  hmin = hmin0
  hmax = hmax0
  e = e0
  tmax = tmax0 
  
  end subroutine init_poinc
  
  subroutine init_errctr(tol0, prsc0)
  ! the error control for termination of Newton Method
  real(kind=dp), intent(in) :: tol0, prsc0
  
  tol = tol0 
  prsc = prsc0
  
  end subroutine init_errctr
  
 
  subroutine init_asympo(asymt, )
! Initialize the index of the control variables and target variables  in the state vector
! The Poincare section is specified as:  y(ind) = sec
! question: symt = ind? NO! for vertical lyapunov orbit, we take y=0 as Poincare section, but use symmetry 2, so cannot join as one parameter

!It is ok to just pass the parameter as input, without assignment by h=h0

  integer, intent(in) :: asymt ! index of symmetry, used to assign the public values accordingly
  real(kind=dp), intent(in) :: sec0, & !poincare section

! -------------- Debug -------------------  
!  debug = 1  ! for debug
  debug = 0 ! normal execution 

 
    
  !! Bound values for RK78, shared by all subroutines that calls rk78 (poinc, plob, fctn)
  print*, 'hmin,hmax, e, tmax,tol,prsc', hmin,hmax, e,tmax,tol,prsc ; read(*,*)
  
  if(debug == 1) then 
! the files to check the iteration of Newton method, matrix norm
! turns out the matrix is in good condition, because gnewnorm is small, that means the error will lead to a small correction(multiplied by the matrix), but pay attention to the numerical error, which is 1.d-13 previously, but I ask the tol<1.d-14, stupid!!!!

    fphi=55;  fgnew=56;  fdx=57
    open(fphi,file='./dat/phinorm.dat',access ='append',status='replace')
    open(fgnew,file='./dat/gnewnorm.dat',access ='append',status='replace')   
    open(fdx,file='./dat/dx.dat',access ='append',status='replace')  
   
    write(fphi,*) '#Phi :    norm_Inf      norm_row (6)' 
    write(fgnew,*) 'g^T(g*g^T)^(-1) :     norm_inf      norm_row(3)' 
    write(fdx,*) 'dxm	 errfm	 dx(3)'
  endif  
     
!  mu = 0.12150586d-1   !mu = 0.012150585609624 --- just to test asym_pomod with halo
! note: to use these symmetries, the invariant transformation in time should be t-->-t

! ***************************  variable specification *************************
! asym   		contol  		target 
! The energy h0 and x0 is fixed, look for the first return to x=x0 plane, we only need to ask four final component in state to be as the initial one, the fifth will be satisfied simultaneously.  And we have 5 initial components to modify(except x0)
! constraint funcion:   ds^2 = ( (yf-y0)^2 + (zf-y0)^2 + (vxf-vx0)^2 + (vyf-vy0)^2 ) to be minimum

! 1   	 (y0,z0, vx0,vy0,vz0) - (2,3, 4,5,6)   (yf, zf, vxf, vyf) 
   
   
  if (symt == 1) then  ! 
    ind = 1;   ntp = 2 ! the Poincare sect y(ind) = sec 
    ctr1 = 2;  ctr2 = 3;  ctr3 = 4;  ! control variables
    tar1 = 5;  tar2 = 6 ! target variables
      
  elseif (symt == 2)  then 
    ind = 2;   ntp = 2   
    ctr1 = 1;  ctr2 = 3;  ctr3 = 5;   
    tar1 = 4;  tar2 = 6  
      
  elseif (symt == 3)  then 
    ind = 3;   ntp = 2   
    ctr1 = 1;  ctr2 = 2;  ctr3 = 6;  
    tar1 = 4;  tar2 = 5 

! axis symmetry      
  elseif (symt == 4)  then  ! x-axis?
    ind =  2;  ntp = 4  
    ctr1 = 1;  ctr2 = 3;  ctr3 = 5;   
    tar1 = 3;  tar2 = 4  
      
! the other cases are to be added, not needed at this moment         
  endif
   
  if (debug == 1) then          
    print*, 'check init_po: symt, ind, y0, control and target var'
    print*, symt, ind,sec, ctr1,ctr2,ctr3, tar1,tar2
    read(*,*)
  endif   
  end subroutine init_asympo
  
  
!***********************************************************************  
  subroutine pofam(yi,npo,imax, dir,ds, fpo, ynew, i, deriv, gr_cj)
! Main subroutine to compute np new p.o.s with one initial guess yi(6)
! To make sure the initial state is periodic, refine it without check
!  the initial condition for all p.o.s are stored in file fpo

! 	Input Varialbles
!   yi(6)	initial condition (ie, x,y=0,z, vx=0, vy, vz=0)
!   npo		number of new p.o.  if np=1, just do differential correction, without continuation
!   imax	number of intersections with Poincare section to get p.o.

!   dir		dirction along the vector field, +1:increase,-1:decrease
!   ds		displacement for new PO
!   fpo		file tag for the initial states of the p.o.s

! 	Output Variables
!   ynew(npo,8)	the initial condition for new p.o.s, npo*8
!     		save the initial p.o in the first row
!   i 		the real number of available p.o.s is (i-1)

!  Routine used: dfcr, fctn, champ, adams
!  Finally revised : 20160218  
! --------------------------------------------------------------------------

  implicit none
  integer, parameter :: n = 42    
!  both the state and the variational equations are desired, so n=42

  integer, intent(in) 	   :: npo, imax, fpo 
  integer, intent(out) 	   :: i
  integer, intent(inout)   :: dir  ! may be reversed
  
  real(kind=dp), intent(in)  :: ds !tol,  prsc, 
  real(kind=dp), intent(out) :: ynew(npo,8)
  real(kind=dp), intent(inout)  :: yi(6) ! initial guess is updated

  external :: deriv, gr_cj  ! vector field  & jacobi constant
  
! Local Variables
  integer :: fout, ispo
  real(kind=dp) :: yf(n), f(2), g(2,3), cham(4,3), tp, yctr(3), &
   		   cj, cj2, dcj, hminim
  
  ynew = 0.d0
  
  if (debug == 1) then 
    print*, 'hmin,hmax,e, tol, prsc ', hmin,hmax,e, tol, prsc 
    read(*,*)
  endif
  
  if (debug == 1) then 
    print*, 'start pofam!' ; !read(*,*)!ck
    call gr_cj(yi, cj) !ck
    print*, 'cj, yi', cj, yi ; !read(*,*) !ck
  endif

  do i = 1, npo
   print* , '***********',i, '-th p.o.','***********'
   print*
   
!  subroutine dfcr(y0,imax, tf,yf,g, cj,dcj,hminim, ispo, deriv, gr_cj)
    call dfcr(yi,imax, tp,yf, g, cj, dcj, hminim, ispo, deriv, gr_cj ) 
    if( ispo == 0 ) return
    
    if(hminim < hmin) then 
      print*, 'hminim', hminim; read(*,*) 
    endif 
    
    tp = ntp * tp 
    ynew(i,:) = (/tp, yi, cj/)
  ! write to file labeled by fpo passed from the main routine, PoInSt.dat  
    write(fpo,'(10d24.16)') tp, yi, cj, dcj, hminim ! fout - PoInSt.dat 
    
    if(debug==1)  then 
      print*, 'finish dfcr! g=', g
      read*
      call gr_cj(yi, cj2)
      
! check to see if there is any point to keep cj, dcj as the output! ckd--cj = cj2 anyway
      print*, 'ck, cj,dcj, cj2, cj+dcj-cj2',cj,dcj, cj2, cj+dcj-cj2 
      write(*,'(10e18.10)') tp, yi, cj, dcj, hminim !  !ck print to screen
      read(*,*)
    endif 
  
! using Adams predictor to get the new initial condition
    if(npo .eq. 1)  return
   
    call champ(i,g, dir, cham) 
    yctr = yi((/ctr1, ctr2, ctr3/))
    
    if (debug == 1) then  
      print*, 'champ! dir', dir !ck
      print*, 'yi, cham', yi, cham(i,:)
      print*, 'before adams, i, ds, yctr, cham', i, ds, yctr, cham
      read(*,*)   
    endif 

    call adams(i, yctr, ds, cham)
    yi((/ctr1, ctr2, ctr3/)) = yctr
!    call adams(i, yi((/ctr1, ctr2, ctr3/)), ds, cham) ! this doesn't work, so introduce another variable yctr

    if (debug == 1) then   
      print*, 'numerical continuation, new state!' !ck
      print*, yi
      print*
      read(*,*) 
    endif 
     
  enddo

  close(fout)

  return
  end subroutine pofam


!**********************************************************************
  subroutine dfcr(y0,imax, tf,yf,g, cj,dcj,hminim, ispo, deriv, gr_cj)
!  Using modified Newton method to refine the p.o, provided with the initial conditon.
!  with 3 free variables instead of 2 variables( inderterminant equations), least-square solution
!  The point is get the derivative of the  target variables with respect to the control variables

! for Newton method, alway check two things:
!   1. the error in target varaibles    --- stop condition, satisfied means finish the correction
!   2. the modulus of the correction    --- if too small, no point to continue ...

!  	Input Variables:
!  y0(6)	initial condition(x,y=0,z, vx=0, vy, vz=0)
!  imax 	the times of the crossing through the section(y=0)

!  	Output Variables:
!  y0, yf	the refined initial and final conditon
!  tf		half period or one period, depends on the number of intsec-imax
!  g(2,3)	JACOBIAN MATRIX OF F(*), 2 by 3,  3 variables and 2 equations
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
 
  implicit none  
  integer, parameter :: n = 42 ! both the state and the variational equations are desired, so n=42
  
  integer, intent(in)  :: imax ! number of crossings of the Poincare section
  integer, intent(out) :: ispo 
  real(kind=dp), intent(out) :: tf, yf(n), g(2,3), cj, dcj, hminim 
  real(kind=dp), intent(inout)  :: y0(6) 
  external :: deriv, gr_cj
  
! Local Varaibles
  integer :: iter, ispc 
  real(kind=dp)  :: f(2), phi(6,6), yi(n),pf(n), dx(3), dxm, fm, cj2
  
  if (debug == 1)   print*, 'start dfcr, y0!', y0   

  iter = 0 
  ispo = 1
  
  do   ! iteration to do the refinement
    iter = iter+1
    
    if(debug == 1) then 
      print*, iter, '-th iteration'   
      read(*,*)
    endif
      
!for Newton method, if the iteration is great than 5 or 6, then there is something wrong, provided that the initial guess is a good one, which is 
! close to the solution
     
    if(iter .ge. 20) then 
      print*, 'Error! Newton method iterates >  20, dxm, fm',  dxm,  fm
      if(fm .gt. tol) ispo = 0 ! tol for fm, and prsc for dxm
      ! consider about the appropriate error control
      return
    endif   
 
    call gr_cj(y0, cj)
    
!  subroutine fctn(x,init,imax, f,g,y,vf,tf, hminim, ispc, deriv, gr_cj)
    call fctn(y0, 0, imax, f, g, yf, pf, tf, hminim, ispc, deriv, gr_cj)  
    
    
    if ( ispc == 0 ) then
      ispo = 0 
      return 
    endif 
   
!   check the variation during one map here, instead of inside the fctn subroutine    
    call gr_cj(yf(1:6), cj2)
    dcj = cj2 - cj  

    if(debug == 1) then
      print*, 'ck cj', dcj, cj, cj2;  !ckd cj=cj2
      print*,'Finish fctn, g=', g
      read(*,*)
    endif 
    
    fm = dsqrt( f(1)*f(1) + f(2)*f(2)) ! error in target varaibles
 ! precision to terminate the correction, 1.d-16 from Gerard, which is of non-sense, because it beyonds the precision of double real
    if ( fm .le. tol)  then 
      print*, 'fm < tol', fm,  tol
!      read*
      return   ! try 1.d-11  >should be no less than the numerical precision we are using: e=1.d-13 (error in rk78)
    endif 
     
    call deltx(f,g, dx)
        
    ! update the initial state for the component as control variables --- check here!
    y0( ctr1 ) = y0( ctr1 ) + dx(1)
    y0( ctr2 ) = y0( ctr2 ) + dx(2)
    y0( ctr3 ) = y0( ctr3 ) + dx(3)
    
    ! the  precision is the modulus of the correction 
    dxm = dsqrt( dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3) )
    
    if(debug == 1) then 
      print*, 'dxm, fm, dx', dxm, fm, dx
      write(fdx, *)  dxm, fm, dx
      print*
      read*
    endif 
    
    if ( dxm .lt. prsc )  then 
      print*, '|dx|<', prsc, 'stop iteration!'
      if(debug == 1) read(*,*)
      
      return
    endif 

  enddo 
   
  return
  end subroutine dfcr


! ***************************************************************************
  subroutine fctn(x,init,imax, f,g,y,vf,tf, hminim, ispc, deriv, gr_cj)
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
!   ind, tar1, tar2, ctr1, ctr2, ctr3  

! period:  symmetric wrt xz plane, 1st crossing at T/2
!          symmetric wrt x-aixs : vertical lyapunov orbit, first crossing in T/4

!  subroutine used: poinc, deriv
!********************************************************************************
  implicit none  
  integer, intent(in)  :: init, imax
  integer, intent(out)  ::  ispc
  real(kind=dp), intent(in)  :: x(6)
  real(kind=dp), intent(out) :: f(2), g(2,3), vf(42), tf, hminim  
  real(kind=dp), intent(inout)  :: y(42)
  
! local varaibles
  real(kind=dp)  ::  phi(6,6), yi(42), f1, f2, & 
  		     phinm(6), phinormi, gnm(2), gnormi, cj
  integer        ::  i 
    
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
  
! the target variables
  f(1) = y(tar1)  ! the first target variale: vx if id = 2, y=0,ctr: x,z,vy
  f(2) = y(tar2)  ! the second target variale: vz
  
  call deriv(0.D0,y, 42, vf)  ! compute the vector field

! the derivative of f  divided by the ind-th component in velocity 
  f1 = -vf(tar1)/vf(ind)  !axf/ vyf   -- ind is the global parameter in module asym_pomod
  f2 = -vf(tar2)/vf(ind)  !azf/ vyf
 
! G is the Jacobi matrix of f(2 target variables) with respect to 3 control variables
! and Phi is the variatioal matrix computed as y(7-42)
  phi = reshape(y(7:42), (/6,6/)) ! better than equivalence declaration.... 
  
  ! there is problem here.....

  if (debug == 1) then 
    do i = 1, 6
      print*, phi(i,:)
    enddo 
    print*, ind,tar1, tar2, ctr1, ctr2, ctr3
    read*
  ! --- check the infinity norm of matrix phi and g
! DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
    phinormi = dlange('I', 6, 6, phi, 6, phinm)
    print*
    print*, 'Phi: norm_inf, norm_row' 
    print*, phinormi, phinm
    write(fphi,*) phinormi, phinm
  endif 
  
  ! to modify to general form ... remember to check every modification
  !  instead of using f1, f2, use f( tar(i) ) instead
  
  g(1,1) = phi(tar1,ctr1) + f1 * phi(ind,ctr1) 
  g(1,2) = phi(tar1,ctr2) + f1 * phi(ind,ctr2) 
  g(1,3) = phi(tar1,ctr3) + f1 * phi(ind,ctr3) 

  g(2,1) = phi(tar2,ctr1) + f2 * phi(ind,ctr1)  
  g(2,2) = phi(tar2,ctr2) + f2 * phi(ind,ctr2)  
  g(2,3) = phi(tar2,ctr3) + f2 * phi(ind,ctr3) 
  

  
  if ( debug == 1) then 
    print*, 'Fctn, g=', g
    read*
  endif 
  
  return
  end subroutine fctn

! **********************************************************************    
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
  
! A little trick to avoid the mis-detection of sign change at the start point, ie, g0 = -1.d-9, g1 = 1.d-9, that is not what we want  
  if(dabs(g) .lt. 1.d-9)  g = 0.d0   ! cancel numberical error, if the compoent of the state is less than 1.e-9, treat as 0 
 
! initial state
  h    = 1.d-3
  t    = 0.d0 ! start time
  hminim = h
  y = yi 
  
  if (debug == 1) then 
    print*, 'Poinc- hmin,hmax,e, h, t, tmax, imax', hmin,hmax,e, h, t, tmax, imax
    read(*,*)
  endif
  
      
  do i = 1, imax
    
    if(dabs(g) .lt. 1.d-9)  g = 0.d0   ! cancel numberical error, if the compoent of the state is less than 1.e-9, treat as 0 
     
    do  !look for the next intersecion across the Poincare section 
      gi = g ! the previous value of g
    
      call gr_rk78(t,y,neq, h,hmin,hmax,e, r,b,f, deriv )
      
      if (debug == 1)   write(*,'(7f10.6)') t, y(1:6)
      
      call sect(y, neq, g, dg)
      
!if it spends too much time to go to the next intersecion across y(ind)=sec plane, seen as failed 
      if (dabs(t) > dabs(tmax) ) then 
        ispc = 0
        print*, 'Maximum time exceeded!, t>tmax', t, tmax
        return
      endif
      
      hminim = dmin1(hminim, h) 
      if (gi*g < 0 ) exit ! terminate of different sign
      
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
  do  
    iter = iter + 1
    if (debug == 1) print*, iter,'-th iteration for intersection with poincare section'
    if (iter > 20 ) then 
     print*, 'Iteration exceeds 20!'; read*
     ispc = 0
     return
    endif 
    
    if(dabs(g) .le. 1.d-13)  exit  ! detect the y(ind)-sec is within the tolerance, values 1.d-13 from Gerard's routine
    
    call deriv(t, y, 42, f)

    dy = 0
    do i = 1,42
      dy = dy + f(i)*dg(i)
    enddo
    dh = - g/dy
    
!    print*, 'before gr_rk78, dh', dh
! to make sure this step works, we need to set hmin to the required stepsize 
! Need to debug here, for possible stuck here....

    call gr_rk78(t,y, neq, dh, dabs(dh),hmax,e, r,b,f, deriv) ! remember that gr_rk78 is only 1 step
    
!    call gr_rk78(t,y, neq, dh, dabs(dh), dabs(dh),e, r,b,f, deriv) ! remember that gr_rk78 is only 1 step
    
!    print*, 'after gr_rk78'
    call sect(y, neq, g, dg)
    
    if(debug == 1) then 
      print*, 'newton, t,g,dh', t, g, dh    
      read(*,*) !ck
    endif
      
  enddo 
 
! we get y=0, t is tp, yi is the intersection	
  yf = y
  tf = t
  
  if(debug == 1) write(*,*)'Poinc finished, g, tf,yf:',g, tf, yf(1:6)  ;  !read(*,*) !ck
  
  return
  end subroutine poinc 
 

! ***************************************************************************
  subroutine deltx(f,g, dx)
!  Nov.26,2014 -- general solution for the Modified Newton method with 3 free variables(control), 2 target equations
!  we have Xk = X(k-1) + deltX(k-1)
!  where we require G( deltX(k-1) ) = - F( X(k-1) )
!  It is better to just do the computation here, and check if we reach to the tolerance in the upper subroutine

!  minimize the norm
!     (deltX)T * Q * deltX
!      Q is a given diagonal positive definite weight matrix,  here,just set it to identity

!  so, the solution of this problem of minima is given by
!     (deltXk-1)= -Q-1 (G)T * (G Q-1 (G)T)-1 * F( X(k-1) ) 
!     (deltXk-1)= -(G)T* (G*(G)T)-1 * F( X(k-1) ) 

!    Input Variables
!  f(*)		target variables to be set 0, (here vx,vz )
!  g(*,*) 	Jacobian matrix of f(*) 

!    Output Variables
!  dx   	the difference for the initial state to meet the final target

!  subroutine used: none
! --------------------------------------------------------------------------
  
  implicit none 
  
  real(kind=dp), intent(in) :: f(2), g(2,3) 
  real(kind=dp), intent(out)  :: dx(3)

! Local Variables  
  integer :: i, j, &
  	     n, info ! for the norm computed by lapack
  real(kind=dp) :: g2(2,2), b(2), det, gnew(3,2), gnm(3), gnormi !, & 
!                   work(4*2), iwork(2), rcond ! for the norm computed by lapack

! g2 = g*g^T
  do  i = 1,2
    do  j = 1,2
      g2(i,j) = g(i,1)*g(j,1)+ g(i,2)*g(j,2)+ g(i,3)*g(j,3)
    enddo
  enddo
  

! In order to compute the norm of G^T*(G*G^T)^(-1), we use matrix to represent it 
  det = g2(1,1)*g2(2,2) - g2(1,2)*g2(2,1)
  
  gnew(1,1) =  g(1,1) * g2(2,2) - g(2,1) * g2(2,1) 
  gnew(1,2) = -g(1,1) * g2(1,2) + g(2,1) * g2(1,1)  
   
  gnew(2,1) =  g(1,2) * g2(2,2) - g(2,2) * g2(2,1) 
  gnew(2,2) = -g(1,2) * g2(1,2) + g(2,2) * g2(1,1)  
  
  gnew(3,1) =  g(1,3) * g2(2,2) - g(2,3) * g2(2,1) 
  gnew(3,2) = -g(1,3) * g2(1,2) + g(2,3) * g2(1,1)  
  
  gnew = gnew/det 
  
  do i = 1, 3
    dx(i) = - ( gnew(i,1) * f(1) + gnew(i,2) * f(2) )
  enddo 
  
  if(debug == 1) then ! check the matrix norm of g^T(g*g^T)^(-1)
    gnormi = dlange('I', 3, 2, gnew, 3, gnm)
    print*, 'g^T(g*g^T)^(-1), : norm_inf,  norm_row(3)'
    print*,  gnormi, gnm
    write(fgnew,*) gnormi, gnm
    print*
    read(*,*)
   
    print*, 'dx', dx  ! ckd, exactly the same
  
    b(1) =  g2(2,2)*f(1) - g2(1,2)*f(2)
    b(2) = -g2(2,1)*f(1) + g2(1,1)*f(2)
 
    do  i = 1,3
      dx(i) = -(g(1,i)*b(1) + g(2,i)*b(2))/det
    enddo 
    print*, 'dx', dx;  read(*,*) 
  endif 
  
  return      
  end subroutine deltx  



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


!***********************************************************************
  SUBROUTINE ADAMS(NP,X,HH,CAM)
!     adams-bashforth formulas

!          input parameters:
!     np         number of computed periodic orbits
!     x(*)       initial conditions of the last p.o.
!     hh         integration step along arc length of solution locus
!     cam(*,*)   (see subroutine camp)

!          output parameters:
!     x(*)       new aproximated initial conditions for p.o.

!  Revised by Yu, 20160219
! --------------------------------------------------------------------------
 
  implicit none 
 
  integer, intent(in) :: np
  real(kind=dp), intent(in) :: hh, cam(4,3) 
  real(kind=dp), intent(inout) :: x(3) !update x(3)
  
! Local Variables 
  integer :: i
  
  do i = 1, 3  
    if (np .ge. 4) then 
      x(i) = x(i)+ hh*(55.*cam(4,i)-59.*cam(3,i)+37.*cam(2,i)-9.*cam(1,i))/24.

    elseif(np == 1) then 
      x(i)=x(i) + hh*cam(1,i)
    
    elseif(np == 2) then 
      x(i)=x(i)+.5*hh*(3.*cam(2,i)-cam(1,i))
    
    elseif(np == 3) then 
      x(i)=x(i)+hh*(23.*cam(3,i)-16.*cam(2,i)+5.*cam(1,i))/12
      
    endif 
  enddo
 
  return
  end subroutine adams

  
!***********************************************************************  
  subroutine champ(np,g,dir, cham)
!     COMPUTATION OF THE VECTOR FIELD GIVING THE CHARACTERISTIC CURVE
!     OF THE FAMILY
!	  INPUT PARAMETERS:
!     NP	     NUMBER OF THE LAST COMPUTED P.O.
!     G(*,*)         JACOBIAN MATRIX OF F(*) (SEE SUBROUTINE TRACA)
!     dir	     SENSE ON THE CHARACTERISTIC CURVE
!	   
! 	  OUTPUT PARAMETERS:
!     cham(*,*)  VECTOR FIELD ON THE CHARACTERISTIC CURVE AT THE LAST 4
!		     POINT, ONLY THE LAST ROW OF cham(*,*) IS COMPUTED
!		     IN THIS SUBROUTINE
! --------------------------------------------------------------------------
 
  implicit none 
 
  integer, intent(in) :: np 
  integer, intent(inout) ::  dir
  real(kind=dp), intent(in) :: g(2,3) 
  real(kind=dp), intent(out) :: cham(4,3)  
  
! Local Variables 
  integer :: i, nnp
  real(kind=dp) :: a(3), sm, c
 
  a(1) =  g(1,2)*g(2,3) - g(1,3)*g(2,2)
  a(2) = -g(1,1)*g(2,3) + g(1,3)*g(2,1)
  a(3) =  g(1,1)*g(2,2) - g(1,2)*g(2,1)
 
  sm = dir*dsqrt(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))

  if(np .gt. 4)  then 
  
    cham(1:3, :) = cham(2:4, :) !checked, assignment ok!
    cham(4, :) = a/sm
    nnp = 4
    
  else 
  
    cham(np, :) = a/sm 
    if (np .eq. 1)  return
    nnp = np
  endif 
 
  c = 0.d0
  do  i = 1,3
    c = c + cham(nnp,i)*cham(nnp-1,i)
  enddo

  if(c .gt. 0.d0) return
  
  !  this part for the reverse of the sense, detect bifurcation,-- need to be understood later
  dir = -dir
  cham(nnp, :) = -cham(nnp, :)
  write(*,*)'reversal in the sense along the characteristic curve'
  return
  end subroutine champ

end module asym_pomod

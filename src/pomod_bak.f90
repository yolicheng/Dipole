module pomod
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
!  3. plob 	plot the orbit 				gr_rk78, deriv, gr_cj 
!         general routine to integrate the orbit in time interval [t0, t0+tf] 

!  4. dfcr   						poinc, deriv, gr_cj , fctn, deltx
!  5. fctn  						poinc, deriv, gr_cj 
!  6. poinc 						sect, gr_rk78, deriv, gr_cj 
!  7. champ, adams, deltx, sect, deriv,gr_cj 	 	none ---- basic subroutines  

! Note: gr_rk78 is outside this module, that is a general integration routine

!   **********  Subroutine  Declaration ************ 
!  subroutine init_po(sym, sec0, hmin, hmax, e, tmax, tol, prsc)
!  subroutine pofam(yi,npo,imax, dir,ds, fpo, ynew, i, deriv, gr_cj)
!  subroutine plob(y0,t0,tf,tdir,ftag, deriv, gr_cj,  y) 

!write(fpo,'(10d24.16)') tp, yi, cj, dcj, hminim ! PoInSt.dat 
!write(ftag,'(7e20.10)') t, y, cj ! po.dat
!    ynew(i,:) = (/tp, yi, cj/)
       
!  subroutine dfcr(y0,imax, tf,yf,g, cj,dcj,hminim, ispo, deriv, gr_cj)
!  subroutine fctn(x,init,imax, f,g,y,vf,tf, hminim, ispc, deriv)
!  subroutine poinc(yi,imax, tf,yf, hminim, ispc, deriv) 
!  subroutine sect(y, g,dg)   
!  subroutine deriv(t,x,n,f)
!  subroutine gr_cj(x, cj)

!   **********   Public  variables   **********  
!    Assigned by subroutine init_lf before call any function from this module 
!  tari, (i=1,2)	target varaibles
!  ctri, (i=1,2,3) 	control varaibles
!  sec, ind	 	the value to specify the Poincare section y(ind) = sec 
!  ntp 			the real period is TP = ntp * tp(to go to the first crossing)
! ************************************************************************************************


implicit none
integer, parameter, private :: dp = kind(1.d0) !to avoid compliance with dp defined in other modules
!( this pofam is supposed for general use, self-contained)

! ******************** Declare Public *****************
! public by default, the most commonly used varaibles which is also accessible by any subroutine that uses this module
real(kind=dp) :: mu  
integer :: ind, ntp, &  !
           tar1, tar2, ctr1, ctr2, ctr3  ! index of the control and target variables


! Declare external for the routines or functions to be called the routines from library lapack  
external dgecon ! subroutine
external gr_rk78

real(kind=dp), external :: dlange ! dlange is a function, declare the type

 ! the files to check the iteration, matrix norm --discard, mostly caused by numerical error
integer, private :: fphi, fgnew, fdx 
  
!  private - only shared by the subroutines, cannot be accessable from outside the module 
!  can avoid passing the shared parameters level-by-level
integer, private      ::  debug
!real(kind=dp), public ::  sec, hmin, hmax, e, & ! Bound values for rk78, for the integration of orbit
!			        tmax, tol, prsc ! tolerance control for the termination of Poincare map and Newton method


real(kind=dp), private, save ::  sec, hmin, hmax, e, & ! Bound values for rk78, for the integration of orbit
			   tmax, tol, prsc ! tolerance control for the termination of Poincare map and Newton method


contains 

!******************** Module-related Parameters ******************************** 
  subroutine init_po(symt, sec0, hmin0, hmax0, e0, tmax0, tol0, prsc0)
! Initialize the index of the control variables and target variables as in the state vector
! The Poincare section is specified as:  y(ind) = sec
! question: symt = ind? NO! for vertical lyapunov orbit, we take y=0 as Poincare section, but use symmetry 2, so cannot join as one parameter

!It is ok to just pass the parameter as input, without assignment by h=h0

  integer, intent(in) :: symt ! index of symmetry, used to assign the public values accordingly
  real(kind=dp), intent(in) :: sec0, & !poincare section
  			       hmax0, hmin0, e0, tmax0, tol0, prsc0 ! Bound values for RK78, shared by all subroutines that calls rk78 (poinc, plob, )
  
        
  debug = 1  ! for debug
!  debug = 0 ! normal execution 

  ! we have to initialze here by assignments, just pass in won't work, neither do declaration of inout attributes 
  sec = sec0
  hmin = hmin0
  hmax = hmax0
  e = e0
  tmax = tmax0 
  tol = tol0 
  prsc = prsc0
!  
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
     
!  mu = 0.12150586d-1   !mu = 0.012150585609624 --- just to test pomod with halo
! note: to use these symmetries, the invariant transformation in time should be t-->-t

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
!    4:x-axis	y=0, xz-plane  		 1,3,3+2  (x, z, vy)  	  T/4    3, 4 (z, vx) (perpendicular to x-axis)
! to be added the other two

! in fact, now we only need to deal with the case 2,  test vt_lyap later----to be modified
! ! symmetry:  index of control var   of target var
! 5: y=0, xz-plane,  1,3,3+2 (x, z, vy)     3+1, 3+3 (vx, vz)      
! 6: z=0, xy-plane,  1,2,3+3 (x, y, vz)     3+1, 3+2 (vx, vy)   
   
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
  elseif (symt == 4)  then 
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
  end subroutine init_po
  
  
!***********************************************************************  
  subroutine pofam(yi,npo,imax, dir,ds, fpo, ynew, i, deriv, gr_cj)
! Main subroutine to compute np new p.o.s with one initial guess yi(6)
! To make sure the initial state is periodic, refine it without check
!  the initial condition for all p.o.s are stored in file fpo

! 	Input Varialbles
!   yi(6)	initial condition (ie, x,y=0,z, vx=0, vy, vz=0)
!   npo		number of new p.o.  if np=0, just do differential correction
!   imax	number of intersections with Poincare section to get p.o.

!   dir		dirction along the vector field, +1:increase,-1:decrease
!   ds		displacement for new PO
!   fpo		file tag for the initial states of the p.o.s

! 	Output Variables
!   ynew	the initial condition for new p.o.s, np*8
!     		save the initial p.o in the first row
!   i 		the number of available p.o.s (i-1)

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
   
!  subroutine dfcr(y0,imax,tol, prsc, tf,yf,g, cj,dcj,hminim, ispc, deriv, gr_cj)
    call dfcr(yi,imax, tp,yf, g, cj, dcj, hminim, ispo, deriv,gr_cj) 
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
      call gr_cj(yi, cj2)
! check to see if there is any point to keep cj, dcj as the output! ckd--cj = cj2 anyway
      print*, 'ck, cj,dcj, cj2, cj+dcj-cj2',cj,dcj, cj2, cj+dcj-cj2 
      write(*,'(10e18.10)') tp, yi, cj, dcj, hminim !  !ck print to screen
      read(*,*)
    endif 
  
! using Adams predictor to get the new initial condition
    if(npo .eq. 1) return
   
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
!    call adams(i, yi((/ctr1, ctr2, ctr3/)), ds, cham)

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
  

!  function used: poinc, deriv, fctn, deltx
! --------------------------------------------------------------------------
 
  implicit none  
  integer, parameter :: n = 42 ! both the state and the variational equations are desired, so n=42
  
  integer, intent(in)  :: imax ! number of crossings of the Poincare section
  integer, intent(out) :: ispo 
  real(kind=dp), intent(out) :: tf, yf(n), g(2,3), cj, dcj, hminim 
  real(kind=dp), intent(inout)  :: y0(6) 
  external :: deriv, gr_cj
  
! local varaibles
  integer :: iter, ispc 
  real(kind=dp)  :: f(2), phi(6,6), yi(n),pf(n), dx(3), dxm, fm, cj2
  
  if (debug == 1)   print*, 'start dfcr, y0!',y0   

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
      if(fm .gt. 1.d-5) ispo = 0
      
!      ispo = 0 ! here is a problem, if it comes to here, that means the there is an intersection, but it is difficult to get the required precision.
!      read(*,*)
      return
    endif   
 
    call gr_cj(y0, cj)
!  subroutine fctn(x,init,imax, f,g,y,vf,tf, hminim, ispc, deriv)
    call fctn(y0, 0, imax, f, g, yf, pf, tf, hminim, ispc, deriv)  
    
    if ( ispc == 0 ) then
      ispo = 0 
      return 
    endif 
   
!   check the variation during one map here, instead of inside the fctn subroutine    
    call gr_cj(yf(1:6), cj2)
    dcj = cj2 - cj  

    if(debug == 1) then
      print*, 'ck cj', dcj, cj, cj2;  !ckd cj=cj2
      read(*,*)
    endif 
    
    fm = dsqrt( f(1)*f(1) + f(2)*f(2)) ! error in target varaibles
 ! precision to terminate the correction, 1.d-16 from Gerard, which is of non-sense, because it beyonds the precision of double real
    if ( fm .le. tol)  return   ! try 1.d-11  < e=1.d-13 (error in rk78)
    
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
      read(*,*)
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
  subroutine fctn(x,init,imax, f,g,y,vf,tf, hminim, ispc, deriv)
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
  		     phinm(6), phinormi, gnm(2), gnormi  
  
  external deriv ! the vector field :  derive(t,x,n,f)
  
  
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

!subroutine poinc(yi,imax, tf,yf, hminim, ispc, deriv) 
  call poinc(yi,imax, tf,y, hminim, ispc, deriv) ! imax determines the number of crossing we are considering as
  
  if (ispc == 0) return !fail to reach to Poincare section, stop  the correction 
  
! the target variables
  f(1) = y(tar1)  ! the first target variale: vx if id = 2, y=0,ctr: x,z,vy
  f(2) = y(tar2)  ! the second target variale: vz
  
  call deriv(0.D0,y, 42, vf)  ! compute the vector field

! the derivative of f  divided by the ind-th component in velocity 
  f1 = -vf(tar1)/vf(ind)  !axf/ vyf   -- ind is the global parameter in module pomod
  f2 = -vf(tar2)/vf(ind)  !azf/ vyf
 
! G is the Jacobi matrix of f(2 target variables) with respect to 3 control variables
! and Phi is the variatioal matrix computed as y(7-42)
  phi = reshape(y(7:42), (/6,6/)) ! better than equivalence declaration.... 


  if (debug == 1) then 
  ! --- check the infinity norm of matrix phi and g
! DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
    phinormi = dlange('I', 6, 6, phi, 6, phinm)
    print*
    print*, 'Phi: norm_inf, norm_row' 
    print*, phinormi, phinm
    write(fphi,*) phinormi, phinm
  endif 
  
  g(1,1) = phi(tar1,ctr1) + f1 * phi(ind,ctr1) 
  g(1,2) = phi(tar1,ctr2) + f1 * phi(ind,ctr2) 
  g(1,3) = phi(tar1,ctr3) + f1 * phi(ind,ctr3) 

  g(2,1) = phi(tar2,ctr1) + f2 * phi(ind,ctr1)  
  g(2,2) = phi(tar2,ctr2) + f2 * phi(ind,ctr2)  
  g(2,3) = phi(tar2,ctr3) + f2 * phi(ind,ctr3) 
  
  gnormi = dlange('I', 2, 3, g, 2, gnm)
  
!  print*
!  print*, 'g: norm_inf, norm_row(2)'
!  print*,  gnormi, gnm
!  write(fg,*) gnormi, gnm
!  print*

  
  return
  end subroutine fctn

! **********************************************************************    
  subroutine poinc(yi,imax, tf,yf, hminim, ispc, deriv) 
! Determination of the imax-th passage of the orbit through the Poincare section defined by subroutine sect(here,y=0) 
!  Note:  only one intersection is record here!
! to make this subroutine more general, use deriv as an argument as  the function to compute the vector field
!   of the common form:  deriv(T,X,N,F)
     
! 	Input Variables
!  yi(*)	initial point 
!  imax		the number of time for intersecting the section

! 	Private Module-based Varaible
!  tmax		maximum time for half revolution, if exceeds, fail to obtain the poincare map -by module

!  	Output Variables
!  tf  	 	time spent by the orbit to go from yi to yf
!  yf(*)  	first cut of the orbit which passes by yi with the surface of section
!  hminim	minimum step used in the integration of the p.o.
!  ispc 	flag of the success to return to the poincare section

! function used: gr_rk78, sect  
! --------------------------------------------------------------------------

  implicit none  
  integer, parameter  ::  neq = 42 ! dimension of the state: 6-state vector, 42-also the variatioal matrix
  
  integer, intent(in)  ::  imax
  integer, intent(out)  ::  ispc
  real(kind=dp), intent(in)  :: yi(neq) 
  real(kind=dp), intent(out) :: tf, yf(neq), hminim
  external deriv  
  
! local varaibles
  integer :: i 
  real(kind=dp)  ::  g, gi,  dg(neq), t, dh, dy, &
  	  	     r(13,neq),b(neq),f(neq), h !, hmin, hmax, e !gr_rk78 
  
  ispc = 1 ! default value is 1, if fails, set to 0
  call sect(yi, g, dg)
  
! A little trick to avoid the mis-detection of sign change at the start point, ie, g0 = -1.d-9, g1 = 1.d-9, that is not what we want  
  if(dabs(g) .lt. 1.d-9)  g = 0.d0   ! cancel numberical error, if the compoent of the state is less than 1.e-9, treat as 0 
 
! initial state
  h    = 1.d-3
  t    = 0.d0 ! start time
  hminim = h
  
  if (debug == 1) then 
    print*, 'Poinc- hmin,hmax,e, h, t, tmax, imax', hmin,hmax,e, h, t, tmax, imax
    read(*,*)
  endif
  
  print*, 'tmax',  tmax 
  read(*,*)
  
  do i = 1, imax
    do  !look for the next intersecion across the Poincare section 
      gi = g ! the previous value of g
    
      print*, 'before rk78- tmax, hmin,hmax,e,ind,sec',tmax, hmin,hmax,e, ind,sec
      print*, 'tar1, tar2, ctr1, ctr2, ctr3 ', tar1, tar2, ctr1, ctr2, ctr3 
      read(*,*)
!SUBROUTINE GR_RK78 (X,Y,N,H,HMI,HMAX,E1,R,B,F,DERIV)
      call gr_rk78(t,yi,neq, h,hmin,hmax,e, r,b,f, deriv)
      
      print*, 'after rk78- tmax, hmin,hmax,e,ind,sec',tmax, hmin,hmax,e, ind,sec
            print*, 'tar1, tar2, ctr1, ctr2, ctr3 ', tar1, tar2, ctr1, ctr2, ctr3 
      read(*,*)
      
      write(*,'(7f10.6)') t, yi(1:6)
      if (debug == 1)   write(*,'(7f10.6)') t, yi(1:6)
      
      call sect(yi,g, dg)
      
      print*, 'tmax',  tmax 
      read(*,*)
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

  do 
    if(dabs(g) .le. 1.d-13)  exit  ! detect the y(ind)-sec is within the tolerance, values 1.d-13 from Gerard's routine
    
    call deriv(t, yi, 42, f)

    dy = 0
    do i = 1,42
      dy = dy + f(i)*dg(i)
    enddo
    dh = - g/dy
    
! to make sure this step works, we need to set hmin to the required stepsize 
    call gr_rk78(t,yi, neq, dh, dabs(dh),hmax,e, r,b,f, deriv) ! remember that gr_rk78 is only 1 step
    call sect(yi, g, dg)
    
    if(debug == 1) then 
      print*, 'newton, t,g,dh', t, g, dh    
      read(*,*) !ck
    endif
      
  enddo 
 
! we get y=0, t is tp, yi is the intersection	
  yf = yi
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
  subroutine sect(y, g, dg)
  ! the surface of section, defined by g =  y(ind) - y0
! 	Input parameters:
! y(*)      	the state vector 

!  	Output parameters:
! g 		funciton that equated to 0 gives the surface of section
! dg(*) 	gradient of function g


!  Global Variables from module  
! ind  		the index of the component of y to be used as Poincare section  
! sec 		the Poincare section to be specified as y(ind) - sec

! 20160218 -by Yu
! --------------------------------------------------------------------------

  implicit none 
  
  real(kind=dp), intent(in) :: y(42) 
  real(kind=dp), intent(out) :: g, dg(42)

  integer :: i

  g = y(ind) - sec
   
  print*, 'Sect- tmax, ind, sec',  tmax,ind,sec
!  print*, 'Sect- tmax, h, hmin, hmax, e',  tmax, h, hmin, hmax, e 
  read(*,*)   
    
  do  i = 1,42
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


!***********************************************************************    
  subroutine plob(y0,t0,tf, tdir, ftag, deriv, gr_cj, y) 
!  this subroutine is to plot a segment of trajectory, with the intial state y0(6), 
!  and time span for integration is [0, tf], always save the energy or jacobi constant in the eight-th column to check the numerical error.
!  And save all the date in file ftag, with 2 blank lines between 2 different orbits

! Note: for backward in time integration, the stepsize h should be negative, and because h=1.d-3 is assigned inside the subroutine
!       we introduce  tdir to control the integration direction    :  1:unstable; -1:stable  
!       h = tdir * h   
   
!  better to do the write loop outside the subroutine--- 
!  in case to change the name of the data file, and write during the do loop

!  Input variables  
!    y0(6) :  inital state of trajectory 
!    t0    :  start time for the integration
!    tf    :  end time for the integration 
!    tdir  :  control the integration direction :  1: forward;  -1: backward  
!    ftag  :  data file to save the result of integration, 1:t, 2-7:y, 8:cj 
!    deriv :  external subroutine to compute the vector field
!   gr_cj  :  external subroutine to compute the conservative quantity: energy or Jacobi Constant

! Output Variables
!    y(6)  :  finial state of trajectory 
  
!  ROUTINE USED: GR_RK78     GR_RTBP(A,B,N,F)
!  revised by Yu  -- 20160219
! --------------------------------------------------------------------------

  implicit none 

  integer, intent(in) :: ftag, tdir   
  real(kind=dp), intent(in) :: y0(6),t0, tf
  real(kind=dp), intent(out) :: y(6) ! As the final state
  external :: deriv, gr_cj ! vector field  && energy or Jacobi constant

! Local Variables
  real(kind=dp) :: r(13,6),b(6),f(6), t, h,  & ! rk78  -- hmin, hmax, e(by module) 
		   cj 

  y = y0    
 
! Integrate the periodic orbits by gr_rk78, write the x, y in txt file with ftag
  t =  t0 !0.d0
  h =  tdir*1.d-3 ! maybe small value is better, for a p.o. with small period, the orbit could be quite corse 

! if the evaluation condition is (t+h .lt. tf), only consider the positive step size, integate forward
! ( dabs(t+h-t0) .lt. dabs(tf-t0) ), both positive and negative cases are considered.

  if(debug == 1) then 
    write(*,*) 'before plob, t0, tf, y0', t0, tf, y0  !ck
    print*, 'hmin, hmax, e, tmax', hmin, hmax, e, tmax; read(*,*)
  endif 

  do while( dabs(t+h-t0) .lt. dabs(tf-t0) ) 
    call gr_cj(y(1:6), cj)
    write(ftag,'(7e20.10)') t, y, cj
    call gr_rk78(t,y,6,h,hmin,hmax,e,r,b,f,  deriv)
    if(debug == 1) print*, 'h, hmin, hmax', h, hmin, hmax !ck
  enddo
     
!if you want the strict periodic orbit, a control of tf must be made
  if(dabs(t-tf) .gt. 1.d-9) then
    h =  tf - t
 ! to make sure the next step can be executed within allowable interval [hmin hmax]
!    hmin should be dabs(h)
    call gr_rk78(t,y,6,h, dabs(h) ,hmax,e,r,b,f, deriv)
  endif 
        
!  write(ftag,'(7e20.10)') t, y
  
  write(ftag,'(7e20.10)') t, y

  if(debug == 1) then 
    write(*,*) 'Finish plob, t0, t, tf, yf', t0, t, tf, y !ck
    read(*,*)   !ck
  endif 
  
  write(ftag,*)  ! better to save as a block than index
  
  end subroutine plob


end module pomod

! 20160406

! For the asymmetric approach, please note that because we are refine w.r.t the initial condition
!   So the Jacobi Matrix w.r.t the correction on the initial state is actually Phi - I  

! Instead of Poincare map approach, we use the variational matrix to do the refinement of the asymmetric p.o., such that we have 
! 5 equations and 6 unkonwns(period included). We are not using the modified Jacobi Matrix in which the differential of the Poincare map w.r.t time t is represent by the other unkonwns in order to decrease the number of unknowns by 1 

! the 7 unkowns are (xi, yi, zi, vxi, vyi, vzi, T) 
! and the 6 constraint equations are satisfied to return to the same initial point 
!  (xfi-xi=0, yf-yi=0, zf-zi=0, vxf-vxi=0, vyf-vyi=0, vzf-vzi=0)

module  asymt_pomod

implicit none
save

integer, parameter, private :: dp = kind(1.d0) !to avoid compliance with dp defined in other modules
!( this pofam is supposed to be  self-contained, and for general use,

! ******************** Declare Public *****************
! public by default, the most commonly used varaibles which are also accessible by any subroutine that uses this module
!real(kind=dp), parameter ::  mu = 0.12150586d-1   !mu = 0.012150585609624  ! this is just to check halo orbit in other polar orbit exploration

integer, private :: debug, dsctr   ! do the assignment from main routine, to aviod the compilation of link before every excution
 
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
real(kind=dp), private ::  hmin, hmax, e, & ! Bound values for rk78, for the integration of orbit
			   tmax, tol, prsc ! tolerance control for the termination of Poincare map and Newton method

contains 


!******************************************************************************* 
! 	Input 
!  y0 	initial state

! 	Output
!  y0 	 	the updated state
!  tf		elapsed time to reach to the next mirror configuration,  which is 1 period for asymmetric p.o. 
!  yf		final state at the epoch of the mirror configuration
!  g 		Jacobi Matrix of the F(the periodicity equation ), will be used for the continuation
!  cj		Jacobi constant
! dcj		difference of the Jacobi constant between the initial and final point along the p.o.
!hminim		the minimum step size during the integration, a reference quantity
! ispo		flag to indicate if the correction succeeds(1) or fails(0)
! niter		number of iteration for the Newton Method, used as criteria to adjust the stepsize for continuation
! deriv		the subroutine to compute vector field of the equation of motion
! gr_cj		the subroutine to compute the Jacobi Constant(or energy)
! 
  subroutine dfcr_asymt(y0, tf,yf,g, cj,dcj,hminim, ispo, niter, deriv, gr_cj)
  
  use po_mod, only : debug
  
  implicit none  
  integer, parameter :: n = 42 ! both the state and the variational equations are desired, so n=42
  
  
  ! number of control variables and target varibales
  integer, parameter :: ntar = 6, nctr = 7
  
  
  integer, intent(out) :: ispo, niter
  real(kind=dp), intent(out) :: tf, yf(n), g(ntar, nctr), cj, dcj, hminim 
  real(kind=dp), intent(inout)  :: y0(6) 
  external :: deriv, gr_cj
  
  
! Local Varaibles
  integer        ::    i, j
  real(kind=dp)  :: f(ntar), phi(6,6), yi(n),pf(n), dx(nctr), dxm, fm, cj2, &
                    yit(6), y0copy(6), &
                    tol, prsc 
                    
  logical :: ok  ! status of the file poit.dat, opened or not, which saves the orbit of all the iterations                  
  
  
  tol = 1.d-10
  prsc = 1.d-10 
  
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

    if (niter > 5) then 
      print*, 'Newton method iterates > 5, dxm, fm',  dxm,  fm
      if(fm .gt. tol) ispo = 0 ! tol for fm, and prsc for dxm
      read*
      return
    endif   
 
    call gr_cj(y0, cj)
    
!  subroutine fctn(x,init, f,g,y,vf,tf, hminim, ispc, deriv, gr_cj)
    call fctn_asymt(y0, 0,  f, g, yf, pf, tf, hminim, ispc, deriv, gr_cj)  
    
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
    
    ! compute the correction 
    call deltx(f,g,  dx)
        
    ! the  precision is the modulus of the correction, or dxm = dnrm2(nctr, dx, 1)
    dxm = dnrm2(nctr, dx, 1)
 
 
    ! if the correction is below the precision, regard as success of the refinement
    if ( dxm .lt. prsc )  then 
      print*, '|dx|<', prsc, 'stop iteration!'
      if(debug == 1)   read(*,*) ! in debug mode, pause to check the result 
      return
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
  end subroutine dfcr_asymt
  
  
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
  subroutine fctn_asymt(x,init, f,g,y,vf,tf, hminim, ispc, deriv, gr_cj)
  
  use po_mod, only : debug
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
 
! G = Phi - vf( tar(i) ) / vf(ind) * Phi, the coefficient matrix of target error to obtain correction
  do i = 1, ntar 
    do j = 1, nctr 
      g(i,j) = phi( tar(i), ctr(j) ) - phi(ind, ctr(j) ) * vf( tar(i) ) / vf(ind) 
    enddo 
  enddo  
 
 ! -----------------------------------------------------
 ! check if the computation of G above is trustable, use 1-3-5  versus 4-6 ---- checked! fine!!!!
! ! 
!  print*, 'Original Jacobi Matrix'
!  do i = 1, ntar
!    print*, g(i,:)
!  enddo 
!  
!  
!  g(1,1) = phi(4,1) - vf(4)/vf(2) * phi(2,1) 
!  g(1,2) = phi(4,3) - vf(4)/vf(2) * phi(2,3) 
!  g(1,3) = phi(4,5) - vf(4)/vf(2) * phi(2,5) 

!  g(2,1) = phi(6,1) - vf(6)/vf(2) * phi(2,1)  
!  g(2,2) = phi(6,3) - vf(6)/vf(2) * phi(2,3)  
!  g(2,3) = phi(6,5) - vf(6)/vf(2) * phi(2,5) 
!  
!  print*, 'Jacobi Matrix by step-by-step computation'
!  do i = 1, ntar
!    print*, g(i,:)
!  enddo 
! 
! read*
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
  

  





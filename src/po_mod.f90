module po_mod 

! ** Idea for computation of symmetric po**   
! --1. use poinc to compute the intersecion with the section y(ind) = p0, 
! --2. use modified Newton method to refine p.o.
!  defaultly we deal with 2 equations(target varaibles) with three unknowns(control variables)


! For the continuation:  
!   Adams predictor + champ(vector field update) 

! 2016-08-07 17:17:28 
!   Add the computation of the differential of the Poincare map 

! 20160411  --- asymmetric periodic orbits 
!  Combine the two approaches: Poincare map + General into this routine 

!    we have 2 possible approaches to deal with the asymmetric p.o. 
!    1. fix one component of the state vector, use this as the poincare section, and call poinc_pomod 
!    2. do the refinement w.r.t X0(6)+TP, altogether 7 unknowns, and try to minimize the difference 
!       between the initial condition and the final condtion, call gr_pomod 

! 20160407 ** NOTE ** 
! For the asymmetric approach, please note that because we are refine w.r.t the initial condition
! So the Jacobi Matrix w.r.t the correction on the initial state is actually Phi - I  


! 3- Error and termination control   -- from Carles' notes 
!    Newton method is so efficient that, after 5 iterations, if it is still not convergent, 
!    we need to check what happens 

!    The tolerance (error in target variables) and presion (of the control variables) is to be careful assigned 

!    -- From Gerard's book, P96, 
!    the percision required for the modified Mewthon's method has been taken equal to 1.d-11 
!    and the bound for local errors in RK78 routine was set to 1.d-13

!    Error in fm = dsqrt( f(1)*f(1) + f(2)*f(2)) ! Gerard uses tol = 1.d-16, quite confusing, should be greater than 1.d-14  

!    tol = 1.d-13; prsc =1.d-11 ! the suggested values
!    but for lf problem, we take  tol = 1.d-11; prsc = 1.d-11  

!    because the error control for rk78 is 1.d-13, all the poincare map-related data has no smaller precision
!    so we cannot ask more precision in Newton method

! 4- Comments: do not use the obsolete syntax goto in f90 

! 5- Fortran has 15 significant digits for  double real numbers 

! 6- Termination control for the poincare map 
!    treat as lose of convergence if it spends too long for the next interation with the poincare section

! 7- -- discard this item!! Check the infinity norm of matrix phi and g - from Lapack!!!!! use LLS solver instead
!     DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
!     phinormi = dlange('I', 6, 6, phi, 6, phinm)
!     subroutine dgecon(norm,n,a,lda,anorm,rcond,work,iwork,info) -- work(4*n), iwork(n)
!     call dgecon('I',6, phi, 6, phinormi,rcond, work,iwork,info)

!  	Contains Subroutines  		  		function used:
!  1. init_po	initialize public varaibles	   	none
!  2. pofam  	main subroutine				dfcr, champ, adams
!                 general routine to integrate the orbit in time interval [t0, t0+tf] 

!  3. dfcr   						poinc, deriv, gr_cj , fctn, deltx
!  4. fctn  						poinc, deriv, gr_cj 
!  5. poinc 						sect, gr_rk78, deriv, gr_cj 
!  6. champ, adams, deltx, sect, deriv,gr_cj 	 	none ---- basic subroutines  

! Note: gr_rk78 is outside this module, that is a general integration routine

!   **********  Subroutine  Declaration ************ 

! ----------------  Module Variable Initialization-------------------
!  subroutine init_errctr(tol0, prsc0)
!  subroutine init_sympo(symt )   ! for symmetric p.o.

!  subroutine pofam(yi,npo,imax, dir,ds, fpo, ynew, i, deriv, gr_cj)
!        write(fpoinst,'(10E24.14)') tp, yi, cj, dcj ! PoInSt.dat 
!        ynew(i,:) = (/tp, yi, cj/)
       
!  subroutine dfcr(y0,imax, tf,yf,g, cj,dcj, ispo, deriv, gr_cj)
!  subroutine fctn(x,init,imax, f,g,y,vf,tf, ispc, deriv, gr_cj)
 

!   *****************   Public  variables   ***************************** 
!  p0, ind              the value to specify the Poincare section y(ind) = p0
!  ntp                  the real period is TP = ntp * tp(to go to the first crossing)
!  tari, (i=1,2)        target varaibles
!  ctri, (i=1,2,3)      control varaibles
! ************************************************************************************** 

use dp_mod
implicit none
save

! ******************** Declare Public *****************
! public by default, the most commonly used varaibles which are also accessible by any subroutine that uses this module
!real(kind=dp), parameter ::  mu = 0.12150586d-1   !mu = 0.012150585609624  ! this is just to check halo orbit in other polar orbit exploration

integer, private :: debug, dsctr, anglectr, isarc, issymt, smh, con_stop  ! do the assignment from main routine, to aviod the compilation of link before every excution
 
integer  :: ntp,  imax, &  ! ntp is the be multiplied to the tf by poinc subroutine to obtain the full period 
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
integer, private       ::  ind, idcor,  tfree
real(kind=dp), private ::  csangle_min,  p0, &    
                           tol, prsc, & ! tolerance control for the termination of Poincare map and Newton method
                           ds_max, ds_min, x0_fxd  ! fctn 

contains  
  include 'init_po.f90'
  
  include 'cor_gr.f90'     ! correction by general approach 
  include 'curv_angle.f90' ! compute the anlge between 3 consecutive points
  
  ! ***********  prediction *****************
  include 'pre_arclen.f90'   ! -- arc length method predictor 
  include 'pre_prd.f90'      ! TODO: un-debugged
!  include 'pre_cj.f90'      ! TODO   haven't done the routine yet 
  
!******************** Module-related Parameters ********************************
! ******************************************************************************  
   ! contol the working mode: To debug, or to execute
! ******************************************************************************  
  subroutine init_debug( debug0, dsctr0, anglectr0, csangle_min0,  isarc0, ds_min0, ds_max0 )
  
  ! TODO: this is too complicated.... disgard this part.... 
  integer, intent(in) :: debug0, dsctr0, isarc0, anglectr0 
  real(kind=dp), intent(in) :: csangle_min0, ds_min0, ds_max0
  
  debug = debug0
   
 ! Automatic control of the stepsize for continuation, with dsctr = 0, we obtain 2 families of p.o.s, which i showed to Gerard.
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
  
  tol  = tol0 
  prsc = prsc0
    
  end subroutine init_errctr
  

! ****************************************************************************** 
!  This  subroutine is to open the file to save the intermediate result, including
!  the files to check the correction,  matrix norm  for every iteration of Newton method 
!  Turns out the matrix is in good condition, because gnewnorm is small, that means the error will lead to a small correction(multiplied by the matrix),
!  but pay attention to the numerical error, which is 1.d-13 previously, but I ask the tol<1.d-14, stupid!!!!
 
! ****************************************************************************** 
  subroutine init_writedx  

!! Bound values for RK78, shared by all subroutines that calls rk78 (poinc, fctn)
  if(debug == 1) then 
    print*, 'tol,prsc',  tol, prsc;  read(*,*)
  endif 
  
  fdx = 32 
  
  if (issymt == 0) then 
    open(fdx,file='./dat/adx.dat',access ='append',status='replace')  
  else if (issymt == 1) then 
    open(fdx,file='./dat/dx.dat',access ='append',status='replace') 
  elseif (issymt == 2) then 
    open(fdx,file='./dat/grdx.dat',access ='append',status='replace') 
  endif 
  
       
  write(fdx,*) '# dxm   errfm       dx    index of correction(', ctr, ')'
    
  end subroutine init_writedx
  
  
!******************************* Subroutines ***********************************
  
  
!************************  Numerical Continuation ****************************** 
! Main subroutine to compute np new p.o.s starting from one initial guess yi(6), 
! and do the numerical continuation along the arc-length parameter.

!  To make sure the initial state is periodic, refine it without check for the first one.
!  the initial condition for all p.o.s are stored in file fpo

! 	Input Varialbles
!  yi(6)    initial condition (ie, x,y=0,z, vx=0, vy, vz=0)
!  tp0      the initial guess for tpo 
!  npo      number of new p.o.  if np=1, just do differential correction, without continuation

!  dir      dirction along the vector field, +1:increase,-1:decrease
!  ds       step size for the continuation (adaptive control), lower and upper bound [ds_min, ds_max] 
!  fpo      file tag for the initial states of the p.o.s

!     Output Variables
!  ynew(npo,8)    the initial condition for new p.o.s, npo*8 ( TP, XPO, CJ )
!                 save the initial p.o in the first row
!  i              the real number of available p.o.s is (i-1)

!  	Private Module-base Varaibles
! tol       tolerance for vx,vz to be zero (ie,1.d-13 or 1.d-16 in gerard's routine)
! prsc      presion to terminate the iteration(ie, 1.d-13)
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
  real(kind=dp), intent(out)    :: ynew(npo, 8)  ! 6+TP+CJ


  external :: deriv, gr_cj  ! vector field  & jacobi constant
  
  ! Local Variables  -- be careful with the dimension
  integer :: ispo, niter, nc, nr, i, & 
             nds, nds_new, incrds,  dir0  
 
             
  real(kind=dp) ::  y0(6), yf(n), g(ntar, nctr), yctr(nctr), &
                cj, cj2, dcj, &
                curv(5, nctr), curv3(3, nctr), cham(4, nctr), cham5(5,nctr), &
                curv_prd(5, nctr+1), curv_cj(5, nctr+1), csangle, vf(nctr), vf6(6) 
 
 
  smh = 0 ! we stop the continuation with too small step size...
   
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
    
    if(con_stop == 1) return 
     
    if( ispo == 1) then 
    ! check the angle between the last 3 consecutive points, if not in the admissable bound, decrease the arc step without update curv
      print*; print*, ' ************** Refined initial state: tp, y0 ******************'
      print*, tp0, y0; print*
      
      ! -------------------- check curve ------------------
      ! if nds > 2, at least 2 available rows in cham and curv , and 3 computed p.o.s
      ! compute the angle between the 3 consecutive points, and decide here if we need to decrease ds or not
      
      nr = min0(nds,5)   ! available rows in curv
      if (tfree == 1 ) nc =  nctr - 1   ! available column in curv 
      if (tfree == 0 ) nc =  nctr  
      
      if (nds .ge. 2 ) then  
        curv3(1:2, :)   = curv( nr-1 : nr,  : )
        curv3(3, 1: nc) = y0(ctr)
        if (tfree == 1)   curv(3, nctr) = tp0  ! period is always the last control variable
       
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
          curv(nds+1, 1: nc)  = y0(ctr) 
          if (tfree == 1) curv(nds+1, nctr) = tp0
    
        else 
          curv(1:4,:) = curv(2:5,:)
          
          curv(5, 1:nc ) = y0(ctr)
          if (tfree == 1) curv(5, nctr) = tp0
        endif 
        
       
        if (isarc == 2) then  ! along the period !!! TODO -- discard at the moment 2017-02-20 16:19:11 
          curv_prd( min0(nds+1, 5), nctr+1) = tp0 

!        elseif(isarc == 3) then ! along the energy...
!          curv_cj( min0(nds+1, 5), nctr+1) = cj
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
    
        yctr(1 : nc) = y0(ctr)  
        if (tfree == 1) yctr(nctr) = tp0
      endif 

!----------------  use arc-length parameter  -----------------------------------
      if ( isarc == 1 .or. ipo < 3) then 
        dir0 = dir 
!  subroutine pre_arclen(ipo, ds, nds, dir, cham, cham5, curv, yctr, c, incrds,  con_stop) 
        call  pre_arclen(ipo, ds, nds, dir, g, cham, cham5, curv, yctr, csangle, incrds, con_stop)  
        if ( isarc /= 1 ) dir = dir0 
      else
       
        call deriv(0.d0,yf(1:6),6,vf6) ! compute the vector field
        
        vf = vf6(ctr)
        
!        print*, 'vf=', vf
!        read*

! TODO: at this moment, haven't make the idea clear...        
        if(isarc == 2) then  ! along the period   --TODO 
          curv_prd(: , 1:nctr) = curv
          call pre_prd( ipo, curv_prd, vf, g,  nds, ds, dir, incrds, csangle, yctr, con_stop)
        
!        elseif(isarc == 3) then! along the energy -- TODO
!          curv_cj(: , 1:nctr) = curv
!          call pre_prd( ipo, curv_cj, vf, g,  nds, ds, dir, incrds, csangle, yctr, con_stop)
        endif 
        
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
      
      ! check if the period is too small 
!      if(tp0 < 1.d-4) then 
!        print*, 'Too small period obtained, not interested!  stop the continuation'
!        con_stop = 1; return 
!      endif 
      
      ! check the component, in order to avoid numerical accumulated errors
      do i = 1,6
        if (dabs( y0(i) ) < 1.d-16) y0(i) = 0.d0
      enddo 
      ynew(ipo,:) = (/tp0, y0, cj/)
      
  ! write to file labeled by fpo passed from the main routine, PoInSt.dat  
      write(fpoinst,'(9e22.12, 1e12.4, i5)')  tp0, y0, cj, dcj, ds, dir  ! fpoinst - PoInSt.dat 
 ! ----------- check the energy ---------------------------   
      if(debug == 1)  then 
        call gr_cj(y0, cj2)
      
! check to see if there is any point to keep cj, dcj as the output! ckd--cj = cj2 anyway
        print*, 'check energy-- cj, dcj, cj2:', cj, dcj, cj2  
        write(*,'(9e22.12, 1e12.4, i5)')  tp0, y0, cj, dcj, ds, dir  !  PoInSt.dat   !  !ck print to screen
        read(*,*)
      endif 
!--------------------------------------------
  
      if(npo .eq. 1)  return  !  refinement for only one p.o.
   
    endif   
    
!  we succeed in orbaining an initial guess for the new p.o, update the value of counter and  state + period(if necessary)     
    ipo = ipo + 1 
    
    if(debug == 1) then 
      print*, 'check ctr', ctr;  print*; read* ! debug y0,  2017-02-21 16:19:54 
      print*, 'size(y0), y0:', size(y0), y0
      print*, 'nc, nctr', nc, nctr      ! not the value assigned in init_asymtpo2
      print*, 'yctr', yctr ; print*; read*
    endif 
   
    print*, 'yctr', yctr
    
    y0(ctr) =  yctr(1:nc)
    if (tfree == 1)  tp0 = yctr(nctr)
    
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

!  We have nctr unkowns and ntar equations(nctr>ntar, inderterminant equations),  
!  using LLS(Linear least square problem) solver and take least-square solution.

!  The main point is get the derivative of the target variables with respect to the control variables


! The cases when this refinement fails to reach convergence
!   1. too many iterations, niter>5.  based on the quadratic convergence of Newton method, 5 iteration is enough!
!   2. For the poincare map approach, too long time to arrive to the next intersection acrossing poincare section,  ispc = 0

! for Newton method, alway check two things:
!   1. the error in target varaibles    --- stop condition, satisfied means finish the correction
!   2. the modulus of the correction    --- if too small, no point to continue ...

!  For the general approach, instead of Poincare map approach, we use the whole variational matrix to do the refinement w.r.t X0 + TP
!  so we 6 equations and 7 unkonwns(period included) 

!  the 7 unkowns are (xi, yi, zi, vxi, vyi, vzi, T) 
!  and the 6 constraint equations are satisfied to return to the same initial point 
!  (xfi-xi=0, yf-yi=0, zf-zi=0, vxf-vxi=0, vyf-vyi=0, vzf-vzi=0)
  
!******************************************************************************* 
!     Input 
!  y0       initial state
!  tp0      inital guess of the period(depend on which approach to apply)

!     Output
!  y0(tp0)  the updated state and period
!  tp0      elapsed time to reach to the next mirror configuration,  which is 1 period for asymmetric p.o. 
!  yf       final state at the epoch of the mirror configuration
!  g        Jacobi Matrix of the F(the periodicity equation ), will be used for the continuation
!  cj       Jacobi constant
!  dcj      difference of the Jacobi constant between the initial and final point along the p.o.
!  ispo     flag to indicate if the correction succeeds(1) or fails(0)
!  niter    number of iteration for the Newton Method, used as criteria to adjust the stepsize for continuation

!  deriv    subroutine to compute vector field of the equation of motion
!  gr_cj    subroutine to compute the Jacobi Constant(or energy)
! --------------------------------------------------------------------------------
  subroutine dfcr(y0, tp0, yf, g, cj,dcj, ispo, niter, deriv, gr_cj)
  
  implicit none  
  integer, parameter :: dp = kind(1.d0), n = 42 ! both the state and the variational equations are desired, so n=42
  
  
  integer, intent(out)       :: ispo, niter
  real(kind=dp), intent(out) :: yf(n), g(ntar, nctr), cj, dcj 
  real(kind=dp), intent(inout)  :: y0(6), tp0
  
  external ::  deriv, gr_cj
  
  
  ! Local Varaibles
  integer        :: i, debug, ispc, info 
  real(kind=dp)  :: f(ntar), pf(n), dx(nctr), dxm, fm, cj2, tp0_copy, pof(6)
                    ! yit(6) !  yi(n), y0copy(6) ! , dx2(nctr) 
!  real(kind=dp)  :: dx2(nctr)                   
  logical :: ok
  
  ! Saved the intermediate orbits. if file poit.dat is opened, close it, use inquire to obtain the status, 
  inquire( unit = 12, opened = ok )
  if(ok)  close(12)
  open(12, file='./dat/poit.dat',access ='append',status='replace')
 
  if (debug == 1)   print*, 'start dfcr, y0!', y0   

  tp0_copy = tp0 
  niter = -1 
  ispo  = 1 ! defaultly set Newton Method as successful
  dx = 0.d0;   dxm = 0.d0
  
  ! initial condition before the refinement
  if (tfree == 1) then 
    write(fdx, *)  'I.C.         ', y0, tp0
    write(*, *)    'I.C.         ', y0, tp0  ! print to screen   
  endif  
!  else
!    tmax = 1.5d0 * tp0 ! the maximum possible time for poincare approach to reach to the imax-th poincare section
!  endif   
  
    
  ! iterate to do the refinement, iter is the counter
  do  
    
    niter = niter+1 ! counter starts from 0
    
    if(debug == 1) then 
      print*, niter, '-th iteration'   
      read(*,*)
    endif
 
      
  ! Check the error in target function before the correction
  ! if niter > 6, too many iterations, bad initital guess, treat as a failure of the refinement, no matter ispo = 1 or not

  ! From Carles' class: For Newton method, which is quadratically convergent, 
  ! if the iteration is great than 5 or 6, then there is something wrong. 
  ! Provided that the initial guess is a good one that is close to the solution.

    if (niter > 6) then 
      print*, 'Newton method iterates > 6, dxm, fm',  dxm,  fm
      ispo = 0
      read*
      return
    endif   
 
    call gr_cj(y0, cj)

    if (issymt == 0) then     
    ! -- General approach: Integrate a whole period for next mirror configuration 
      
      ! subroutine gr_fctn(y0, tp0, f, g, yf, deriv, gr_cj )
      call gr_fctn(y0,  tp0,  f, g, yf, deriv, gr_cj)

    else                      
    ! -- Poincare map approach   

      call fctn(y0, 0, tp0_copy, f, g, tp0, yf, pf, ispc, deriv, gr_cj) 
      if(ispc == 0) then   ! we miss the case that the orbit escapes
        ispo = 0;  con_stop = 1;  return;
      endif 
       
    endif 
    
    call plob(y0, 0.d0,  tp0, 6, 6, 1, 1, 12, deriv, gr_cj, pof) 
    close(12)  ! close the file poit.dat to be able to plot for check !! discard.... 
    
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
    write(fdx,*) 'I.C.         ', y0(1:6)!finial condition of the poincare map point
    write(fdx,*) 'F.C.         ', yf(1:6)!finial condition of the poincare map point
    ! difference between the final and the initial condition, check only after a full revolution
    if(ntp==1) write(fdx,*) 'F.C.-I.C.    ', yf(1:6)-y0(1:6)
    
    !print to screen    
    write(*,*) 'I.C.         ', y0(1:6)!finial condition of the poincare map point
    write(*,*) 'F.C.         ', yf(1:6)
    if(ntp==1)  write(*,*) 'F.C.-I.C.    ', yf(1:6)-y0(1:6) 
    
    if(debug == 1) read*
   
 ! precision to terminate the correction, 1.d-16 from Gerard, which is of non-sense, because it beyonds the precision of double real
    if ( fm .le. tol)  then 
      print*, 'fm < tol', fm,  tol
      return   ! try 1.d-11  >should be no less than the numerical precision we are using: e=1.d-13 (error in rk78)
    endif 
    
    !  solve the LLS problem to obtain the correction     
    !  subroutine deltx_lns(f, g,  dx)
    !    call deltx_lns(f, g,  dx) ! discard this one, and use more secure one, deltx-keep-updated one 
    
    ! subroutine deltx( nr, nc, nrhs, a, b, x, info) 
    call deltx( ntar, nctr, 1, g, -f, dx, info) 

!    print*,  'Check deltx_lns and deltx, dx = :'
!    print*,  dx 
!    print*,  dx 
!    print*; read* 
    
    
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
    if(issymt == 0 ) then 
   ! general case: 7 control varaibles
      y0(ctr) = y0(ctr) + dx(1: nctr-1)
      tp0 = tp0 + dx(nctr)
      
    else 
      y0(ctr) = y0(ctr) + dx 
    endif 
    
    
    write(fdx,*) ! add a blank line to seperate iterations
    write(fdx,*) niter+1, '-th iteration'
    write(fdx,*) '|DX0|,|f|,DX0', dxm, fm, dx 
    
    if(issymt == 0 ) then 
      write(fdx,*) 'U.I.C.       ', y0, tp0 ! the updated initial conditon
    else 
      write(fdx,*) 'U.I.C.       ', y0 ! the updated initial conditon   
    endif 
    
    ! print to screen 
    write(*,*)
    write(*,*) niter+1, '-th iteration'
    write(*, *)  '|DX0|,|f|,DX0', dxm, fm, dx
    if(issymt == 0 ) then 
      write(fdx,*) 'U.I.C.       ', y0, tp0 ! the updated initial conditon
    else 
      write(fdx,*) 'U.I.C.       ', y0 ! the updated initial conditon   
    endif 
    
    if(debug == 1)  read*
    
  enddo 
   
  return
  end subroutine dfcr
  

!********************************************************************************
! Compute the Jacobi matrix g  and the error f  for Newton method 
!     Input 
!  x        initial state, (x,y,z,vx,vy,vz)
!  init     flag of the assignment of current state, y 
!  y(42)    state vector + STM, at the initial and final configuration 

!     Output 
!  f        dimension ntar, error to be corrected    
!  g        dimension ntar-by-nctr, Jacobi matix
!  tf       elapsed time to reach next mirror configuration 
!  vf       vector field at next mirror configuration
!  ispc     flag of status of the poinc, 1:successful, 0: failure

! deriv, gr_cj     external routine for vector field + energy 
! -------------------------------------------------------------------

  subroutine fctn(x,init, tp0, f,g,  tf,y,vf, ispc, deriv, gr_cj)
  
  implicit none  
  integer, intent(in)   :: init  ! flag of the assignment of initial condition x(6)
  integer, intent(out)  ::  ispc
  real(kind=dp), intent(in)  :: x(6), tp0  
  real(kind=dp), intent(out) :: f(ntar), g(ntar, nctr), vf(42), tf 
  real(kind=dp), intent(inout)  :: y(42)
  
  
  ! local varaibles
  integer        ::  i, j 
  real(kind=dp)  ::  phi(6,6), yi(42), hminim, tf0  
    
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

! subroutine poinc(sti, ndim, nvar, tdir, tf, stf, hminim, ispc, deriv, gr_cj)  --from  poinc.f90
  call poinc(yi, 6, 42, 1, tf, y, hminim, ispc, deriv, gr_cj)  !  --new one 
    
  ! how  to deal with the failure (ispc = 0) while trying to obtain small h in gr_rk78  
   
  ! if the time spend for the imax-th crossing is too small (<tp0/ntp/2), there is something wrong.
  ! we continue the integration from yi 
   if(debug == 1) then 
     print*, 'tf, tp0, yf by poinc: ', tf, tp0, y(1:6); print*; read*
   endif 
   
   if(ispc == 0 .and.  hminim < 1.d-14)  then 
      con_stop = 1; return
   endif 
   
!  do while (dabs(tf) < dabs(tp0)/ntp/2)   
!    print*, 'Too small tf, continue!' 
!    print*, 'tf, yf by poinc: ', tf, y(1:6); print*; read* 
!    tf0 = tf 
!    yi  = y 
!    call poinc(yi, 6, 42, 1, tf, y, hminim, ispc, deriv, gr_cj)  !  --new one 
!    tf = tf0 + tf  
!  enddo 
       

  ! Do we need to check hminim? for correction to the crossing, it is reasonable to 
  ! use step less than 1.d-6     
  if(hminim < 1.d-6 ) then 
     print*, 'Toop small step: hminim < 1.d-6', hminim;  read(*,*) 
  endif 
    
  if (ispc == 0)  return !fail to reach to Poincare section, stop  the correction 
  
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
  ! remember not to forget to abstract the identity matrix  ! ckd!!! 
  if (issymt == 0 ) then 
    print*, 'Asymmteric case! Phi= Phi-I' !ckd!
    if(debug == 1) read*
    
    do i = 1, 6
      phi(i,i) = phi(i,i) - 1.d0
    enddo
  endif 
 
  ! this is only for Poincare-map approach....
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
 
    print*, 'Jacobi Matrix of target variable w.r.t control variables' 
    do i = 1, ntar 
      write(*,*)  phi(tar(i), ctr)
    enddo 
    read*

  endif 
  
  return
  end subroutine fctn
  
  
   
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

!     Input 
!  ntar     number of rows in G 
!  nctr     number of columns in G
!  g        Jacobi Matrix of F w.r.t (X0 + T) 
!  f        Target function
!   
!     Output 
!  dx       the least-norm correction by sloving G * DX = -F 

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
  RCOND = 1.d-14
 
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

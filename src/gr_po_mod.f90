! Compute periodic orbit for a spacial problem, the phase space is 6d 

module gr_po_mod 

use dp_mod
use pi_mod 

implicit none

! ******************** Declare Public *****************
! public by default, the most commonly used varaibles which are also accessible by any subroutine that uses this module

real(kind=dp) :: tol, tol_err


integer, private :: debug  
 
integer  :: ntp,  &     ! ntp is the be multiplied to the tf by poinc subroutine to obtain the full period 
            ntar, nctr, nc, tfree  ! number of target variables(equations) and control variables(free variables)


! declare  allocatable arrays to store the index of target and control variables  -- to check
integer, allocatable :: tar(:), ctr(:)


! -------------- declaration of subroutines from lapack or nag( external package)
real(kind=dp), external          :: dlange ! compute the norm of array,  dlange is a function, declare the type
real(kind=dp), external, private :: dnrm2 ! compute the norm of array, from package NAG Fortran Library Routine

external  deriv,  gr_cj ! i don't know if this helps

    
!  private - only shared by the subroutines, cannot be accessable from outside the module 
integer, private       ::  ind,  idcor
real(kind=dp), private ::  x0_fxd  

 ! the files to check the iteration, matrix norm --discard, mostly caused by numerical error
integer, private :: fdx 
    
contains 
  
  
! ******************************************************************************  
  ! the error control for termination of Newton Method
! ******************************************************************************  
  subroutine init_errctr(debug0, tol0, tol_err0)
    integer, intent(in)     ::  debug0 
    real(kind=dp), intent(in) :: tol0, tol_err0  
    debug   = debug0
    tol     = tol0 
    tol_err = tol_err0
  end subroutine init_errctr
  
! ****************************************************************************** 
  ! ****************************************************************************** 
! The initializaiton for asymmetric p.o.:  idcor = (5, 6, 7) 

! ****************************************************************************** 
  subroutine init_gr_po(idcor0, ind0, x0_fxd0)
  
  ! this is try to add a constraint to fix a component of the state vector, y(ind0)= x0_fxd0
  ! so the control variables are still all the six components, but the number of target function increases by 1 
  
  integer,intent(in) :: idcor0, ind0
  real(kind=dp),intent(in) :: x0_fxd0
  
  integer ::  allocatestatus, nrow, ncol, i, j


!  module-based variables  
  idcor = idcor0 ! supposed to be 5 
  tfree  = 0 ! by default, fixed-time shooting
  
  ind    = ind0
  x0_fxd = x0_fxd
  
  if(idcor == 1) then  ! fix time shooting + fix one coordinate
    tfree = 0 
    nctr = 5
    ntar = 5
    
    nrow = ntar 
    ncol = nctr 
  
  elseif( idcor == 5 ) then   !-- only finish this one 
  !-- doesn't work very well, the continuation only proceed for a few orbits 
  !  free time shooting + additional constraint
    nctr = 7
    ntar = 6
    
    tfree = 1
    nrow = ntar - 1
    ncol = nctr - 1 
  
  elseif (idcor == 7) then  ! improve idcor=5 strategy
  ! free time shooting, deliminate ind from tar to satisfy x(ind) = x0_fxd 
    nctr = 6  ! five state components + TP
    ntar = 5
    
    tfree = 1 
    nrow = ntar 
    ncol = nctr-1   
    
    
  elseif( idcor == 6 ) then 
    ! free time shooting -- no additional constraint....
    tfree = 1
    nctr = 7
    ntar = 6 
     
    nrow = ntar
    ncol = nctr - 1 
    
  end if
  
 ! allocate the memory for tar and ctr 
  allocate ( tar(nrow), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"

  allocate ( ctr(ncol), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"
  
  
 ! assign the vaules for tar and ctr 
    
  if ( idcor == 5  ) then ! the current used one, additional constraint: y(ind) = x0_fxd
    tar = (/1,2,3,4,5/)   ! x_6^f is satisfied by conservation energy law
    ctr = (/1,2,3,4,5,6/) 
  
  elseif( idcor == 6  ) then ! --TODO with full state vector as the control variables    
    tar = (/1,2,3,4,5,6/) 
    ctr = (/1,2,3,4,5,6/) 
  
  elseif( idcor == 7 .or. idcor == 1)  then  ! -- to approve idcor = 5 
    tar   = (/1, 2, 3, 4, 5 /)
     
    ! instead of an additional equation, fix x(ind)=x0_fxd as the initial condition,
    ! and eliminate ind-th component from ctr, as well as ind-th column in Phi
    j = 0
    do i= 1, 6, 1
      print*, 'ind=', ind, 'i=', i 
      read*
      
      if( i /= ind ) then 
        j = j + 1
        ctr(j) = i
      endif 
    end do
      
  end if  
  
  nc = nctr 
  if(tfree == 1) nc = nctr - 1 
  print*,' ind, x0_fxd', ind, x0_fxd
  read*
  
  print*, 'idcor=', idcor, 'nctr=', nctr, 'ntar=', ntar
  print*, 'ctr=', ctr, 'tar=', tar 
  read* 
  end subroutine init_gr_po
  

!********************************************************************
!  Using modified Newton method to refine the p.o,   using Lapack linear equation solver

! The cases when this refinement fails to reach convergence
!   1. too many iterations, niter>6.  based on the quadratic convergence of Newton method, 6 iteration is enough!
!   2. too long time to arrive to the next intersection acrossing poincare section, ispc = 0

! for Newton method, alway check two things:
!   1. the error in target varaibles    --- stop condition, satisfied means finish the correction
!   2. the modulus of the correction    --- if too small, no point to continue ...

!  	Input Variables:
!  y0(n)    initial condition 

!  	Output Variables:
!  y0, yf   the refined initial and final conditon
!  tf	      half period or one period, depends on the number of intsec-imax
!  g(ntar, nctr)	JACOBIAN MATRIX OF F(*), ntar by  nctr,  for symmetric p.o. :  3 variables and 2 equations
  
! --------------------------------------------------------------------------
  subroutine refn_po(y0, tp0, yf, dg, niter, ispo, deriv, gr_cj)  
   
  implicit none  
  real(kind=dp), intent(out)    :: yf(42), dg(ntar, nctr)
  real(kind=dp), intent(inout)  :: y0(6), tp0
  integer, intent(out) :: ispo 
  external :: deriv, gr_cj
  
! Local Varaibles
  integer ::  i,   info, niter 
  real(kind=dp)  :: ferr(ntar), dx(nctr), dxm, fm 
     
  if (debug == 1)   then 
     print*, 'start refn_p0, tp0, y0!'
     print*, tp0,  y0; read*
  endif    
  
  niter = -1 
  ispo = 1 ! defaultly set Newton Method as successful
  dx = 0.d0;  dxm = 0.d0
  
  
! iteration to do the refinement, iter is the counter
  do  
    niter = niter+1 ! counter starts from 0
         
   ! Newton method  is quadratically convergent. 
   ! if the iteration is great than 5 or 6, then there is something wrong, 
   ! provided that the initial guess is a good one, which is  close to the solution.
    if (niter > 10) then 
      print*, 'Newton method iterates > 10, dxm, fm',  dxm,  fm
      if(fm .gt. tol_err)  ispo = 0  
      read*;    return
    endif 
    
    write(fdx,*) '***********', niter, 't-th iteration ***************'
  
    !  subroutine gr_fctn(y0, tp0, f, g, yf, deriv, gr_cj )
    call gr_fctn(y0,  tp0,  ferr, dg,  yf, deriv, gr_cj)
    
    call deltx( ntar, nctr, 1, dg, -ferr, dx, info)

    if(debug == 1) then
      do i = 1, ntar, 1
        print*,  dg(i,:)
      enddo  
      print*;   read(*,*)
    endif 
    
    fm  = dnrm2(ntar, ferr, 1) / dsqrt(dble(ntar))
    dxm = dnrm2(nctr, dx, 1) / dsqrt(dble(nctr))
    
    print*, 'tol, tol_err', tol, tol_err; print*; read*
    if ( fm .le. tol_err)  then 
      print*, 'fm < tol', fm,  tol_err
      return    
    endif 
    
    ! if the correction is below the precision, regard as success of the refinement
    if ( dxm .lt. tol )  then 
      print*, '|dx|<', tol, 'stop iteration!'
      read*; return
    endif 
    
    ! Update the initial state, + TP(tfree==1)
    y0(ctr) = y0(ctr) + dx
    if(tfree == 1) tp0 = tp0 + dx(nctr)
    
    write(fdx,*) ! add a blank line to seperate iterations
    write(fdx,*) niter+1, '-th iteration'
    write(fdx,*) '|DX0|,|f|,DX0', dxm, fm, dx 
    write(fdx,*) 'U.I.C.       ',tp0,  y0 ! the updated 
    
  enddo 
   
  return
  end subroutine refn_po

!********************************************************************************  
! For the asymmetric approach, please note that because we are refine w.r.t the initial condition
! So the Jacobi Matrix w.r.t the correction on the initial state is actually Phi - I  

! Instead of Poincare map approach, we use the variational matrix to do the refinement of the asymmetric p.o.,
! such that we have 5 equations and 6 unkonwns(period included). 

! We are not using the modified Jacobi Matrix in which the differential 
! of the Poincare map w.r.t time t is represent by the other unkonwns 
! in order to decrease the number of unknowns by 1.  

! the 7 unkowns are (xi, yi, zi, vxi, vyi, vzi, T) 
! and the 6 constraint equations are satisfied to return to the same initial point 
! (xfi-xi=0, yf-yi=0, zf-zi=0, vxf-vxi=0, vyf-vyi=0, vzf-vzi=0)

! If we fix one component in the control variable, it is better to remove the corresponding 
! column in the Jacobi Matrix, to avoid a zero determinant for the characteristic curve computation. 
 
subroutine gr_fctn(y0, tp0, f, g, yf, deriv, gr_cj )
  
  implicit none  
  real(kind=dp), intent(in)   ::  y0(6), tp0
  real(kind=dp), intent(out)  ::  g(ntar, nctr), f(ntar),  yf(42) 

  external ::  deriv, gr_cj
  
  ! Local Varaibles
  integer        ::  i, j, ispl, ftag, nc
  real(kind=dp)  ::  phi(6, 6), yi(42), vf(6)  

  
  ! number of control variable as state variable, if we change the period,
  ! it should be nctr-1
  nc = nctr
  if(tfree == 1) nc = nctr - 1
  
  print*, 'tfree= ', tfree, ', nctr=',nctr, ', ntar=', ntar 
  read*
  
  ! Integrate for an approximated period of tp0 
  yi = 0.d0
  yi(1:6) = y0
  yi(7:42:7) = 1.d0 


  ! In debug mode, save the integration data of the orbit, other wise, do not save the data of the orbit     
  if(debug == 1) then 
    ispl = 1;    ftag = 12  ! 12 is opened before calling the dfcn 

  else 
    ispl = 0;    ftag = 0 
  endif  
  
!  subroutine plob(y0, t0, tf, ndim, nvar, tdir, ispl, ftag, deriv, gr_cj,  y) 
  call plob(yi, 0.d0, tp0, 6, 42, 1, ispl, ftag, deriv, gr_cj,  yf)  
    
  if(debug == 1) then 
    print*, 'tf, yf';    print*,  tp0, yf(1:6);     print*;  read*
  endif   
  
  ! compute the Jacobi Matrix of F w.r.t (X0,T)
  ! G = [Phi, VF]
  phi = reshape( yf(7:42), (/6,6/) ) ! state transition matrix
   
  if (debug == 1) then  
    print*,'Phi'
    do i = 1, 6
      print*, phi(i,:)
    enddo
    read*
  endif  
  
  g(:, 1:nc) = phi(tar, ctr)
  
  ! for asymmetric p.o., refine to the initial condition 
!  if( issymt == 0 ) then
    do i = 1, ntar, 1
      do j = 1, nc, 1
        if( tar(i) == ctr(j) ) then
          g(i,j) = g(i,j) - 1.d0
        end if
      end do
    end do
!  end if  
  
  
  ! in free-time shooting, the last colomn is the velocity 
  if( tfree == 1) then 
    call deriv(tp0, yf(1:6),  6, vf)
    g(:, nctr) = vf(tar) 
  endif

     
  if (debug == 1) then    
    print*,'Jacobi Matrix: G'
    do i = 1, ntar
      print*, g(i,:)
    enddo
    read*
  endif  
 
  ! refine to the initial condtion 
  f = yf(tar) - y0(tar)
  
  if (debug == 1) then   
    print*,'Constraint function: F'
    print*, f
    read*
  endif 
  
end subroutine gr_fctn   
  
!***********************************************************************
!     ****   pred_h   ****
! Prediction for the new p.o. along the energy, based on the two previous refined p.o.s 
! and use linear extrapolation

!       Input Variables 
!  x0_pre     dimension 2-by-nctr, the free parameter of the two previous p.o.s 
!  h_pre      dimension 2, the parameter (energy or period) of the two previous p.o.s
!  dh         step size of energy to continue  
!  ipo        number of current refined p.o.s 


!!!!  ispo       flag to show if we succeed the refinement with the last guess
!             if not, we replace the first one by the middle point of the previous good ones.

!       Output Variables 
!  x0, h       prediction for the new p.o.            

!  Routine Used:
!     

!  Finally Revised by Yu -- 20160929
!***********************************************************************
subroutine pred_h( x0_pre, h_pre, dh, ipo, x0, h)

use dp_mod
implicit none
 
! Input  and Output Declaration  
integer, intent(in)            ::  ipo  
real(kind=dp), intent(in)      ::  x0_pre(2, nctr), h_pre(2), dh
real(kind=dp), intent(out)     ::  x0(nctr), h 
  
integer :: i  
!  if(debug == 1 .or. dabs(dh) < 1.d-5) then 
    print*,'We have already', ipo, 'refined p.o.s. step =', dh
    print*, 'The previous 2 p.o.s:'
    print*, 'x0_pre'
    do i = 1, 2
      print*, x0_pre(i, :)
    enddo  
    
    print*, 'h_pre',  h_pre; print*; read*
!  endif   
  
  
  if( ipo <= 1 ) then
    h  = h_pre(1) + dh
    x0 = x0_pre(1, :) 
    return
  
  else 
  
    ! discard the linear extrapolation for a middle point, since the slope will be the same
    ! we take an extrapolation   
    x0 = x0_pre(2, :) + ( x0_pre(2, :) - x0_pre(1, :) ) / ( h_pre(2) - h_pre(1) ) * dh
    h = h_pre(2) + dh
  endif  
  
  print*,    
  return  
end subroutine pred_h


!***********************************************************************
! Adams-bashforth predictor method, Linear multistep method is to computationally 
! solve the ordinary differential equations

! 20160329  Multistep methods refer to several previous points and derivative values. 
!           Adams method uses fixed step size for all the previous points.

!    So be careful with the input cam 

!          input parameters:
!     np         number of computed periodic orbits
!     x(*)       initial conditions of the last p.o.
!     hh         integration step along arc length of solution locus
!     cam(*,*)   (see subroutine camp)

!          output parameters:
!     x(*)       new aproximated initial conditions for p.o.

!  Revised by Yu, 20160219
! --------------------------------------------------------------------------
  SUBROUTINE ADAMS(NP,X,HH,CAM)
  implicit none 
 
  integer, intent(in) :: np
  real(kind=dp), intent(in) :: hh, cam(4, nctr) 
  real(kind=dp), intent(inout) :: x(nctr) !update x(3)
  
! Local Variables 
  integer :: i
  
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
 
  return
  end subroutine adams

  
!***********************************************************************  
  subroutine champ(np,g,dir, cham)
! this one is hard to modify for the asymmetric p.o. case, we have to redefine every component of the derivative of the arc length parameter with respect to control variables

! add the theoretical part......

! 20160316 -- todo
! forget about it at this moment, because we are not going to do the continuation now. 
! the point is to compute one asymmetric p.o., and then explore if we could use other attributes
  
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
  real(kind=dp), intent(in) :: g(ntar, nctr) 
  real(kind=dp), intent(out) :: cham(4, nctr)  
  
! Local Variables 
  integer :: i,   nnp
  real(kind=dp) :: a(nctr),  c, gi(ntar, ntar), det 
                    
 
  ! --- this maybe the key part to be modified
  
  ! the vector field of the characteristic curve. which is its integral curve. 
  ! the characteristic curve is obtained by integrating the vector field a using Adams method
  
  ! the following is the computation of the vector field a, the components of a are the derivative of the control variable w.r.t the arc-length paramter s.
  
!------------------------  to be modified --------------------------------------
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
  
  !
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
 
! ------------------------------------------------------------------------------
  if (debug == 1) then 
    print*, 'check the kernel by detmat' ! --- ckd, match the one by direct computation 
    print*, a;  print*
    
    print*, 'check g in champ'
    do i = 1, ntar
      print*, g(i, :)
    enddo
    
    print*;  read*
  endif 

  if(np .gt. 4)  then 
    cham(1:3, :) = cham(2:4, :) !checked, assignment ok!
    cham(4, :)   = a 
    nnp = 4
    
  else 
    cham(np, :) = a  
    if (np .eq. 1)   return
    nnp = np
    
  endif 
 
! detect whether the new point obtained is a bifurcating point. ---- Modify this part 
! If the vector product of the new vector field and the previous one is < 0, change the sign of the vector field
  
  c = 0.d0
  do  i = 1,3
    c = c + cham(nnp,i)*cham(nnp-1,i)
  enddo

  if( c .gt. 0.d0 )  return
  
!  this part for the reverse of the sense, detect bifurcation? -- need to be understood later
!  comment this part... to see what happens, the comparision of the result shows what???? 
!  Add the conclusion here! 

  dir = -dir
  cham(nnp, :) = - cham(nnp, :)
  write(*,*)  'Reversal in the Sense along the Characteristic Curve!'
  
  return
  end subroutine champ
  
end 

! Compute periodic orbit for a planar problem, which means the phase space is 4d 
! If we take advantage of the symmetry of the system, we will reduce the problem to 2d,
! for instance, we have  two unknows (x, vx) when the system is symmetric w.r.t. x-axis in opposite time direction 

! Furthermore, we want the initial and final velocity  to be perpendicular to the symmetric object (x axis in this case )
! that is, vx_i = 0, 
! we have only ndim/2-1 = 4/2-1 = 1 equations : vx_f = 0
! with only 1 unknown : x0

! this simplicify the problem a lot! 

! use PoincMap_XXXX to compute the differential 
!module symt_pomod 

module po_plf_mod 

use dp_mod
use plf_mod


implicit none

! ******************** Declare Public *****************
! public by default, the most commonly used varaibles which are also accessible by any subroutine that uses this module

real(kind=dp) :: tol, tol_err
integer, parameter, private  :: ndim = 4, n = 2 


integer, private :: debug  
 
integer  :: ntp,  &     ! ntp is the be multiplied to the tf by poinc subroutine to obtain the full period 
            ntar, nctr  ! number of target variables(equations) and control variables(free variables)


! declare  allocatable arrays to store the index of target and control variables  -- to check
integer, allocatable :: tar(:), ctr(:)


! -------------- declaration of subroutines from lapack or nag( external package)
real(kind=dp), external          :: dlange ! compute the norm of array,  dlange is a function, declare the type
real(kind=dp), external, private :: dnrm2 ! compute the norm of array, from package NAG Fortran Library Routine


 ! the files to check the iteration, matrix norm --discard, mostly caused by numerical error
integer, private :: fdx 
    
contains 
  
  
! ******************************************************************************  
  ! the error control for termination of Newton Method
! ******************************************************************************  
  subroutine init_errctr 
 
  tol     = 1.d-10 
  tol_err = 1.d-10
  
  
  end subroutine init_errctr
  
! ****************************************************************************** 
! ind_zero is the index of target varaibles in the reduced state 
! for example, we are in R4 : (x, y, vx, vy), and take poincare section y=0
! then we have 2 freedoms (x, vx), and we want vx=0, so nzero = 1, and ind_zero = 2

! For R6 case, a little more complicated 

! ****************************************************************************** 
  subroutine init_po_plf( ind_zero )
!  use poinc_mod, only : ind_fun 
  
  integer, intent(in) ::  ind_zero(1)  ! index of symmetry, used to assign the public values accordingly

  integer ::  allocatestatus, i, kc, kt
  
  
  debug = 0
  print*, 'check if n=2, if not, something is wrong, n=', n, 'ndim=', ndim
  print*; read*
  
  call init_errctr
  
  
! --- allocate memory for ctr and tar ----
!  if(h0_fxd == 1) then 
    ! if we fix the energy level, we only have 1 free parameter 
    nctr = ndim/2 - 1
    ntar = ndim/2 - 1
!  else 
       
!  endif 
  
  
  allocate ( tar(ntar), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"

  allocate ( ctr(nctr), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"
   
  ntp  = 2  
  
  ! TODO:  higher dimension, 
  ! for the moment, we only deal with planar problem, do not make it too difficult 
  
  ! assign the index of control parameters and target parameters, they together 
  ! form ind_fun( the free parameters in Poincare map) 
  kc = 1; kt = 1
  do i = 1, n
    if( i  == ind_zero(kt) ) then 
      tar(kt) = i ! ind_fun(i)
      kt = kt + 1
    else  
      
      ctr(kc) = i ! ind_fun(i)
      kc = kc + 1   
    endif    
  enddo 
  
  ! TODO: consider the other kind of symmetry? for the moment, no!
  if (debug == 1) then          
    print*, 'check:   control and target variables'
    print*,  ctr, ',', tar 
    read*
  endif   
  
  end subroutine init_po_plf


!  !***********************************************************************
!  !     ****  pofam    ****
!  !
!  
!  !       Input Variables 
!  !  a            
!  
!  !       Output Variables 
!  !  a            
!  
!  !  Routine Used:
!  !     
!  
!  !  Finally Revised by Yu -- 20160
!  !***********************************************************************
!  subroutine ( in, out  optional)
!  
!  use dp_mod
!  implicit none
!   
!  ! Input  and Output Declaration   
!  real(kind=dp), intent(in)      ::    
!  real(kind=dp), intent(out)     ::    
!   
!  ! Local Variable
!  real(kind=dp)  ::
!    
!      
!    source
!  
!    return  
!  end subroutine 
  
!**************************   Differential Correction  ************************* 
!  Using modified Newton method to refine the p.o,   using Lapack linear equation solver

! The cases when this refinement fails to reach convergence
!   1. too many iterations, niter>6.  based on the quadratic convergence of Newton method, 6 iteration is enough!
!   2. too long time to arrive to the next intersection acrossing poincare section, ispc = 0

! for Newton method, alway check two things:
!   1. the error in target varaibles    --- stop condition, satisfied means finish the correction
!   2. the modulus of the correction    --- if too small, no point to continue ...

!  	Input Variables:
!  y0(n)    initial condition on poincare section, dimension n
!  imax     the times of the crossing through the section(y=0)

!  	Output Variables:
!  y0, yf   the refined initial and final conditon
!  tf	      half period or one period, depends on the number of intsec-imax
!  g(ntar, nctr)	JACOBIAN MATRIX OF F(*), ntar by  nctr,  for symmetric p.o. :  3 variables and 2 equations
  
!  function used: poinc, deriv, gr_cj, fctn, deltx
! --------------------------------------------------------------------------
  subroutine refn_po(y0, tf, yf, dg, niter, ispo) !, deriv, gr_cj)
   
  implicit none  
  real(kind=dp), intent(out)    :: tf, yf(ndim), dg(ntar, nctr) 
  real(kind=dp), intent(inout)  :: y0(n) 
  integer, intent(out) :: ispo 
  external :: PoincMap_tf_plf
  
  
! Local Varaibles
  integer ::  i,   info, niter, ispc
  real(kind=dp)  :: ferr(ntar), dx(nctr), dxm, fm 
     
  if (debug == 1)   then 
     print*, 'start dfcr, y0!', y0; read*
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

! subroutine gdg_poinc_lf( pv0, tf, pvf, dg, ferr, gr_tf_poinc)
    call gdg_poinc_lf(y0, tf, yf, dg, ferr, ispc, PoincMap_tf_plf) 
    
    if(ispc == 0) then 
      ispo = 0; return
    endif 
    
    call deltx( ntar, nctr, 1, dg, -ferr, dx, info)

    if(debug == 1) then
      do i = 1,  ntar
        print*,  dg(i,:)
      enddo  
      print*;       read(*,*)
    endif 
    
    fm = dnrm2(ntar,  ferr, 1) / dsqrt(dble(ntar))
    dxm = dnrm2(nctr, dx, 1) / dsqrt(dble(nctr))
    
    
    if ( fm .le. tol_err)  then 
      print*, 'fm < tol', fm,  tol_err
      return    
    endif 
    
    ! if the correction is below the precision, regard as success of the refinement
    if ( dxm .lt. tol )  then 
      print*, '|dx|<', tol, 'stop iteration!'
      read*; return
    endif 
    
    ! Update the initial state, for all the control variables 
    y0(ctr) = y0(ctr) + dx
    
    write(fdx,*) ! add a blank line to seperate iterations
    write(fdx,*) niter+1, '-th iteration'
    write(fdx,*) '|DX0|,|f|,DX0', dxm, fm, dx 
    write(fdx,*) 'U.I.C.       ', y0 ! the updated 
    
    ! print to screen 
!    write(*,*)
!    write(*,*) niter+1, '-th iteration'
!    write(*, *)  '|DX0|,|f|,DX0', dxm, fm, dx
!    write(*,*)   'U.I.C.       ', y0 
!    
  enddo 
   
  return
  end subroutine refn_po

!***********************************************************************
!     ****   gdg_poinc_lf   ****
! Given the initial point, Compute the Jacobi matrix for for the linear equations, and the error to refined
! since this is for the computation of periodic orbit, so we need to return the return time 

!  by Newton Method, as input of deltx 
!            dg *  dX = - ferr

!       Input Variables 
!  pv0      the initial state of P.O.            

!       Output Variables 
!  tf       the return time 
!  dg      dimension ntar X nctr, the Jacobi matrix of the linear equation  
! ferr     dimension ntar, the error to cancel

! Module- based Varaibles 
!   ndim, ntar, nctr, n(dimension of poincare map)

!  Routine Used:
!     PoincMap_tf_plf

!  Finally Revised by Yu -- 20160921
!***********************************************************************
subroutine gdg_poinc_lf( pv0, tf, pvf, dg, ferr, ispc, gr_tf_poinc)
use poinc_mod, only : ind_fun

implicit none
 
! Input  and Output Declaration   
real(kind=dp), intent(in)      :: pv0(n)  
integer, intent(out)           :: ispc  
real(kind=dp), intent(out)     :: tf, pvf(ndim),  dg(ntar, nctr), ferr(nctr)   
 
external :: gr_tf_poinc 
 
! Local Variable
integer  :: i, ispl
real(kind=dp)  :: dpdx(n, n), pv_pc(n)
logical :: ok  ! status of the file poit.dat, opened or not, which saves the orbit of all the iterations      
   
  ! Saved the intermediate orbits. if file poit.dat is opened, close it, use inquire to obtain the status, 
  inquire( unit = 444, opened = ok )
  if(ok)  then 
    close(444)
    open(444, file='./dat/lf_poit.dat',access ='append', status='old')
  else 
    open(444, file='./dat/lf_poit.dat',access ='append', status='replace')
  endif 
  
  ispl = 1
!  print*, 'Do u want to save the intermediate p.o.? Yes(1), No(0)'
!  read(*,*) ispl 
    
!  subroutine PoincMap_plf( pvin, tf,  pvf,  isdp, dp_dx)
  call  gr_tf_poinc( pv0, tf, pvf, 1, dpdx, ispc)
  
  ! we have to deal with the case, that we fail to reach to the 
  if(ispc == 0) return 
  
  pv_pc =  pvf(ind_fun) 
  ferr =  pv_pc(tar) 
  dg = dpdx(tar, ctr) 
    
  if(debug == 1) then 
    print*, 'After half revolution '
    print* ,' tf, pvf ', tf, pvf;  print*;  
    print*, 'Differential of Poincare map:  ', shape(dpdx)
    
    do i = 1, n, 1
      print*, dpdx(i, :)
    enddo 
    
    print*; read*
  
    ! the error 
    print*, 'tar = ', tar 
    print*, 'ferr = ', ferr ; print*; read*
  
    print*, 'Jacobi matrix: '
    print*,  dg 
    print*; read*
  endif  
  
  return  
end subroutine gdg_poinc_lf

!***********************************************************************
!     ****   pred_h   ****
! Prediction for the new p.o. along the energy, based on the two previous refined p.o.s 
! and use linear extrapolation

!       Input Variables 
!  x0_pre     dimension 2, the free parameter of the two previous p.o.s 
!  h_pre      dimension 2, the energy level of the two previous p.o.s
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
subroutine pred_h( x0_pre, h_pre, dh, ipo , x0, h)

use dp_mod
implicit none
 
! Input  and Output Declaration  
integer, intent(in)     ::  ipo !, ispo 
real(kind=dp), intent(in)      ::  x0_pre(2), h_pre(2), dh
!real(kind=dp), intent(inout)   ::  dh  
real(kind=dp), intent(out)     ::  x0, h 
 
  if(debug == 1 .or. dabs(dh) < 1.d-5) then 
    print*,'We have already', ipo, 'refined p.o.s. step =', dh
    print*, 'The previous 2 p.o.s:'
    print*, 'x0', x0_pre
    print*, 'h', h_pre ; print*; read*
  endif   
  
  
  if( ipo <= 1 ) then
    h  = h_pre(1) + dh
    x0 = x0_pre(1) 
    return
  else 
  
  ! -- discard, if this is linear extrapolation, there is no point in adding a middle point, the slope will be the same
!    if(ispo == 0) then 
!      x0_pre(1) = ( x0_pre(2) + x0_pre(1) )/ 2
!      h_pre(1) = ( x0_pre(2) + x0_pre(1) )/ 2
!    endif 
      
    x0 = x0_pre(2) + ( x0_pre(2) - x0_pre(1) ) / ( h_pre(2) - h_pre(1) ) * dh
    h = h_pre(2) + dh
  endif  
    
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

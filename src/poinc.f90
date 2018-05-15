!  The general routine to do poincare map for any given vector field (+ energy gr_cj),
!  with the poincare section specified by the subroutine 'sect'. 

! 2016-10-13 12:07:02 
!  remove the plot within poinc, we can do it seperately after the this routine,
!  do not try to mix too many things together

!  2016-12-02 16:54:46 -- CKD
!  the  substitude of the Bi-Section in the failure of Newton method
!  If the Newton method fails when the curve close to the refinement point 
!  is actually quadratic as a function of the control variable 
!  Use bi-section method instead, introduce dt_mid as the distance to the middle point,
!  take the sign of dh calculated by Newton method.  


! *******************************************************************************    
!  Determination of the imax-th return map of an orbit  to the Poincare section 
! 
!  **Note:  only one intersection is record here! 

!    To make this subroutine more general, we 
!    1. use deriv as the vector field, with a common format:  deriv(T,X,N,F)
!    2. specify the Poincare section by the subroutine 'sect', pv(ind) = p0 
!    3. TODO: compute the differential of Poincare map w.r.t. (X,T)... 
!       the initial state and arrivial time of the imax-time return map to the Poincare section 
!    3. save the orbit in file ./dat/ob_poinc.dat in debug mode(debug == 1) 
!    4. control the direction of velocity that passes the Poincare section 
!       to avoid 2 intersections within 1 revolution if we only want to keep the one with certain direction 

!    4. TODO: note it's hard to do a general routine, too many factors to consider
!    5. TODO: we put the planar RTBP also as a possible application of this routine 

!       for example:
!          1- too small the period
!          2- too much time for the next intersection
!          3- wrong detection of the intersecion at the very first steps
!          4- appropriate initial step size to avoid the above situations.

!     Proposed solution: 
!          introduce paramter tmax = 2*prd(for p.o.s)  for the time control, take 
!          the initial step size to be min(1.d-3, tmax/100)


!   TODO: should we detect also if the orbit (position component) is out of some domain??? or only control the elasped time...

! If we only want to obtain the Poincare map, n=6 is enough.  but if we want to compute the differential of the Poincare map, n has to be equal to 42  

! *****  Input Variables ******
!  sti(*)        dimension n, the initial state to start with 
!                ** NOTE **, the variational matrix should be initialized before the call of poinc 

!  ndim          dimension of the phase space, 6: spatial case; 4: planar case 
!  nvar          = ndim, if we only need to compute the phase space ; 
!                  ndim*(ndim+1), if we also want to compute the variational matrix 
!  tdir          time sence for integration, 1: forward (WU), -1: backward(WS) 


!  ****** Output Variables ******
!  tf            time spent by the orbit to go from sti to stf
!  stf(*)         imax-th cut of the orbit which passes by sti with the surface of section

!  hminim        minimum step used in the integration of the p.o., useful?
!                TODO: check if hminim is smaller than the lower bound of the integration stepsize hmin  
!  ispc          flag of the success to return to the poincare section

!  ***** Pass by module poinc_mod ******
!  ind, p0       Parameters that define the Poincare section pv(ind) = p0
!                for example, to get z=0 as the poincare section, we take ind=3, p0=0.d0
!  dir           the direction of the velocity considered when cross Poincare section 
!                1: v(ind)>0; -1: v(ind)<0; 0: keep both 
!  imax          the number of time of the intersecion with the Poinare section to consider 
!  tmax          maximal time for the imax-th intersecion, if exceeds this value, fail to get the Poincare map 


! function used: 
!       deriv, gr_cj,  sect(self-included) 
 
! Version 0.1                  -  2016-08-05 19:50:46 
! Finally revised by Yu        -  2016-10-13 12:09:11 

! TODO 2016-12-14 12:19:39   
! -- add the integration direction: tdir: 1-unstable, -1-stable 
! --------------------------------------------------------------------------

subroutine poinc(sti, ndim, nvar, tdir, tf, stf, hminim, ispc, deriv, gr_cj) 

use dp_mod 
use poinc_mod, only : ind, p0, dir, imax, tmax, xmax, hmin, hmax, e
implicit none  

! Input and Output Declaration  
integer, intent(in)   ::  ndim, nvar, tdir  
integer, intent(out)  ::  ispc

real(kind=dp), intent(in)  :: sti(nvar) 
real(kind=dp), intent(out) :: tf, stf(nvar), hminim

!  External subroutines
external ::  deriv, gr_cj
  
! Local varaibles
integer        ::  i, iter, debug 
real(kind=dp)  ::  g, gi, ti, dt_mid, dg(nvar), t, dh, tol, xm_max, tfbi, gfbi, &     ! sect + Newton method 
                   y(nvar), yvf(nvar), r(13,nvar),b(nvar), f(nvar), h, cj, hmin0   ! gr_rk78 
  
  ! TODO: set back to 0 after debug 
  debug = 0 ! Debug mode, for unchecked routines 
  
  if(debug == 1)  then  
    print*, 'check tdir: ', tdir; print*; read*
  endif 
    
  ispc = 1 ! default value is 1, if fails update to 0
  call sect(sti, nvar, ind, p0, g, dg)
  
  ! Tried different values: 1.d-3, 1.d-4, take 1.d-3 since the accumulated error
  ! is relatively smaller.
  ! h inlcude the time sense with sign \pm
  h = tdir*dmin1(1.d-3, tmax/1.d2)
  
  ! since hmin in gr_rk78 may be modified, use hmin0 to keep the original value
  hmin0 = hmin 
  
  ! the tolerance of the intersection, 1.d-13 from Gerard's routine
  ! TODO: small? big??
  
  tol = 1.d-13  
  
  ! initial time, state and step size to start 
  t = 0.d0 
  y = sti 
  
  ! the minimal step size 
  hminim = dabs(h) 
  
  if (debug == 1) then 
    print*, 'Check Poinc: hmin, hmax, e, h, t, tmax, xmax', &
            hmin, hmax,e, h, tmax, xmax, imax 
    print*, 'ndim, nvar', ndim, nvar             
    print*; read*
  endif
  
      
  do i = 1, imax
  
    ! A trick to avoid the wrong decect of crossing at the first step,
    ! treat as 0 if st(ind) is less than 1.d-9
    ! TODO: what is the appropriate value for this lower bound?? how many first steps to avoid??? 
    
    if(dabs(g) .lt. 1.d-9)  g = 0.d0 
    
    ! --- look for the next cut through the Poincare section  --- 
    do 
      ti = t;  gi = g ! the previous t and g
    
      if(debug == 1)      print*, t, h, y(1:ndim) 
          
     ! before gr_rk78, so do not save an extra point on the other side of the Poincare section 
      call gr_cj(y(1:ndim), cj)
      call gr_rk78(t, y, nvar, h, hmin,hmax,e, r,b,f, deriv )
      
      
      if(dabs(h) < hmin) then 
        print*, 'Poinc: Too small h need! Failure!'
        ispc = 0; hminim = 1.d-15 
        return 
      endif   

      hminim = dmin1(hminim, h) 
      
      if(dabs(h) < hmin0) then 
        print*, 'Too small step size needed for integration! Stop poinc' 
        print*, 'hmin = ', hmin0, 'h=', h; read*
        ispc = 0
        return
      endif 
      
      xm_max  = maxval( dabs(y(1:ndim/2) )  )
      if( xm_max > xmax ) then
        ispc = 0
        print*, 'Orbit escape!', y(1:ndim/2), xmax
        return
      endif
      
      call sect(y, nvar, ind, p0, g, dg)
      
! The detection of the next intersection fails if it spends too much time, 
! use as practical constraint.

! TODO: On the other hand, if the elapsed time is too short, 
!  we will get very small period or it is a wrong detection, which is not of
!  practical meaning, discard as well. How to handle this?

! Do not consider p.o. with period less than hmin... not too much pratical meaning.... 
      if (dabs(t) > dabs(tmax) ) then 
        ispc = 0
        print*, 'Maximum time exceeded!, t>tmax', t, tmax
        return
      endif
    
      ! terminate when detect different signs of g, cross the Poincare section in the prescribed direction 
      
      if (gi*g < 0.d0 ) then 
        
        if( dir == 0) exit  ! no requirement in the direction 
        
         ! when the velocity of the crossing point is along the same direction as dir, 
         ! start refinement, otherwise continue the integration
         
         ! the direction of the velocity has the same sign as g, which means
         ! if g>0, we goes from gi<0 --> g>0, with positive velocity
         ! otherwise, we goes from gi>0 --> g<0 with negative velocity
          
        if(  g * dir > 0.d0)  exit
        
      endif 
        
    enddo
    
  enddo  
 
  if(debug == 1) then 
    print*, 'Check the cut of the Poincare section at t, gi, g', t,  gi, g      
    read* !ckd
  endif  
 
 
! --- REFINEMENT OF THE INTERSECTION POINT stf(*) USING THE NEWTON'S METHOD --
!  If we detect two consecutive points with different signs, for sure there exists 
!  a crossing in between, so dh is of same sign with (ti-t), and with smaller modulus.
!  we have |dh| < |ti-t| and dh*(ti-t) > 0 simultaneously

!  TO GET A ZERO OF THE FUNCTION G (SEE SUBROUTINE SECT), precision is 1.d-13 ?

!  Newton method may fail when the refined variable is quadratic function of the correction 
!  Two possible cases might appear
!  -- if the current point (t, g) is on the right half of the quadratic curve, 
!     we predict along the wrong direction, dh*(ti-t) < 0
!  -- if (t,g) is close to the critical point of the quadratic curve, we might get 
!     dh in the right direction but with too big modulus.

! -- either of the above cases occurs, we use bisection method instead 
!    so we have to update the two ending points ti-t, no matter Newton method works or not.
  
  
  iter = 0
 
  do  
    iter = iter + 1
 
    if (debug == 1)  print*, iter,'-th iteration for refinement of Poincare map'
   
    ! Sometimes we need to use Bi-section method as substitude, which is slower in 
    ! convergency, so allow more iterations, 12 as default
    if (iter > 12 ) then 
      print*, 'Iteration exceeds 12  for Poincare map! Exit! ';  read*
      ispc = 0
      return
    endif 
   
 !  if (dabs(g) .le. 1.d-13) exit   !  1.d-13 is from Gerard's routine
    if (dabs(g) .le. tol)  exit
    
    call deriv(t, y, nvar, yvf)
    call dt_poinc( g, dg, yvf, nvar, dh)  ! ckd 
    
    ! ti is updated as one ending point that is on the other side of current g
    ! ti might be used without being initialized
    dt_mid = (ti-t) / 2.d0
      
    if( dh*dt_mid < 0.d0 .or. dabs(dh ) >  dabs(dt_mid) * 2.d0) then 
      ! ** NOTE** Newton method fails, use bisection method instead 
      !           take middle point, we need to update ti and t carefully every step
      if(debug==1) print*, 'Newton method fails! Use bisection method'
      dh   = dt_mid
    endif  
      
    tfbi = t;  gfbi = g   ! the other ending point  
    call gr_rk78(t, y, nvar, dh,  dabs(dh/10), dabs(dh), e,r,b,f, deriv )
    call sect(y, nvar, ind, p0, g, dg)
    
    ! take ti as the point that is on the other side of current point (t,g) 
    if( gfbi * g < 0.d0) then 
      ti = tfbi; gi = gfbi  
    endif   
     
    if(debug == 1) then 
      print*, iter, 't-g-dh-ti-t-dt_mid:', t, g, dh, ti, t, dt_mid;  print*; read*   
    endif
      
  enddo 
 
! we get y=0, t is tp, sti is the intersection      
  stf = y
  tf = t
  
  call gr_cj(y(1:ndim), cj)
  
  if(debug == 1) then 
    print*,  'Poinc finished, g, tf, stf:',  g, tf
    print*,  stf(1:ndim);  print* 
  endif 

  return
end subroutine poinc 
 
!*****************************************************************************
!  the surface of section, defined by g =  y(ind) - p0
!       Input parameters:
! y(*)        the state vector 
! n        the dimension of the state y, 6: position+velocity; 42: + variational equation
! ind, p0     Poincare section: y(ind) = p0 

!        Output parameters:
! g           funciton that is equated to 0 gives the surface of section
! dg(*)       gradient of function g


!  Global Variables from module  
! ind              theg index of the component of y to be used as Poincare section  
! sec              the Poincare section to be specified as y(ind) - sec

! 20160218 -by Yu
! --------------------------------------------------------------------------
subroutine sect(y, n, ind, p0, g, dg)

implicit none 
integer, parameter :: dp = kind(1.d0)
  
integer, intent(in)        :: n, ind
real(kind=dp), intent(in)  :: y(n), p0 
real(kind=dp), intent(out) :: g, dg(n)

integer :: i

  g = y(ind) - p0
   
  do  i = 1, n
    dg(i) = 0
  enddo 

  dg(ind) = 1.d0

  return       
end subroutine sect

!***********************************************************************
!     ****   dt_poinc   ****
! The change of arrival time needed for the required return to the Poincare section. 
! See Carles' note:  On the Analytical and Numerical Approximation of Invariant Manifolds, P13
!    detect the sign change of the flow \varphi(t) along the integration, 
!    the intersection point is refined by Newton's method

!    dt = -g / ( \Delta g( \varphi(t, X) ), f( \varphi(t, X) )
!    where \varphi is the flow with the initial conditon X,  and t is the elapsed time 
!      f is the vector field at current point (\varphi(t, X) 
!      \Delta g( \varphi(t, X) ) is the differential of g w.r.t. the state of current point, obtained from 'sect'

!    the above equation is obtained by the linear approximation of g( \varphi( t+dt, x1 ) ) w.r.t the time t.

!       Input Variables 
! g        function that equals to 0 gives the surface of Poincare section       
! dg       the gradient of g w.r.t. the state 
! f        vector filed at the evaluated point 
! n     dimension of the state, 6: pos+vel; 42: pv + variatioal matrix 

!  Comment: g, dg are computed  by subroutine 'sect'

!       Output Variables 
!  dt      time change that is needed for the refinement of the intersection        

!  Routine Used:
!    None  

!  Finally Revised by Yu -- 20160807 
!***********************************************************************
subroutine dt_poinc( g, dg, f, n, dt)
implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration  
integer, intent(in)     ::  n
real(kind=dp), intent(in)      :: g, dg(n), f(n)
real(kind=dp), intent(out)     :: dt 
 
! Local Variable
integer :: i 
real(kind=dp)  :: dy 
  
  dy = 0
  do i = 1, n
    dy = dy + f(i)*dg(i)
  enddo
    
  dt = - g/dy

  return  
end subroutine dt_poinc
  
 
!!***********************************************************************
!!     ****   diffpoinc   ****
!! The general routine to numercially compute the differential of Poincare map   
!! See Carles' note:  On the Analytical and Numerical Approximation of Invariant Manifolds, P14

!! We consider a simple situation: the standard first(imax-th) return map to a surface of section  Sigma

!! Starting at x0 \in Sigma (the Poincare section), one defines P:  Sigma --> Sigma by 
!! P( x0 ) = xf = \varphi( t(x0), x0 ) such that \psi(t) = g( \varphi( t(x0), x0 ) ) = 0
!! the arrival time depends, of course, on the initial point. The refinement of the
!! intersection is done by routine 'dt_poinc'
! 
!! let xf \in Sigma be the arrival point. Then, by differentiation w.r.t. x0 one has 

!!  d P / d x0 = f(xf) d t / d x0 + d \varphi( t(x0), x0 )  / d x0 

!!  where d t / d x0  is computed by the differentiation of g( \varphi( t(x0), x0 ) ) = 0
!!  and we have 
!!     D g(xf) * ( f(xf) dt / d x0 + d \varphi / d x0 ) = 0

!!     where f is the vector field at current point (\varphi(t, x0),  
!!     D g(xf) is the differential of g w.r.t. the state of the arrival point 

!!d P / d x0 = - f(xf) / ( \nabla g(xf), f(xf) ) D g(xf) * d \varphi / d x0 + d \varphi / d x0

!! if we take g(xf) = xf(ind) = 0, then the donominator ( \nabla g(xf), f(xf) ) = f(xf)_ind
!! We assume the matrix of the first variational equations is given by d \varphi / d x0  = a_ij

!! then the (n-1)-by-(n-1) matrix d P / d x0 | (x0)  restricted to the n-1 components expect ind-th one is given by 

!!  d P / d x0 |k,j =  a_k,j - f(xf)_k / f(xf)_ind * a_n,j ,  1<=j,k<=n && j,k /= ind 

!! In order to make it more general, we only compute the necessary components of d P / d x0 


!!       Input Variables 
!! phi       variational matrix at time t, (State transition matrix)   
!! f         vector filed at the arrival point  xf, dimension 6, velocity+acceleration
!! ind       the index of the component of xf set to be zero to define 
!!           the Poincare sectionxf(ind) =0 
!! nr, nc    dimension of dpdx, nr-by-nc 
!! para      index of the free parameter, dimension nc  
!! fun       index of target fcuntions, dimension nr
!!            
!! ***Comment***: for asymmetric p.o., para = fun = [1,2,3,4,5] if ind=6
!!                for symmetric p.o., it depends, the values may vary a lot
! 
!!       Output Variables 
!!  dpdx     differential of Poincare map associated with target components specified by fun,
!!           w.r.t. the parameters specified by para.  dimension nr-by-nc
! 
! 
!!  Routine Used:
!!    None  

!!  Finally Revised by Yu -- 20160808 
!!***********************************************************************
!subroutine diffpoinc( phi, f, ind, nr, nc, para, fun, dpdx)
!implicit none
!integer, parameter :: dp = kind(1.d0)   

!! Input  and Output Declaration  
!integer, intent(in)     ::  ind, nr, nc, para, fun 
!real(kind=dp), intent(in)      ::  phi(6,6), f(6)  
!real(kind=dp), intent(out)     ::  dpdx(nr, nc)
! 
!! Local Variable
!integer :: i, j
!  
!  dpdx = 0.d0
!  
!  ! compute the component one-by-one, the component of dpdx with index i,j corresponds to 
!  ! the one specified by fun(i) and para(j) in phi 
!  do i = 1, nr 
!    do j = 1, nc 
!      dpdx(i,j) = phi( fun(i), para(j) ) -  f( fun(i) ) / f(ind) * phi( ind, para(j) ) 
!    enddo 
!  enddo  

!  return  
!end subroutine diffpoinc  

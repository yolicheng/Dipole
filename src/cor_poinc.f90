! 2017-02-20 09:23:34  
! take only subroutine fctn, and discard poinc and sect within this file, instead use the ones from poinc_mod, more sophiscated ones with less bugs. 

! 20160411 
!  Final update:  This file contains is to refine p.o. with poincare map approach, for either symmetric one or asymmetric one. 
!                 The general approach in which the correction is all the whole state vector + the peroid is done is gr_pomod

!********************************************************************************
  subroutine fctn(x,init, tmax, f,g, tf,y,vf, ispc, deriv, gr_cj)
  
  implicit none  
  integer, intent(in)  :: init 
  integer, intent(out)  ::  ispc
  real(kind=dp), intent(in)  :: x(6), tmax
  real(kind=dp), intent(out) :: f(ntar), g(ntar, nctr), vf(42), tf 
  real(kind=dp), intent(inout)  :: y(42)
  
! local varaibles
  integer        ::  i, j 
  real(kind=dp)  ::  phi(6,6), yi(42), hminim !, f1, f2 !, & 
!                     phinm(6), phinormi, gnm(ntar), gnormi, cj, &
!                     vftar(ntar,1)  
    
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

!subroutine poinc(sti, ndim, nvar, tdir, tf, stf, hminim, ispc, deriv, gr_cj)  --from  poinc.f90

!  subroutine poinc(yi,imax, tmax, neq, tf,yf, hminim, ispc, deriv, gr_cj) 
  call poinc(yi,imax, tmax, 42, tf,y,  hminim, ispc, deriv, gr_cj) ! imax determines the number of crossing we are considering as
      
  if(hminim < hmin ) then 
     print*, 'hminim < hmin', hminim, hmin;  read(*,*) 
  endif 
    
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
 ! remember not to forget to abstract the identity matrix  ! ckd!!! 
  if (issymt == 0 ) then 
    print*, 'Asymmteric case! Phi= Phi-I' !ckd!
    if(debug == 1) read*
    
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
 
    print*, 'Jacobi Matrix of target variable w.r.t control variables' 
    do i = 1, ntar 
      write(*,*)  phi(tar(i), ctr)
    enddo 
    read*

  endif 
  
  return
  end subroutine fctn
  

!*******************************************************************************    
  subroutine poinc(yi,imax, tmax, neq, tf,yf, hminim, ispc, deriv, gr_cj) 
! Determination of the imax-th passage of the orbit through the Poincare section defined by subroutine sect(here,y=0) 
!  Note:  only one intersection is record here!
! to make this subroutine more general, use deriv as an argument as  the function to compute the vector field
!   of the common form:  deriv(T,X,N,F)
     
!       Input Variables
!  yi(*)      initial point  
!  imax            the number of time for intersecting the section
!  neq             dimension of the state: 6--position+velocity  42: ...+ variational matrix

!       Private Module-based Varaible
!  tmax, hmax, hmin, e, dp

!        Output Variables
!  tf               time spent by the orbit to go from yi to yf
!  yf(*)        first cut of the orbit which passes by yi with the surface of section
!  hminim      minimum step used in the integration of the p.o.
!  ispc       flag of the success to return to the poincare section

! function used: gr_rk78, sect  
! --------------------------------------------------------------------------

  implicit none  
!  integer, parameter  ::  neq = 42 ! dimension of the state: 6-state vector, 42-also the variatioal matrix
  
  integer, intent(in)  ::  imax, neq
  integer, intent(out)  ::  ispc
  real(kind=dp), intent(in)  :: yi(neq), tmax 
  real(kind=dp), intent(out) :: tf, yf(neq), hminim
  external deriv, gr_cj
  
! Local varaibles
  integer :: i, iter
  real(kind=dp)  ::  g, gi,  dg(neq), t, dh, dy, &
                     y(neq), r(13,neq),b(neq),f(neq), h  ! gr_rk78 
  
  ispc = 1 ! default value is 1, if fails, set to 0
  call sect(yi, neq, g, dg)
  
! initial state
  h = 1.d-3
  t = 0.d0 ! initial time
  hminim = h
  y = yi 
  
  if (debug == 1) then 
    print*, 'Poinc- hmin,hmax,e, h, t, tmax, imax', hmin,hmax,e, h, t, tmax, imax
    print*;    read*
  endif
  
      
  do i = 1, imax
  
    ! A little trick to cancel numberical error, if the compoent of the state is less than 1.e-9, treat as 0 
    if(dabs(g) .lt. 1.d-9)  g = 0.d0 
      
    do  !look for the next intersecion across the Poincare section 
      gi = g ! the previous value of g
    
      call gr_rk78(t,y,neq, h,hmin,hmax,e, r,b,f, deriv )
     
      if (debug == 1)   write(12, '(8f18.14)') t, y(1:6), h ! print the correction process on the screen
      
      call sect(y, neq, g, dg)
      
! if it spends too much time to go to the next intersecion across y(ind)=sec plane, seen as failed 
! otherwise, if the elapsed time is too short, we will get very small period, which is not of too much point, discard as well.
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
    print*, 'ck crossing the Poincare section at t', t,  gi, g      
    read(*,*) !ck
  endif  

!   REFINEMENT OF THE INTERSECTION POINT YF(*) USING THE NEWTON'S METHOD
!   TO GET A ZERO OF THE FUNCTION G (SEE SUBROUTINE SECCIO), precision is 1.d-13 ?
! pay attention here, there is possibility that it falls into endless iteration here
  iter = 0
  if (debug == 1) print*, 'Iteration for poincare map! t - g - dh'
 
 
  do  
  
    iter = iter + 1
 
    if (debug == 1)     print*, iter,'-th iteration for intersection with poincare section'
   
    if (iter > 6 ) then 
     print*, 'Iteration exceeds 6 !'; read*
     ispc = 0
     return
    endif 
   
    
!    if (dabs(g) .le. 1.d-13) exit  ! detect the y(ind)-sec is within the tolerance, values 1.d-13 from Gerard's routine
    if (dabs(g) .le. tol_poinc) exit
    
    call deriv(t, y, neq, f)

    dy = 0
    do i = 1, neq
      dy = dy + f(i)*dg(i)
    enddo
    dh = - g/dy
    
    if (dabs(dh ) > hmax) then 
      print*, 'Too big step dh needed, dh=', dh, '>', 'hmax=',hmax
      ispc = 0 ; read* 
      return 
    endif 
    
   if(debug == 1)  print*, 'before gr_rk78, dh', dh
! to make sure this step works, we need to set hmin to the required stepsize 
! Need to debug here, for possible stuck here....
   
    call gr_rk78(t,y, neq, dh, dabs(dh),hmax,e, r,b,f, deriv) ! remember that gr_rk78 is only 1 step
    
!    call gr_rk78(t,y, neq, dh, dabs(dh), dabs(dh),e, r,b,f, deriv) ! remember that gr_rk78 is only 1 step
    
    call sect(y, neq, g, dg)
    
    if(debug == 1) then 
      print*, iter, ':', t, g, dh    
      read*   !ckd
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
! the surface of section, defined by g =  y(ind) - y0
!       Input parameters:
! y(*)            the state vector 
! neq             the dimension of the state y, 6: position+velocity; 42: + variational equation

!        Output parameters:
! g             funciton that equated to 0 gives the surface of section
! dg(*)       gradient of function g


!  Global Variables from module  
! ind              theg index of the component of y to be used as Poincare section  
! p0              the Poincare section to be specified as y(ind) - sec

! 20160218 -by Yu
! --------------------------------------------------------------------------
  subroutine sect(y, neq, g, dg)
  implicit none 
  
  integer, intent(in) :: neq
  real(kind=dp), intent(in) :: y(neq) 
  real(kind=dp), intent(out) :: g, dg(neq)

  integer :: i

  g = y(ind) - p0
   
  do  i = 1, neq
    dg(i) = 0
  enddo 

  dg(ind) = 1.d0

  return       
  end subroutine sect
  
  

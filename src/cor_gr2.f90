! the problem with this approach lies in the continuation 
! for the computation of the kernel  v = DF, if we keep the ind-th component as one of the control variables, 
!  then the last row of DF = [1 0 0 0 0 0 0 ], question: the last element shoule be dx(ind) or 0? with 0 fails to compute the det, because it is equal to
!  try dx(ind) 


! ****************************************************************************** 
! The initializaiton for asymmetric p.o.:  idcor = (4,5,6) 
!  Initialize the index of the control variables and target variables as in the state vector

! Module based varaibles
!  tfree, idcor (4,5,6)

! TODO  idcor = 4,6,7 ,  only finish  5
! ****************************************************************************** 
  subroutine init_asymtpo2(idcor0, ind0, x0_fxd0)
  
  ! this is try to add a constraint to fix a component of the state vector, y(ind0)= x0_fxd0
  ! so the control variables are still all the six components, but the number of target function increases by 1 
  
  integer,intent(in) :: idcor0, ind0
  real(kind=dp),intent(in) :: x0_fxd0
 
  
  integer ::  allocatestatus, nrow, ncol, i, j


!  module-based variables  
  idcor = idcor0 ! supposed to be 5 
  issymt = 0 
  tfree  = 0 ! by default, fixed-time shooting
 
 
! ! --- allocate memory for ctr and tar  for asymmetric p.o.  
  if( idcor == 4 ) then ! TODO
  ! poincare map approach, by default, we search for the first crossing through the poincare section
    imax = 1
    ind = ind0
    sec = x0_fxd0 
    
    ntar = 4
    nctr = 5
    
    nrow = ntar
    ncol = nctr 
  elseif( idcor == 5 ) then 
  !  free time shooting + additional constraint
    ind = ind0 
    x0_fxd = x0_fxd0
    nctr = 7
    ntar = 6
    
    tfree = 1
    nrow = ntar - 1
    ncol = nctr - 1 
    
    
  elseif( idcor == 6 ) then 
  !  free time shooting -- no additional constraint....
    nctr = 7
    ntar = 6 
     
    tfree = 1    
    nrow = ntar
    ncol = nctr - 1 
    
  elseif (idcor == 7) then 
  ! free time shooting, deliminate ind from tar to satisfy x(ind) = x0_fxd 
    nctr = 6
    ntar = 5
    
    tfree = 1 
    nrow = ntar 
    ncol = nctr-1 
  end if
  
    
  ntp = 1
 
 
 ! allocate the memory for tar and ctr 
  allocate ( tar(nrow), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"

  allocate ( ctr(ncol), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"
  
  
 ! assign the vaules for tar and ctr 
  if ( idcor == 5  ) then ! the current used one, additional constraint: y(ind) = x0_fxd
    tar = (/1,2,3,4,5/) 
    ctr = (/1,2,3,4,5,6/) 
  
  elseif( idcor == 6  ) then ! --TODO with full state vector as the control variables    
    tar = (/1,2,3,4,5,6/) 
    ctr = (/1,2,3,4,5,6/) 
! 
  elseif( idcor == 7)  then  ! --TODO, to approve idcor = 5 
    tar = (/1, 2, 3, 4, 5 /)
     
!   instead of an additional equation, fix x(ind)=x0_fxd as the initial condition, as eliminate ind from ctr 
    j = 0
    do i= 1, 6, 1
      print*, 'ind=', ind, 'i=', i 
      read*
      
      if( i /= ind ) then 
        
        j = j + 1
        ctr(j) = i
      endif 
      
    end do
     
!  elseif( idcor == 4) ! --TODO, poincare approach.... 
      
  end if  
  
  
  print*,' ind, x0_fxd', ind, x0_fxd
  read*
  
  print*, 'idcor=', idcor, 'nctr=', nctr, 'ntar=', ntar
  print*, 'ctr=', ctr, 'tar=', tar 
  read* 
  
  end subroutine init_asymtpo2
  
! This file contains the subroutines that are realted to general approach 
! 1.  gr_fctn2

! this is to test with only five constraint equations, how to add another constraint, 
! fix a component? energy? period? 


! For the asymmetric approach, please note that because we are refine w.r.t the initial condition
!   So the Jacobi Matrix w.r.t the correction on the initial state is actually Phi - I  

! Instead of Poincare map approach, we use the variational matrix to do the refinement of the asymmetric p.o., such that we have 
! 5 equations and 6 unkonwns(period included). We are not using the modified Jacobi Matrix in which the differential of the Poincare map w.r.t time t is represent by the other unkonwns in order to decrease the number of unknowns by 1 

! the 7 unkowns are (xi, yi, zi, vxi, vyi, vzi, T) 
! and the 6 constraint equations are satisfied to return to the same initial point 
!  (xfi-xi=0, yf-yi=0, zf-zi=0, vxf-vxi=0, vyf-vyi=0, vzf-vzi=0)
 
! appropriate for issymt = 3
subroutine gr_fctn2(y0, tp0, f, g, yf, deriv, gr_cj )

! module based Varaibles   
! ntar, nctr, tar, debug 
  
  implicit none  
  integer, parameter :: dp = kind(1.d0)   
  real(kind=dp), intent(in)   ::  y0(6), tp0
  real(kind=dp), intent(out)  ::  g(ntar, nctr), f(ntar),  yf(42) 

  external ::  deriv, gr_cj
  
! Local Varaibles
  integer        ::  i, ispl, ftag 
  real(kind=dp)  ::  phi(6,6), yi(42), vf(6)  

! Integrate for an approximate period of tp0 
  yi = 0.d0
  yi(1:6) = y0
  yi(7:42:7) = 1.d0 

! In debug mode, save the integration data of the orbit, other wise, do not save the data of the orbit     
  if(debug == 1) then 
    ispl = 1
    ftag = 12  ! 12 is opened before calling the dfcn 

  else 
    ispl = 0
    ftag = 0 
  endif  
  
  
!   subroutine plob(y0,t0,tf, n, tdir, ftag, ispl,  deriv, gr_cj,  y)   
  call plob(yi, 0.d0, tp0, 42, 1, ftag, ispl, deriv, gr_cj,  yf)  
    
  if(debug == 1) then 
    print*, 'tf, yf'
    print*, tp0, yf(1:6)
    print*;  read*
  endif   
  
! compute the Jacobi Matrix of F w.r.t (X0,T)
!   G = [Phi, VF]
  phi = reshape( yf(7:42), (/6,6/) ) ! state transition matrix
   
  if (debug == 1) then  
    print*,'Phi'
    do i = 1, 6
      print*, phi(i,:)
    enddo
    read*
  endif  
  
  do i = 1, 6
    phi(i,i) = phi(i,i) - 1.d0
  enddo 
   
 ! the first 5 rows are component of the state,  the first 6 columns
  g(1:5, 1:nctr-1) = phi(tar, ctr)
! this is to modify, pay attention how many free variables to keep at last.... 
 
 ! the sixth row = x(ind) - x0_fxd   
   g(6, ind) = 1.d0
    
    
! subroutine gr_lf(t, y, neq, f)
  call deriv(tp0, yf(1:6),  6, vf)
  
  ! the first 5 rows   
  g(1:5, nctr) = vf(tar) ! the 7-th columns is the vector field at epoch tp0
  
  ! the sixth row 
  g(6, nctr) = 0.d0 
  
  if (debug == 1) then    
    print*,'Jacobi Matrix: G'
    do i = 1, ntar
      print*, g(i,:)
    enddo
    read*
  endif  
 
 ! the target function,  y(ind) is fixed to be  the initial value ---- x0_fxd
  f(1:5) = yf(tar) - y0(tar)
  f(6)   = y0(ind) - x0_fxd  ! the point is how to set the fixed value of y0(ind) here? 
  
  if (debug == 1) then   
    print*, 'ind=', ind, 'x0_fxd=', x0_fxd 
    print*,'Constraint function: F'
    print*, f
    read*
  endif 
  
end subroutine gr_fctn2   
  
  
   

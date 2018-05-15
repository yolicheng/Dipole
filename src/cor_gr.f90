! ****************************************************************************** 
! The initializaiton for general approach to refine p.o.: 
! Initialize the index of the control variables and target variables as in the state vector
!  ** At this moment, I didn't use this, use init_asympo instead

! ****************************************************************************** 
  subroutine init_grpo(idcor0) 
  
  integer, intent(in) :: idcor0 
  
  integer ::  allocatestatus
  
  idcor = idcor0 ! supposed to be 6 
  
  issymt = 0 ! asymmetric
  tfree = 1  ! free-time shooting
  
  ntp = 1
! --- allocate memory for ctr and tar 
  nctr = 7
  ntar = 6 
  
  allocate ( tar(ntar), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"

  allocate ( ctr(nctr-1), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"
  
  tar = (/1,2,3,4,5,6/)
    
  ! the control variables in fact have 7 components, X0(6) + TP,
  ! here we only need the index of the ones that are the component of the state  
  ctr = (/1,2,3,4,5,6/) 
  
  end subroutine init_grpo 
  

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

! module based Varaibles   
! ntar, nctr, tar, tfree,  debug 

! approapriate for all cases, except poincare map approach 
  
  implicit none  
  real(kind=dp), intent(in)   ::  y0(6), tp0
  real(kind=dp), intent(out)  ::  g(ntar, nctr), f(ntar),  yf(42) 

  external ::  deriv, gr_cj
  
! Local Varaibles
  integer        ::  i, j, ispl, ftag, ncol, nrow 
  real(kind=dp)  ::  phi(6, 6), yi(42), vf(6)  

  print*, 'check the symmetry and approach, issymt=', issymt, 'idcor=', idcor, 'tfree=', tfree
  if(debug == 1)  read*
  
! number of column in G w.r.t the state equals to the size of ctr 
! ncol = nctr for fixed-time shooting, while ncol = nctr -1 for free-time shooting
  ncol = size(ctr)

! number of rows in G w.r.t the state equals to the size of tar 
  nrow = size(tar) 
  print*, 'tfree=', tfree, ' ,nctr=',nctr, ' ,ncol=', ncol, 'nrow=', nrow
!  read*
  
  
! Integrate for an approximated period of tp0 
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
  
!  subroutine plob(y0, t0, tf, ndim, nvar, tdir, ispl, ftag, deriv, gr_cj,  y) 
  call plob(yi, 0.d0, tp0, 6, 42, 1, ispl, ftag, deriv, gr_cj,  yf)  
    
  if(debug == 1) then 
    print*, 'tf, yf'
    print*, tp0, yf(1:6)
    print*;  read*
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
  
  ! If idcor == 7, we remove the ind-th column, 
  g(1:nrow, 1:ncol) = phi(tar, ctr)
  
! Pay attention to the difference between symmetric and asymmetric p.o., and if it is free-time shooting 

! for asymmetric case, refine to the initial condition 
  if( issymt == 0 ) then
    do i = 1, nrow, 1
      do j = 1, ncol, 1
        if( tar(i) == ctr(j) ) then
          g(i,j) = g(i,j) - 1.d0
        end if
      end do
    end do
  end if  
  
  
! in free-time shooting, the last colomn is the velocity 
  if( tfree == 1) then 
   
    ! subroutine gr_lf(t, y, neq, f)
    call deriv(tp0, yf(1:6),  6, vf)
  
    ! the last columns is the vector field at epoch tp0   
    g(1:nrow, nctr) = vf(tar) 
    
  endif

  ! for additional constraint of fixing one of the component in state vector
  ! as the last row in constraint matrix equation
  if( idcor == 5 ) then
    g(ntar, ind) = 1.d0
    f(ntar) = y0(ind) - x0_fxd
  end if 
     
  if (debug == 1) then    
    print*,'Jacobi Matrix: G'
    do i = 1, ntar
      print*, g(i,:)
    enddo
    read*
  endif  
 
  ! the target functions that are related to the state 
  if( issymt == 1 ) then
    ! for symmetric p.o., the target varaibles should be zeros
    f(1:nrow) = yf(tar) !- 0.d0
  
  else
    ! for asymmetric p.o., refine to the initial condtion 
    f(1:nrow) = yf(tar) - y0(tar)
  end if
 
  
  if (debug == 1) then   
    print*,'Constraint function: F'
    print*, f
    read*
  endif 
  
end subroutine gr_fctn   
  
  
   

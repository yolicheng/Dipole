!  **** My Poincare map Module*******
! TODO-- check!!!!  
!  init_poinc is duplicated with the one inside po_mod.f90 

!  put poincare map related routines here, include 
!   poinc, poinc_n,  diffpoinc 
!    commonly used parameters: ind, p0,  dir, imax, tmax, ndim 
                         
!  the Poincare section: X(ind) = p0  

!  Assume we have a 2d map that maps R^2(y, vy) to R^2 (y, vy) 
!  and vx is a function of vx = f(x0, y, vy, h0)

!  so vx is not free, the differential of the Poincare map should be w.r.t. the two unknowns (y, vy) 
!  taking the  vx = f(x0, y, vy, h0) into consideration 

!  It's suggested to write the equations in a compact way, and use matrix to express the Poincare map  
!  such that we will not miss anything

! Assume we have \hat P: (x,y,vx, vy)_0 --> (x,y,vx, vy)_f, actually the flow \phi
! introduce a map \gamma:  (y, vy)_0 --> (x,y,vx, vy)_0 
! A projection Pi (the inverse of gamma): (x,y,vx, vy)_f -->  (y, vy)_f 


! Poincare map is P:  (y, vy)_0 --> (y, vy)_f, can be expressed as a series of compositions:
  
!   Pi(   \hat P  ( \gamma( (y, vy)_0 ) )  )   )

!  so d P / d (y, vy)_0  = d Pi / d \hat P *  d \hat P / d \gamma *  d \gamma / d (y, vy)_0
! where the first term: 
!   d Pi / d \hat P (x,y,vx, vy)_f  = d Pi /  [0 1 0 0; 0 0 0 1]   ! the second and fourth rows of STM

! the second term 
! d \hat P / d \gamma = d (x,y,vx, vy)_f / d (x,y,vx, vy)_0 = STM 

! the third term --- ** NOTE ** be careful with this one, vx= f(x0, y, vy, h0) 
! d \gamma / d (y, vy)_0 = d (x,y,vx, vy)_0 / d (y, vy)_0   
!                        = [0    0; 1    0; d vx/ d y     d vx / d vy;  0    1 ]

module poinc_mod

use dp_mod   
use pi_mod 
implicit none
save 

! we are going to use these values in other routines, mainly PoincMap_XXXX
! so declare them public 

! --- public ---
real(kind=dp) :: p0, tmax, xmax, h0,  hmin, hmax, e   
integer       :: n,  ind_vel, ind, dir, imax,  isv 


! 2017-03-17 09:21:16  add flag isv, to avoid meaning computation of extra points during the curve refinement 

! we do not know the size of ind_fun in the beginning, declare allocatable                         
integer, allocatable  ::  ind_fun(:) ! poinc       



!  --- private ---                 
integer, private       :: ndim, debug 


contains

! TODO: put these two subroutines in this module! for more convenient use... 
!  include  'poinc.f90'
!  include  'poinc_n.f90'
    
!***********************************************************************
!     ****   init_poinc   ****
!  the dimension of related arrays 
!  for Poincare map: ind, p0, tmax, dir, imax  
!***********************************************************************
subroutine init_poinc(ind0, p00, dir0, imax0, tmax0, xmax0, ndim0) 
! Assignment of shared variables for poincare map-related subroutine  + rk78 (for integration)

! imax is explicitly assigned here, which is also implicitly assigned in init_asymtpo and init_symtpo
! so if we want to modify a little bit how the routine works, remember to call this one after these two subroutines init_asymtpo and init_symtpo
 
implicit none
integer, intent(in)          :: ind0, dir0,  imax0, ndim0
real(kind=dp), intent(in)    :: p00,  tmax0, xmax0 

integer :: i, k

  debug = 0
!  print*, 'For poinc, Debug (1) or not? Put 2 to check dp_dx'
!  read(*,*)  debug

!! debug = 2 
  
  ! By default, we set the flag of cj2v to successful
  isv = 1 
  
 ! Initialize for the integration by gr_rk78 (the most commonly used value) 
  hmax = 1.d-1 ! Suggested by Alex, not to be too big. 1.d0 is too big 
  hmin = 1.d-6  
  e    = 1.d-14
  
  ! -- dimension declaration
  ndim  = ndim0
!  npvar = ndim*(ndim+1) ! for variational matrix computation 
  
  n     = ndim - 2 ! the dimension of Poincare map 
 
!  ! -- energy 
!  h0 = h00 
  
  ! -- For Poincare map 
  ind  = ind0
  p0   = p00
  
  ! the velocity in the direction of pv(ind) 
  ind_vel = ind + ndim / 2  
    
  dir  = dir0
  imax = imax0 
  tmax = tmax0
  xmax = xmax0
  
!  assign the index of the component in  free parameters and functions 
!  allocate(ind_para(ndim)) 
  allocate(ind_fun(n)) 
  
!  ind_para =  (/ (i, i = 1, ndim) / ! implied do loop

  k = 1
  do i = 1, ndim, 1
!    ind_para(i) = i 
    if(i == ind .or. i == ind_vel) cycle 
    ind_fun(k) = i
    k = k+1
  end do
  
  
!if(debug == 1) then 
  print*, 'check ind, p0, tmax, dir, imax',  ind, p0, tmax, dir, imax
  print*, 'ind_fun', ind_fun
  print*, 'finish init_poinc!'
  print*;  read*  
!endif    
  
  return  
end subroutine init_poinc


!***********************************************************************
!     ****  init_h0    ****
!  Set the value of energy level seperately in case that we do the continuation 
!  by varing the energy level
!***********************************************************************
subroutine init_poinc_h0( h00)
implicit none
 
! Input and Output Declaration   
real(kind=dp), intent(in)      ::   h00  
 
 ! -- energy 
  h0 = h00 
  
  print*, 'poinc_mod, cj0= ', h0 
  print*; read*
  
  return  
end subroutine init_poinc_h0


end 



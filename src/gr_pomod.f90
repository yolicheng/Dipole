! This file contains the subroutines that are realted to general approach 
! 1.  gr_fctn

! For the asymmetric approach, please note that because we are refine w.r.t the initial condition
!   So the Jacobi Matrix w.r.t the correction on the initial state is actually Phi - I  

! Instead of Poincare map approach, we use the variational matrix to do the refinement of the asymmetric p.o., such that we have 
! 5 equations and 6 unkonwns(period included). We are not using the modified Jacobi Matrix in which the differential of the Poincare map w.r.t time t is represent by the other unkonwns in order to decrease the number of unknowns by 1 

! the 7 unkowns are (xi, yi, zi, vxi, vyi, vzi, T) 
! and the 6 constraint equations are satisfied to return to the same initial point 
!  (xfi-xi=0, yf-yi=0, zf-zi=0, vxf-vxi=0, vyf-vyi=0, vzf-vzi=0)
 
 
subroutine gr_fctn(y0, tp0, f, g, yf, deriv, gr_cj )

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

! to save the integration data of the orbit or not     
  if(debug == 1) then 
    ispl = 1
    ftag = 12 
!   
!    print*, 'Initial Phi:' !ckd-- fine!
!    phi = reshape( yi(7:42), (/6,6/) ) ! state transition matrix
!    do i = 1, 6
!      print*, phi(i,:)
!    enddo
!    read*
   
  else 
  
    ispl = 0
    ftag = 0 
  endif  
  
  
!   subroutine plob(y0,t0,tf, n, tdir, ftag, ispl,  deriv, gr_cj,  y)   
  call plob(yi, 0.d0, tp0, 42, 1, ftag, ispl, deriv, gr_cj,  yf)  
    
  print*, 'tf, yf'
  print*, tp0, yf(1:6)
  read*
  print*
  
  
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
  
  g(:, 1:6) = phi ! the first 6 columns
    
  do i = 1, 6
    g(i,i) = g(i,i) - 1.d0
  enddo 
    
    
! subroutine gr_lf(t, y, neq, f)
  call deriv(tp0, yf(1:6),  6, vf)
    
  g(:, 7) = vf ! the 7-th columns is the vector field at epoch tp0
  
  if (debug == 1) then    
    print*,'Jacobi Matrix: G'
    do i = 1, ntar
      print*, g(i,:)
    enddo
    read*
  endif  
 
 ! the target function
  f = yf(tar) - y0(tar)
    
end subroutine gr_fctn   
  
  
   

!  **** My curve Module*******
!  Compute and refine the invariant curve obtained by Poincare map method,  with the energy and rotation number are both fixedd 

!  the idea is different from the Time t map, 
!  so do a module seperately, although they have some similar part, but do not mix the two methods

! following several steps : 

!  -- 1. take a Poincare map, the section is specified by pv(ind) = p0, for example, x = x_L4 

!  -- 2. In order to get equispaced points to do Fourier analysis, 
!        -- 2.1  compute the rotation number  within [0,1] by Carles' note, multiply by 2pi to obtain rho \in [0, 2pi]
!        -- 2.2  sort all the points by the value \rho, and then do linear interpolation within [0, 2pi] for y and vy
!                the last component vx is obtained by the prescribed energy vx = f(H0, x,y,vy)
!        -- 2.3  using the general Fourier.f routine to compute the Fourier coefficients, which provide the initial seed for refinement
!                here, we need to check the approximated Fourier representation and the interpolated points 

!  -- 3. do the refinement using Newton Method by imposing the invariance equation 
!        ** NOTE **
!        3.1  there is only one indetermination in Poincare map method, the phase shift indetermination, 
!             it can be avoided by fixing C1_1 = 0, this is done by remove C1_1 from the unknows, although in the system equations
!             it is expressed by adding another equation  C1_1 = 0. 
!             since if we keep C1_1 as one unknown, then the Jacobi Matrix will be singular, which makes it difficult for Newton Method 
!        3.2  by fixing the energy H0 and the rotation number \rho, we fix a particular invariant curve 
!        3.3  the unknowns are only the Fourier coefficients: C0, CK, SK 

!  so, we have the number of unknowns : n * (2*nf+1) - 1 
!     where n is the dimension of the phase space to be evaluated,  n = 3 for PRTBP, nf is the number of Fourier modes. 
!  the number of equations 
!     n * (2*nf+1) + 1, since we have one more equation h=H0 to get specify the energy level  

 
! ** Note ** 
!    1. we have 2 more equations than unknows, but Josep Maria claims the kernel dimension is zero... 
!    2. use the general name for vector field and energy (Jacobi costant)
!    3. subroutines from libraries, need to declare external attribute **
!         external  ::  deltx, trigrec_c, trigrec_s   ! -- libgr.a 

module curve_poinc_mod

! TODO: what if dp and pi, pi2 conflict with the ones in the main routine?
use dp_mod ! dp and pi, pi2 
use pi_mod 

! dimension declaration and termination control of Newton method
!use tori_mod, only: nf, nf2, nf2p1, ndim, nitmax, tol, tol_err 

implicit none
save 

! subroutines from libraries, declare external attribute, 
! this is only necessary for module to call external subroutines or functions 
external  ::  deltx,  trigrec_c, trigrec_s   ! -- libgr.a 

real(kind=dp), private ::  tmax, p0, &  ! poinc
                           tol, tol_err ! for tori_mod TODO
                           
! --- For declaration of array ---- 
integer, private   ::  n, nc_dg, nr_dg, nc_dgall, nr_dgall, nc_dgxi1, & ! dimension
                       ind, dir, imax, & ! poinc 
                       nf, nf2, nf2p1, ndim, nitmax ! for tori_mod TODO

! --- For save -- 
! to write the output of refinment process
integer, private :: fwork, fout 


contains

!***********************************************************************
!     ****   init_curve_poinc   ****
!  the dimension of related arrays 
!  From tori_mod:  nf, nf2, nf2p1, ndim 
!  for Poincare map: ind, p0, tmax 
!***********************************************************************
subroutine init_curve_poinc(nf0, ndim0, nitmax0, tol0, tol_err0, ind0, p00, dir0, imax0, tmax0 ) 

implicit none
integer, intent(in)          ::  ind0, dir0, imax0, &
                                 nf0, ndim0, nitmax0  
real(kind=dp), intent(in)    ::  p00, tmax0, &
                                 tol0, tol_err0   

! for the tori_mod... TODO: put this part into tori_mod later
  ndim  = ndim0
  nf    = nf0
  nf2   = 2*nf  
  nf2p1 = nf2 + 1
  
  ! termination control for Newton method 
  tol     = tol0
  tol_err = tol_err0
  nitmax  = nitmax0 
  
!  print*, 'check the initialization of curve_poinc_mod module'
!  print*, 'nf, ndim, tol, tol_err, nitmax', nf, ndim, tol, tol_err, nitmax
!  read*
  
  
! ---- commonly used value for dimension declaration  
  ! number of free components of the phase 
  ! since one component is specified by Poincare map, we only need to compute the rest ones 
    
  ! number of columns and rows in dg, which is of dimension ndim * (nf2p1)+1 - by -  ndim * (nf2p1)-1
  nr_dg = ndim * nf2p1 + 1 ! one is or energy 
  nc_dg = nr_dg - 2        ! fix c1_x = 0 to avoid the phase shift indetermination
  
  ! dimension of dgall, the number of row is the same as nr_dg 
  ! number of the column is 1 less, since we need to fix c1_(ind)
  
  ! dimension of dg for one value of xi 
!  nr_dgxi1 = nf2p1      ! withou the energy equation
!  nc_dgxi1 = ndim*nf2p1 ! number of unknowns, all the coefficients for (T,y,vx,vy), each is of dimension nf2p1 

  ! For Poincare map 
  ind = ind0
  p0  = p00
  tmax = tmax0  
  dir = dir0
  imax = imax0 
  
!  print*, 'check ind, p0, tmax, dir, imax',  ind, p0, tmax, dir, imax
!  read*
  
!  for save - TODO: for the moment, not used
   fwork = 6
!  fwork = 222 
!   open(fwork, file='curve_poinc.work')
  
  print*, 'finish init_curve_poinc!'
  read* 
  
  return  
end subroutine init_curve_poinc

  
! **********************************************************************
! This routine is to refine the invariant curve, 
! with the initial guess of the coefficients computed by gr_four   
 

!**NOTE: Alex's note about the Poincare map method 
!  The tori we are looking for are embedded in 2–parametric families, which can be paramtrized by 2 parameters among h, rho   
!  The curve on the torus is fixed by thePoincare section, so we do not have curve indetermination

!  Therefore, in order to compute a torus, we  need to 
!    1. set a coordinate of C1 (S1) equal to zero and eliminate it, in order to get rid of the phase shift indetermination. 
!       Here we take C1_ind(S1_ind), ind is specified by routine  'fix_index', according to the strategy in Josep-Maria's thesis

!    2. eliminate two unknonws h and rho, in order to fix a particular curve,
!       by looking at the Poincare map representation of the torus, we can see we have bunch of 
!       invariant curve, which can be specified by the rotation number rho 
!       so, we eliminate h and rho       

!  So we end up with a   nr ( = 1+  ndim*nf2p1 ) - by -  nc (= 6*nf2p1 - 1) system of linear equations, 

!  with unique solution but which has more equations than uknowns. 
!  This is not a problem, as long as we use the general rutine we have described, specifying the kernel dimension to be zero

! TODO: why the kernel dimension is zero? 

!   we have to 
!   The linear equations to be solved:
!                  dg *  dX = - ferr
! 
!  1.  the whole diffrential dg: D F/ d A (C_0, C's, S's) 
!  2.  the error: ferr,  to be refined by Newton method 

!  -- the Linear equation F to slove for Newton Method : 
!     F =    [ H( \varphi(xi_0)) - h = 0;                                   ! energy 
!         F1: \varphi(xi_i + rho) - P (xi_i) = 0 ],  i = 1,..., Nf    ! the invariance condition 

!       Input Variables 
!  h        energy level to fix
!  rho      the prescribed rotation number           

!       Input-Output Variables 
! c0, ck, sk:  All the Fourier coefficients to be refined 
      
!       Output Variables 
!  isref    flag to show if the refinement is successful or not           

! Finally Revised by Yu -- 20160826
!----------------------------------------
subroutine refine_curve_poinc( h, rho, c0, ck, sk, isref, deriv, gr_cj, deriv_cj)

implicit none

! Input  and Output Declaration   
real(kind=dp), intent(in)       ::  h  ! the prescribed energy 
real(kind=dp), intent(inout)    ::  rho, c0(ndim),  ck(ndim,nf), sk(ndim,nf) 
integer, intent(out)     ::  isref  
  
! Local Variable
integer        :: ind_cs1,  is_ck, & ! the index of unknowns C1 to eliminate 
                  iter, i, j, info, num_incr, nrhs,  &
                  ncol, k, kcol, kcol_end  
                       
real(kind=dp)  :: dg(nr_dg, nc_dg), ferr(nr_dg), dx(nc_dg), & ! Netwon method 
                  dc0(ndim), dck(ndim,nf), dsk(ndim,nf),  & ! correction  
                  dxm_pre, dxm, dferrm, &   ! termination control
                  dg_all(nr_dg, nc_dg+1)
                    
real(kind=dp), external :: dnrm2 
external :: deriv, gr_cj, deriv_cj, poinc, deltx   

  ! it seems we never enter this routine .... why???
  
!  print*, 'in refine_curve_poinc, ndim=', ndim 
!  read* 
  
  nrhs = 1  ! only one column in right-hand side 
   
! refine iteratively 
  iter = 0          ! index of iteration 
  num_incr = 0      ! number of times when the modulus correction increases

  !  TODO: note that to call gr_lf and gr_cjlf, we need to call lf_mod, and initialize in the main routine 
  
  ! we do not need to compute ind_cs1 for every iteration, so compute it before the do loop   
  ! subroutine fix_index( ck, sk, ind_cs1, is_ck)
  call fix_index(ck, sk, ind_cs1, is_ck) 
  
!  print*, 'ind=', ind,  'ind_cs1=', ind_cs1
!  print*, 'if ind=ind_cs1  there is something wrong!'
!  read*
!  
  ! fix the corresponding coefficient to 0, and do not update during the Newton method 
  if(is_ck == 1) then 
    ck(ind_cs1, 1) = 0.d0
  else 
    sk(ind_cs1, 1) = 0.d0
  endif 
  
  do 
    ! initialize all the correcitons to be zeros  
    dc0 = 0.d0
    dck = 0.d0
    dsk = 0.d0
    
    iter = iter + 1 
    
    if(iter > nitmax) then 
      write(fwork,*) 'Too many iterations! ', iter
      isref = 0 
      return
    endif 
    
    write(fwork,*) '***********', iter, 't-th iteration ***************'
    
    ! -- differential + error - 
    ! subroutine  gdg_curve_poinc( h, rho, c0, ck, sk, ferr_all, dg_all, deriv, gr_cj)
    call  gdg_curve_poinc( h, rho,  c0, ck, sk,  ferr, dg_all, deriv, gr_cj, deriv_cj)
    
    ! -- extract the corresponding sub-matrix from ferr_all, and dg_all 
    
    !    we fix the coordinates of  c1_ind (is_ck=1) or s1_ind (is_ck=0)
    !    the actual Jacobi Matrix should be dg_all eliminate : 
    
    !     1.  if is_ck==1,  the 2 + (ind-1)*nf2p1      -th column  -- set to zero
    !            is_ck==0,  the 2 + nf + (ind-1)*nf2p1 -th column  
    
    !     where ind is the component of pv for poincare section pv(ind) = p0, 
    !     the block (c0,ck,sk) starts from  1+(ind-1)*nf2p1, and has nf2p1 elements
    !     and c1 (s1) is the second (third) component in each block, index is 1+(ind-1)*nf2p1 + 1 (1+nf)
    ! 
    !  -- eliminate the column related to c1_{ind_cs1} if is_ck = 1 from the unknowns 
    !     and fix the value to be 0 
    
    ! when we update the Fourier coefficients, we need to do the opposite for the corresponding column
    ! if is_ck = 1, keep s1_{ind_cs1};  elseif is_ck = 0,  keep c1_{ind_cs1} 
    
    ! we need to copy the all the columns from  1 to  1+(ind_cs1-1)*nf2p1 
    ! if is_ck = 1, we eliminate c1_{ind_cs1},  copy the columns before and after this one, 
    ! just to skip this very column 
    ! 
  
    ! index of the first column for ind-th block - 1
    ncol = (ind_cs1 - 1) * nf2p1 + 1
    
    ! copy all the column before c1_{ind_cs1}
    if(is_ck == 1) then 
      dg(:, 1 : ncol)        = dg_all(:, 1 : ncol)
      dg(:, 1+ncol : nc_dg)  = dg_all(:, 2+ncol : nc_dg+1)
         
    else ! (is_ck == 0), eliminate the column numbered ncol+nf+1 
      dg(:, 1 : ncol+nf)       = dg_all(:, 1 : ncol+nf)
      dg(:, 1+ncol+nf : nc_dg) = dg_all(:, 2+ncol+nf : nc_dg+1)
    endif 
    
    ! use dxm_pre to detect the modulus of the error, in principle, it should decrease 
    ! we allow it increases for two times, if it exceeds, the refinement fails 
    if( iter > 1)  dxm_pre = dxm
    print*, 'nr_dg, nc_dg, ndim, nf', nr_dg, nc_dg, ndim, nf
    read*
!    subroutine deltx( nr, nc, nrhs, a, b, x, info)
     print*, 'dg '
    do i = 1, nr_dg, 1
      write(*,'(12f18.8)') dg(i, :)
    end do
    print*; read*
    
    call deltx( nr_dg, nc_dg, nrhs, dg, ferr, dx, info)

    write(*, '(12f18.10)') dx 
    read* 
    
    if( info /= 0 ) then 
      write(fwork,*)  'SSL Solver Failed to converge!';  read*
      isref = 0
      return
      
    else   
 
    !  Terminate the iteration, if the change trend of the modulus of the dx fluctuates, stop
    !  in principle, the modulus of dx is decreasing, we allow it increases twice as suggested by Gerard 

      dxm    = dnrm2(nc_dg, dx,   1)
      dferrm = dnrm2(nr_dg, ferr, 1)
      
      print*, 'dxm=',dxm, 'derr=', dferrm, 'tol=', tol;  read*; read* !ck
      
      !TODO: if ferr is small enough, stop the iteration and return also.....
      if( dferrm < tol_err ) then
        print*, 'Error less than tol_err', dferrm, tol_err
        isref = 1
        return  
      end if
      
      if(iter > 1 .and. dxm > dxm_pre) then 
      
        num_incr = num_incr + 1 
        
        if(num_incr >= 3 ) then ! allow the correction  dx to increase no more than twice... 
          isref = 0
          print*, 'Terminate refinement! |dx| starts increasing!'; read*
          return 
        endif 
          
      endif  

      if(dxm < tol) then 
        write(fwork,*) 'Modulus of correction less than tol! Succeed!', dxm 
        return
      endif  
    
      ! -- print the old values of the unknowns to check
      !    note that c0_x is fixed.... 
      write(fwork,*) 'Old: X0 -- C0  \\  ck, | sk  for \\ (k = 1, ', nf , ')'
      
      ! write the coefficients C's and S's into file fcs.dat 
      write(fwork, *)  c0 
      
      do k = 1, nf, 1
        write(fwork, *)  ck(:, k), sk(:, k) 
      enddo
      
      ! -- print the correction to check
      write(fwork,*) 'Correction: dX --  dc0 \\  dck, | dsk \\ (k = 1, ', nf , ')'
      
      
      ! -- Update C's and S's
      ! ** Note: one coordinate of C1 is fixed (of ind ) -- the same rule as dg 
      !          so do not update the corresponding coefficient, but should we put it to 0? when? 
      ! if is_ck = 1, keep s1_ind; elseif is_ck = 0,  keep c1_ind 
    
      do k = 1, ndim, 1
        ! the start index of k-th block (component)
        ncol = 1 + (k - 1) * nf2p1 
        
        print*, 'ncol=', ncol, 'for', k, '-th component'
        
        if(k == ind_cs1) then  ! we skip the c1 or s1
          
          if(is_ck == 1) then 
            dck(k, 1) = 0.d0
            dck(k, 2:nf) = dx(1+ncol  : ncol+nf-1)
            
            dsk(k, :) = dx(ncol+nf : ncol+2*nf-1)
          else  
          
            dck(k, :) = dx(1+ncol : ncol+nf)
            
            dsk(k, 1) = 0.d0
            dsk(k, 2:nf) = dx(1+ncol+nf : ncol+2*nf-1)
          endif   
               
        else
          ! for the components after ind_cs1-th one, copy the whole block
          if(k > ind_cs1)  ncol = ncol -1  
          
          ! this part is the same for the blocks not equal to ind_cs1
          dck(k, :) = dx(ncol+1:  ncol+nf)
          dsk(k, :) = dx(ncol+1+nf:  ncol+2*nf)
        endif 
        
        ! update C0  
        dc0(k) =  dx(ncol)
        
      end do

      ! --- update the unknowns: rho + t2 + c0 + ck + sk
      c0  = c0 + dc0
      ck  = ck + dck
      sk  = sk + dsk 
        
      ! --- write the correcitons in the working file to check ?? TODO 
               ! write the coefficients C's and S's into file fcs.dat 
      write(fwork, *)   dc0 
      do k = 1, nf, 1
        write(fwork, *)  dck(:,k), dsk(:, k) 
      enddo
        
    endif   ! info == 1
   
  enddo  ! iteration   

  return  
end subroutine refine_curve_poinc

!***********************************************************************
!     ****   fix_index   ****
!  Specify the index of unknowns to be fixed to avoid indeterminations 
  
!  There are two indeterminations:    2--phase shift,
!  As suggested by Josep-Maria, we could fix one coordinate of C0 and C1, 
!  with index as ind0 and ind1, respectively 
    
! -- as a consequence, we have  ndim*(2*nf+1) - 1 =  ndim*nf2p1- 1 unknowns ...
!    but we still have   1 + ndim*(2*nf+1) equations.... 
    
!** Comment by Josep-Maria: Advanced course on long time integrations, P59, Title: Computation of a torus
!   This is not a problem, as long as we use the general rutine we have described, specifying the kernel dimension to be zero.

!-- PhD thesis, chapter 6, Methodolody, Section 6.2.5,   P100 

!-- Strategy to eliminate unknowns to avoid indeterminations

!   - 2. To avoid the indetermination in phase shift, mostly for continuation cosideration
!       2.1. Eliminate C1_J or C1_J, where J is chosen in order to have || (C1_J, S1_J) ||_2  =  max_{i=1,..., ndim} || (C1_i, S1_i) ||_2
!            The reason for this choice is heuristic and tries to prevent asking conditions on harmonics of low amplitudes.

!       2.2  Between C1_J and B1_J, we select the coordinate for which min(|C1_J|, |S1_J| ) is achieved. This avoids the situation:
!            assume before the correction, C1_J=0 and |S1_J| = alpha. If we eliminate S1_J, we are forcing  
!            || (C1_J, S1_J) ||_2 > alpha and therefore we will not be able to continue if || (C1_J, S1_J) ||_2 > decreases
!            along the family. If |C1_J| is different from zero but small, we will be able to continue but the continuation step
!            will be artificially reduced.  

!** Another strategy is: Set C(k)_J = 0 as long as S(k)_J remains different from zero, that is  ( C(k)_J, S(k)_J) \= 0 
!   Usually , we take k=1 (the x component) and J has been chosen in order to maximize || ( C(k)_J, S(k)_J) ||_2. 

!   Comment: the first strategy is better, because it doesn't depend on a Fourier coefficients being different from zero 
!   and C1_1 = C1_x = 0, we need to introduce 2 variables as the indice of unknowns to fix the value 
    

!       Input Variables 
!  ck,sk    see gdg_curve 

!       Output Variables 
!  ind      the index of component in C1 to be fixed as zero.             

!       Module-based Varaibles
!  ndim, nf   for purpose of size declaration

!  routine Used:
!    none  

!  Finally Revised by Yu -- 20160731
!***********************************************************************
subroutine fix_index(ck, sk, ind_cs1, is_ck)

implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration 
!integer, intent(in)            ::  nf
integer, intent(out)           ::  ind_cs1, is_ck   
real(kind=dp), intent(in)      ::  ck(ndim, nf), sk(ndim, nf)  
 
! Local Variable
integer :: i, ind_loc(1)
real(kind=dp)  :: nm(ndim)
  
  ! evaluate the coefficients of the x componet || (C1_J, S1_J) ||_2 for J = 1, ..., nf 
  do i = 1, ndim, 1
    nm(i) = dsqrt( ck(i, 1)**2 + sk(i, 1)**2 )
  end do  
  write(*,*) nm
  print*, 'ind=', ind 
  
  nm(ind) = 0
  
  ! take J that achieves  max { || (C1_J, S1_J) ||_2 } 
  ind_loc = maxloc(nm)
  
  ! to extract the value of ind_loc 
  ind_cs1 = ind_loc(1)
  
  
  print*, 'ind_cs1 =',  ind_cs1
  read* 
  
  ! between C1_J and S1_J, select the one that achieves min ( |C1_J|, |S1_J|)
  if( dabs( ck(ind_cs1, 1) ) .le. dabs( sk(ind_cs1, 1) ) ) then 
    is_ck = 1
  else 
    is_ck = 0
  endif     
  
  return  
end subroutine fix_index

!*********************************************************************** 
!  This routine computes the coefficient matrix and error for the linear equations to be solved
!  by Newton Method, as input of deltx_lnx 
!                  dg *  dX = - ferr
!                 
!  For the nf values of xi_i, i = 0, ..., 2*nf 
!  
!  1.  the whole diffrential dg: D F/ d   A(C_0, C's, S's)  
!  2.  the error: ferr,  to be refined by Newton method 

! where the matrix A include components T,y,vx,vy
! the Fourier coefficients of the return time T also needs to be refined
! to save memory, we store them in the column ind (the Poincare section) 


!  -- the Linear equation F to slove for Newton Method : 
!     F = [ H( \varphi(xi_0)) - h = 0;                                   ! energy 
!       F1: \varphi(xi_i + rho) - P (xi_i) = 0 ],  i = 1,..., Nf   ! the invariance condition  
!          Note that the invariance conditon is for all the  ndim components  
!          the one chosen to be Poincare section is replace by the return time T 
!          For instantce, if ind =1, we have (T, y, vx, vy)
  
!  the unknowns are written in column vector, of dimension 3 + ndim *(2*Nf+1) :
!  ---  X = [  h,  A_T (A_x),   A_y, ... , A_vz ]  if ind=1... 
  
! the component of dg =     
!!     -- the first row,  d H /  d C_0    !, where d H/ dh = -1  -- discard for the moment
  
!     -- the rest row,   d F1 / d A(C_0, C's, S's) 

!        for each xi, compute by subroutine gdg_inv_xi1
!        do not forget:  -d P( \varphi(xi_i) ) / d A  
  
! -- 2. TO make everything clear, we compute the full matrix to aviod confusion
!       Do the substraction of corresponding submatrix in the routine refine_curve_poinc

!       Input Variables 
!  h            the presribed vaule of energy 
! rho           rotation number: 2pi / omega2 * omega1
! c0, ck, sk   the initial unkowns to be refined, in total is of dimension 6 * (2*nf+1)


!       Output Variables 
! ferr_all    the error in F to be refined, dimension:   1 + 6*(2*nf+1) 
! dg_all       the differential of F w.r.t. the unkowns, 
!             dimension:  [ 1 + 6*(2*nf+1) ] - by -  [ 3 + 6*(2*nf+1) ]        
! 
!        Module-Based Varaibles
! ndim, nf, nf, nf2p1,  pi2 
! nr_dg = 6 * ( 2*nf+1 ) + 1, nc_dg = 6 * ( 2*nf+1 ) 


!  routine Used: -- TODO 
!     None

!  Finally Revised by Yu -- 20160731
! ---- TODO : check 20160826
!------------------------------------------------------------------------
subroutine  gdg_curve_poinc( h, rho, c0, ck, sk, ferr_all, dg_all, deriv, gr_cj, deriv_cj)

implicit none

! Input  and Output Declaration   
real(kind=dp), intent(in)    ::  h, rho,  c0(ndim), ck(ndim, nf), sk(ndim, nf) 
real(kind=dp), intent(out)   ::  ferr_all(nr_dg),  dg_all(nr_dg, nc_dg+1) 
 
! Local Variable
integer        :: i, k, irow  
real(kind=dp)  :: pv0(ndim), dxi, xi, dg_xi1(ndim, ndim*nf2p1),  ferr_xi1(ndim), &
                  h0, dh(ndim) 

external :: deriv, gr_cj, deriv_cj
  
  ! Initialize as zeros, then we only need to deal with nonzero components 
  dg_all = 0.d0 
  
  ! ---  The Jacobi Matrix of the general function F ----
 
  ! -- the first row of F is  H( \varphi(xi_0)) - h = 0 
  ! -- We should put the energy equation outside this routine 
  !    here, since it is only dependent on the constant term C0, different values of xi has no effect. 
  !    as well as the differential of H w.r.t. C0 
  !    where xi_0   = C_0 ([x,y,z,vx,vy,vz])
  
  ! -- d H / d C_0 ^(x,y,z, vx,vy,vz)
  !    the index of C_0^(x,y,z, vx,vy,vz) :    4 + i* (2*nf+1), where i = 0, 5 
  
  !    C_0 =  H( \varphi(0)) ) is obtained  by substituding xi=0 into the truncated Fourier Series \varphi
  !    which is the actual value of (x,y,z, vx, vy,vz)_0 to compute the energy h
  !    however, c0(ind) is replaced by c0_T, so we need to introduce a new array 
  pv0 = c0 
  pv0(ind) = p0
  
  call deriv_cj(pv0, dh)
  
  ! replace C0 of the ind-th block by 0.d0, since h is independent on the return time T 
  dh(ind) = 0.d0
  ! the index of each component of C0, it depends on how we store the unknows, 1, 1+nf2p1, 1+2*nf2p1
  dg_all(1, 1: nc_dg : nf2p1)  =  dh  
  
  print*,  'dh=', dh 
  read* 
  
  !---- H0 - h  as the firt equation  
  call gr_cj(pv0, h0)
  ferr_all(ndim+1) = h0 - h
  
!  ----- deal with argument xi_i one by one  ---------   

! time step size for discretisize the parameter space, dxi = 2*pi/(1+2*nf), 
! and we need  to evaluate xi_i = i * dxi, where i = 0, ..., 2*nf 

  ! step size for xi 
  dxi   = pi2 / nf2p1
  
  do i = 1, nf2p1, 1
  
    ! the current xi to be evaluated
    xi = (i-1) * dxi
    
    ! the starting row for each block of xi_i, i = 0, ..., 2*nf 
    irow = 2 + (i-1)*ndim
    
    print*, 'i=', i, 'xi=', xi, 'starting row:', irow 
    print*, 'check rho=', rho 
    read*
        
    !  subroutine gdg_inv_xi1( xi, rho, c0, ck, sk, dg_xi1, ferr_xi1, deriv, gr_cj)
    call gdg_inv_xi1( xi, rho, c0, ck, sk, dg_xi1, ferr_xi1, deriv, gr_cj)
    
    print*, 'ferr_xi1 = ', ferr_xi1 
    print*, 'dg_xi1 = '
    do k =  1, ndim, 1
      write(*, '(12f18.8)') dg_xi1(k,:)
    end do  
      
    read* 
    
    ! the sub-block in dg associated with one value of xi, ndim rows starting from irow
    dg_all(irow : irow+ndim-1, :) = dg_xi1
    
    ! error in F (the invariance equations)
    ferr_all(irow : irow+ndim-1 ) =  ferr_xi1 
    
  end do 
  
  return 
end subroutine gdg_curve_poinc


!***********************************************************************
!     ****   dvphi_da   ****
!   for a certain value of xi, the Jacobi Matrix of the invariant curve \varphi 
!   w.r.t. the  Fourier coefficients A.  
!   
!   so dvphi_da is a matrix with ndim same block in diagonal line 
!   we store by Fourier modes, k = 0, 1, ..., nf 

!   each block is a row vector of dimension 1 + 2Nf: 
!   1  cos(xi) sin(xi) cos(2*xi) sin(2*x) ... 

!   we denote A = [ C0, C's, S's] as the Fourier coefficients, then A is of dimenson 6 * (2*Nf+1), write as components A = [A_x, A_y, A_z, A_vx, A_vy, A_vz] 
  
!  where \varphi as a truncated Fourier series.
!        \varphi (xi) = C_0 + Sigma_ k from 1 to Nf of { ( C_k cos(k*xi) + S_k cos(k*xi)) }
  
!       Input Variables 
!  xi       the argument for the invariant curve \varphi 
                    
!       Output Variables 
!  dgda      the differential (Jacobi Matrix) of \varphi w.r.t. the Fourier coefficients A            
 
!       Module-based Varaibles
! nf, nf2, nf2p1

!  routine Used:
!   trigrec_c, trigrec_s 

!  Finally Revised by Yu -- 20160729
!***********************************************************************
subroutine dvphi_da_cal( xi, dgda)
implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration   
real(kind=dp), intent(in)      ::  xi  
real(kind=dp), intent(out)     ::  dgda(ndim, ndim*nf2p1) 
 
! Local Variable
integer        :: i, j, k, jcol
real(kind=dp)  :: cn(nf), sn(nf)

  print*;  print*, '-------- check dvphi/da ----------'
  
  print*, 'xi = ', xi
  dgda = 0.d0 ! initialization
  
 ! compute the trigometric function  xi 
  call trigrec_c( xi, nf, cn)  
  call trigrec_s( xi, nf, sn)
   
  
  ! The Jacobi Matrix of the general function F (the ndim coordinates)  w.r.t. the  Foureier coefficients A 
  !  each is of dimension ndim * ( 2*Nf+1 ) and are of the same block in the diagonal line, the rest components are zeros
  
  do i = 1, ndim, 1
  
     ! the starting  index of i-th block  for i-th component
     jcol = (i-1)*nf2p1 + 1  
     
     ! C_0
     dgda( i,  jcol) = 1.d0       
     
     ! the 2*nf components are cos(k*xi) for k =1, ...,  nf ,  sin(k*xi) for k =1, ...,  nf 
     dgda( i,  jcol+1 : jcol+nf)       =  cn  ! cos( k*xi ) 
     dgda( i,  jcol+nf+1 : jcol+2*nf)  =  sn  ! sin( k*xi )
     
     print*, 
     print*, 'cn = '
     write(*, '(12f18.8)') cn
     
     print*, 'sn='
     write(*, '(12f18.8)') sn 
     
     print*, 'dgda = '
     write(*, '(12f18.8)') dgda(i, :)
     read*
  end do
 
  return  
end subroutine dvphi_da_cal


!***********************************************************************
!     ****   gdg_inv_xi1  ****
!  for a certain value of xi, the Jacobi Matrix of the invariance equation 
!  we fix the energy level here, so h is eliminated from the unknowns 

!  **NOTE** 
!    F1:  
!    \varphi(xi + rho) - P ( \varphi(xi) ) = 0,  i = 1,..., Nf, for components (y,vx,vy)
!                  T(xi)  - P_T            = 0 
!-- where P_T is the return time , 
!   the C_T and S_T are computed in the first iteration   
!   so  the refinement in the first iteration only to improve the precision of the Fourier anay


!  the unknowns are written as a column vector, of dimension ndim *(2*Nf+1):
!  ---  X = [  A_T,  A_y, ... , A_vy ] (ind=1 as example) 

!   we denote A = [ C0, C's, S's] as the Fourier coefficients, then A is of dimenson ndim * (2*Nf+1), 
!   A = [A_x, A_y, A_z, A_vx, A_vy, A_vz], and ind-th column is replaced by the derivative w.r.t. the return time T in Poincare map method 
    
!  the component of dg =
!  -- the ind-th  column:   d F1 / d A_T ,     where A_T = ( C0_T, CK_T,  SK_T)  
!     since    d \varphi( xi+ rho ) / d T  = 0, independent on T, applying chain rule, 
!     we have  d F1 / d A_T = - d P / d T  * d T / d A_T
!     where  d T / d A_T is one of the ndim same blocks in dvphi_da, so take dvphi_da(1:nf2p1) instead
!               
  
!      d P / d T is the  vector field at  point P( \varphi(xi_i)  ) 
!       (since the flow is the solution of the vector field, and the vector field is the derivative of the flow w.r.t. time )
!    
!   -- d F1 / d A_(y,vx,vy) ---  ( as a block of dimension 2*Nf + 1) 
!      we also apply the chain rule to compute  d P / d A_(y,vx,vy) 

! Finally, we have: 
!    d \varphi(xi+rho) / d A_(y,vx,vy) - [ d P / d \varphi(xi) ] * d \varphi(xi) / d A_(y,vx,vy)
!            
!
! where   \varphi is a truncated Fourier series evaluated at xi 
!         \varphi (xi+rho) is the same Fourier series evaluated at xi+rho 
  
!       Input Variables  (initial guess)
!  xi           the argument for the invariant curve  \varphi 
!  rho          rotation number:  2pi / omega2 * omega1
!  c0, ck, sk   dimension ndim-by-nf,  the unknowns C_k and S_k (k=1, ..., nf), for each component
!               (x,y,z,vx,vy,vz), the values are different. 

                    
!       Output Variables 
!  dg_xi1      dimension: ndim(6 or 4), ndim * ( 2*nf+1 ) 
!              the differential (Jacobi Matrix) of invariance equations w.r.t. A_{T, y, vx, vy}
 
! ferr_xi1     dimesion ndim, the error to be refined by Netwon method 
!              F1      = \varphi(xi_i + rho) - P (xi_i)  
!              F1(ind) = T(xi) - P_T

!     Module based Varaibles
! ndim, nf      dimension of the phase space, number of fourier modes
! for poinc:  ind, p0, dir, imax

!  routine Used:
!   trigrec_c, trigrec_s 

!  Finally Revised by Yu -- 20160826
!***********************************************************************
subroutine gdg_inv_xi1( xi, rho, c0, ck, sk, dg_xi1, ferr_xi1, deriv, gr_cj)
implicit none

! Input  and Output Declaration   
!integer, intent(in)            ::  nf 
real(kind=dp), intent(in)      ::  xi, rho, c0(ndim), ck(ndim, nf), sk(ndim, nf) 
                                   
real(kind=dp), intent(out)     ::  dg_xi1(ndim, ndim*nf2p1), ferr_xi1(ndim) 
 
! Local Variable
integer        :: i, k, npvar, jcol_st, jcol_end, para(ndim), ispc   
real(kind=dp)  :: cn(nf2), sn(nf2), cn_rho(nf2), sn_rho(nf2), &
                  pv0(ndim), vf0(ndim), pv0_rho(ndim),  &  
                  pci( ndim*(ndim+1) ), pcf( ndim*(ndim+1) ), pv_pc(ndim), stm(ndim,ndim), & !  poinc 
                  dp_dx(ndim, ndim), dp_da(ndim, ndim*nf2p1),  hminim, tf,  &  ! poinc
                  dvphi_da(ndim, ndim*nf2p1), dvphi_rho_da(ndim, ndim*nf2p1), vf_pc(ndim)  ! for dg_xi1

!  External subroutines                  
external :: deriv, gr_cj 
  
  dg_xi1 = 0.d0 ! initialization
  
 ! the trigometric function at k* (xi+rho) for k = 1_2Nf
  call trigrec_c( xi+rho, nf2, cn_rho)  
  call trigrec_s( xi+rho, nf2, sn_rho)
 
  ! ---- The Jacobi Matrix of the general function F (the ndim coordinates) w.r.t. A---
  
  ! Evaluate the invariant curve \varphi(xi), with the initial guess of coefficients (c0, ck, sk)
  ! TODO: \varphi(xi) is given value or not... I think it's not, we need to recompute for each pv0 
  call varphi( xi, c0, ck, sk, pv0) 
  
  ! start from pv0, compute the first return to the Poincare section 
  ! we need to compute the variational matrix, so the size of initial state is  n = ndim*(ndim+1)
  
! initialize also the variational matrix
  npvar = ndim*(ndim+1)
  
  print*, 'npvar =', npvar, 'ndim=', ndim  !ck
  read*
  
  pci = 0.d0
  pci(1:ndim) = pv0 
  pci(ind) = p0 
  pci(ndim+1 : npvar : ndim+1) = 1.d0 
  
  print*, 'initial state befor poinc '
  write(*, '(4f18.8)') pci 
  read* 
  
!subroutine poinc(sti, ndim, n, ind, p0, dir, imax, tmax, ispl, fob, tf,stf, hminim, ispc,  deriv, gr_cj) 
  ! without plotting the orbit, ispl=0
  call poinc(pci, ndim, npvar, ind, p0,  dir, imax, tmax, 0, 6, tf, pcf, hminim, ispc, deriv, gr_cj)
  
  print*, 'final state after poinc '
  print*, 'tf = ', tf 
  write(*, '(4f18.8)') pcf 
  read* 
  
  pv_pc = pcf(1:ndim) ! the final state 
  call deriv(0.d0, pv_pc,  ndim, vf_pc)  
  
  print*, 'pv_pc = ', pv_pc 
  print*, 'vf_pc = ', vf_pc 
  read*  
  
  stm = reshape( pcf(ndim+1:npvar), (/ndim,ndim/) ) ! the STM 
  
  ! differential of Poincare map w.r.t. the initial state 
  ! in order to make it more general, compute all the components of dpdx0
  ! and the replace the ind-th column by d p / d A_T
  
  !  subroutine diffpoinc( phi, f,  ndim, ind, nr, nc, para, fun, dpdx)
  para(1:ndim) = (/ (i, i = 1, ndim) /) ! implied do loop 
  call diffpoinc( stm, vf_pc,  ndim, ind, ndim, ndim, para, para, dp_dx)
  
  ! since A0(ind) is fixed, we should assign the ind-th column of dp_dx to zeros 
  dp_dx(:, ind) = 0.d0 
  if(ispc == 0) then 
      print*, 'Poinc failed!'  
      read*  
  endif 
                 
  ! \varphi( xi+rho ) = pv0_rho
  call varphi( xi+rho, c0, ck, sk, pv0_rho) ! the truncated Fourier series 
  
  ! the error to be corrected by Newton Method:  \varphi(xi+rho) - P( \varphi(xi) ) = 0
  ! for the components (y, vx, vy)
  ! first deal with the columns that are related to A_T for the phase components 
  
  print*, 'pv0 =', pv0
  print*, 'pv0_rho = ', pv0_rho 

  
  ferr_xi1 = pv0_rho - pv_pc  
  
  ! then deal with the rows and columns related to T(xi) - P_T = 0  
  ferr_xi1(ind) = pv0(ind) - tf 
   
  print*, 'ferr_xi1 = ', ferr_xi1 
  read* 
    
  !  d \varphi( xi+rho ) / d A 
  call dvphi_da_cal( xi+rho, dvphi_rho_da)
  
  !  --- compute d \phi / dA  ---------
  !  d \varphi(xi) / d A 
  call dvphi_da_cal( xi, dvphi_da)
  
  ! the ind-th block that related to T should be set to zeros 
  ! d P / dA  =  d pv_pc{y,vx,vy} / d \varphi(xi) * d \varphi(xi) / d A
  !           = STM * dvphi_da
  
  dp_da = matmul(dp_dx, dvphi_da) 
  
  ! here, we need to replace the ind-th column of dp_da, which is related to d T
!  ncol = (ind-1)*nf2p1 + 1 ! starting column for T
!  dp_da(:, ncol : ncol + nf) =   
!  
  dg_xi1 = dvphi_rho_da - dp_da 
  
  
  print*, 'dp_da = dp_dx X  dvphi_da' 
  do i = 1, ndim, 1
    write(*, '(12f18.8)')  dp_da(i, :)
    print*; read*
  enddo 
  
  print*; read*
  
  print*, 'dg_dx1' 
  do i = 1, ndim, 1
    write(*, '(12f18.8)')  dg_xi1(i, :)
    print*; read*
  enddo 
  
  print*; read*
  
      
  ! --- ind-th column of d F1 / d A be replaced by d F1 / d A_T = - d P ( \varphi(xi) ) / d A_T
  !     which is :  - d P / d T * d T / d A_T  
  !     where d P / d T is  the vector field at pv_pc 
 
  ! ndim components (x,y,z, vx, vy, vz), where ind-th component is replaced by T(xi)
  do i = 1, ndim, 1  
     
   ! d F1 / d A(x,y,z,vx,vy,vz):   ndim rows
   ! each is of dimension ndim * ( 2*Nf+1 ), the same block in the diagonal line, the rest are zeros
   
    ! the start index and end index of the column for ind-th block 
    jcol_st  = (ind-1)* nf2p1 + 1 
    jcol_end = jcol_st + nf2
   
    print*, 'i,  jcol_st, jcol_end ', i, jcol_st, jcol_end   ! ck 
    read*  !ck
     
    if( i /= ind) then 
   
    ! -- replace  the A_T- related  elements by d F1 / d A_T ---- 

    ! we apply the chain rule to get:  d F1 / d A_T =  0 -  d P / d T * d T / d A_T  
  
    ! d P / d A_T =  -vf_pc *  d T / d A 
    ! where  d T / d A  is the ind-th block in dvphi_da, 
    ! since all the block in dvphi_da are the same,  so we use the first block, 
    ! which is the first nf2p1 elements in the first row  
    
      dg_xi1( i, jcol_st : jcol_end ) = -vf_pc(i) *  dvphi_da(1, 1:nf2p1)       
      print*, 'i /= ind,  d G / d A_T',  dg_xi1( i, jcol_st : jcol_end )
      read* 
      
    else 
      
      ! for the ind-row, it is the differential of equation T(xi) - P_T = 0
      ! P_T is computed by poinc, and T(xi) is represented by Truncated Fourier series
      ! with the coefficients stored in c0(ind), ck(ind,:), sk(ind,:)
      
      ! the non-zero d T / d A_T  is one of the block in dvphi_da, so this part is easy 
      
      print*, 'i = ind, deal with the return time T'; read* 
      dg_xi1(i, :) = 0.d0
      dg_xi1(i, jcol_st : jcol_end) = dvphi_da(1, 1:nf2p1)      
    endif
     
     ! -- ck --
     print*, 'i = ', i 
     print*, 'dvphi_rho_da' 
     write(*, '(12f18.8)')  dvphi_rho_da(i, :)
     
     print*; print*, 'dg_xi1' 
     write(*, '(12f18.8)')   dg_xi1(i, :) 
     print*; read* 
       
  end do
 
  return  
end subroutine gdg_inv_xi1 


!***********************************************************************
!     ****   varphi   ****
!  Compute the truncated Fourier series of the invariant curve \varphi(xi) 
!  given the coefficients and initial phases

!       Input Variables 
!  xi           the argument for the invariant curve \varphi 
!  c0, ck, sk   dimension 6-by-nf,  the unknowns C_k and S_k (k=1, ..., nf), for each component
!               (x,y,z,vx,vy,vz), the values are different. 

!       Output Variables 
!  pv           the state vector: position + velocity            

!      Module-based Varaibles
!  ndim, nf  
! 
!  routine Used:
!     trigrec_c, trigrec_s

! TODO: probably we will not need to compute pv using a standalone subroutine 
!       because at the point when pv is needed, we have already computed  cn, sn 

!  Finally Revised by Yu -- 20160730
!***********************************************************************
subroutine varphi( xi, c0, ck, sk, pv)
implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration 
real(kind=dp), intent(in)    :: xi, c0(ndim), ck(ndim, nf), sk(ndim, nf)   
real(kind=dp), intent(out)   :: pv(ndim) 
 
! Local Variable
integer :: i, k 
real(kind=dp)  :: cn(nf), sn(nf)
  
 ! compute the trigometric function at k*xi 
  call trigrec_c( xi, nf, cn)  
  call trigrec_s( xi, nf, sn)
  
  do i = 1, ndim, 1
    pv(i) = c0(i)
    do k = 1, nf, 1
       pv(i) = pv(i) + ck(i,k) * cn(k) + sk(i,k) * cn(k)
    end do
  end do  

  return  
end subroutine varphi


end module curve_poinc_mod
